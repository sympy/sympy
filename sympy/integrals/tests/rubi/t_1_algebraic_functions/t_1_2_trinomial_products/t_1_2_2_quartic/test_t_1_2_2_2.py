"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.2 (d x)^m (a+b x^2+c x^4)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, m, p = symbols('a b c d m p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/4)
    F = 3*sqrt(a)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/4)*asinh(sqrt(b)*x/sqrt(a))/(8*sqrt(b)*(1 + b*x**2/a)**(sympy.S(3)/2)) + 3*a*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/4)/(8*a + 8*b*x**2) + x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/4)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_2():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/4)
    F = sqrt(a)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/4)*asinh(sqrt(b)*x/sqrt(a))/(2*sqrt(b)*sqrt(1 + b*x**2/a)) + x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/4)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_3():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-1)/4)
    F = sqrt(a)*sqrt(1 + b*x**2/a)*asinh(sqrt(b)*x/sqrt(a))/(sqrt(b)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_4():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-3)/4)
    F = x*(a + b*x**2)/(a*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_5():
    f = 1/(a**2 + 2*a*x**2 + b + x**4)
    F = -sqrt(2)*atan((-sqrt(2)*x + sqrt(-a + sqrt(a**2 + b)))/sqrt(a + sqrt(a**2 + b)))/(4*sqrt(a + sqrt(a**2 + b))*sqrt(a**2 + b)) + sqrt(2)*atan((sqrt(2)*x + sqrt(-a + sqrt(a**2 + b)))/sqrt(a + sqrt(a**2 + b)))/(4*sqrt(a + sqrt(a**2 + b))*sqrt(a**2 + b)) - sqrt(2)*log(x**2 - sqrt(2)*x*sqrt(-a + sqrt(a**2 + b)) + sqrt(a**2 + b))/(8*sqrt(-a + sqrt(a**2 + b))*sqrt(a**2 + b)) + sqrt(2)*log(x**2 + sqrt(2)*x*sqrt(-a + sqrt(a**2 + b)) + sqrt(a**2 + b))/(8*sqrt(-a + sqrt(a**2 + b))*sqrt(a**2 + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_6():
    f = 1/(a**2 + 2*a*x**2 + x**4 - 1)
    F = -atan(x/sqrt(a + 1))/(2*sqrt(a + 1)) - atanh(x/sqrt(1 - a))/(2*sqrt(1 - a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_7():
    f = 1/(a**2 + 2*a*x**2 + x**4 + 1)
    F = -sqrt(2)*atan((-sqrt(2)*x + sqrt(-a + sqrt(a**2 + 1)))/sqrt(a + sqrt(a**2 + 1)))/(4*sqrt(a + sqrt(a**2 + 1))*sqrt(a**2 + 1)) + sqrt(2)*atan((sqrt(2)*x + sqrt(-a + sqrt(a**2 + 1)))/sqrt(a + sqrt(a**2 + 1)))/(4*sqrt(a + sqrt(a**2 + 1))*sqrt(a**2 + 1)) - sqrt(2)*log(x**2 - sqrt(2)*x*sqrt(-a + sqrt(a**2 + 1)) + sqrt(a**2 + 1))/(8*sqrt(-a + sqrt(a**2 + 1))*sqrt(a**2 + 1)) + sqrt(2)*log(x**2 + sqrt(2)*x*sqrt(-a + sqrt(a**2 + 1)) + sqrt(a**2 + 1))/(8*sqrt(-a + sqrt(a**2 + 1))*sqrt(a**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_8():
    f = 1/(x**4 - 5*x**2 + 4)
    F = -atanh(x/2)/6 + atanh(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_9():
    f = 1/(x**4 + 4*x**2 + 3)
    F = atan(x)/2 - sqrt(3)*atan(sqrt(3)*x/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_10():
    f = 1/(x**4 + 5*x**2 + 9)
    F = -log(x**2 - x + 3)/12 + log(x**2 + x + 3)/12 - sqrt(11)*atan(sqrt(11)*(1 - 2*x)/11)/66 + sqrt(11)*atan(sqrt(11)*(2*x + 1)/11)/66
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_11():
    f = 1/(x**4 - x**2 + 1)
    F = -sqrt(3)*log(x**2 - sqrt(3)*x + 1)/12 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/12 + atan(2*x - sqrt(3))/2 + atan(2*x + sqrt(3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_12():
    f = 1/(x**4 + 2*x**2 + 2)
    F = -log(x**2 - x*sqrt(-2 + 2*sqrt(2)) + sqrt(2))/(8*sqrt(-1 + sqrt(2))) + log(x**2 + x*sqrt(-2 + 2*sqrt(2)) + sqrt(2))/(8*sqrt(-1 + sqrt(2))) - sqrt(-1 + sqrt(2))*atan((-2*x + sqrt(-2 + 2*sqrt(2)))/sqrt(2 + 2*sqrt(2)))/4 + sqrt(-1 + sqrt(2))*atan((2*x + sqrt(-2 + 2*sqrt(2)))/sqrt(2 + 2*sqrt(2)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_13():
    f = 1/sqrt(-3*x**4 + 5*x**2 + 2)
    F = elliptic_f(asin(sqrt(2)*x/2), -6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_14():
    f = 1/sqrt(-3*x**4 + 4*x**2 + 2)
    F = sqrt(sympy.S(1)/3 + sqrt(10)/6)*elliptic_f(asin(x*sqrt(-1 + sqrt(10)/2)), sympy.S(-7)/3 - 2*sqrt(10)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_15():
    f = 1/sqrt(-3*x**4 + 3*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(sqrt(6)*x/sqrt(3 + sqrt(33))), sympy.S(-7)/4 - sqrt(33)/4)/sqrt(-3 + sqrt(33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_16():
    f = 1/sqrt(-3*x**4 + 2*x**2 + 2)
    F = elliptic_f(asin(sqrt(3)*x/sqrt(1 + sqrt(7))), sympy.S(-4)/3 - sqrt(7)/3)/sqrt(-1 + sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_17():
    f = 1/sqrt(-3*x**4 + x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(-3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_18():
    f = 1/sqrt(2 - 3*x**4)
    F = 6**(sympy.S(3)/4)*elliptic_f(asin(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), -1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_19():
    f = 1/sqrt(-3*x**4 - x**2 + 2)
    F = sqrt(3)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_20():
    f = 1/sqrt(-3*x**4 - 2*x**2 + 2)
    F = elliptic_f(asin(sqrt(3)*x/sqrt(-1 + sqrt(7))), sympy.S(-4)/3 + sqrt(7)/3)/sqrt(1 + sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_21():
    f = 1/sqrt(-3*x**4 - 3*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(sqrt(6)*x/sqrt(-3 + sqrt(33))), sympy.S(-7)/4 + sqrt(33)/4)/sqrt(3 + sqrt(33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_22():
    f = 1/sqrt(-3*x**4 - 4*x**2 + 2)
    F = sqrt(sympy.S(-1)/3 + sqrt(10)/6)*elliptic_f(asin(x*sqrt(1 + sqrt(10)/2)), sympy.S(-7)/3 + 2*sqrt(10)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_23():
    f = 1/sqrt(-3*x**4 - 5*x**2 + 2)
    F = sqrt(6)*elliptic_f(asin(sqrt(3)*x), sympy.S(-1)/6)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_24():
    f = 1/sqrt(-2*x**4 + 7*x**2 + 3)
    F = sqrt(2)*elliptic_f(asin(2*x/sqrt(7 + sqrt(73))), sympy.S(-61)/12 - 7*sqrt(73)/12)/sqrt(-7 + sqrt(73))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_25():
    f = 1/sqrt(-2*x**4 + 6*x**2 + 3)
    F = sqrt(sympy.S.Half + sqrt(15)/6)*elliptic_f(asin(x*sqrt(-1 + sqrt(15)/3)), -4 - sqrt(15))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_26():
    f = 1/sqrt(-2*x**4 + 5*x**2 + 3)
    F = elliptic_f(asin(sqrt(3)*x/3), -6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_27():
    f = 1/sqrt(-2*x**4 + 4*x**2 + 3)
    F = elliptic_f(asin(sqrt(2)*x/sqrt(2 + sqrt(10))), sympy.S(-7)/3 - 2*sqrt(10)/3)/sqrt(-2 + sqrt(10))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_28():
    f = 1/sqrt(-2*x**4 + 3*x**2 + 3)
    F = sqrt(2)*elliptic_f(asin(2*x/sqrt(3 + sqrt(33))), sympy.S(-7)/4 - sqrt(33)/4)/sqrt(-3 + sqrt(33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_29():
    f = 1/sqrt(-2*x**4 + 2*x**2 + 3)
    F = elliptic_f(asin(sqrt(2)*x/sqrt(1 + sqrt(7))), sympy.S(-4)/3 - sqrt(7)/3)/sqrt(-1 + sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_30():
    f = 1/sqrt(-2*x**4 + x**2 + 3)
    F = sqrt(2)*elliptic_f(asin(sqrt(6)*x/3), sympy.S(-3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_31():
    f = 1/sqrt(3 - 2*x**4)
    F = 6**(sympy.S(3)/4)*elliptic_f(asin(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), -1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_32():
    f = 1/sqrt(-2*x**4 - x**2 + 3)
    F = sqrt(3)*elliptic_f(asin(x), sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_33():
    f = 1/sqrt(-2*x**4 - 2*x**2 + 3)
    F = elliptic_f(asin(sqrt(2)*x/sqrt(-1 + sqrt(7))), sympy.S(-4)/3 + sqrt(7)/3)/sqrt(1 + sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_34():
    f = 1/sqrt(-2*x**4 - 3*x**2 + 3)
    F = sqrt(2)*elliptic_f(asin(2*x/sqrt(-3 + sqrt(33))), sympy.S(-7)/4 + sqrt(33)/4)/sqrt(3 + sqrt(33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_35():
    f = 1/sqrt(-2*x**4 - 4*x**2 + 3)
    F = elliptic_f(asin(sqrt(2)*x/sqrt(-2 + sqrt(10))), sympy.S(-7)/3 + 2*sqrt(10)/3)/sqrt(2 + sqrt(10))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_36():
    f = 1/sqrt(-2*x**4 - 5*x**2 + 3)
    F = sqrt(6)*elliptic_f(asin(sqrt(2)*x), sympy.S(-1)/6)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_37():
    f = 1/sqrt(-2*x**4 - 6*x**2 + 3)
    F = sqrt(sympy.S(-1)/2 + sqrt(15)/6)*elliptic_f(asin(x*sqrt(1 + sqrt(15)/3)), -4 + sqrt(15))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_38():
    f = 1/sqrt(-2*x**4 - 7*x**2 + 3)
    F = sqrt(2)*elliptic_f(asin(2*x/sqrt(-7 + sqrt(73))), sympy.S(-61)/12 + 7*sqrt(73)/12)/sqrt(7 + sqrt(73))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_39():
    f = 1/sqrt(3*x**4 + 5*x**2 - 2)
    F = sqrt(7)*sqrt(x**2 + 2)*sqrt(3*x**2 - 1)*elliptic_f(asin(sqrt(14)*x/(2*sqrt(3*x**2 - 1))), sympy.S(6)/7)/(7*sqrt(3*x**4 + 5*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_40():
    f = 1/sqrt(3*x**4 + 4*x**2 - 2)
    F = 10**(sympy.S(3)/4)*sqrt((-x**2*(2 - sqrt(10)) + 2)/(-x**2*(2 + sqrt(10)) + 2))*sqrt(x**2*(2 + sqrt(10)) - 2)*elliptic_f(asin(2**(sympy.S(3)/4)*5**(sympy.S(1)/4)*x/sqrt(x**2*(2 + sqrt(10)) - 2)), sqrt(10)/10 + sympy.S.Half)/(20*sqrt(3*x**4 + 4*x**2 - 2)*sqrt(1/(-x**2*(2 + sqrt(10)) + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_41():
    f = 1/sqrt(3*x**4 + 3*x**2 - 2)
    F = sqrt(2)*33**(sympy.S(3)/4)*sqrt((-x**2*(3 - sqrt(33)) + 4)/(-x**2*(3 + sqrt(33)) + 4))*sqrt(x**2*(3 + sqrt(33)) - 4)*elliptic_f(asin(sqrt(2)*33**(sympy.S(1)/4)*x/sqrt(x**2*(3 + sqrt(33)) - 4)), sqrt(33)/22 + sympy.S.Half)/(132*sqrt(3*x**4 + 3*x**2 - 2)*sqrt(1/(-x**2*(3 + sqrt(33)) + 4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_42():
    f = 1/sqrt(3*x**4 + 2*x**2 - 2)
    F = 7**(sympy.S(3)/4)*sqrt((-x**2*(1 - sqrt(7)) + 2)/(-x**2*(1 + sqrt(7)) + 2))*sqrt(x**2*(1 + sqrt(7)) - 2)*elliptic_f(asin(sqrt(2)*7**(sympy.S(1)/4)*x/sqrt(x**2*(1 + sqrt(7)) - 2)), sqrt(7)/14 + sympy.S.Half)/(14*sqrt(3*x**4 + 2*x**2 - 2)*sqrt(1/(-x**2*(1 + sqrt(7)) + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_43():
    f = 1/sqrt(3*x**4 + x**2 - 2)
    F = sqrt(5)*sqrt(x**2 + 1)*sqrt(3*x**2 - 2)*elliptic_f(asin(sqrt(5)*x/sqrt(3*x**2 - 2)), sympy.S(3)/5)/(5*sqrt(3*x**4 + x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_44():
    f = 1/sqrt(3*x**4 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((sqrt(6)*x**2 + 2)/(-sqrt(6)*x**2 + 2))*sqrt(sqrt(6)*x**2 - 2)*elliptic_f(asin(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/sqrt(sqrt(6)*x**2 - 2)), sympy.S.Half)/(12*sqrt(3*x**4 - 2)*sqrt(1/(-sqrt(6)*x**2 + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_45():
    f = 1/sqrt(3*x**4 - x**2 - 2)
    F = sqrt(5)*sqrt(x**2 - 1)*sqrt(3*x**2 + 2)*elliptic_f(asin(sqrt(10)*x/(2*sqrt(x**2 - 1))), sympy.S(2)/5)/(5*sqrt(3*x**4 - x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_46():
    f = 1/sqrt(3*x**4 - 2*x**2 - 2)
    F = 7**(sympy.S(3)/4)*sqrt((x**2*(1 + sqrt(7)) + 2)/(x**2*(1 - sqrt(7)) + 2))*sqrt(-x**2*(1 - sqrt(7)) - 2)*elliptic_f(asin(sqrt(2)*7**(sympy.S(1)/4)*x/sqrt(-x**2*(1 - sqrt(7)) - 2)), sympy.S.Half - sqrt(7)/14)/(14*sqrt(3*x**4 - 2*x**2 - 2)*sqrt(1/(x**2*(1 - sqrt(7)) + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_47():
    f = 1/sqrt(3*x**4 - 3*x**2 - 2)
    F = sqrt(2)*33**(sympy.S(3)/4)*sqrt((x**2*(3 + sqrt(33)) + 4)/(x**2*(3 - sqrt(33)) + 4))*sqrt(-x**2*(3 - sqrt(33)) - 4)*elliptic_f(asin(sqrt(2)*33**(sympy.S(1)/4)*x/sqrt(-x**2*(3 - sqrt(33)) - 4)), sympy.S.Half - sqrt(33)/22)/(132*sqrt(3*x**4 - 3*x**2 - 2)*sqrt(1/(x**2*(3 - sqrt(33)) + 4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_48():
    f = 1/sqrt(3*x**4 - 4*x**2 - 2)
    F = 10**(sympy.S(3)/4)*sqrt((x**2*(2 + sqrt(10)) + 2)/(x**2*(2 - sqrt(10)) + 2))*sqrt(-x**2*(2 - sqrt(10)) - 2)*elliptic_f(asin(2**(sympy.S(3)/4)*5**(sympy.S(1)/4)*x/sqrt(-x**2*(2 - sqrt(10)) - 2)), sympy.S.Half - sqrt(10)/10)/(20*sqrt(3*x**4 - 4*x**2 - 2)*sqrt(1/(x**2*(2 - sqrt(10)) + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_49():
    f = 1/sqrt(3*x**4 - 5*x**2 - 2)
    F = sqrt(7)*sqrt(x**2 - 2)*sqrt(3*x**2 + 1)*elliptic_f(asin(sqrt(7)*x/sqrt(x**2 - 2)), sympy.S(1)/7)/(7*sqrt(3*x**4 - 5*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_50():
    f = 1/sqrt(2*x**4 + 7*x**2 - 3)
    F = sqrt(3)*73**(sympy.S(3)/4)*sqrt((-x**2*(7 - sqrt(73)) + 6)/(-x**2*(7 + sqrt(73)) + 6))*sqrt(x**2*(7 + sqrt(73)) - 6)*elliptic_f(asin(sqrt(2)*73**(sympy.S(1)/4)*x/sqrt(x**2*(7 + sqrt(73)) - 6)), 7*sqrt(73)/146 + sympy.S.Half)/(438*sqrt(2*x**4 + 7*x**2 - 3)*sqrt(1/(-x**2*(7 + sqrt(73)) + 6)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_51():
    f = 1/sqrt(2*x**4 + 6*x**2 - 3)
    F = sqrt(2)*3**(sympy.S(1)/4)*5**(sympy.S(3)/4)*sqrt((-x**2*(3 - sqrt(15)) + 3)/(-x**2*(3 + sqrt(15)) + 3))*sqrt(x**2*(3 + sqrt(15)) - 3)*elliptic_f(asin(15**(sympy.S(1)/4)*sqrt(2)*x/sqrt(x**2*(3 + sqrt(15)) - 3)), sqrt(15)/10 + sympy.S.Half)/(30*sqrt(2*x**4 + 6*x**2 - 3)*sqrt(1/(-x**2*(3 + sqrt(15)) + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_52():
    f = 1/sqrt(2*x**4 + 5*x**2 - 3)
    F = sqrt(7)*sqrt(x**2 + 3)*sqrt(2*x**2 - 1)*elliptic_f(asin(sqrt(21)*x/(3*sqrt(2*x**2 - 1))), sympy.S(6)/7)/(7*sqrt(2*x**4 + 5*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_53():
    f = 1/sqrt(2*x**4 + 4*x**2 - 3)
    F = 2**(sympy.S(1)/4)*sqrt(3)*5**(sympy.S(3)/4)*sqrt((-x**2*(2 - sqrt(10)) + 3)/(-x**2*(2 + sqrt(10)) + 3))*sqrt(x**2*(2 + sqrt(10)) - 3)*elliptic_f(asin(2**(sympy.S(3)/4)*5**(sympy.S(1)/4)*x/sqrt(x**2*(2 + sqrt(10)) - 3)), sqrt(10)/10 + sympy.S.Half)/(30*sqrt(2*x**4 + 4*x**2 - 3)*sqrt(1/(-x**2*(2 + sqrt(10)) + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_54():
    f = 1/sqrt(2*x**4 + 3*x**2 - 3)
    F = 11**(sympy.S(3)/4)*3**(sympy.S(1)/4)*sqrt((-x**2*(3 - sqrt(33)) + 6)/(-x**2*(3 + sqrt(33)) + 6))*sqrt(x**2*(3 + sqrt(33)) - 6)*elliptic_f(asin(sqrt(2)*33**(sympy.S(1)/4)*x/sqrt(x**2*(3 + sqrt(33)) - 6)), sqrt(33)/22 + sympy.S.Half)/(66*sqrt(2*x**4 + 3*x**2 - 3)*sqrt(1/(-x**2*(3 + sqrt(33)) + 6)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_55():
    f = 1/sqrt(2*x**4 + 2*x**2 - 3)
    F = sqrt(6)*7**(sympy.S(3)/4)*sqrt((-x**2*(1 - sqrt(7)) + 3)/(-x**2*(1 + sqrt(7)) + 3))*sqrt(x**2*(1 + sqrt(7)) - 3)*elliptic_f(asin(sqrt(2)*7**(sympy.S(1)/4)*x/sqrt(x**2*(1 + sqrt(7)) - 3)), sqrt(7)/14 + sympy.S.Half)/(42*sqrt(2*x**4 + 2*x**2 - 3)*sqrt(1/(-x**2*(1 + sqrt(7)) + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_56():
    f = 1/sqrt(2*x**4 + x**2 - 3)
    F = sqrt(5)*sqrt(x**2 - 1)*sqrt(2*x**2 + 3)*elliptic_f(asin(sqrt(15)*x/(3*sqrt(x**2 - 1))), sympy.S(3)/5)/(5*sqrt(2*x**4 + x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_57():
    f = 1/sqrt(2*x**4 - 3)
    F = 6**(sympy.S(1)/4)*sqrt((sqrt(6)*x**2 + 3)/(-sqrt(6)*x**2 + 3))*sqrt(sqrt(6)*x**2 - 3)*elliptic_f(asin(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/sqrt(sqrt(6)*x**2 - 3)), sympy.S.Half)/(6*sqrt(2*x**4 - 3)*sqrt(1/(-sqrt(6)*x**2 + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_58():
    f = 1/sqrt(2*x**4 - x**2 - 3)
    F = sqrt(5)*sqrt(x**2 + 1)*sqrt(2*x**2 - 3)*elliptic_f(asin(sqrt(5)*x/sqrt(2*x**2 - 3)), sympy.S(2)/5)/(5*sqrt(2*x**4 - x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_59():
    f = 1/sqrt(2*x**4 - 2*x**2 - 3)
    F = sqrt(6)*7**(sympy.S(3)/4)*sqrt((x**2*(1 + sqrt(7)) + 3)/(x**2*(1 - sqrt(7)) + 3))*sqrt(-x**2*(1 - sqrt(7)) - 3)*elliptic_f(asin(sqrt(2)*7**(sympy.S(1)/4)*x/sqrt(-x**2*(1 - sqrt(7)) - 3)), sympy.S.Half - sqrt(7)/14)/(42*sqrt(2*x**4 - 2*x**2 - 3)*sqrt(1/(x**2*(1 - sqrt(7)) + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_60():
    f = 1/sqrt(2*x**4 - 3*x**2 - 3)
    F = 11**(sympy.S(3)/4)*3**(sympy.S(1)/4)*sqrt((x**2*(3 + sqrt(33)) + 6)/(x**2*(3 - sqrt(33)) + 6))*sqrt(-x**2*(3 - sqrt(33)) - 6)*elliptic_f(asin(sqrt(2)*33**(sympy.S(1)/4)*x/sqrt(-x**2*(3 - sqrt(33)) - 6)), sympy.S.Half - sqrt(33)/22)/(66*sqrt(2*x**4 - 3*x**2 - 3)*sqrt(1/(x**2*(3 - sqrt(33)) + 6)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_61():
    f = 1/sqrt(2*x**4 - 4*x**2 - 3)
    F = 2**(sympy.S(1)/4)*sqrt(3)*5**(sympy.S(3)/4)*sqrt((x**2*(2 + sqrt(10)) + 3)/(x**2*(2 - sqrt(10)) + 3))*sqrt(-x**2*(2 - sqrt(10)) - 3)*elliptic_f(asin(2**(sympy.S(3)/4)*5**(sympy.S(1)/4)*x/sqrt(-x**2*(2 - sqrt(10)) - 3)), sympy.S.Half - sqrt(10)/10)/(30*sqrt(2*x**4 - 4*x**2 - 3)*sqrt(1/(x**2*(2 - sqrt(10)) + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_62():
    f = 1/sqrt(2*x**4 - 5*x**2 - 3)
    F = sqrt(7)*sqrt(x**2 - 3)*sqrt(2*x**2 + 1)*elliptic_f(asin(sqrt(7)*x/sqrt(x**2 - 3)), sympy.S(1)/7)/(7*sqrt(2*x**4 - 5*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_63():
    f = 1/sqrt(3*x**4 + 5*x**2 + 2)
    F = sqrt(2)*sqrt((3*x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt(3*x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_64():
    f = 1/sqrt(3*x**4 + 4*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 4*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/6)/(12*sqrt(3*x**4 + 4*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_65():
    f = 1/sqrt(3*x**4 + 3*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 3*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/8)/(12*sqrt(3*x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_66():
    f = 1/sqrt(3*x**4 + 2*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 2*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/12)/(12*sqrt(3*x**4 + 2*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_67():
    f = 1/sqrt(3*x**4 + x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/24)/(12*sqrt(3*x**4 + x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_68():
    f = 1/sqrt(3*x**4 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half)/(12*sqrt(3*x**4 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_69():
    f = 1/sqrt(3*x**4 - x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/24 + sympy.S.Half)/(12*sqrt(3*x**4 - x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_70():
    f = 1/sqrt(3*x**4 - 2*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 2*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/12 + sympy.S.Half)/(12*sqrt(3*x**4 - 2*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_71():
    f = 1/sqrt(3*x**4 - 3*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 3*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/8 + sympy.S.Half)/(12*sqrt(3*x**4 - 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_72():
    f = 1/sqrt(3*x**4 - 4*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 4*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/6 + sympy.S.Half)/(12*sqrt(3*x**4 - 4*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_73():
    f = 1/sqrt(3*x**4 - 5*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 5*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half + 5*sqrt(6)/24)/(12*sqrt(3*x**4 - 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_74():
    f = 1/sqrt(3*x**4 - 6*x**2 + 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 6*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half + sqrt(6)/4)/(12*sqrt(3*x**4 - 6*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_75():
    f = 1/sqrt(2*x**4 + 9*x**2 + 3)
    F = sqrt((x**2*(9 - sqrt(57)) + 6)/(x**2*(sqrt(57) + 9) + 6))*(x**2*(sqrt(57) + 9) + 6)*elliptic_f(atan(x*sqrt(sqrt(57)/6 + sympy.S(3)/2)), sympy.S(-19)/4 + 3*sqrt(57)/4)/(sqrt(6*sqrt(57) + 54)*sqrt(2*x**4 + 9*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_76():
    f = 1/sqrt(2*x**4 + 8*x**2 + 3)
    F = sqrt((x**2*(4 - sqrt(10)) + 3)/(x**2*(sqrt(10) + 4) + 3))*(x**2*(sqrt(10) + 4) + 3)*elliptic_f(atan(x*sqrt(sqrt(10)/3 + sympy.S(4)/3)), sympy.S(-10)/3 + 4*sqrt(10)/3)/(sqrt(3*sqrt(10) + 12)*sqrt(2*x**4 + 8*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_77():
    f = 1/sqrt(2*x**4 + 7*x**2 + 3)
    F = sqrt(6)*sqrt((x**2 + 3)/(2*x**2 + 1))*(2*x**2 + 1)*elliptic_f(atan(sqrt(2)*x), sympy.S(5)/6)/(6*sqrt(2*x**4 + 7*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_78():
    f = 1/sqrt(2*x**4 + 6*x**2 + 3)
    F = sqrt((x**2*(3 - sqrt(3)) + 3)/(x**2*(sqrt(3) + 3) + 3))*(x**2*(sqrt(3) + 3) + 3)*elliptic_f(atan(x*sqrt(sqrt(3)/3 + 1)), -1 + sqrt(3))/(sqrt(3*sqrt(3) + 9)*sqrt(2*x**4 + 6*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_79():
    f = 1/sqrt(2*x**4 + 5*x**2 + 3)
    F = sqrt(3)*sqrt((2*x**2 + 3)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S(1)/3)/(3*sqrt(2*x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_80():
    f = 1/sqrt(2*x**4 + 4*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 4*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/6)/(12*sqrt(2*x**4 + 4*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_81():
    f = 1/sqrt(2*x**4 + 3*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 3*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/8)/(12*sqrt(2*x**4 + 3*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_82():
    f = 1/sqrt(2*x**4 + 2*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 2*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/12)/(12*sqrt(2*x**4 + 2*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_83():
    f = 1/sqrt(2*x**4 + x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/24)/(12*sqrt(2*x**4 + x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_84():
    f = 1/sqrt(2*x**4 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half)/(12*sqrt(2*x**4 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_85():
    f = 1/sqrt(2*x**4 - x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/24 + sympy.S.Half)/(12*sqrt(2*x**4 - x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_86():
    f = 1/sqrt(2*x**4 - 2*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 2*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/12 + sympy.S.Half)/(12*sqrt(2*x**4 - 2*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_87():
    f = 1/sqrt(2*x**4 - 3*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 3*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/8 + sympy.S.Half)/(12*sqrt(2*x**4 - 3*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_88():
    f = 1/sqrt(2*x**4 - 4*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 4*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/6 + sympy.S.Half)/(12*sqrt(2*x**4 - 4*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_89():
    f = 1/sqrt(2*x**4 - 5*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 5*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half + 5*sqrt(6)/24)/(12*sqrt(2*x**4 - 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_90():
    f = 1/sqrt(2*x**4 - 6*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 6*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half + sqrt(6)/4)/(12*sqrt(2*x**4 - 6*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_91():
    f = 1/sqrt(2*x**4 - 7*x**2 + 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 7*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half + 7*sqrt(6)/24)/(12*sqrt(2*x**4 - 7*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_92():
    f = 1/sqrt(-2*x**4 + 7*x**2 - 3)
    F = -sqrt(5)*elliptic_f(acos(sqrt(3)*x/3), sympy.S(6)/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_93():
    f = 1/sqrt(-2*x**4 + 6*x**2 - 3)
    F = -sqrt(2)*3**(sympy.S(3)/4)*elliptic_f(acos(x*sqrt(1 - sqrt(3)/3)), sympy.S.Half + sqrt(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_94():
    f = 1/sqrt(-2*x**4 + 5*x**2 - 3)
    F = -elliptic_f(acos(sqrt(6)*x/3), 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_95():
    f = 1/sqrt(-2*x**4 + 4*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 4*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/6 + sympy.S.Half)/(12*sqrt(-2*x**4 + 4*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_96():
    f = 1/sqrt(-2*x**4 + 3*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 3*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/8 + sympy.S.Half)/(12*sqrt(-2*x**4 + 3*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_97():
    f = 1/sqrt(-2*x**4 + 2*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - 2*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/12 + sympy.S.Half)/(12*sqrt(-2*x**4 + 2*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_98():
    f = 1/sqrt(-2*x**4 + x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 - x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sqrt(6)/24 + sympy.S.Half)/(12*sqrt(-2*x**4 + x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_99():
    f = 1/sqrt(-2*x**4 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half)/(12*sqrt(-2*x**4 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_100():
    f = 1/sqrt(-2*x**4 - x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/24)/(12*sqrt(-2*x**4 - x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_101():
    f = 1/sqrt(-2*x**4 - 2*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 2*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/12)/(12*sqrt(-2*x**4 - 2*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_102():
    f = 1/sqrt(-2*x**4 - 3*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 3*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/8)/(12*sqrt(-2*x**4 - 3*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_103():
    f = 1/sqrt(-2*x**4 - 4*x**2 - 3)
    F = 6**(sympy.S(3)/4)*sqrt((2*x**4 + 4*x**2 + 3)/(sqrt(6)*x**2 + 3)**2)*(sqrt(6)*x**2 + 3)*elliptic_f(2*atan(2**(sympy.S(1)/4)*3**(sympy.S(3)/4)*x/3), sympy.S.Half - sqrt(6)/6)/(12*sqrt(-2*x**4 - 4*x**2 - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_104():
    f = 1/sqrt(-2*x**4 - 5*x**2 - 3)
    F = sqrt(3)*sqrt(2*x**2 + 3)*elliptic_f(atan(x), sympy.S(1)/3)/(3*sqrt((2*x**2 + 3)/(x**2 + 1))*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_105():
    f = 1/sqrt(-3*x**4 + 6*x**2 - 2)
    F = -sqrt(2)*3**(sympy.S(3)/4)*elliptic_f(acos(sqrt(3)*x/sqrt(sqrt(3) + 3)), sympy.S.Half + sqrt(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_106():
    f = 1/sqrt(-3*x**4 + 5*x**2 - 2)
    F = -elliptic_f(acos(x), 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_107():
    f = 1/sqrt(-3*x**4 + 4*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 4*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/6 + sympy.S.Half)/(12*sqrt(-3*x**4 + 4*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_108():
    f = 1/sqrt(-3*x**4 + 3*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 3*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/8 + sympy.S.Half)/(12*sqrt(-3*x**4 + 3*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_109():
    f = 1/sqrt(-3*x**4 + 2*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - 2*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/12 + sympy.S.Half)/(12*sqrt(-3*x**4 + 2*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_110():
    f = 1/sqrt(-3*x**4 + x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 - x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sqrt(6)/24 + sympy.S.Half)/(12*sqrt(-3*x**4 + x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_111():
    f = 1/sqrt(-3*x**4 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half)/(12*sqrt(-3*x**4 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_112():
    f = 1/sqrt(-3*x**4 - x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/24)/(12*sqrt(-3*x**4 - x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_113():
    f = 1/sqrt(-3*x**4 - 2*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 2*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/12)/(12*sqrt(-3*x**4 - 2*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_114():
    f = 1/sqrt(-3*x**4 - 3*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 3*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/8)/(12*sqrt(-3*x**4 - 3*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_115():
    f = 1/sqrt(-3*x**4 - 4*x**2 - 2)
    F = 6**(sympy.S(3)/4)*sqrt((3*x**4 + 4*x**2 + 2)/(sqrt(6)*x**2 + 2)**2)*(sqrt(6)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*3**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(6)/6)/(12*sqrt(-3*x**4 - 4*x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_116():
    f = 1/sqrt(-3*x**4 - 5*x**2 - 2)
    F = -sqrt(2)*sqrt(-3*x**2 - 2)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_117():
    f = 1/sqrt(5*x**4 + 5*x**2 + 2)
    F = 10**(sympy.S(3)/4)*sqrt((5*x**4 + 5*x**2 + 2)/(sqrt(10)*x**2 + 2)**2)*(sqrt(10)*x**2 + 2)*elliptic_f(2*atan(2**(sympy.S(3)/4)*5**(sympy.S(1)/4)*x/2), sympy.S.Half - sqrt(10)/8)/(20*sqrt(5*x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_118():
    f = 1/sqrt(4*x**4 + 5*x**2 + 2)
    F = 2**(sympy.S(1)/4)*sqrt((4*x**4 + 5*x**2 + 2)/(sqrt(2)*x**2 + 1)**2)*(sqrt(2)*x**2 + 1)*elliptic_f(2*atan(2**(sympy.S(1)/4)*x), sympy.S.Half - 5*sqrt(2)/16)/(4*sqrt(4*x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_119():
    f = 1/sqrt(3*x**4 + 5*x**2 + 2)
    F = sqrt(2)*sqrt((3*x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt(3*x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_120():
    f = 1/sqrt(2*x**4 + 5*x**2 + 2)
    F = sqrt((x**2 + 2)/(2*x**2 + 1))*(2*x**2 + 1)*elliptic_f(atan(sqrt(2)*x), sympy.S(3)/4)/(2*sqrt(2*x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_121():
    f = 1/sqrt(x**4 + 5*x**2 + 2)
    F = sqrt((x**2*(5 - sqrt(17)) + 4)/(x**2*(sqrt(17) + 5) + 4))*(x**2*(sqrt(17) + 5) + 4)*elliptic_f(atan(x*sqrt(sqrt(17) + 5)/2), sympy.S(-17)/4 + 5*sqrt(17)/4)/(2*sqrt(sqrt(17) + 5)*sqrt(x**4 + 5*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_122():
    f = 1/sqrt(-x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(sqrt(2)*x/sqrt(5 + sqrt(33))), sympy.S(-29)/4 - 5*sqrt(33)/4)/sqrt(-5 + sqrt(33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_123():
    f = 1/sqrt(-2*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(2*x/sqrt(5 + sqrt(41))), sympy.S(-33)/8 - 5*sqrt(41)/8)/sqrt(-5 + sqrt(41))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_124():
    f = 1/sqrt(-3*x**4 + 5*x**2 + 2)
    F = elliptic_f(asin(sqrt(2)*x/2), -6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_125():
    f = 1/sqrt(-4*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(2*sqrt(2)*x/sqrt(5 + sqrt(57))), sympy.S(-41)/16 - 5*sqrt(57)/16)/sqrt(-5 + sqrt(57))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_126():
    f = 1/sqrt(-5*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(sqrt(10)*x/sqrt(5 + sqrt(65))), sympy.S(-9)/4 - sqrt(65)/4)/sqrt(-5 + sqrt(65))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_127():
    f = 1/sqrt(-6*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(2*sqrt(3)*x/sqrt(5 + sqrt(73))), sympy.S(-49)/24 - 5*sqrt(73)/24)/sqrt(-5 + sqrt(73))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_128():
    f = 1/sqrt(-7*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(-7)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_129():
    f = 1/sqrt(-8*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(4*x/sqrt(5 + sqrt(89))), sympy.S(-57)/32 - 5*sqrt(89)/32)/sqrt(-5 + sqrt(89))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_130():
    f = 1/sqrt(-9*x**4 + 5*x**2 + 2)
    F = sqrt(2)*elliptic_f(asin(3*sqrt(2)*x/sqrt(5 + sqrt(97))), sympy.S(-61)/36 - 5*sqrt(97)/36)/sqrt(-5 + sqrt(97))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_131():
    f = x**2*(b*x**2 + c*x**4)
    F = b*x**5/5 + c*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_132():
    f = x*(b*x**2 + c*x**4)
    F = b*x**4/4 + c*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_133():
    f = b*x**2 + c*x**4
    F = b*x**3/3 + c*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_134():
    f = (b*x**2 + c*x**4)/x
    F = b*x**2/2 + c*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_135():
    f = (b*x**2 + c*x**4)/x**2
    F = b*x + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_136():
    f = (b*x**2 + c*x**4)/x**3
    F = b*log(x) + c*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_137():
    f = (b*x**2 + c*x**4)/x**4
    F = -b/x + c*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_138():
    f = (b*x**2 + c*x**4)/x**5
    F = -b/(2*x**2) + c*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_139():
    f = (b*x**2 + c*x**4)/x**6
    F = -b/(3*x**3) - c/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_140():
    f = (b*x**2 + c*x**4)/x**7
    F = -b/(4*x**4) - c/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_141():
    f = (b*x**2 + c*x**4)/x**8
    F = -b/(5*x**5) - c/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_142():
    f = (b*x**2 + c*x**4)**2
    F = b**2*x**5/5 + 2*b*c*x**7/7 + c**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_143():
    f = (b*x**2 + c*x**4)**2/x
    F = b**2*x**4/4 + b*c*x**6/3 + c**2*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_144():
    f = (b*x**2 + c*x**4)**2/x**2
    F = b**2*x**3/3 + 2*b*c*x**5/5 + c**2*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_145():
    f = (b*x**2 + c*x**4)**2/x**3
    F = (b + c*x**2)**3/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_146():
    f = (b*x**2 + c*x**4)**2/x**4
    F = b**2*x + 2*b*c*x**3/3 + c**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_147():
    f = (b*x**2 + c*x**4)**2/x**5
    F = b**2*log(x) + b*c*x**2 + c**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_148():
    f = (b*x**2 + c*x**4)**2/x**6
    F = -b**2/x + 2*b*c*x + c**2*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_149():
    f = (b*x**2 + c*x**4)**2/x**7
    F = -b**2/(2*x**2) + 2*b*c*log(x) + c**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_150():
    f = (b*x**2 + c*x**4)**2/x**8
    F = -b**2/(3*x**3) - 2*b*c/x + c**2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_151():
    f = (b*x**2 + c*x**4)**2/x**9
    F = -b**2/(4*x**4) - b*c/x**2 + c**2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_152():
    f = (b*x**2 + c*x**4)**2/x**10
    F = -b**2/(5*x**5) - 2*b*c/(3*x**3) - c**2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_153():
    f = (b*x**2 + c*x**4)**2/x**11
    F = -(b + c*x**2)**3/(6*b*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_154():
    f = (b*x**2 + c*x**4)**2/x**12
    F = -b**2/(7*x**7) - 2*b*c/(5*x**5) - c**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_155():
    f = (b*x**2 + c*x**4)**3/x**2
    F = b**3*x**5/5 + 3*b**2*c*x**7/7 + b*c**2*x**9/3 + c**3*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_156():
    f = (b*x**2 + c*x**4)**3/x**3
    F = -b*(b + c*x**2)**4/(8*c**2) + (b + c*x**2)**5/(10*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_157():
    f = (b*x**2 + c*x**4)**3/x**4
    F = b**3*x**3/3 + 3*b**2*c*x**5/5 + 3*b*c**2*x**7/7 + c**3*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_158():
    f = (b*x**2 + c*x**4)**3/x**5
    F = (b + c*x**2)**4/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_159():
    f = (b*x**2 + c*x**4)**3/x**6
    F = b**3*x + b**2*c*x**3 + 3*b*c**2*x**5/5 + c**3*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_160():
    f = (b*x**2 + c*x**4)**3/x**7
    F = b**3*log(x) + 3*b**2*c*x**2/2 + 3*b*c**2*x**4/4 + c**3*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_161():
    f = (b*x**2 + c*x**4)**3/x**8
    F = -b**3/x + 3*b**2*c*x + b*c**2*x**3 + c**3*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_162():
    f = (b*x**2 + c*x**4)**3/x**9
    F = -b**3/(2*x**2) + 3*b**2*c*log(x) + 3*b*c**2*x**2/2 + c**3*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_163():
    f = (b*x**2 + c*x**4)**3/x**10
    F = -b**3/(3*x**3) - 3*b**2*c/x + 3*b*c**2*x + c**3*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_164():
    f = (b*x**2 + c*x**4)**3/x**11
    F = -b**3/(4*x**4) - 3*b**2*c/(2*x**2) + 3*b*c**2*log(x) + c**3*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_165():
    f = (b*x**2 + c*x**4)**3/x**12
    F = -b**3/(5*x**5) - b**2*c/x**3 - 3*b*c**2/x + c**3*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_166():
    f = (b*x**2 + c*x**4)**3/x**13
    F = -b**3/(6*x**6) - 3*b**2*c/(4*x**4) - 3*b*c**2/(2*x**2) + c**3*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_167():
    f = (b*x**2 + c*x**4)**3/x**14
    F = -b**3/(7*x**7) - 3*b**2*c/(5*x**5) - b*c**2/x**3 - c**3/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_168():
    f = (b*x**2 + c*x**4)**3/x**15
    F = -(b + c*x**2)**4/(8*b*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_169():
    f = (b*x**2 + c*x**4)**3/x**16
    F = -b**3/(9*x**9) - 3*b**2*c/(7*x**7) - 3*b*c**2/(5*x**5) - c**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_170():
    f = (b*x**2 + c*x**4)**3/x**17
    F = -(b + c*x**2)**4/(10*b*x**10) + c*(b + c*x**2)**4/(40*b**2*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_171():
    f = x**10/(b*x**2 + c*x**4)
    F = b**(sympy.S(7)/2)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(9)/2) - b**3*x/c**4 + b**2*x**3/(3*c**3) - b*x**5/(5*c**2) + x**7/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_172():
    f = x**9/(b*x**2 + c*x**4)
    F = -b**3*log(b + c*x**2)/(2*c**4) + b**2*x**2/(2*c**3) - b*x**4/(4*c**2) + x**6/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_173():
    f = x**8/(b*x**2 + c*x**4)
    F = -b**(sympy.S(5)/2)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(7)/2) + b**2*x/c**3 - b*x**3/(3*c**2) + x**5/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_174():
    f = x**7/(b*x**2 + c*x**4)
    F = b**2*log(b + c*x**2)/(2*c**3) - b*x**2/(2*c**2) + x**4/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_175():
    f = x**6/(b*x**2 + c*x**4)
    F = b**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(5)/2) - b*x/c**2 + x**3/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_176():
    f = x**5/(b*x**2 + c*x**4)
    F = -b*log(b + c*x**2)/(2*c**2) + x**2/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_177():
    f = x**4/(b*x**2 + c*x**4)
    F = -sqrt(b)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(3)/2) + x/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_178():
    f = x**3/(b*x**2 + c*x**4)
    F = log(b + c*x**2)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_179():
    f = x**2/(b*x**2 + c*x**4)
    F = atan(sqrt(c)*x/sqrt(b))/(sqrt(b)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_180():
    f = x/(b*x**2 + c*x**4)
    F = log(x)/b - log(b + c*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_181():
    f = 1/(b*x**2 + c*x**4)
    F = -1/(b*x) - sqrt(c)*atan(sqrt(c)*x/sqrt(b))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_182():
    f = 1/(x*(b*x**2 + c*x**4))
    F = -1/(2*b*x**2) - c*log(x)/b**2 + c*log(b + c*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_183():
    f = 1/(x**2*(b*x**2 + c*x**4))
    F = -1/(3*b*x**3) + c/(b**2*x) + c**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_184():
    f = 1/(x**3*(b*x**2 + c*x**4))
    F = -1/(4*b*x**4) + c/(2*b**2*x**2) + c**2*log(x)/b**3 - c**2*log(b + c*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_185():
    f = 1/(x**4*(b*x**2 + c*x**4))
    F = -1/(5*b*x**5) + c/(3*b**2*x**3) - c**2/(b**3*x) - c**(sympy.S(5)/2)*atan(sqrt(c)*x/sqrt(b))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_186():
    f = 1/(x**5*(b*x**2 + c*x**4))
    F = -1/(6*b*x**6) + c/(4*b**2*x**4) - c**2/(2*b**3*x**2) - c**3*log(x)/b**4 + c**3*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_187():
    f = x**12/(b*x**2 + c*x**4)**2
    F = -7*b**(sympy.S(5)/2)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(9)/2)) + 7*b**2*x/(2*c**4) - 7*b*x**3/(6*c**3) - x**7/(2*c*(b + c*x**2)) + 7*x**5/(10*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_188():
    f = x**11/(b*x**2 + c*x**4)**2
    F = b**3/(2*c**4*(b + c*x**2)) + 3*b**2*log(b + c*x**2)/(2*c**4) - b*x**2/c**3 + x**4/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_189():
    f = x**10/(b*x**2 + c*x**4)**2
    F = 5*b**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(7)/2)) - 5*b*x/(2*c**3) - x**5/(2*c*(b + c*x**2)) + 5*x**3/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_190():
    f = x**9/(b*x**2 + c*x**4)**2
    F = -b**2/(2*c**3*(b + c*x**2)) - b*log(b + c*x**2)/c**3 + x**2/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_191():
    f = x**8/(b*x**2 + c*x**4)**2
    F = -3*sqrt(b)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(5)/2)) - x**3/(2*c*(b + c*x**2)) + 3*x/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_192():
    f = x**7/(b*x**2 + c*x**4)**2
    F = b/(2*c**2*(b + c*x**2)) + log(b + c*x**2)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_193():
    f = x**6/(b*x**2 + c*x**4)**2
    F = -x/(2*c*(b + c*x**2)) + atan(sqrt(c)*x/sqrt(b))/(2*sqrt(b)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_194():
    f = x**5/(b*x**2 + c*x**4)**2
    F = -1/(2*c*(b + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_195():
    f = x**4/(b*x**2 + c*x**4)**2
    F = x/(2*b*(b + c*x**2)) + atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(3)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_196():
    f = x**3/(b*x**2 + c*x**4)**2
    F = 1/(2*b*(b + c*x**2)) + log(x)/b**2 - log(b + c*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_197():
    f = x**2/(b*x**2 + c*x**4)**2
    F = 1/(2*b*x*(b + c*x**2)) - 3/(2*b**2*x) - 3*sqrt(c)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_198():
    f = x/(b*x**2 + c*x**4)**2
    F = -c/(2*b**2*(b + c*x**2)) - 1/(2*b**2*x**2) - 2*c*log(x)/b**3 + c*log(b + c*x**2)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_199():
    f = (b*x**2 + c*x**4)**(-2)
    F = 1/(2*b*x**3*(b + c*x**2)) - 5/(6*b**2*x**3) + 5*c/(2*b**3*x) + 5*c**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_200():
    f = 1/(x*(b*x**2 + c*x**4)**2)
    F = -1/(4*b**2*x**4) + c**2/(2*b**3*(b + c*x**2)) + c/(b**3*x**2) + 3*c**2*log(x)/b**4 - 3*c**2*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_201():
    f = 1/(x**2*(b*x**2 + c*x**4)**2)
    F = 1/(2*b*x**5*(b + c*x**2)) - 7/(10*b**2*x**5) + 7*c/(6*b**3*x**3) - 7*c**2/(2*b**4*x) - 7*c**(sympy.S(5)/2)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_202():
    f = x**14/(b*x**2 + c*x**4)**3
    F = 35*b**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/(8*c**(sympy.S(9)/2)) - 35*b*x/(8*c**4) - x**7/(4*c*(b + c*x**2)**2) - 7*x**5/(8*c**2*(b + c*x**2)) + 35*x**3/(24*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_203():
    f = x**13/(b*x**2 + c*x**4)**3
    F = b**3/(4*c**4*(b + c*x**2)**2) - 3*b**2/(2*c**4*(b + c*x**2)) - 3*b*log(b + c*x**2)/(2*c**4) + x**2/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_204():
    f = x**12/(b*x**2 + c*x**4)**3
    F = -15*sqrt(b)*atan(sqrt(c)*x/sqrt(b))/(8*c**(sympy.S(7)/2)) - x**5/(4*c*(b + c*x**2)**2) - 5*x**3/(8*c**2*(b + c*x**2)) + 15*x/(8*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_205():
    f = x**11/(b*x**2 + c*x**4)**3
    F = -b**2/(4*c**3*(b + c*x**2)**2) + b/(c**3*(b + c*x**2)) + log(b + c*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_206():
    f = x**10/(b*x**2 + c*x**4)**3
    F = -x**3/(4*c*(b + c*x**2)**2) - 3*x/(8*c**2*(b + c*x**2)) + 3*atan(sqrt(c)*x/sqrt(b))/(8*sqrt(b)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_207():
    f = x**9/(b*x**2 + c*x**4)**3
    F = x**4/(4*b*(b + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_208():
    f = x**8/(b*x**2 + c*x**4)**3
    F = -x/(4*c*(b + c*x**2)**2) + x/(8*b*c*(b + c*x**2)) + atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_209():
    f = x**7/(b*x**2 + c*x**4)**3
    F = -1/(4*c*(b + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_210():
    f = x**6/(b*x**2 + c*x**4)**3
    F = x/(4*b*(b + c*x**2)**2) + 3*x/(8*b**2*(b + c*x**2)) + 3*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(5)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_211():
    f = x**5/(b*x**2 + c*x**4)**3
    F = 1/(4*b*(b + c*x**2)**2) + 1/(2*b**2*(b + c*x**2)) + log(x)/b**3 - log(b + c*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_212():
    f = x**4/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x*(b + c*x**2)**2) + 5/(8*b**2*x*(b + c*x**2)) - 15/(8*b**3*x) - 15*sqrt(c)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_213():
    f = x**3/(b*x**2 + c*x**4)**3
    F = -c/(4*b**2*(b + c*x**2)**2) - c/(b**3*(b + c*x**2)) - 1/(2*b**3*x**2) - 3*c*log(x)/b**4 + 3*c*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_214():
    f = x**2/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x**3*(b + c*x**2)**2) + 7/(8*b**2*x**3*(b + c*x**2)) - 35/(24*b**3*x**3) + 35*c/(8*b**4*x) + 35*c**(sympy.S(3)/2)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_215():
    f = x/(b*x**2 + c*x**4)**3
    F = c**2/(4*b**3*(b + c*x**2)**2) - 1/(4*b**3*x**4) + 3*c**2/(2*b**4*(b + c*x**2)) + 3*c/(2*b**4*x**2) + 6*c**2*log(x)/b**5 - 3*c**2*log(b + c*x**2)/b**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_216():
    f = (b*x**2 + c*x**4)**(-3)
    F = 1/(4*b*x**5*(b + c*x**2)**2) + 9/(8*b**2*x**5*(b + c*x**2)) - 63/(40*b**3*x**5) + 21*c/(8*b**4*x**3) - 63*c**2/(8*b**5*x) - 63*c**(sympy.S(5)/2)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_217():
    f = 1/(x*(b*x**2 + c*x**4)**3)
    F = -1/(6*b**3*x**6) - c**3/(4*b**4*(b + c*x**2)**2) + 3*c/(4*b**4*x**4) - 2*c**3/(b**5*(b + c*x**2)) - 3*c**2/(b**5*x**2) - 10*c**3*log(x)/b**6 + 5*c**3*log(b + c*x**2)/b**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_218():
    f = x**5*sqrt(b*x**2 + c*x**4)
    F = -5*b**4*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(128*c**(sympy.S(7)/2)) + 5*b**2*(b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(128*c**3) - 5*b*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(48*c**2) + x**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_219():
    f = x**3*sqrt(b*x**2 + c*x**4)
    F = b**3*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(5)/2)) - b*(b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(16*c**2) + (b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_220():
    f = x*sqrt(b*x**2 + c*x**4)
    F = -b**2*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(3)/2)) + (b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_221():
    f = sqrt(b*x**2 + c*x**4)/x
    F = b*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*sqrt(c)) + sqrt(b*x**2 + c*x**4)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_222():
    f = sqrt(b*x**2 + c*x**4)/x**3
    F = sqrt(c)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4)) - sqrt(b*x**2 + c*x**4)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_223():
    f = sqrt(b*x**2 + c*x**4)/x**5
    F = -(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*b*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_224():
    f = sqrt(b*x**2 + c*x**4)/x**7
    F = -(b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*b*x**8) + 2*c*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*b**2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_225():
    f = sqrt(b*x**2 + c*x**4)/x**9
    F = -(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*b*x**10) + 4*c*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(35*b**2*x**8) - 8*c**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*b**3*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_226():
    f = sqrt(b*x**2 + c*x**4)/x**11
    F = -(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*b*x**12) + 2*c*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(21*b**2*x**10) - 8*c**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*b**3*x**8) + 16*c**3*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(315*b**4*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_227():
    f = sqrt(b*x**2 + c*x**4)/x**13
    F = -(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*b*x**14) + 8*c*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(99*b**2*x**12) - 16*c**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(231*b**3*x**10) + 64*c**3*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(1155*b**4*x**8) - 128*c**4*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3465*b**5*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_228():
    f = x**4*sqrt(b*x**2 + c*x**4)
    F = 8*b**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*c**3*x**3) - 4*b*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(35*c**2*x) + x*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_229():
    f = x**2*sqrt(b*x**2 + c*x**4)
    F = -2*b*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*c**2*x**3) + (b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_230():
    f = sqrt(b*x**2 + c*x**4)
    F = (b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_231():
    f = sqrt(b*x**2 + c*x**4)/x**2
    F = -sqrt(b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4)) + sqrt(b*x**2 + c*x**4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_232():
    f = sqrt(b*x**2 + c*x**4)/x**4
    F = -sqrt(b*x**2 + c*x**4)/(2*x**3) - c*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_233():
    f = sqrt(b*x**2 + c*x**4)/x**6
    F = -sqrt(b*x**2 + c*x**4)/(4*x**5) - c*sqrt(b*x**2 + c*x**4)/(8*b*x**3) + c**2*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_234():
    f = sqrt(b*x**2 + c*x**4)/x**8
    F = -sqrt(b*x**2 + c*x**4)/(6*x**7) - c*sqrt(b*x**2 + c*x**4)/(24*b*x**5) + c**2*sqrt(b*x**2 + c*x**4)/(16*b**2*x**3) - c**3*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(16*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_235():
    f = x**3*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -3*b**5*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(256*c**(sympy.S(7)/2)) + 3*b**3*(b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(256*c**3) - b*(b + 2*c*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(32*c**2) + (b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_236():
    f = x*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 3*b**4*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(128*c**(sympy.S(5)/2)) - 3*b**2*(b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(128*c**2) + (b + 2*c*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(16*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_237():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x
    F = -b**3*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(3)/2)) + b*(b + 2*c*x**2)*sqrt(b*x**2 + c*x**4)/(16*c) + (b*x**2 + c*x**4)**(sympy.S(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_238():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**3
    F = 3*b**2*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*sqrt(c)) + 3*b*sqrt(b*x**2 + c*x**4)/8 + (b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_239():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**5
    F = 3*b*sqrt(c)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/2 + 3*c*sqrt(b*x**2 + c*x**4)/2 - (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_240():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**7
    F = c**(sympy.S(3)/2)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4)) - c*sqrt(b*x**2 + c*x**4)/x**2 - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_241():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**9
    F = -(b*x**2 + c*x**4)**(sympy.S(5)/2)/(5*b*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_242():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**11
    F = -(b*x**2 + c*x**4)**(sympy.S(5)/2)/(7*b*x**12) + 2*c*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(35*b**2*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_243():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**13
    F = -(b*x**2 + c*x**4)**(sympy.S(5)/2)/(9*b*x**14) + 4*c*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(63*b**2*x**12) - 8*c**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(315*b**3*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_244():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**15
    F = -(b*x**2 + c*x**4)**(sympy.S(5)/2)/(11*b*x**16) + 2*c*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(33*b**2*x**14) - 8*c**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(231*b**3*x**12) + 16*c**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(1155*b**4*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_245():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**17
    F = -(b*x**2 + c*x**4)**(sympy.S(5)/2)/(13*b*x**18) + 8*c*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(143*b**2*x**16) - 16*c**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(429*b**3*x**14) + 64*c**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3003*b**4*x**12) - 128*c**4*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15015*b**5*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_246():
    f = x**6*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 128*b**4*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15015*c**5*x**5) - 64*b**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3003*c**4*x**3) + 16*b**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(429*c**3*x) - 8*b*x*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(143*c**2) + x**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(13*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_247():
    f = x**4*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -16*b**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(1155*c**4*x**5) + 8*b**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(231*c**3*x**3) - 2*b*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(33*c**2*x) + x*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(11*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_248():
    f = x**2*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 8*b**2*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(315*c**3*x**5) - 4*b*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(63*c**2*x**3) + (b*x**2 + c*x**4)**(sympy.S(5)/2)/(9*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_249():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -2*b*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(35*c**2*x**5) + (b*x**2 + c*x**4)**(sympy.S(5)/2)/(7*c*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_250():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**2
    F = (b*x**2 + c*x**4)**(sympy.S(5)/2)/(5*c*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_251():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**4
    F = -b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4)) + b*sqrt(b*x**2 + c*x**4)/x + (b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_252():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**6
    F = -3*sqrt(b)*c*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/2 + 3*c*sqrt(b*x**2 + c*x**4)/(2*x) - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_253():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**8
    F = -3*c*sqrt(b*x**2 + c*x**4)/(8*x**3) - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*x**7) - 3*c**2*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_254():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**10
    F = -c*sqrt(b*x**2 + c*x**4)/(8*x**5) - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*x**9) - c**2*sqrt(b*x**2 + c*x**4)/(16*b*x**3) + c**3*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(16*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_255():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**12
    F = -c*sqrt(b*x**2 + c*x**4)/(16*x**7) - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(8*x**11) - c**2*sqrt(b*x**2 + c*x**4)/(64*b*x**5) + 3*c**3*sqrt(b*x**2 + c*x**4)/(128*b**2*x**3) - 3*c**4*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(128*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_256():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**14
    F = -3*c*sqrt(b*x**2 + c*x**4)/(80*x**9) - (b*x**2 + c*x**4)**(sympy.S(3)/2)/(10*x**13) - c**2*sqrt(b*x**2 + c*x**4)/(160*b*x**7) + c**3*sqrt(b*x**2 + c*x**4)/(128*b**2*x**5) - 3*c**4*sqrt(b*x**2 + c*x**4)/(256*b**3*x**3) + 3*c**5*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(256*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_257():
    f = x**7/sqrt(b*x**2 + c*x**4)
    F = -5*b**3*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(7)/2)) + 5*b**2*sqrt(b*x**2 + c*x**4)/(16*c**3) - 5*b*x**2*sqrt(b*x**2 + c*x**4)/(24*c**2) + x**4*sqrt(b*x**2 + c*x**4)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_258():
    f = x**5/sqrt(b*x**2 + c*x**4)
    F = 3*b**2*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(5)/2)) - 3*b*sqrt(b*x**2 + c*x**4)/(8*c**2) + x**2*sqrt(b*x**2 + c*x**4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_259():
    f = x**3/sqrt(b*x**2 + c*x**4)
    F = -b*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*c**(sympy.S(3)/2)) + sqrt(b*x**2 + c*x**4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_260():
    f = x/sqrt(b*x**2 + c*x**4)
    F = atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_261():
    f = 1/(x*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_262():
    f = 1/(x**3*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(3*b*x**4) + 2*c*sqrt(b*x**2 + c*x**4)/(3*b**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_263():
    f = 1/(x**5*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(5*b*x**6) + 4*c*sqrt(b*x**2 + c*x**4)/(15*b**2*x**4) - 8*c**2*sqrt(b*x**2 + c*x**4)/(15*b**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_264():
    f = 1/(x**7*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(7*b*x**8) + 6*c*sqrt(b*x**2 + c*x**4)/(35*b**2*x**6) - 8*c**2*sqrt(b*x**2 + c*x**4)/(35*b**3*x**4) + 16*c**3*sqrt(b*x**2 + c*x**4)/(35*b**4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_265():
    f = x**4/sqrt(b*x**2 + c*x**4)
    F = -2*b*sqrt(b*x**2 + c*x**4)/(3*c**2*x) + x*sqrt(b*x**2 + c*x**4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_266():
    f = x**2/sqrt(b*x**2 + c*x**4)
    F = sqrt(b*x**2 + c*x**4)/(c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_267():
    f = 1/sqrt(b*x**2 + c*x**4)
    F = -atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_268():
    f = 1/(x**2*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(2*b*x**3) + c*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_269():
    f = 1/(x**4*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(4*b*x**5) + 3*c*sqrt(b*x**2 + c*x**4)/(8*b**2*x**3) - 3*c**2*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_270():
    f = x**9/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 15*b**2*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(7)/2)) - 15*b*sqrt(b*x**2 + c*x**4)/(8*c**3) - x**6/(c*sqrt(b*x**2 + c*x**4)) + 5*x**2*sqrt(b*x**2 + c*x**4)/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_271():
    f = x**7/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -3*b*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*c**(sympy.S(5)/2)) - x**4/(c*sqrt(b*x**2 + c*x**4)) + 3*sqrt(b*x**2 + c*x**4)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_272():
    f = x**5/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**2/(c*sqrt(b*x**2 + c*x**4)) + atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/c**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_273():
    f = x**3/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = x**2/(b*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_274():
    f = x/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(b + 2*c*x**2)/(b**2*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_275():
    f = 1/(x*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**2*sqrt(b*x**2 + c*x**4)) - 4*sqrt(b*x**2 + c*x**4)/(3*b**2*x**4) + 8*c*sqrt(b*x**2 + c*x**4)/(3*b**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_276():
    f = 1/(x**3*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**4*sqrt(b*x**2 + c*x**4)) - 6*sqrt(b*x**2 + c*x**4)/(5*b**2*x**6) + 8*c*sqrt(b*x**2 + c*x**4)/(5*b**3*x**4) - 16*c**2*sqrt(b*x**2 + c*x**4)/(5*b**4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_277():
    f = 1/(x**5*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**6*sqrt(b*x**2 + c*x**4)) - 8*sqrt(b*x**2 + c*x**4)/(7*b**2*x**8) + 48*c*sqrt(b*x**2 + c*x**4)/(35*b**3*x**6) - 64*c**2*sqrt(b*x**2 + c*x**4)/(35*b**4*x**4) + 128*c**3*sqrt(b*x**2 + c*x**4)/(35*b**5*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_278():
    f = x**6/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**3/(c*sqrt(b*x**2 + c*x**4)) + 2*sqrt(b*x**2 + c*x**4)/(c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_279():
    f = x**4/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x/(c*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_280():
    f = x**2/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = x/(b*sqrt(b*x**2 + c*x**4)) - atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_281():
    f = (b*x**2 + c*x**4)**(sympy.S(-3)/2)
    F = 1/(b*x*sqrt(b*x**2 + c*x**4)) - 3*sqrt(b*x**2 + c*x**4)/(2*b**2*x**3) + 3*c*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_282():
    f = 1/(x**2*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**3*sqrt(b*x**2 + c*x**4)) - 5*sqrt(b*x**2 + c*x**4)/(4*b**2*x**5) + 15*c*sqrt(b*x**2 + c*x**4)/(8*b**3*x**3) - 15*c**2*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_283():
    f = x**3/sqrt(-4*x**4 + 3*x**2)
    F = -sqrt(-4*x**4 + 3*x**2)/8 + 3*asin(8*x**2/3 - 1)/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_284():
    f = x**3/sqrt(-4*x**4 - 3*x**2)
    F = -sqrt(-4*x**4 - 3*x**2)/8 - 3*asin(8*x**2/3 + 1)/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_285():
    f = x**3/sqrt(4*x**4 + 3*x**2)
    F = sqrt(4*x**4 + 3*x**2)/8 - 3*atanh(2*x**2/sqrt(4*x**4 + 3*x**2))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_286():
    f = x**3/sqrt(4*x**4 - 3*x**2)
    F = sqrt(4*x**4 - 3*x**2)/8 + 3*atanh(2*x**2/sqrt(4*x**4 - 3*x**2))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_287():
    f = x**3/sqrt(a*x**2 + b*x**4)
    F = -a*atanh(sqrt(b)*x**2/sqrt(a*x**2 + b*x**4))/(2*b**(sympy.S(3)/2)) + sqrt(a*x**2 + b*x**4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_288():
    f = x**3/sqrt(a*x**2 - b*x**4)
    F = a*atan(sqrt(b)*x**2/sqrt(a*x**2 - b*x**4))/(2*b**(sympy.S(3)/2)) - sqrt(a*x**2 - b*x**4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_289():
    f = x**(sympy.S(7)/2)*(b*x**2 + c*x**4)
    F = 2*b*x**(sympy.S(13)/2)/13 + 2*c*x**(sympy.S(17)/2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_290():
    f = x**(sympy.S(5)/2)*(b*x**2 + c*x**4)
    F = 2*b*x**(sympy.S(11)/2)/11 + 2*c*x**(sympy.S(15)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_291():
    f = x**(sympy.S(3)/2)*(b*x**2 + c*x**4)
    F = 2*b*x**(sympy.S(9)/2)/9 + 2*c*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_292():
    f = sqrt(x)*(b*x**2 + c*x**4)
    F = 2*b*x**(sympy.S(7)/2)/7 + 2*c*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_293():
    f = (b*x**2 + c*x**4)/sqrt(x)
    F = 2*b*x**(sympy.S(5)/2)/5 + 2*c*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_294():
    f = (b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    F = 2*b*x**(sympy.S(3)/2)/3 + 2*c*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_295():
    f = (b*x**2 + c*x**4)/x**(sympy.S(5)/2)
    F = 2*b*sqrt(x) + 2*c*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_296():
    f = (b*x**2 + c*x**4)/x**(sympy.S(7)/2)
    F = -2*b/sqrt(x) + 2*c*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_297():
    f = x**(sympy.S(7)/2)*(b*x**2 + c*x**4)**2
    F = 2*b**2*x**(sympy.S(17)/2)/17 + 4*b*c*x**(sympy.S(21)/2)/21 + 2*c**2*x**(sympy.S(25)/2)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_298():
    f = x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**2
    F = 2*b**2*x**(sympy.S(15)/2)/15 + 4*b*c*x**(sympy.S(19)/2)/19 + 2*c**2*x**(sympy.S(23)/2)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_299():
    f = x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**2
    F = 2*b**2*x**(sympy.S(13)/2)/13 + 4*b*c*x**(sympy.S(17)/2)/17 + 2*c**2*x**(sympy.S(21)/2)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_300():
    f = sqrt(x)*(b*x**2 + c*x**4)**2
    F = 2*b**2*x**(sympy.S(11)/2)/11 + 4*b*c*x**(sympy.S(15)/2)/15 + 2*c**2*x**(sympy.S(19)/2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_301():
    f = (b*x**2 + c*x**4)**2/sqrt(x)
    F = 2*b**2*x**(sympy.S(9)/2)/9 + 4*b*c*x**(sympy.S(13)/2)/13 + 2*c**2*x**(sympy.S(17)/2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_302():
    f = (b*x**2 + c*x**4)**2/x**(sympy.S(3)/2)
    F = 2*b**2*x**(sympy.S(7)/2)/7 + 4*b*c*x**(sympy.S(11)/2)/11 + 2*c**2*x**(sympy.S(15)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_303():
    f = (b*x**2 + c*x**4)**2/x**(sympy.S(5)/2)
    F = 2*b**2*x**(sympy.S(5)/2)/5 + 4*b*c*x**(sympy.S(9)/2)/9 + 2*c**2*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_304():
    f = (b*x**2 + c*x**4)**2/x**(sympy.S(7)/2)
    F = 2*b**2*x**(sympy.S(3)/2)/3 + 4*b*c*x**(sympy.S(7)/2)/7 + 2*c**2*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_305():
    f = x**(sympy.S(7)/2)*(b*x**2 + c*x**4)**3
    F = 2*b**3*x**(sympy.S(21)/2)/21 + 6*b**2*c*x**(sympy.S(25)/2)/25 + 6*b*c**2*x**(sympy.S(29)/2)/29 + 2*c**3*x**(sympy.S(33)/2)/33
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_306():
    f = x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**3
    F = 2*b**3*x**(sympy.S(19)/2)/19 + 6*b**2*c*x**(sympy.S(23)/2)/23 + 2*b*c**2*x**(sympy.S(27)/2)/9 + 2*c**3*x**(sympy.S(31)/2)/31
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_307():
    f = x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**3
    F = 2*b**3*x**(sympy.S(17)/2)/17 + 2*b**2*c*x**(sympy.S(21)/2)/7 + 6*b*c**2*x**(sympy.S(25)/2)/25 + 2*c**3*x**(sympy.S(29)/2)/29
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_308():
    f = sqrt(x)*(b*x**2 + c*x**4)**3
    F = 2*b**3*x**(sympy.S(15)/2)/15 + 6*b**2*c*x**(sympy.S(19)/2)/19 + 6*b*c**2*x**(sympy.S(23)/2)/23 + 2*c**3*x**(sympy.S(27)/2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_309():
    f = (b*x**2 + c*x**4)**3/sqrt(x)
    F = 2*b**3*x**(sympy.S(13)/2)/13 + 6*b**2*c*x**(sympy.S(17)/2)/17 + 2*b*c**2*x**(sympy.S(21)/2)/7 + 2*c**3*x**(sympy.S(25)/2)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_310():
    f = (b*x**2 + c*x**4)**3/x**(sympy.S(3)/2)
    F = 2*b**3*x**(sympy.S(11)/2)/11 + 2*b**2*c*x**(sympy.S(15)/2)/5 + 6*b*c**2*x**(sympy.S(19)/2)/19 + 2*c**3*x**(sympy.S(23)/2)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_311():
    f = (b*x**2 + c*x**4)**3/x**(sympy.S(5)/2)
    F = 2*b**3*x**(sympy.S(9)/2)/9 + 6*b**2*c*x**(sympy.S(13)/2)/13 + 6*b*c**2*x**(sympy.S(17)/2)/17 + 2*c**3*x**(sympy.S(21)/2)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_312():
    f = (b*x**2 + c*x**4)**3/x**(sympy.S(7)/2)
    F = 2*b**3*x**(sympy.S(7)/2)/7 + 6*b**2*c*x**(sympy.S(11)/2)/11 + 2*b*c**2*x**(sympy.S(15)/2)/5 + 2*c**3*x**(sympy.S(19)/2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_313():
    f = x**(sympy.S(13)/2)/(b*x**2 + c*x**4)
    F = sqrt(2)*b**(sympy.S(7)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(11)/4)) - sqrt(2)*b**(sympy.S(7)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(11)/4)) - sqrt(2)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)) + sqrt(2)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)) - 2*b*x**(sympy.S(3)/2)/(3*c**2) + 2*x**(sympy.S(7)/2)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_314():
    f = x**(sympy.S(11)/2)/(b*x**2 + c*x**4)
    F = -sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)) - 2*b*sqrt(x)/c**2 + 2*x**(sympy.S(5)/2)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_315():
    f = x**(sympy.S(9)/2)/(b*x**2 + c*x**4)
    F = -sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(7)/4)) + sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(7)/4)) + sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)) - sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)) + 2*x**(sympy.S(3)/2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_316():
    f = x**(sympy.S(7)/2)/(b*x**2 + c*x**4)
    F = sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(5)/4)) - sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(5)/4)) + sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)) - sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)) + 2*sqrt(x)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_317():
    f = x**(sympy.S(5)/2)/(b*x**2 + c*x**4)
    F = sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(1)/4)*c**(sympy.S(3)/4)) - sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(1)/4)*c**(sympy.S(3)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*c**(sympy.S(3)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_318():
    f = x**(sympy.S(3)/2)/(b*x**2 + c*x**4)
    F = -sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(3)/4)*c**(sympy.S(1)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_319():
    f = sqrt(x)/(b*x**2 + c*x**4)
    F = -2/(b*sqrt(x)) - sqrt(2)*c**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(5)/4)) + sqrt(2)*c**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(5)/4)) + sqrt(2)*c**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)) - sqrt(2)*c**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_320():
    f = 1/(sqrt(x)*(b*x**2 + c*x**4))
    F = -2/(3*b*x**(sympy.S(3)/2)) + sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(7)/4)) - sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(7)/4)) + sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)) - sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_321():
    f = 1/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4))
    F = -2/(5*b*x**(sympy.S(5)/2)) + 2*c/(b**2*sqrt(x)) + sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(9)/4)) - sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(9)/4)) - sqrt(2)*c**(sympy.S(5)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) + sqrt(2)*c**(sympy.S(5)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_322():
    f = 1/(x**(sympy.S(5)/2)*(b*x**2 + c*x**4))
    F = -2/(7*b*x**(sympy.S(7)/2)) + 2*c/(3*b**2*x**(sympy.S(3)/2)) - sqrt(2)*c**(sympy.S(7)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(11)/4)) + sqrt(2)*c**(sympy.S(7)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(11)/4)) - sqrt(2)*c**(sympy.S(7)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4)) + sqrt(2)*c**(sympy.S(7)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_323():
    f = 1/(x**(sympy.S(7)/2)*(b*x**2 + c*x**4))
    F = -2/(9*b*x**(sympy.S(9)/2)) + 2*c/(5*b**2*x**(sympy.S(5)/2)) - 2*c**2/(b**3*sqrt(x)) - sqrt(2)*c**(sympy.S(9)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(9)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(9)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(9)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_324():
    f = x**(sympy.S(19)/2)/(b*x**2 + c*x**4)**2
    F = -9*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(13)/4)) + 9*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(13)/4)) - 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)) + 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)) - 9*b*sqrt(x)/(2*c**3) - x**(sympy.S(9)/2)/(2*c*(b + c*x**2)) + 9*x**(sympy.S(5)/2)/(10*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_325():
    f = x**(sympy.S(17)/2)/(b*x**2 + c*x**4)**2
    F = -7*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(11)/4)) + 7*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(11)/4)) + 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)) - 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)) - x**(sympy.S(7)/2)/(2*c*(b + c*x**2)) + 7*x**(sympy.S(3)/2)/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_326():
    f = x**(sympy.S(15)/2)/(b*x**2 + c*x**4)**2
    F = 5*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(9)/4)) - 5*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(9)/4)) + 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)) - 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)) - x**(sympy.S(5)/2)/(2*c*(b + c*x**2)) + 5*sqrt(x)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_327():
    f = x**(sympy.S(13)/2)/(b*x**2 + c*x**4)**2
    F = -x**(sympy.S(3)/2)/(2*c*(b + c*x**2)) + 3*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) - 3*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_328():
    f = x**(sympy.S(11)/2)/(b*x**2 + c*x**4)**2
    F = -sqrt(x)/(2*c*(b + c*x**2)) - sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) + sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(3)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_329():
    f = x**(sympy.S(9)/2)/(b*x**2 + c*x**4)**2
    F = x**(sympy.S(3)/2)/(2*b*(b + c*x**2)) + sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) - sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(5)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_330():
    f = x**(sympy.S(7)/2)/(b*x**2 + c*x**4)**2
    F = sqrt(x)/(2*b*(b + c*x**2)) - 3*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(7)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_331():
    f = x**(sympy.S(5)/2)/(b*x**2 + c*x**4)**2
    F = 1/(2*b*sqrt(x)*(b + c*x**2)) - 5/(2*b**2*sqrt(x)) - 5*sqrt(2)*c**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(9)/4)) + 5*sqrt(2)*c**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(9)/4)) + 5*sqrt(2)*c**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4)) - 5*sqrt(2)*c**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_332():
    f = x**(sympy.S(3)/2)/(b*x**2 + c*x**4)**2
    F = 1/(2*b*x**(sympy.S(3)/2)*(b + c*x**2)) - 7/(6*b**2*x**(sympy.S(3)/2)) + 7*sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(11)/4)) - 7*sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(11)/4)) + 7*sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(11)/4)) - 7*sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_333():
    f = sqrt(x)/(b*x**2 + c*x**4)**2
    F = 1/(2*b*x**(sympy.S(5)/2)*(b + c*x**2)) - 9/(10*b**2*x**(sympy.S(5)/2)) + 9*c/(2*b**3*sqrt(x)) + 9*sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(13)/4)) - 9*sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(13)/4)) - 9*sqrt(2)*c**(sympy.S(5)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4)) + 9*sqrt(2)*c**(sympy.S(5)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_334():
    f = 1/(sqrt(x)*(b*x**2 + c*x**4)**2)
    F = 1/(2*b*x**(sympy.S(7)/2)*(b + c*x**2)) - 11/(14*b**2*x**(sympy.S(7)/2)) + 11*c/(6*b**3*x**(sympy.S(3)/2)) - 11*sqrt(2)*c**(sympy.S(7)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(15)/4)) + 11*sqrt(2)*c**(sympy.S(7)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(15)/4)) - 11*sqrt(2)*c**(sympy.S(7)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(15)/4)) + 11*sqrt(2)*c**(sympy.S(7)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_335():
    f = 1/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**2)
    F = 1/(2*b*x**(sympy.S(9)/2)*(b + c*x**2)) - 13/(18*b**2*x**(sympy.S(9)/2)) + 13*c/(10*b**3*x**(sympy.S(5)/2)) - 13*c**2/(2*b**4*sqrt(x)) - 13*sqrt(2)*c**(sympy.S(9)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(17)/4)) + 13*sqrt(2)*c**(sympy.S(9)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(17)/4)) + 13*sqrt(2)*c**(sympy.S(9)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(17)/4)) - 13*sqrt(2)*c**(sympy.S(9)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_336():
    f = x**(sympy.S(23)/2)/(b*x**2 + c*x**4)**3
    F = 45*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*c**(sympy.S(13)/4)) - 45*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*c**(sympy.S(13)/4)) + 45*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)) - 45*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)) - x**(sympy.S(9)/2)/(4*c*(b + c*x**2)**2) - 9*x**(sympy.S(5)/2)/(16*c**2*(b + c*x**2)) + 45*sqrt(x)/(16*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_337():
    f = x**(sympy.S(21)/2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(7)/2)/(4*c*(b + c*x**2)**2) - 7*x**(sympy.S(3)/2)/(16*c**2*(b + c*x**2)) + 21*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) - 21*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) - 21*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) + 21*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(1)/4)*c**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_338():
    f = x**(sympy.S(19)/2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(5)/2)/(4*c*(b + c*x**2)**2) - 5*sqrt(x)/(16*c**2*(b + c*x**2)) - 5*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) + 5*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) - 5*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) + 5*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(3)/4)*c**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_339():
    f = x**(sympy.S(17)/2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(3)/2)/(4*c*(b + c*x**2)**2) + 3*x**(sympy.S(3)/2)/(16*b*c*(b + c*x**2)) + 3*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) - 3*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(5)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_340():
    f = x**(sympy.S(15)/2)/(b*x**2 + c*x**4)**3
    F = -sqrt(x)/(4*c*(b + c*x**2)**2) + sqrt(x)/(16*b*c*(b + c*x**2)) - 3*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + 3*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(7)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_341():
    f = x**(sympy.S(13)/2)/(b*x**2 + c*x**4)**3
    F = x**(sympy.S(3)/2)/(4*b*(b + c*x**2)**2) + 5*x**(sympy.S(3)/2)/(16*b**2*(b + c*x**2)) + 5*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) - 5*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) - 5*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) + 5*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(9)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_342():
    f = x**(sympy.S(11)/2)/(b*x**2 + c*x**4)**3
    F = sqrt(x)/(4*b*(b + c*x**2)**2) + 7*sqrt(x)/(16*b**2*(b + c*x**2)) - 21*sqrt(2)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) - 21*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(11)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_343():
    f = x**(sympy.S(9)/2)/(b*x**2 + c*x**4)**3
    F = 1/(4*b*sqrt(x)*(b + c*x**2)**2) + 9/(16*b**2*sqrt(x)*(b + c*x**2)) - 45/(16*b**3*sqrt(x)) - 45*sqrt(2)*c**(sympy.S(1)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(13)/4)) + 45*sqrt(2)*c**(sympy.S(1)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(13)/4)) + 45*sqrt(2)*c**(sympy.S(1)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(13)/4)) - 45*sqrt(2)*c**(sympy.S(1)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_344():
    f = x**(sympy.S(7)/2)/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x**(sympy.S(3)/2)*(b + c*x**2)**2) + 11/(16*b**2*x**(sympy.S(3)/2)*(b + c*x**2)) - 77/(48*b**3*x**(sympy.S(3)/2)) + 77*sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(15)/4)) - 77*sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(15)/4)) + 77*sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(15)/4)) - 77*sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_345():
    f = x**(sympy.S(5)/2)/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x**(sympy.S(5)/2)*(b + c*x**2)**2) + 13/(16*b**2*x**(sympy.S(5)/2)*(b + c*x**2)) - 117/(80*b**3*x**(sympy.S(5)/2)) + 117*c/(16*b**4*sqrt(x)) + 117*sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(17)/4)) - 117*sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(17)/4)) - 117*sqrt(2)*c**(sympy.S(5)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(17)/4)) + 117*sqrt(2)*c**(sympy.S(5)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_346():
    f = x**(sympy.S(3)/2)/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x**(sympy.S(7)/2)*(b + c*x**2)**2) + 15/(16*b**2*x**(sympy.S(7)/2)*(b + c*x**2)) - 165/(112*b**3*x**(sympy.S(7)/2)) + 55*c/(16*b**4*x**(sympy.S(3)/2)) - 165*sqrt(2)*c**(sympy.S(7)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(19)/4)) + 165*sqrt(2)*c**(sympy.S(7)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(19)/4)) - 165*sqrt(2)*c**(sympy.S(7)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(19)/4)) + 165*sqrt(2)*c**(sympy.S(7)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(19)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_347():
    f = sqrt(x)/(b*x**2 + c*x**4)**3
    F = 1/(4*b*x**(sympy.S(9)/2)*(b + c*x**2)**2) + 17/(16*b**2*x**(sympy.S(9)/2)*(b + c*x**2)) - 221/(144*b**3*x**(sympy.S(9)/2)) + 221*c/(80*b**4*x**(sympy.S(5)/2)) - 221*c**2/(16*b**5*sqrt(x)) - 221*sqrt(2)*c**(sympy.S(9)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(21)/4)) + 221*sqrt(2)*c**(sympy.S(9)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(21)/4)) + 221*sqrt(2)*c**(sympy.S(9)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(21)/4)) - 221*sqrt(2)*c**(sympy.S(9)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(21)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_348():
    f = 1/(sqrt(x)*(b*x**2 + c*x**4)**3)
    F = 1/(4*b*x**(sympy.S(11)/2)*(b + c*x**2)**2) + 19/(16*b**2*x**(sympy.S(11)/2)*(b + c*x**2)) - 285/(176*b**3*x**(sympy.S(11)/2)) + 285*c/(112*b**4*x**(sympy.S(7)/2)) - 95*c**2/(16*b**5*x**(sympy.S(3)/2)) + 285*sqrt(2)*c**(sympy.S(11)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(23)/4)) - 285*sqrt(2)*c**(sympy.S(11)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(23)/4)) + 285*sqrt(2)*c**(sympy.S(11)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(23)/4)) - 285*sqrt(2)*c**(sympy.S(11)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(23)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_349():
    f = x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)
    F = -28*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 14*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 28*b**3*x**(sympy.S(3)/2)*(b + c*x**2)/(195*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 28*b**2*sqrt(x)*sqrt(b*x**2 + c*x**4)/(585*c**2) + 4*b*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/(117*c) + 2*x**(sympy.S(9)/2)*sqrt(b*x**2 + c*x**4)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_350():
    f = x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)
    F = 10*b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - 20*b**2*sqrt(b*x**2 + c*x**4)/(231*c**2*sqrt(x)) + 4*b*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(77*c) + 2*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_351():
    f = x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)
    F = 4*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 2*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 4*b**2*x**(sympy.S(3)/2)*(b + c*x**2)/(15*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 4*b*sqrt(x)*sqrt(b*x**2 + c*x**4)/(45*c) + 2*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_352():
    f = sqrt(x)*sqrt(b*x**2 + c*x**4)
    F = -2*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4)) + 4*b*sqrt(b*x**2 + c*x**4)/(21*c*sqrt(x)) + 2*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_353():
    f = sqrt(b*x**2 + c*x**4)/sqrt(x)
    F = -4*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 2*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 4*b*x**(sympy.S(3)/2)*(b + c*x**2)/(5*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 2*sqrt(x)*sqrt(b*x**2 + c*x**4)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_354():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    F = 2*b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4)) + 2*sqrt(b*x**2 + c*x**4)/(3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_355():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(5)/2)
    F = -4*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/sqrt(b*x**2 + c*x**4) + 2*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/sqrt(b*x**2 + c*x**4) + 4*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)/((sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 2*sqrt(b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_356():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(7)/2)
    F = -2*sqrt(b*x**2 + c*x**4)/(3*x**(sympy.S(5)/2)) + 2*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_357():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(9)/2)
    F = -2*sqrt(b*x**2 + c*x**4)/(5*x**(sympy.S(7)/2)) + 4*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(5*b*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 4*c*sqrt(b*x**2 + c*x**4)/(5*b*x**(sympy.S(3)/2)) - 4*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 2*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_358():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(11)/2)
    F = -2*sqrt(b*x**2 + c*x**4)/(7*x**(sympy.S(9)/2)) - 4*c*sqrt(b*x**2 + c*x**4)/(21*b*x**(sympy.S(5)/2)) - 2*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_359():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(13)/2)
    F = -2*sqrt(b*x**2 + c*x**4)/(9*x**(sympy.S(11)/2)) - 4*c*sqrt(b*x**2 + c*x**4)/(45*b*x**(sympy.S(7)/2)) - 4*c**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(15*b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 4*c**2*sqrt(b*x**2 + c*x**4)/(15*b**2*x**(sympy.S(3)/2)) + 4*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 2*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_360():
    f = sqrt(b*x**2 + c*x**4)/x**(sympy.S(15)/2)
    F = -2*sqrt(b*x**2 + c*x**4)/(11*x**(sympy.S(13)/2)) - 4*c*sqrt(b*x**2 + c*x**4)/(77*b*x**(sympy.S(9)/2)) + 20*c**2*sqrt(b*x**2 + c*x**4)/(231*b**2*x**(sympy.S(5)/2)) + 10*c**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_361():
    f = x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -56*b**(sympy.S(17)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 28*b**(sympy.S(17)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 56*b**4*x**(sympy.S(3)/2)*(b + c*x**2)/(1105*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 56*b**3*sqrt(x)*sqrt(b*x**2 + c*x**4)/(3315*c**2) + 8*b**2*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/(663*c) + 12*b*x**(sympy.S(9)/2)*sqrt(b*x**2 + c*x**4)/221 + 2*x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_362():
    f = sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 4*b**(sympy.S(15)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - 8*b**3*sqrt(b*x**2 + c*x**4)/(231*c**2*sqrt(x)) + 8*b**2*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(385*c) + 4*b*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)/55 + 2*x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_363():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/sqrt(x)
    F = 8*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 4*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 8*b**3*x**(sympy.S(3)/2)*(b + c*x**2)/(65*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 8*b**2*sqrt(x)*sqrt(b*x**2 + c*x**4)/(195*c) + 4*b*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/39 + 2*sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_364():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = -4*b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4)) + 8*b**2*sqrt(b*x**2 + c*x**4)/(77*c*sqrt(x)) + 12*b*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/77 + 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_365():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = -8*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 4*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 8*b**2*x**(sympy.S(3)/2)*(b + c*x**2)/(15*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 4*b*sqrt(x)*sqrt(b*x**2 + c*x**4)/15 + 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_366():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(7)/2)
    F = 4*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(7*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4)) + 4*b*sqrt(b*x**2 + c*x**4)/(7*sqrt(x)) + 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_367():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(9)/2)
    F = -24*b**(sympy.S(5)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 12*b**(sympy.S(5)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 24*b*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)/((5*sqrt(b) + 5*sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 12*c*sqrt(x)*sqrt(b*x**2 + c*x**4)/5 - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_368():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(11)/2)
    F = 4*b**(sympy.S(3)/4)*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*sqrt(b*x**2 + c*x**4)) + 4*c*sqrt(b*x**2 + c*x**4)/(3*sqrt(x)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_369():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(13)/2)
    F = -24*b**(sympy.S(1)/4)*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 12*b**(sympy.S(1)/4)*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 24*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/((5*sqrt(b) + 5*sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 12*c*sqrt(b*x**2 + c*x**4)/(5*x**(sympy.S(3)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*x**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_370():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(15)/2)
    F = -4*c*sqrt(b*x**2 + c*x**4)/(7*x**(sympy.S(5)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*x**(sympy.S(13)/2)) + 4*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(7*b**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_371():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(17)/2)
    F = -4*c*sqrt(b*x**2 + c*x**4)/(15*x**(sympy.S(7)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*x**(sympy.S(15)/2)) + 8*c**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(15*b*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 8*c**2*sqrt(b*x**2 + c*x**4)/(15*b*x**(sympy.S(3)/2)) - 8*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 4*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_372():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(19)/2)
    F = -12*c*sqrt(b*x**2 + c*x**4)/(77*x**(sympy.S(9)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*x**(sympy.S(17)/2)) - 8*c**2*sqrt(b*x**2 + c*x**4)/(77*b*x**(sympy.S(5)/2)) - 4*c**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*b**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_373():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(21)/2)
    F = -4*c*sqrt(b*x**2 + c*x**4)/(39*x**(sympy.S(11)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(13*x**(sympy.S(19)/2)) - 8*c**2*sqrt(b*x**2 + c*x**4)/(195*b*x**(sympy.S(7)/2)) - 8*c**(sympy.S(7)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(65*b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 8*c**3*sqrt(b*x**2 + c*x**4)/(65*b**2*x**(sympy.S(3)/2)) + 8*c**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 4*c**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_374():
    f = (b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(23)/2)
    F = -4*c*sqrt(b*x**2 + c*x**4)/(55*x**(sympy.S(13)/2)) - 2*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*x**(sympy.S(21)/2)) - 8*c**2*sqrt(b*x**2 + c*x**4)/(385*b*x**(sympy.S(9)/2)) + 8*c**3*sqrt(b*x**2 + c*x**4)/(231*b**2*x**(sympy.S(5)/2)) + 4*c**(sympy.S(15)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_375():
    f = x**(sympy.S(13)/2)/sqrt(b*x**2 + c*x**4)
    F = -15*b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) + 30*b**2*sqrt(b*x**2 + c*x**4)/(77*c**3*sqrt(x)) - 18*b*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(77*c**2) + 2*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)/(11*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_376():
    f = x**(sympy.S(11)/2)/sqrt(b*x**2 + c*x**4)
    F = -14*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 7*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 14*b**2*x**(sympy.S(3)/2)*(b + c*x**2)/(15*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 14*b*sqrt(x)*sqrt(b*x**2 + c*x**4)/(45*c**2) + 2*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/(9*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_377():
    f = x**(sympy.S(9)/2)/sqrt(b*x**2 + c*x**4)
    F = 5*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - 10*b*sqrt(b*x**2 + c*x**4)/(21*c**2*sqrt(x)) + 2*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_378():
    f = x**(sympy.S(7)/2)/sqrt(b*x**2 + c*x**4)
    F = 6*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 3*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 6*b*x**(sympy.S(3)/2)*(b + c*x**2)/(5*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 2*sqrt(x)*sqrt(b*x**2 + c*x**4)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_379():
    f = x**(sympy.S(5)/2)/sqrt(b*x**2 + c*x**4)
    F = -b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4)) + 2*sqrt(b*x**2 + c*x**4)/(3*c*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_380():
    f = x**(sympy.S(3)/2)/sqrt(b*x**2 + c*x**4)
    F = -2*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 2*x**(sympy.S(3)/2)*(b + c*x**2)/(sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_381():
    f = sqrt(x)/sqrt(b*x**2 + c*x**4)
    F = x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_382():
    f = 1/(sqrt(x)*sqrt(b*x**2 + c*x**4))
    F = 2*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)/(b*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 2*sqrt(b*x**2 + c*x**4)/(b*x**(sympy.S(3)/2)) - 2*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_383():
    f = 1/(x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*sqrt(b*x**2 + c*x**4)/(3*b*x**(sympy.S(5)/2)) - c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_384():
    f = 1/(x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*sqrt(b*x**2 + c*x**4)/(5*b*x**(sympy.S(7)/2)) - 6*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(5*b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 6*c*sqrt(b*x**2 + c*x**4)/(5*b**2*x**(sympy.S(3)/2)) + 6*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 3*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_385():
    f = 1/(x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*sqrt(b*x**2 + c*x**4)/(7*b*x**(sympy.S(9)/2)) + 10*c*sqrt(b*x**2 + c*x**4)/(21*b**2*x**(sympy.S(5)/2)) + 5*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_386():
    f = 1/(x**(sympy.S(9)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*sqrt(b*x**2 + c*x**4)/(9*b*x**(sympy.S(11)/2)) + 14*c*sqrt(b*x**2 + c*x**4)/(45*b**2*x**(sympy.S(7)/2)) + 14*c**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(15*b**3*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 14*c**2*sqrt(b*x**2 + c*x**4)/(15*b**3*x**(sympy.S(3)/2)) - 14*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 7*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_387():
    f = 1/(x**(sympy.S(11)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*sqrt(b*x**2 + c*x**4)/(11*b*x**(sympy.S(13)/2)) + 18*c*sqrt(b*x**2 + c*x**4)/(77*b**2*x**(sympy.S(9)/2)) - 30*c**2*sqrt(b*x**2 + c*x**4)/(77*b**3*x**(sympy.S(5)/2)) - 15*c**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*b**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_388():
    f = x**(sympy.S(17)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 15*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(14*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) - 15*b*sqrt(b*x**2 + c*x**4)/(7*c**3*sqrt(x)) - x**(sympy.S(11)/2)/(c*sqrt(b*x**2 + c*x**4)) + 9*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(7*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_389():
    f = x**(sympy.S(15)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 21*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - 21*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - 21*b*x**(sympy.S(3)/2)*(b + c*x**2)/(5*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(9)/2)/(c*sqrt(b*x**2 + c*x**4)) + 7*sqrt(x)*sqrt(b*x**2 + c*x**4)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_390():
    f = x**(sympy.S(13)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -5*b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(6*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(7)/2)/(c*sqrt(b*x**2 + c*x**4)) + 5*sqrt(b*x**2 + c*x**4)/(3*c**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_391():
    f = x**(sympy.S(11)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -3*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) + 3*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(5)/2)/(c*sqrt(b*x**2 + c*x**4)) + 3*x**(sympy.S(3)/2)*(b + c*x**2)/(c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_392():
    f = x**(sympy.S(9)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**(sympy.S(3)/2)/(c*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(1)/4)*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_393():
    f = x**(sympy.S(7)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = x**(sympy.S(5)/2)/(b*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(3)/2)*(b + c*x**2)/(b*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_394():
    f = x**(sympy.S(5)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = x**(sympy.S(3)/2)/(b*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(5)/4)*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_395():
    f = x**(sympy.S(3)/2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = sqrt(x)/(b*sqrt(b*x**2 + c*x**4)) + 3*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)/(b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 3*sqrt(b*x**2 + c*x**4)/(b**2*x**(sympy.S(3)/2)) - 3*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) + 3*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_396():
    f = sqrt(x)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 1/(b*sqrt(x)*sqrt(b*x**2 + c*x**4)) - 5*sqrt(b*x**2 + c*x**4)/(3*b**2*x**(sympy.S(5)/2)) - 5*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_397():
    f = 1/(sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)) - 7*sqrt(b*x**2 + c*x**4)/(5*b**2*x**(sympy.S(7)/2)) - 21*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(5*b**3*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 21*c*sqrt(b*x**2 + c*x**4)/(5*b**3*x**(sympy.S(3)/2)) + 21*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - 21*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_398():
    f = 1/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)) - 9*sqrt(b*x**2 + c*x**4)/(7*b**2*x**(sympy.S(9)/2)) + 15*c*sqrt(b*x**2 + c*x**4)/(7*b**3*x**(sympy.S(5)/2)) + 15*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(14*b**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_399():
    f = 1/(x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 1/(b*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)) - 11*sqrt(b*x**2 + c*x**4)/(9*b**2*x**(sympy.S(11)/2)) + 77*c*sqrt(b*x**2 + c*x**4)/(45*b**3*x**(sympy.S(7)/2)) + 77*c**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(b + c*x**2)/(15*b**4*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 77*c**2*sqrt(b*x**2 + c*x**4)/(15*b**4*x**(sympy.S(3)/2)) - 77*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) + 77*c**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(30*b**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_400():
    f = (c*x)**m*(b*x**2 + c*x**4)**3
    F = b**3*x**7*(c*x)**m/(m + 7) + 3*b**2*c*x**9*(c*x)**m/(m + 9) + 3*b*c**2*x**11*(c*x)**m/(m + 11) + c**3*x**13*(c*x)**m/(m + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_401():
    f = (c*x)**m*(b*x**2 + c*x**4)**2
    F = b**2*x**5*(c*x)**m/(m + 5) + 2*b*c*x**7*(c*x)**m/(m + 7) + c**2*x**9*(c*x)**m/(m + 9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_402():
    f = (c*x)**m*(b*x**2 + c*x**4)
    F = b*(c*x)**(m + 3)/(c**3*(m + 3)) + (c*x)**(m + 5)/(c**4*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_403():
    f = (c*x)**m/(b*x**2 + c*x**4)
    F = -(c*x)**m*hyper((1, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), -c*x**2/b)/(b*x*(1 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_404():
    f = (c*x)**m/(b*x**2 + c*x**4)**2
    F = -(c*x)**m*hyper((2, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), -c*x**2/b)/(b**2*x**3*(3 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_405():
    f = (c*x)**m/(b*x**2 + c*x**4)**3
    F = -(c*x)**m*hyper((3, m/2 + sympy.S(-5)/2), (m/2 + sympy.S(-3)/2,), -c*x**2/b)/(b**3*x**5*(5 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_406():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**2*x**4/4 + a*b*x**6/3 + b**2*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_407():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**2*x**3/3 + 2*a*b*x**5/5 + b**2*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_408():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**2*x**2/2 + a*b*x**4/2 + b**2*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_409():
    f = a**2 + 2*a*b*x**2 + b**2*x**4
    F = a**2*x + 2*a*b*x**3/3 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_410():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x
    F = a**2*log(x) + a*b*x**2 + b**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_411():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**2
    F = -a**2/x + 2*a*b*x + b**2*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_412():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**3
    F = -a**2/(2*x**2) + 2*a*b*log(x) + b**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_413():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**4
    F = -a**2/(3*x**3) - 2*a*b/x + b**2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_414():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**5
    F = -a**2/(4*x**4) - a*b/x**2 + b**2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_415():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**6
    F = -a**2/(5*x**5) - 2*a*b/(3*x**3) - b**2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_416():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**7
    F = -a**2/(6*x**6) - a*b/(2*x**4) - b**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_417():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/x**8
    F = -a**2/(7*x**7) - 2*a*b/(5*x**5) - b**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_418():
    f = x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**4*x**7/7 + 4*a**3*b*x**9/9 + 6*a**2*b**2*x**11/11 + 4*a*b**3*x**13/13 + b**4*x**15/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_419():
    f = x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**2*(a + b*x**2)**5/(10*b**3) - a*(a + b*x**2)**6/(6*b**3) + (a + b*x**2)**7/(14*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_420():
    f = x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**4*x**5/5 + 4*a**3*b*x**7/7 + 2*a**2*b**2*x**9/3 + 4*a*b**3*x**11/11 + b**4*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_421():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -a*(a + b*x**2)**5/(10*b**2) + (a + b*x**2)**6/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_422():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**4*x**3/3 + 4*a**3*b*x**5/5 + 6*a**2*b**2*x**7/7 + 4*a*b**3*x**9/9 + b**4*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_423():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = (a + b*x**2)**5/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_424():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**4*x + 4*a**3*b*x**3/3 + 6*a**2*b**2*x**5/5 + 4*a*b**3*x**7/7 + b**4*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_425():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x
    F = a**4*log(x) + 2*a**3*b*x**2 + 3*a**2*b**2*x**4/2 + 2*a*b**3*x**6/3 + b**4*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_426():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**2
    F = -a**4/x + 4*a**3*b*x + 2*a**2*b**2*x**3 + 4*a*b**3*x**5/5 + b**4*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_427():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**3
    F = -a**4/(2*x**2) + 4*a**3*b*log(x) + 3*a**2*b**2*x**2 + a*b**3*x**4 + b**4*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_428():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**4
    F = -a**4/(3*x**3) - 4*a**3*b/x + 6*a**2*b**2*x + 4*a*b**3*x**3/3 + b**4*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_429():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**5
    F = -a**4/(4*x**4) - 2*a**3*b/x**2 + 6*a**2*b**2*log(x) + 2*a*b**3*x**2 + b**4*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_430():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**6
    F = -a**4/(5*x**5) - 4*a**3*b/(3*x**3) - 6*a**2*b**2/x + 4*a*b**3*x + b**4*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_431():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**7
    F = -a**4/(6*x**6) - a**3*b/x**4 - 3*a**2*b**2/x**2 + 4*a*b**3*log(x) + b**4*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_432():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**8
    F = -a**4/(7*x**7) - 4*a**3*b/(5*x**5) - 2*a**2*b**2/x**3 - 4*a*b**3/x + b**4*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_433():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**9
    F = -a**4/(8*x**8) - 2*a**3*b/(3*x**6) - 3*a**2*b**2/(2*x**4) - 2*a*b**3/x**2 + b**4*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_434():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**10
    F = -a**4/(9*x**9) - 4*a**3*b/(7*x**7) - 6*a**2*b**2/(5*x**5) - 4*a*b**3/(3*x**3) - b**4/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_435():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**11
    F = -(a + b*x**2)**5/(10*a*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_436():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**12
    F = -a**4/(11*x**11) - 4*a**3*b/(9*x**9) - 6*a**2*b**2/(7*x**7) - 4*a*b**3/(5*x**5) - b**4/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_437():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**13
    F = -(a + b*x**2)**5/(12*a*x**12) + b*(a + b*x**2)**5/(60*a**2*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_438():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**14
    F = -a**4/(13*x**13) - 4*a**3*b/(11*x**11) - 2*a**2*b**2/(3*x**9) - 4*a*b**3/(7*x**7) - b**4/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_439():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**15
    F = -a**4/(14*x**14) - a**3*b/(3*x**12) - 3*a**2*b**2/(5*x**10) - a*b**3/(2*x**8) - b**4/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_440():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/x**16
    F = -a**4/(15*x**15) - 4*a**3*b/(13*x**13) - 6*a**2*b**2/(11*x**11) - 4*a*b**3/(9*x**9) - b**4/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_441():
    f = x**8*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*x**9/9 + 6*a**5*b*x**11/11 + 15*a**4*b**2*x**13/13 + 4*a**3*b**3*x**15/3 + 15*a**2*b**4*x**17/17 + 6*a*b**5*x**19/19 + b**6*x**21/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_442():
    f = x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -a**3*(a + b*x**2)**7/(14*b**4) + 3*a**2*(a + b*x**2)**8/(16*b**4) - a*(a + b*x**2)**9/(6*b**4) + (a + b*x**2)**10/(20*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_443():
    f = x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*x**7/7 + 2*a**5*b*x**9/3 + 15*a**4*b**2*x**11/11 + 20*a**3*b**3*x**13/13 + a**2*b**4*x**15 + 6*a*b**5*x**17/17 + b**6*x**19/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_444():
    f = x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**2*(a + b*x**2)**7/(14*b**3) - a*(a + b*x**2)**8/(8*b**3) + (a + b*x**2)**9/(18*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_445():
    f = x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*x**5/5 + 6*a**5*b*x**7/7 + 5*a**4*b**2*x**9/3 + 20*a**3*b**3*x**11/11 + 15*a**2*b**4*x**13/13 + 2*a*b**5*x**15/5 + b**6*x**17/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_446():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -a*(a + b*x**2)**7/(14*b**2) + (a + b*x**2)**8/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_447():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*x**3/3 + 6*a**5*b*x**5/5 + 15*a**4*b**2*x**7/7 + 20*a**3*b**3*x**9/9 + 15*a**2*b**4*x**11/11 + 6*a*b**5*x**13/13 + b**6*x**15/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_448():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = (a + b*x**2)**7/(14*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_449():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*x + 2*a**5*b*x**3 + 3*a**4*b**2*x**5 + 20*a**3*b**3*x**7/7 + 5*a**2*b**4*x**9/3 + 6*a*b**5*x**11/11 + b**6*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_450():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x
    F = a**6*log(x) + 3*a**5*b*x**2 + 15*a**4*b**2*x**4/4 + 10*a**3*b**3*x**6/3 + 15*a**2*b**4*x**8/8 + 3*a*b**5*x**10/5 + b**6*x**12/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_451():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**2
    F = -a**6/x + 6*a**5*b*x + 5*a**4*b**2*x**3 + 4*a**3*b**3*x**5 + 15*a**2*b**4*x**7/7 + 2*a*b**5*x**9/3 + b**6*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_452():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**3
    F = -a**6/(2*x**2) + 6*a**5*b*log(x) + 15*a**4*b**2*x**2/2 + 5*a**3*b**3*x**4 + 5*a**2*b**4*x**6/2 + 3*a*b**5*x**8/4 + b**6*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_453():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**4
    F = -a**6/(3*x**3) - 6*a**5*b/x + 15*a**4*b**2*x + 20*a**3*b**3*x**3/3 + 3*a**2*b**4*x**5 + 6*a*b**5*x**7/7 + b**6*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_454():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**5
    F = -a**6/(4*x**4) - 3*a**5*b/x**2 + 15*a**4*b**2*log(x) + 10*a**3*b**3*x**2 + 15*a**2*b**4*x**4/4 + a*b**5*x**6 + b**6*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_455():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**6
    F = -a**6/(5*x**5) - 2*a**5*b/x**3 - 15*a**4*b**2/x + 20*a**3*b**3*x + 5*a**2*b**4*x**3 + 6*a*b**5*x**5/5 + b**6*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_456():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**7
    F = -a**6/(6*x**6) - 3*a**5*b/(2*x**4) - 15*a**4*b**2/(2*x**2) + 20*a**3*b**3*log(x) + 15*a**2*b**4*x**2/2 + 3*a*b**5*x**4/2 + b**6*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_457():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**8
    F = -a**6/(7*x**7) - 6*a**5*b/(5*x**5) - 5*a**4*b**2/x**3 - 20*a**3*b**3/x + 15*a**2*b**4*x + 2*a*b**5*x**3 + b**6*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_458():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**9
    F = -a**6/(8*x**8) - a**5*b/x**6 - 15*a**4*b**2/(4*x**4) - 10*a**3*b**3/x**2 + 15*a**2*b**4*log(x) + 3*a*b**5*x**2 + b**6*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_459():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**10
    F = -a**6/(9*x**9) - 6*a**5*b/(7*x**7) - 3*a**4*b**2/x**5 - 20*a**3*b**3/(3*x**3) - 15*a**2*b**4/x + 6*a*b**5*x + b**6*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_460():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**11
    F = -a**6/(10*x**10) - 3*a**5*b/(4*x**8) - 5*a**4*b**2/(2*x**6) - 5*a**3*b**3/x**4 - 15*a**2*b**4/(2*x**2) + 6*a*b**5*log(x) + b**6*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_461():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**12
    F = -a**6/(11*x**11) - 2*a**5*b/(3*x**9) - 15*a**4*b**2/(7*x**7) - 4*a**3*b**3/x**5 - 5*a**2*b**4/x**3 - 6*a*b**5/x + b**6*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_462():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**13
    F = -a**6/(12*x**12) - 3*a**5*b/(5*x**10) - 15*a**4*b**2/(8*x**8) - 10*a**3*b**3/(3*x**6) - 15*a**2*b**4/(4*x**4) - 3*a*b**5/x**2 + b**6*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_463():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**14
    F = -a**6/(13*x**13) - 6*a**5*b/(11*x**11) - 5*a**4*b**2/(3*x**9) - 20*a**3*b**3/(7*x**7) - 3*a**2*b**4/x**5 - 2*a*b**5/x**3 - b**6/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_464():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**15
    F = -(a + b*x**2)**7/(14*a*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_465():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**16
    F = -a**6/(15*x**15) - 6*a**5*b/(13*x**13) - 15*a**4*b**2/(11*x**11) - 20*a**3*b**3/(9*x**9) - 15*a**2*b**4/(7*x**7) - 6*a*b**5/(5*x**5) - b**6/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_466():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**17
    F = -(a + b*x**2)**7/(16*a*x**16) + b*(a + b*x**2)**7/(112*a**2*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_467():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**18
    F = -a**6/(17*x**17) - 2*a**5*b/(5*x**15) - 15*a**4*b**2/(13*x**13) - 20*a**3*b**3/(11*x**11) - 5*a**2*b**4/(3*x**9) - 6*a*b**5/(7*x**7) - b**6/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_468():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**19
    F = -(a + b*x**2)**7/(18*a*x**18) + b*(a + b*x**2)**7/(72*a**2*x**16) - b**2*(a + b*x**2)**7/(504*a**3*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_469():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**20
    F = -a**6/(19*x**19) - 6*a**5*b/(17*x**17) - a**4*b**2/x**15 - 20*a**3*b**3/(13*x**13) - 15*a**2*b**4/(11*x**11) - 2*a*b**5/(3*x**9) - b**6/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_470():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**21
    F = -(a + b*x**2)**7/(20*a*x**20) + b*(a + b*x**2)**7/(60*a**2*x**18) - b**2*(a + b*x**2)**7/(240*a**3*x**16) + b**3*(a + b*x**2)**7/(1680*a**4*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_471():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/x**22
    F = -a**6/(21*x**21) - 6*a**5*b/(19*x**19) - 15*a**4*b**2/(17*x**17) - 4*a**3*b**3/(3*x**15) - 15*a**2*b**4/(13*x**13) - 6*a*b**5/(11*x**11) - b**6/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_472():
    f = x**11/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**5/(2*b**6*(a + b*x**2)) + 5*a**4*log(a + b*x**2)/(2*b**6) - 2*a**3*x**2/b**5 + 3*a**2*x**4/(4*b**4) - a*x**6/(3*b**3) + x**8/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_473():
    f = x**9/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -a**4/(2*b**5*(a + b*x**2)) - 2*a**3*log(a + b*x**2)/b**5 + 3*a**2*x**2/(2*b**4) - a*x**4/(2*b**3) + x**6/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_474():
    f = x**7/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**3/(2*b**4*(a + b*x**2)) + 3*a**2*log(a + b*x**2)/(2*b**4) - a*x**2/b**3 + x**4/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_475():
    f = x**5/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -a**2/(2*b**3*(a + b*x**2)) - a*log(a + b*x**2)/b**3 + x**2/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_476():
    f = x**3/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a/(2*b**2*(a + b*x**2)) + log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_477():
    f = x/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -1/(2*b*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_478():
    f = 1/(x*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*(a + b*x**2)) + log(x)/a**2 - log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_479():
    f = 1/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -b/(2*a**2*(a + b*x**2)) - 1/(2*a**2*x**2) - 2*b*log(x)/a**3 + b*log(a + b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_480():
    f = 1/(x**5*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -1/(4*a**2*x**4) + b**2/(2*a**3*(a + b*x**2)) + b/(a**3*x**2) + 3*b**2*log(x)/a**4 - 3*b**2*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_481():
    f = x**10/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 9*a**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(11)/2)) - 9*a**3*x/(2*b**5) + 3*a**2*x**3/(2*b**4) - 9*a*x**5/(10*b**3) - x**9/(2*b*(a + b*x**2)) + 9*x**7/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_482():
    f = x**8/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -7*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(9)/2)) + 7*a**2*x/(2*b**4) - 7*a*x**3/(6*b**3) - x**7/(2*b*(a + b*x**2)) + 7*x**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_483():
    f = x**6/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 5*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) - 5*a*x/(2*b**3) - x**5/(2*b*(a + b*x**2)) + 5*x**3/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_484():
    f = x**4/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -3*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(5)/2)) - x**3/(2*b*(a + b*x**2)) + 3*x/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_485():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -x/(2*b*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_486():
    f = 1/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = x/(2*a*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_487():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*x*(a + b*x**2)) - 3/(2*a**2*x) - 3*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_488():
    f = 1/(x**4*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*x**3*(a + b*x**2)) - 5/(6*a**2*x**3) + 5*b/(2*a**3*x) + 5*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_489():
    f = 1/(x**6*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*x**5*(a + b*x**2)) - 7/(10*a**2*x**5) + 7*b/(6*a**3*x**3) - 7*b**2/(2*a**4*x) - 7*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_490():
    f = x**11/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**5/(6*b**6*(a + b*x**2)**3) - 5*a**4/(4*b**6*(a + b*x**2)**2) + 5*a**3/(b**6*(a + b*x**2)) + 5*a**2*log(a + b*x**2)/b**6 - 2*a*x**2/b**5 + x**4/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_491():
    f = x**9/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -a**4/(6*b**5*(a + b*x**2)**3) + a**3/(b**5*(a + b*x**2)**2) - 3*a**2/(b**5*(a + b*x**2)) - 2*a*log(a + b*x**2)/b**5 + x**2/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_492():
    f = x**7/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**3/(6*b**4*(a + b*x**2)**3) - 3*a**2/(4*b**4*(a + b*x**2)**2) + 3*a/(2*b**4*(a + b*x**2)) + log(a + b*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_493():
    f = x**5/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = x**6/(6*a*(a + b*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_494():
    f = x**3/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a/(6*b**2*(a + b*x**2)**3) - 1/(4*b**2*(a + b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_495():
    f = x/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -1/(6*b*(a + b*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_496():
    f = 1/(x*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*(a + b*x**2)**3) + 1/(4*a**2*(a + b*x**2)**2) + 1/(2*a**3*(a + b*x**2)) + log(x)/a**4 - log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_497():
    f = 1/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = -b/(6*a**2*(a + b*x**2)**3) - b/(2*a**3*(a + b*x**2)**2) - 3*b/(2*a**4*(a + b*x**2)) - 1/(2*a**4*x**2) - 4*b*log(x)/a**5 + 2*b*log(a + b*x**2)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_498():
    f = 1/(x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = b**2/(6*a**3*(a + b*x**2)**3) + 3*b**2/(4*a**4*(a + b*x**2)**2) - 1/(4*a**4*x**4) + 3*b**2/(a**5*(a + b*x**2)) + 2*b/(a**5*x**2) + 10*b**2*log(x)/a**6 - 5*b**2*log(a + b*x**2)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_499():
    f = x**12/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -231*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(16*b**(sympy.S(13)/2)) + 231*a**2*x/(16*b**6) - 77*a*x**3/(16*b**5) - x**11/(6*b*(a + b*x**2)**3) - 11*x**9/(24*b**2*(a + b*x**2)**2) - 33*x**7/(16*b**3*(a + b*x**2)) + 231*x**5/(80*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_500():
    f = x**10/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = 105*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(16*b**(sympy.S(11)/2)) - 105*a*x/(16*b**5) - x**9/(6*b*(a + b*x**2)**3) - 3*x**7/(8*b**2*(a + b*x**2)**2) - 21*x**5/(16*b**3*(a + b*x**2)) + 35*x**3/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_501():
    f = x**8/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -35*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(16*b**(sympy.S(9)/2)) - x**7/(6*b*(a + b*x**2)**3) - 7*x**5/(24*b**2*(a + b*x**2)**2) - 35*x**3/(48*b**3*(a + b*x**2)) + 35*x/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_502():
    f = x**6/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -x**5/(6*b*(a + b*x**2)**3) - 5*x**3/(24*b**2*(a + b*x**2)**2) - 5*x/(16*b**3*(a + b*x**2)) + 5*atan(sqrt(b)*x/sqrt(a))/(16*sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_503():
    f = x**4/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -x**3/(6*b*(a + b*x**2)**3) - x/(8*b**2*(a + b*x**2)**2) + x/(16*a*b**2*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_504():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -x/(6*b*(a + b*x**2)**3) + x/(24*a*b*(a + b*x**2)**2) + x/(16*a**2*b*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_505():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(-2)
    F = x/(6*a*(a + b*x**2)**3) + 5*x/(24*a**2*(a + b*x**2)**2) + 5*x/(16*a**3*(a + b*x**2)) + 5*atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_506():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*x*(a + b*x**2)**3) + 7/(24*a**2*x*(a + b*x**2)**2) + 35/(48*a**3*x*(a + b*x**2)) - 35/(16*a**4*x) - 35*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_507():
    f = 1/(x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*x**3*(a + b*x**2)**3) + 3/(8*a**2*x**3*(a + b*x**2)**2) + 21/(16*a**3*x**3*(a + b*x**2)) - 35/(16*a**4*x**3) + 105*b/(16*a**5*x) + 105*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_508():
    f = 1/(x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*x**5*(a + b*x**2)**3) + 11/(24*a**2*x**5*(a + b*x**2)**2) + 33/(16*a**3*x**5*(a + b*x**2)) - 231/(80*a**4*x**5) + 77*b/(16*a**5*x**3) - 231*b**2/(16*a**6*x) - 231*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(16*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_509():
    f = x**15/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**7/(10*b**8*(a + b*x**2)**5) - 7*a**6/(8*b**8*(a + b*x**2)**4) + 7*a**5/(2*b**8*(a + b*x**2)**3) - 35*a**4/(4*b**8*(a + b*x**2)**2) + 35*a**3/(2*b**8*(a + b*x**2)) + 21*a**2*log(a + b*x**2)/(2*b**8) - 3*a*x**2/b**7 + x**4/(4*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_510():
    f = x**13/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -a**6/(10*b**7*(a + b*x**2)**5) + 3*a**5/(4*b**7*(a + b*x**2)**4) - 5*a**4/(2*b**7*(a + b*x**2)**3) + 5*a**3/(b**7*(a + b*x**2)**2) - 15*a**2/(2*b**7*(a + b*x**2)) - 3*a*log(a + b*x**2)/b**7 + x**2/(2*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_511():
    f = x**11/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**5/(10*b**6*(a + b*x**2)**5) - 5*a**4/(8*b**6*(a + b*x**2)**4) + 5*a**3/(3*b**6*(a + b*x**2)**3) - 5*a**2/(2*b**6*(a + b*x**2)**2) + 5*a/(2*b**6*(a + b*x**2)) + log(a + b*x**2)/(2*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_512():
    f = x**9/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = x**10/(10*a*(a + b*x**2)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_513():
    f = x**7/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = x**8/(10*a*(a + b*x**2)**5) + x**8/(40*a**2*(a + b*x**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_514():
    f = x**5/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -a**2/(10*b**3*(a + b*x**2)**5) + a/(4*b**3*(a + b*x**2)**4) - 1/(6*b**3*(a + b*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_515():
    f = x**3/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a/(10*b**2*(a + b*x**2)**5) - 1/(8*b**2*(a + b*x**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_516():
    f = x/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -1/(10*b*(a + b*x**2)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_517():
    f = 1/(x*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*(a + b*x**2)**5) + 1/(8*a**2*(a + b*x**2)**4) + 1/(6*a**3*(a + b*x**2)**3) + 1/(4*a**4*(a + b*x**2)**2) + 1/(2*a**5*(a + b*x**2)) + log(x)/a**6 - log(a + b*x**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_518():
    f = 1/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = -b/(10*a**2*(a + b*x**2)**5) - b/(4*a**3*(a + b*x**2)**4) - b/(2*a**4*(a + b*x**2)**3) - b/(a**5*(a + b*x**2)**2) - 5*b/(2*a**6*(a + b*x**2)) - 1/(2*a**6*x**2) - 6*b*log(x)/a**7 + 3*b*log(a + b*x**2)/a**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_519():
    f = 1/(x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = b**2/(10*a**3*(a + b*x**2)**5) + 3*b**2/(8*a**4*(a + b*x**2)**4) + b**2/(a**5*(a + b*x**2)**3) + 5*b**2/(2*a**6*(a + b*x**2)**2) - 1/(4*a**6*x**4) + 15*b**2/(2*a**7*(a + b*x**2)) + 3*b/(a**7*x**2) + 21*b**2*log(x)/a**8 - 21*b**2*log(a + b*x**2)/(2*a**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_520():
    f = x**16/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -9009*a**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(256*b**(sympy.S(17)/2)) + 9009*a**2*x/(256*b**8) - 3003*a*x**3/(256*b**7) - x**15/(10*b*(a + b*x**2)**5) - 3*x**13/(16*b**2*(a + b*x**2)**4) - 13*x**11/(32*b**3*(a + b*x**2)**3) - 143*x**9/(128*b**4*(a + b*x**2)**2) - 1287*x**7/(256*b**5*(a + b*x**2)) + 9009*x**5/(1280*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_521():
    f = x**14/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = 3003*a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(256*b**(sympy.S(15)/2)) - 3003*a*x/(256*b**7) - x**13/(10*b*(a + b*x**2)**5) - 13*x**11/(80*b**2*(a + b*x**2)**4) - 143*x**9/(480*b**3*(a + b*x**2)**3) - 429*x**7/(640*b**4*(a + b*x**2)**2) - 3003*x**5/(1280*b**5*(a + b*x**2)) + 1001*x**3/(256*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_522():
    f = x**12/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -693*sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(256*b**(sympy.S(13)/2)) - x**11/(10*b*(a + b*x**2)**5) - 11*x**9/(80*b**2*(a + b*x**2)**4) - 33*x**7/(160*b**3*(a + b*x**2)**3) - 231*x**5/(640*b**4*(a + b*x**2)**2) - 231*x**3/(256*b**5*(a + b*x**2)) + 693*x/(256*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_523():
    f = x**10/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -x**9/(10*b*(a + b*x**2)**5) - 9*x**7/(80*b**2*(a + b*x**2)**4) - 21*x**5/(160*b**3*(a + b*x**2)**3) - 21*x**3/(128*b**4*(a + b*x**2)**2) - 63*x/(256*b**5*(a + b*x**2)) + 63*atan(sqrt(b)*x/sqrt(a))/(256*sqrt(a)*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_524():
    f = x**8/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -x**7/(10*b*(a + b*x**2)**5) - 7*x**5/(80*b**2*(a + b*x**2)**4) - 7*x**3/(96*b**3*(a + b*x**2)**3) - 7*x/(128*b**4*(a + b*x**2)**2) + 7*x/(256*a*b**4*(a + b*x**2)) + 7*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(3)/2)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_525():
    f = x**6/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -x**5/(10*b*(a + b*x**2)**5) - x**3/(16*b**2*(a + b*x**2)**4) - x/(32*b**3*(a + b*x**2)**3) + x/(128*a*b**3*(a + b*x**2)**2) + 3*x/(256*a**2*b**3*(a + b*x**2)) + 3*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(5)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_526():
    f = x**4/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -x**3/(10*b*(a + b*x**2)**5) - 3*x/(80*b**2*(a + b*x**2)**4) + x/(160*a*b**2*(a + b*x**2)**3) + x/(128*a**2*b**2*(a + b*x**2)**2) + 3*x/(256*a**3*b**2*(a + b*x**2)) + 3*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(7)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_527():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -x/(10*b*(a + b*x**2)**5) + x/(80*a*b*(a + b*x**2)**4) + 7*x/(480*a**2*b*(a + b*x**2)**3) + 7*x/(384*a**3*b*(a + b*x**2)**2) + 7*x/(256*a**4*b*(a + b*x**2)) + 7*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(9)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_528():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(-3)
    F = x/(10*a*(a + b*x**2)**5) + 9*x/(80*a**2*(a + b*x**2)**4) + 21*x/(160*a**3*(a + b*x**2)**3) + 21*x/(128*a**4*(a + b*x**2)**2) + 63*x/(256*a**5*(a + b*x**2)) + 63*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(11)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_529():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*x*(a + b*x**2)**5) + 11/(80*a**2*x*(a + b*x**2)**4) + 33/(160*a**3*x*(a + b*x**2)**3) + 231/(640*a**4*x*(a + b*x**2)**2) + 231/(256*a**5*x*(a + b*x**2)) - 693/(256*a**6*x) - 693*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_530():
    f = 1/(x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*x**3*(a + b*x**2)**5) + 13/(80*a**2*x**3*(a + b*x**2)**4) + 143/(480*a**3*x**3*(a + b*x**2)**3) + 429/(640*a**4*x**3*(a + b*x**2)**2) + 3003/(1280*a**5*x**3*(a + b*x**2)) - 1001/(256*a**6*x**3) + 3003*b/(256*a**7*x) + 3003*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_531():
    f = 1/(x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*x**5*(a + b*x**2)**5) + 3/(16*a**2*x**5*(a + b*x**2)**4) + 13/(32*a**3*x**5*(a + b*x**2)**3) + 143/(128*a**4*x**5*(a + b*x**2)**2) + 1287/(256*a**5*x**5*(a + b*x**2)) - 9009/(1280*a**6*x**5) + 3003*b/(256*a**7*x**3) - 9009*b**2/(256*a**8*x) - 9009*b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_532():
    f = 1/(x**4 + 2*x**2 + 1)
    F = x/(2*x**2 + 2) + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_533():
    f = x/(x**4 + 2*x**2 + 1)
    F = -1/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_534():
    f = x**2/(x**4 + 2*x**2 + 1)
    F = -x/(2*x**2 + 2) + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_535():
    f = x**3/(x**4 + 2*x**2 + 1)
    F = log(x**2 + 1)/2 + 1/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_536():
    f = x/(x**4 - 18*x**2 + 81)
    F = 1/(18 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_537():
    f = x**3/(x**4 - 8*x**2 + 16)
    F = log(4 - x**2)/2 + 2/(4 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_538():
    f = x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*x**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*a + 6*b*x**2) + b*x**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a + 8*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_539():
    f = x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -a*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*b**2) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_540():
    f = x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_541():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x
    F = a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + b*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_542():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**3
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_543():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**5
    F = -(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_544():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**7
    F = -(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a*x**6) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(12*a**2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_545():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**9
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**6*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_546():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**11
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_547():
    f = x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2) + b*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_548():
    f = x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + b*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_549():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_550():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**2
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + b*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_551():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**4
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_552():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**6
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_553():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**8
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_554():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**10
    F = -a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_555():
    f = x**9*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**10*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*a + 10*b*x**2) + a**2*b*x**12*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2) + 3*a*b**2*x**14*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*a + 14*b*x**2) + b**3*x**16*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*a + 16*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_556():
    f = x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a + 8*b*x**2) + 3*a**2*b*x**10*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*a + 10*b*x**2) + a*b**2*x**12*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2) + b**3*x**14*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*a + 14*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_557():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -a*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(8*b**2) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_558():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = (a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_559():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x
    F = a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 3*a**2*b*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + 3*a*b**2*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2) + b**3*x**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*a + 6*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_560():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**3
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 3*a*b**2*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + b**3*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_561():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**5
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**4*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + b**3*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_562():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**7
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**6*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**4*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_563():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**9
    F = -(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_564():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**11
    F = -(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(8*a*x**10) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(40*a**2*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_565():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**13
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*x**12*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**6*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_566():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**15
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*x**14*(a + b*x**2)) - a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**12*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_567():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**17
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*x**16*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*x**14*(a + b*x**2)) - a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**12*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_568():
    f = x**8*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + 3*a**2*b*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + 3*a*b**2*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + b**3*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*a + 15*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_569():
    f = x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + a**2*b*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 3*a*b**2*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + b**3*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_570():
    f = x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2) + 3*a**2*b*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + a*b**2*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + b**3*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_571():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 3*a**2*b*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2) + 3*a*b**2*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + b**3*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_572():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(a + b*x**2)**3 + a**2*b*x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(a + b*x**2)**3 + 3*a*b**2*x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(5*(a + b*x**2)**3) + b**3*x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(7*(a + b*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_573():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**2
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 3*a**2*b*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + a*b**2*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b**3*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_574():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**4
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 3*a*b**2*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b**3*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_575():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**6
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2)) - a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**3*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + b**3*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_576():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**8
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2)) - a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**3*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_577():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**10
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_578():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**12
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**9*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_579():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**14
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**9*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_580():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/x**16
    F = -a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*x**15*(a + b*x**2)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_581():
    f = x**13*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**14*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*a + 14*b*x**2) + 5*a**4*b*x**16*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*a + 16*b*x**2) + 5*a**3*b**2*x**18*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + a**2*b**3*x**20*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + 5*a*b**4*x**22*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(22*a + 22*b*x**2) + b**5*x**24*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(24*a + 24*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_582():
    f = x**11*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**12*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*a + 12*b*x**2) + 5*a**4*b*x**14*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*a + 14*b*x**2) + 5*a**3*b**2*x**16*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a + 8*b*x**2) + 5*a**2*b**3*x**18*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + a*b**4*x**20*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2) + b**5*x**22*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(22*a + 22*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_583():
    f = x**9*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**4*(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*b**5) - 2*a**3*(a + b*x**2)**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*b**5) + 3*a**2*(a + b*x**2)**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*b**5) - 2*a*(a + b*x**2)**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*b**5) + (a + b*x**2)**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(20*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_584():
    f = x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -a**3*(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*b**4) + 3*a**2*(a + b*x**2)**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*b**4) - 3*a*(a + b*x**2)**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*b**4) + (a + b*x**2)**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(18*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_585():
    f = x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**2*(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*b**3) - a*(a + b*x**2)**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*b**3) + (a + b*x**2)**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_586():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -a*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(12*b**2) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(7)/2)/(14*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_587():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = (a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_588():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x
    F = a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 5*a**4*b*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + 5*a**3*b**2*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + 5*a**2*b**3*x**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 5*a*b**4*x**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a + 8*b*x**2) + b**5*x**10*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*a + 10*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_589():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**3
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 5*a**3*b**2*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 5*a**2*b**3*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + 5*a*b**4*x**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*a + 6*b*x**2) + b**5*x**8*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*a + 8*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_590():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**5
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**4*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 5*a**2*b**3*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 5*a*b**4*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2) + b**5*x**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*a + 6*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_591():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**7
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**6*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**4*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**2*(a + b*x**2)) + 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + 5*a*b**4*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2) + b**5*x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*a + 4*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_592():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**9
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**6*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**4*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**2*(a + b*x**2)) + 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2) + b**5*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*a + 2*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_593():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**11
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**6*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**4*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**2*(a + b*x**2)) + b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*log(x)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_594():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**13
    F = -(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*a*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_595():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**15
    F = -(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(12*a*x**14) + (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(7)/2)/(84*a**2*x**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_596():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**17
    F = -(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*a*x**16) + b*(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(56*a**2*x**14) - b**2*(a + b*x**2)**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(336*a**3*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_597():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**19
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(18*x**18*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*x**16*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**14*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*x**12*(a + b*x**2)) - a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**10*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**8*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_598():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**21
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(20*x**20*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(18*x**18*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**16*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**14*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*x**12*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(10*x**10*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_599():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**23
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(22*x**22*(a + b*x**2)) - a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*x**20*(a + b*x**2)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**18*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*x**16*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*x**14*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(12*x**12*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_600():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**25
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(24*x**24*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(22*x**22*(a + b*x**2)) - a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*x**20*(a + b*x**2)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**18*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*x**16*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(14*x**14*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_601():
    f = x**12*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + a**4*b*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 10*a**3*b**2*x**17*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*a + 17*b*x**2) + 10*a**2*b**3*x**19*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*a + 19*b*x**2) + 5*a*b**4*x**21*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*a + 21*b*x**2) + b**5*x**23*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(23*a + 23*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_602():
    f = x**10*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + 5*a**4*b*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + 2*a**3*b**2*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 10*a**2*b**3*x**17*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*a + 17*b*x**2) + 5*a*b**4*x**19*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*a + 19*b*x**2) + b**5*x**21*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*a + 21*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_603():
    f = x**8*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + 5*a**4*b*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + 10*a**3*b**2*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + 2*a**2*b**3*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 5*a*b**4*x**17*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*a + 17*b*x**2) + b**5*x**19*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*a + 19*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_604():
    f = x**6*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + 5*a**4*b*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + 10*a**3*b**2*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + 10*a**2*b**3*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + a*b**4*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + b**5*x**17*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*a + 17*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_605():
    f = x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2) + 5*a**4*b*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + 10*a**3*b**2*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + 10*a**2*b**3*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + 5*a*b**4*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2) + b**5*x**15*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*a + 15*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_606():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + a**4*b*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 10*a**3*b**2*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + 10*a**2*b**3*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2) + 5*a*b**4*x**11*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*a + 11*b*x**2) + b**5*x**13*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*a + 13*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_607():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(a + b*x**2)**5 + 5*a**4*b*x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(3*(a + b*x**2)**5) + 2*a**3*b**2*x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(a + b*x**2)**5 + 10*a**2*b**3*x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(7*(a + b*x**2)**5) + 5*a*b**4*x**9*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(9*(a + b*x**2)**5) + b**5*x**11*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(11*(a + b*x**2)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_608():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**2
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 5*a**4*b*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 10*a**3*b**2*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + 2*a**2*b**3*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 5*a*b**4*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2) + b**5*x**9*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*a + 9*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_609():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**4
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 10*a**3*b**2*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 10*a**2*b**3*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + a*b**4*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b**5*x**7*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*a + 7*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_610():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**6
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 10*a**2*b**3*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + 5*a*b**4*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2) + b**5*x**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*a + 5*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_611():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**8
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**5*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + 5*a*b**4*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b**5*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*a + 3*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_612():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**10
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - 2*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**5*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2)) + b**5*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_613():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**12
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - 2*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**5*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_614():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**14
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(x**5*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_615():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**16
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*x**15*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*x**5*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_616():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**18
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*x**17*(a + b*x**2)) - a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**15*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*x**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_617():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**20
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*x**19*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*x**17*(a + b*x**2)) - 2*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**15*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*x**9*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_618():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**22
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*x**21*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*x**19*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*x**17*(a + b*x**2)) - 2*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**15*(a + b*x**2)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*x**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_619():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/x**24
    F = -a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(23*x**23*(a + b*x**2)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*x**21*(a + b*x**2)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*x**19*(a + b*x**2)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*x**17*(a + b*x**2)) - a*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*x**15*(a + b*x**2)) - b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*x**13*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_620():
    f = x**5/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**2*(a + b*x**2)*log(a + b*x**2)/(2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - a*x**2*(a + b*x**2)/(2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x**4*(a + b*x**2)/(4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_621():
    f = x**3/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -a*(a + b*x**2)*log(a + b*x**2)/(2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_622():
    f = x/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (a + b*x**2)*log(a + b*x**2)/(2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_623():
    f = 1/(x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (a + b*x**2)*log(x)/(a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*log(a + b*x**2)/(2*a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_624():
    f = 1/(x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (-a - b*x**2)/(2*a*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - b*(a + b*x**2)*log(x)/(a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + b*(a + b*x**2)*log(a + b*x**2)/(2*a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_625():
    f = x**4/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**(sympy.S(3)/2)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - a*x*(a + b*x**2)/(b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x**3*(a + b*x**2)/(3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_626():
    f = x**2/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -sqrt(a)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x*(a + b*x**2)/(b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_627():
    f = 1/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_628():
    f = 1/(x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -(a + b*x**2)/(a*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(b)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_629():
    f = 1/(x**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (-a - b*x**2)/(3*a*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + b*(a + b*x**2)/(a**2*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + b**(sympy.S(3)/2)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_630():
    f = x**7/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3/(4*b**4*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*a**2/(2*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*a*(a + b*x**2)*log(a + b*x**2)/(2*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x**2*(a + b*x**2)/(2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_631():
    f = x**5/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -a**2/(4*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + a/(b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*log(a + b*x**2)/(2*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_632():
    f = x/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -1/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_633():
    f = 1/(x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1/(2*a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*log(x)/(a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*log(a + b*x**2)/(2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_634():
    f = 1/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = -b/(4*a**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - b/(a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)/(2*a**3*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*b*(a + b*x**2)*log(x)/(a**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*b*(a + b*x**2)*log(a + b*x**2)/(2*a**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_635():
    f = x**4/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -x**3/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*x/(8*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (3*a + 3*b*x**2)*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_636():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -x/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x/(8*a*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_637():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-3)/2)
    F = x*(a + b*x**2)/(4*a*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)) + 3*x*(a + b*x**2)**2/(8*a**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)) + 3*(a + b*x**2)**3*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_638():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*x*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5/(8*a**2*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (15*a + 15*b*x**2)/(8*a**3*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 15*sqrt(b)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_639():
    f = 1/(x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*x**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7/(8*a**2*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (35*a + 35*b*x**2)/(24*a**3*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*b*(a + b*x**2)/(8*a**4*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*b**(sympy.S(3)/2)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_640():
    f = x**11/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5/(8*b**6*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*a**4/(6*b**6*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*a**3/(2*b**6*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*a**2/(b**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*a*(a + b*x**2)*log(a + b*x**2)/(2*b**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x**2*(a + b*x**2)/(2*b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_641():
    f = x**9/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -a**4/(8*b**5*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*a**3/(3*b**5*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*a**2/(2*b**5*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*a/(b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*log(a + b*x**2)/(2*b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_642():
    f = x**7/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = x**8/(8*a*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_643():
    f = x**5/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = x**6/(8*a*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)) + x**6/(24*a**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_644():
    f = x**3/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a/(8*b**2*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)) - 1/(6*b**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_645():
    f = x/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -1/(8*b*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_646():
    f = 1/(x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1/(6*a**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1/(4*a**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1/(2*a**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*log(x)/(a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*log(a + b*x**2)/(2*a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_647():
    f = 1/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = -b/(8*a**2*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - b/(3*a**3*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*b/(4*a**4*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 2*b/(a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)/(2*a**5*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*b*(a + b*x**2)*log(x)/(a**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*b*(a + b*x**2)*log(a + b*x**2)/(2*a**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_648():
    f = x**6/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -x**5/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*x**3/(48*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*x/(64*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*x/(128*a*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (5*a + 5*b*x**2)*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(3)/2)*b**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_649():
    f = x**4/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -x**3/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - x/(16*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x/(64*a*b**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*x/(128*a**2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (3*a + 3*b*x**2)*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(5)/2)*b**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_650():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -x/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x/(48*a*b*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*x/(192*a**2*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*x/(128*a**3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (5*a + 5*b*x**2)*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_651():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-5)/2)
    F = x*(a + b*x**2)/(8*a*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)) + 7*x*(a + b*x**2)**2/(48*a**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)) + 35*x*(a + b*x**2)**3/(192*a**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)) + 35*x*(a + b*x**2)**4/(128*a**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)) + 35*(a + b*x**2)**5*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(9)/2)*sqrt(b)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_652():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*x*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3/(16*a**2*x*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 21/(64*a**3*x*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 105/(128*a**4*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (315*a + 315*b*x**2)/(128*a**5*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 315*sqrt(b)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_653():
    f = 1/(x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*x**3*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 11/(48*a**2*x**3*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 33/(64*a**3*x**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 231/(128*a**4*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (385*a + 385*b*x**2)/(128*a**5*x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1155*b*(a + b*x**2)/(128*a**6*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1155*b**(sympy.S(3)/2)*(a + b*x**2)*atan(sqrt(b)*x/sqrt(a))/(128*a**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_654():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3)
    F = 3*3**(sympy.S(3)/4)*a**2*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(2)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(5*b**2*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3)) + 3*x*(a + b*x**2)/(5*b*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_655():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-1)/3)
    F = -3**(sympy.S(3)/4)*a*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(2)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(b*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_656():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3))
    F = 3**(sympy.S(3)/4)*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(2)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3)) - (a + b*x**2)/(a*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_657():
    f = x**2/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)
    F = 9*3**(sympy.S(1)/4)*a**2*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b**2*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) - 3*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*b**2*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) - 9*a*x*(1 + b*x**2/a)**(sympy.S(4)/3)/(2*b*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)*(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)) - 3*x*(a + b*x**2)/(2*b*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_658():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(-2)/3)
    F = -3*3**(sympy.S(1)/4)*a*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*b*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*b*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) + 3*x*(1 + b*x**2/a)**(sympy.S(4)/3)/(2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)*(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)) + 3*x*(a + b*x**2)/(2*a*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_659():
    f = 1/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3))
    F = 5*3**(sympy.S(1)/4)*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) - 5*sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 + b*x**2/a)**(sympy.S(2)/3) + (1 + b*x**2/a)**(sympy.S(1)/3) + 1)/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 + b*x**2/a)**(sympy.S(4)/3)*(1 - (1 + b*x**2/a)**(sympy.S(1)/3))*elliptic_f(asin((-(1 + b*x**2/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(6*x*sqrt(-(1 - (1 + b*x**2/a)**(sympy.S(1)/3))/(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) - 5*b*x*(1 + b*x**2/a)**(sympy.S(4)/3)/(2*a*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)*(-(1 + b*x**2/a)**(sympy.S(1)/3) - sqrt(3) + 1)) + (3*a + 3*b*x**2)/(2*a*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3)) - 5*(a + b*x**2)**2/(2*a**2*x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_660():
    f = (d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a**2*(d*x)**(sympy.S(7)/2)/(7*d) + 4*a*b*(d*x)**(sympy.S(11)/2)/(11*d**3) + 2*b**2*(d*x)**(sympy.S(15)/2)/(15*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_661():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a**2*(d*x)**(sympy.S(5)/2)/(5*d) + 4*a*b*(d*x)**(sympy.S(9)/2)/(9*d**3) + 2*b**2*(d*x)**(sympy.S(13)/2)/(13*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_662():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a**2*(d*x)**(sympy.S(3)/2)/(3*d) + 4*a*b*(d*x)**(sympy.S(7)/2)/(7*d**3) + 2*b**2*(d*x)**(sympy.S(11)/2)/(11*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_663():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/sqrt(d*x)
    F = 2*a**2*sqrt(d*x)/d + 4*a*b*(d*x)**(sympy.S(5)/2)/(5*d**3) + 2*b**2*(d*x)**(sympy.S(9)/2)/(9*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_664():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(3)/2)
    F = -2*a**2/(d*sqrt(d*x)) + 4*a*b*(d*x)**(sympy.S(3)/2)/(3*d**3) + 2*b**2*(d*x)**(sympy.S(7)/2)/(7*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_665():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(5)/2)
    F = -2*a**2/(3*d*(d*x)**(sympy.S(3)/2)) + 4*a*b*sqrt(d*x)/d**3 + 2*b**2*(d*x)**(sympy.S(5)/2)/(5*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_666():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(7)/2)
    F = -2*a**2/(5*d*(d*x)**(sympy.S(5)/2)) - 4*a*b/(d**3*sqrt(d*x)) + 2*b**2*(d*x)**(sympy.S(3)/2)/(3*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_667():
    f = (d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = 2*a**4*(d*x)**(sympy.S(7)/2)/(7*d) + 8*a**3*b*(d*x)**(sympy.S(11)/2)/(11*d**3) + 4*a**2*b**2*(d*x)**(sympy.S(15)/2)/(5*d**5) + 8*a*b**3*(d*x)**(sympy.S(19)/2)/(19*d**7) + 2*b**4*(d*x)**(sympy.S(23)/2)/(23*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_668():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = 2*a**4*(d*x)**(sympy.S(5)/2)/(5*d) + 8*a**3*b*(d*x)**(sympy.S(9)/2)/(9*d**3) + 12*a**2*b**2*(d*x)**(sympy.S(13)/2)/(13*d**5) + 8*a*b**3*(d*x)**(sympy.S(17)/2)/(17*d**7) + 2*b**4*(d*x)**(sympy.S(21)/2)/(21*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_669():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = 2*a**4*(d*x)**(sympy.S(3)/2)/(3*d) + 8*a**3*b*(d*x)**(sympy.S(7)/2)/(7*d**3) + 12*a**2*b**2*(d*x)**(sympy.S(11)/2)/(11*d**5) + 8*a*b**3*(d*x)**(sympy.S(15)/2)/(15*d**7) + 2*b**4*(d*x)**(sympy.S(19)/2)/(19*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_670():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/sqrt(d*x)
    F = 2*a**4*sqrt(d*x)/d + 8*a**3*b*(d*x)**(sympy.S(5)/2)/(5*d**3) + 4*a**2*b**2*(d*x)**(sympy.S(9)/2)/(3*d**5) + 8*a*b**3*(d*x)**(sympy.S(13)/2)/(13*d**7) + 2*b**4*(d*x)**(sympy.S(17)/2)/(17*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_671():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/(d*x)**(sympy.S(3)/2)
    F = -2*a**4/(d*sqrt(d*x)) + 8*a**3*b*(d*x)**(sympy.S(3)/2)/(3*d**3) + 12*a**2*b**2*(d*x)**(sympy.S(7)/2)/(7*d**5) + 8*a*b**3*(d*x)**(sympy.S(11)/2)/(11*d**7) + 2*b**4*(d*x)**(sympy.S(15)/2)/(15*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_672():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/(d*x)**(sympy.S(5)/2)
    F = -2*a**4/(3*d*(d*x)**(sympy.S(3)/2)) + 8*a**3*b*sqrt(d*x)/d**3 + 12*a**2*b**2*(d*x)**(sympy.S(5)/2)/(5*d**5) + 8*a*b**3*(d*x)**(sympy.S(9)/2)/(9*d**7) + 2*b**4*(d*x)**(sympy.S(13)/2)/(13*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_673():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**2/(d*x)**(sympy.S(7)/2)
    F = -2*a**4/(5*d*(d*x)**(sympy.S(5)/2)) - 8*a**3*b/(d**3*sqrt(d*x)) + 4*a**2*b**2*(d*x)**(sympy.S(3)/2)/d**5 + 8*a*b**3*(d*x)**(sympy.S(7)/2)/(7*d**7) + 2*b**4*(d*x)**(sympy.S(11)/2)/(11*d**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_674():
    f = (d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = 2*a**6*(d*x)**(sympy.S(7)/2)/(7*d) + 12*a**5*b*(d*x)**(sympy.S(11)/2)/(11*d**3) + 2*a**4*b**2*(d*x)**(sympy.S(15)/2)/d**5 + 40*a**3*b**3*(d*x)**(sympy.S(19)/2)/(19*d**7) + 30*a**2*b**4*(d*x)**(sympy.S(23)/2)/(23*d**9) + 4*a*b**5*(d*x)**(sympy.S(27)/2)/(9*d**11) + 2*b**6*(d*x)**(sympy.S(31)/2)/(31*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_675():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = 2*a**6*(d*x)**(sympy.S(5)/2)/(5*d) + 4*a**5*b*(d*x)**(sympy.S(9)/2)/(3*d**3) + 30*a**4*b**2*(d*x)**(sympy.S(13)/2)/(13*d**5) + 40*a**3*b**3*(d*x)**(sympy.S(17)/2)/(17*d**7) + 10*a**2*b**4*(d*x)**(sympy.S(21)/2)/(7*d**9) + 12*a*b**5*(d*x)**(sympy.S(25)/2)/(25*d**11) + 2*b**6*(d*x)**(sympy.S(29)/2)/(29*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_676():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = 2*a**6*(d*x)**(sympy.S(3)/2)/(3*d) + 12*a**5*b*(d*x)**(sympy.S(7)/2)/(7*d**3) + 30*a**4*b**2*(d*x)**(sympy.S(11)/2)/(11*d**5) + 8*a**3*b**3*(d*x)**(sympy.S(15)/2)/(3*d**7) + 30*a**2*b**4*(d*x)**(sympy.S(19)/2)/(19*d**9) + 12*a*b**5*(d*x)**(sympy.S(23)/2)/(23*d**11) + 2*b**6*(d*x)**(sympy.S(27)/2)/(27*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_677():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/sqrt(d*x)
    F = 2*a**6*sqrt(d*x)/d + 12*a**5*b*(d*x)**(sympy.S(5)/2)/(5*d**3) + 10*a**4*b**2*(d*x)**(sympy.S(9)/2)/(3*d**5) + 40*a**3*b**3*(d*x)**(sympy.S(13)/2)/(13*d**7) + 30*a**2*b**4*(d*x)**(sympy.S(17)/2)/(17*d**9) + 4*a*b**5*(d*x)**(sympy.S(21)/2)/(7*d**11) + 2*b**6*(d*x)**(sympy.S(25)/2)/(25*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_678():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/(d*x)**(sympy.S(3)/2)
    F = -2*a**6/(d*sqrt(d*x)) + 4*a**5*b*(d*x)**(sympy.S(3)/2)/d**3 + 30*a**4*b**2*(d*x)**(sympy.S(7)/2)/(7*d**5) + 40*a**3*b**3*(d*x)**(sympy.S(11)/2)/(11*d**7) + 2*a**2*b**4*(d*x)**(sympy.S(15)/2)/d**9 + 12*a*b**5*(d*x)**(sympy.S(19)/2)/(19*d**11) + 2*b**6*(d*x)**(sympy.S(23)/2)/(23*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_679():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/(d*x)**(sympy.S(5)/2)
    F = -2*a**6/(3*d*(d*x)**(sympy.S(3)/2)) + 12*a**5*b*sqrt(d*x)/d**3 + 6*a**4*b**2*(d*x)**(sympy.S(5)/2)/d**5 + 40*a**3*b**3*(d*x)**(sympy.S(9)/2)/(9*d**7) + 30*a**2*b**4*(d*x)**(sympy.S(13)/2)/(13*d**9) + 12*a*b**5*(d*x)**(sympy.S(17)/2)/(17*d**11) + 2*b**6*(d*x)**(sympy.S(21)/2)/(21*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_680():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**3/(d*x)**(sympy.S(7)/2)
    F = -2*a**6/(5*d*(d*x)**(sympy.S(5)/2)) - 12*a**5*b/(d**3*sqrt(d*x)) + 10*a**4*b**2*(d*x)**(sympy.S(3)/2)/d**5 + 40*a**3*b**3*(d*x)**(sympy.S(7)/2)/(7*d**7) + 30*a**2*b**4*(d*x)**(sympy.S(11)/2)/(11*d**9) + 4*a*b**5*(d*x)**(sympy.S(15)/2)/(5*d**11) + 2*b**6*(d*x)**(sympy.S(19)/2)/(19*d**13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_681():
    f = (d*x)**(sympy.S(11)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -9*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(11)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(13)/4)) + 9*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(11)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(13)/4)) - 9*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(11)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(13)/4)) + 9*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(11)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(13)/4)) - 9*a*d**5*sqrt(d*x)/(2*b**3) - d*(d*x)**(sympy.S(9)/2)/(2*b*(a + b*x**2)) + 9*d**3*(d*x)**(sympy.S(5)/2)/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_682():
    f = (d*x)**(sympy.S(9)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -7*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(9)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(11)/4)) + 7*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(9)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(11)/4)) + 7*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(9)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(11)/4)) - 7*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(9)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(11)/4)) - d*(d*x)**(sympy.S(7)/2)/(2*b*(a + b*x**2)) + 7*d**3*(d*x)**(sympy.S(3)/2)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_683():
    f = (d*x)**(sympy.S(7)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 5*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(7)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(9)/4)) - 5*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(7)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*b**(sympy.S(9)/4)) + 5*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(7)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(9)/4)) - 5*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(7)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*b**(sympy.S(9)/4)) - d*(d*x)**(sympy.S(5)/2)/(2*b*(a + b*x**2)) + 5*d**3*sqrt(d*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_684():
    f = (d*x)**(sympy.S(5)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -d*(d*x)**(sympy.S(3)/2)/(2*b*(a + b*x**2)) + 3*sqrt(2)*d**(sympy.S(5)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) + 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(1)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_685():
    f = (d*x)**(sympy.S(3)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -d*sqrt(d*x)/(2*b*(a + b*x**2)) - sqrt(2)*d**(sympy.S(3)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_686():
    f = sqrt(d*x)/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (d*x)**(sympy.S(3)/2)/(2*a*d*(a + b*x**2)) + sqrt(2)*sqrt(d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) - sqrt(2)*sqrt(d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(5)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_687():
    f = 1/(sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = sqrt(d*x)/(2*a*d*(a + b*x**2)) - 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(d)) - 3*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 3*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_688():
    f = 1/((d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*d*sqrt(d*x)*(a + b*x**2)) - 5/(2*a**2*d*sqrt(d*x)) - 5*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(9)/4)*d**(sympy.S(3)/2)) + 5*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(9)/4)*d**(sympy.S(3)/2)) + 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(9)/4)*d**(sympy.S(3)/2)) - 5*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(9)/4)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_689():
    f = 1/((d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) - 7/(6*a**2*d*(d*x)**(sympy.S(3)/2)) + 7*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(11)/4)*d**(sympy.S(5)/2)) - 7*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(11)/4)*d**(sympy.S(5)/2)) + 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(11)/4)*d**(sympy.S(5)/2)) - 7*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(11)/4)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_690():
    f = 1/((d*x)**(sympy.S(7)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = 1/(2*a*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 9/(10*a**2*d*(d*x)**(sympy.S(5)/2)) + 9*b/(2*a**3*d**3*sqrt(d*x)) + 9*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(13)/4)*d**(sympy.S(7)/2)) - 9*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(16*a**(sympy.S(13)/4)*d**(sympy.S(7)/2)) - 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(13)/4)*d**(sympy.S(7)/2)) + 9*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(8*a**(sympy.S(13)/4)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_691():
    f = (d*x)**(sympy.S(19)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -663*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(19)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(21)/4)) + 663*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(19)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(21)/4)) - 663*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(19)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(21)/4)) + 663*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(19)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(21)/4)) - 663*a*d**9*sqrt(d*x)/(64*b**5) - d*(d*x)**(sympy.S(17)/2)/(6*b*(a + b*x**2)**3) - 17*d**3*(d*x)**(sympy.S(13)/2)/(48*b**2*(a + b*x**2)**2) - 221*d**5*(d*x)**(sympy.S(9)/2)/(192*b**3*(a + b*x**2)) + 663*d**7*(d*x)**(sympy.S(5)/2)/(320*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_692():
    f = (d*x)**(sympy.S(17)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -385*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(17)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(19)/4)) + 385*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(17)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(19)/4)) + 385*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(17)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(19)/4)) - 385*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(17)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(19)/4)) - d*(d*x)**(sympy.S(15)/2)/(6*b*(a + b*x**2)**3) - 5*d**3*(d*x)**(sympy.S(11)/2)/(16*b**2*(a + b*x**2)**2) - 55*d**5*(d*x)**(sympy.S(7)/2)/(64*b**3*(a + b*x**2)) + 385*d**7*(d*x)**(sympy.S(3)/2)/(192*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_693():
    f = (d*x)**(sympy.S(15)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = 195*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(15)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(17)/4)) - 195*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(15)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*b**(sympy.S(17)/4)) + 195*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(15)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(17)/4)) - 195*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(15)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*b**(sympy.S(17)/4)) - d*(d*x)**(sympy.S(13)/2)/(6*b*(a + b*x**2)**3) - 13*d**3*(d*x)**(sympy.S(9)/2)/(48*b**2*(a + b*x**2)**2) - 39*d**5*(d*x)**(sympy.S(5)/2)/(64*b**3*(a + b*x**2)) + 195*d**7*sqrt(d*x)/(64*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_694():
    f = (d*x)**(sympy.S(13)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*(d*x)**(sympy.S(11)/2)/(6*b*(a + b*x**2)**3) - 11*d**3*(d*x)**(sympy.S(7)/2)/(48*b**2*(a + b*x**2)**2) - 77*d**5*(d*x)**(sympy.S(3)/2)/(192*b**3*(a + b*x**2)) + 77*sqrt(2)*d**(sympy.S(13)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) + 77*sqrt(2)*d**(sympy.S(13)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(1)/4)*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_695():
    f = (d*x)**(sympy.S(11)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*(d*x)**(sympy.S(9)/2)/(6*b*(a + b*x**2)**3) - 3*d**3*(d*x)**(sympy.S(5)/2)/(16*b**2*(a + b*x**2)**2) - 15*d**5*sqrt(d*x)/(64*b**3*(a + b*x**2)) - 15*sqrt(2)*d**(sympy.S(11)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + 15*sqrt(2)*d**(sympy.S(11)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) - 15*sqrt(2)*d**(sympy.S(11)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + 15*sqrt(2)*d**(sympy.S(11)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(3)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_696():
    f = (d*x)**(sympy.S(9)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*(d*x)**(sympy.S(7)/2)/(6*b*(a + b*x**2)**3) - 7*d**3*(d*x)**(sympy.S(3)/2)/(48*b**2*(a + b*x**2)**2) + 7*d**3*(d*x)**(sympy.S(3)/2)/(64*a*b**2*(a + b*x**2)) + 7*sqrt(2)*d**(sympy.S(9)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) - 7*sqrt(2)*d**(sympy.S(9)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) - 7*sqrt(2)*d**(sympy.S(9)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) + 7*sqrt(2)*d**(sympy.S(9)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(5)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_697():
    f = (d*x)**(sympy.S(7)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*(d*x)**(sympy.S(5)/2)/(6*b*(a + b*x**2)**3) - 5*d**3*sqrt(d*x)/(48*b**2*(a + b*x**2)**2) + 5*d**3*sqrt(d*x)/(192*a*b**2*(a + b*x**2)) - 5*sqrt(2)*d**(sympy.S(7)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + 5*sqrt(2)*d**(sympy.S(7)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) - 5*sqrt(2)*d**(sympy.S(7)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + 5*sqrt(2)*d**(sympy.S(7)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(7)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_698():
    f = (d*x)**(sympy.S(5)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*(d*x)**(sympy.S(3)/2)/(6*b*(a + b*x**2)**3) + d*(d*x)**(sympy.S(3)/2)/(16*a*b*(a + b*x**2)**2) + 5*d*(d*x)**(sympy.S(3)/2)/(64*a**2*b*(a + b*x**2)) + 5*sqrt(2)*d**(sympy.S(5)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - 5*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - 5*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) + 5*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(9)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_699():
    f = (d*x)**(sympy.S(3)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = -d*sqrt(d*x)/(6*b*(a + b*x**2)**3) + d*sqrt(d*x)/(48*a*b*(a + b*x**2)**2) + 7*d*sqrt(d*x)/(192*a**2*b*(a + b*x**2)) - 7*sqrt(2)*d**(sympy.S(3)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + 7*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) - 7*sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + 7*sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_700():
    f = sqrt(d*x)/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = (d*x)**(sympy.S(3)/2)/(6*a*d*(a + b*x**2)**3) + 3*(d*x)**(sympy.S(3)/2)/(16*a**2*d*(a + b*x**2)**2) + 15*(d*x)**(sympy.S(3)/2)/(64*a**3*d*(a + b*x**2)) + 15*sqrt(2)*sqrt(d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) - 15*sqrt(2)*sqrt(d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) - 15*sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) + 15*sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(13)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_701():
    f = 1/(sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = sqrt(d*x)/(6*a*d*(a + b*x**2)**3) + 11*sqrt(d*x)/(48*a**2*d*(a + b*x**2)**2) + 77*sqrt(d*x)/(192*a**3*d*(a + b*x**2)) - 77*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 77*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)*sqrt(d)) - 77*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 77*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_702():
    f = 1/((d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*d*sqrt(d*x)*(a + b*x**2)**3) + 13/(48*a**2*d*sqrt(d*x)*(a + b*x**2)**2) + 39/(64*a**3*d*sqrt(d*x)*(a + b*x**2)) - 195/(64*a**4*d*sqrt(d*x)) - 195*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(17)/4)*d**(sympy.S(3)/2)) + 195*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(17)/4)*d**(sympy.S(3)/2)) + 195*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(17)/4)*d**(sympy.S(3)/2)) - 195*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(17)/4)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_703():
    f = 1/((d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**3) + 5/(16*a**2*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**2) + 55/(64*a**3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) - 385/(192*a**4*d*(d*x)**(sympy.S(3)/2)) + 385*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(19)/4)*d**(sympy.S(5)/2)) - 385*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(19)/4)*d**(sympy.S(5)/2)) + 385*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(19)/4)*d**(sympy.S(5)/2)) - 385*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(19)/4)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_704():
    f = 1/((d*x)**(sympy.S(7)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**2)
    F = 1/(6*a*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**3) + 17/(48*a**2*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**2) + 221/(192*a**3*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 663/(320*a**4*d*(d*x)**(sympy.S(5)/2)) + 663*b/(64*a**5*d**3*sqrt(d*x)) + 663*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(21)/4)*d**(sympy.S(7)/2)) - 663*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(512*a**(sympy.S(21)/4)*d**(sympy.S(7)/2)) - 663*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(21)/4)*d**(sympy.S(7)/2)) + 663*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(256*a**(sympy.S(21)/4)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_705():
    f = (d*x)**(sympy.S(27)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -69615*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(27)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(29)/4)) + 69615*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(27)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(29)/4)) - 69615*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(27)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(29)/4)) + 69615*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(27)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(29)/4)) - 69615*a*d**13*sqrt(d*x)/(4096*b**7) - d*(d*x)**(sympy.S(25)/2)/(10*b*(a + b*x**2)**5) - 5*d**3*(d*x)**(sympy.S(21)/2)/(32*b**2*(a + b*x**2)**4) - 35*d**5*(d*x)**(sympy.S(17)/2)/(128*b**3*(a + b*x**2)**3) - 595*d**7*(d*x)**(sympy.S(13)/2)/(1024*b**4*(a + b*x**2)**2) - 7735*d**9*(d*x)**(sympy.S(9)/2)/(4096*b**5*(a + b*x**2)) + 13923*d**11*(d*x)**(sympy.S(5)/2)/(4096*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_706():
    f = (d*x)**(sympy.S(25)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -33649*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(25)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(27)/4)) + 33649*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(25)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(27)/4)) + 33649*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(25)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(27)/4)) - 33649*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(25)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(27)/4)) - d*(d*x)**(sympy.S(23)/2)/(10*b*(a + b*x**2)**5) - 23*d**3*(d*x)**(sympy.S(19)/2)/(160*b**2*(a + b*x**2)**4) - 437*d**5*(d*x)**(sympy.S(15)/2)/(1920*b**3*(a + b*x**2)**3) - 437*d**7*(d*x)**(sympy.S(11)/2)/(1024*b**4*(a + b*x**2)**2) - 4807*d**9*(d*x)**(sympy.S(7)/2)/(4096*b**5*(a + b*x**2)) + 33649*d**11*(d*x)**(sympy.S(3)/2)/(12288*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_707():
    f = (d*x)**(sympy.S(23)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = 13923*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(23)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(25)/4)) - 13923*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(23)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*b**(sympy.S(25)/4)) + 13923*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(23)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(25)/4)) - 13923*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(23)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*b**(sympy.S(25)/4)) - d*(d*x)**(sympy.S(21)/2)/(10*b*(a + b*x**2)**5) - 21*d**3*(d*x)**(sympy.S(17)/2)/(160*b**2*(a + b*x**2)**4) - 119*d**5*(d*x)**(sympy.S(13)/2)/(640*b**3*(a + b*x**2)**3) - 1547*d**7*(d*x)**(sympy.S(9)/2)/(5120*b**4*(a + b*x**2)**2) - 13923*d**9*(d*x)**(sympy.S(5)/2)/(20480*b**5*(a + b*x**2)) + 13923*d**11*sqrt(d*x)/(4096*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_708():
    f = (d*x)**(sympy.S(21)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(19)/2)/(10*b*(a + b*x**2)**5) - 19*d**3*(d*x)**(sympy.S(15)/2)/(160*b**2*(a + b*x**2)**4) - 19*d**5*(d*x)**(sympy.S(11)/2)/(128*b**3*(a + b*x**2)**3) - 209*d**7*(d*x)**(sympy.S(7)/2)/(1024*b**4*(a + b*x**2)**2) - 1463*d**9*(d*x)**(sympy.S(3)/2)/(4096*b**5*(a + b*x**2)) + 4389*sqrt(2)*d**(sympy.S(21)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(1)/4)*b**(sympy.S(23)/4)) - 4389*sqrt(2)*d**(sympy.S(21)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(1)/4)*b**(sympy.S(23)/4)) - 4389*sqrt(2)*d**(sympy.S(21)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(1)/4)*b**(sympy.S(23)/4)) + 4389*sqrt(2)*d**(sympy.S(21)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(1)/4)*b**(sympy.S(23)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_709():
    f = (d*x)**(sympy.S(19)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(17)/2)/(10*b*(a + b*x**2)**5) - 17*d**3*(d*x)**(sympy.S(13)/2)/(160*b**2*(a + b*x**2)**4) - 221*d**5*(d*x)**(sympy.S(9)/2)/(1920*b**3*(a + b*x**2)**3) - 663*d**7*(d*x)**(sympy.S(5)/2)/(5120*b**4*(a + b*x**2)**2) - 663*d**9*sqrt(d*x)/(4096*b**5*(a + b*x**2)) - 663*sqrt(2)*d**(sympy.S(19)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(3)/4)*b**(sympy.S(21)/4)) + 663*sqrt(2)*d**(sympy.S(19)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(3)/4)*b**(sympy.S(21)/4)) - 663*sqrt(2)*d**(sympy.S(19)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(3)/4)*b**(sympy.S(21)/4)) + 663*sqrt(2)*d**(sympy.S(19)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(3)/4)*b**(sympy.S(21)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_710():
    f = (d*x)**(sympy.S(17)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(15)/2)/(10*b*(a + b*x**2)**5) - 3*d**3*(d*x)**(sympy.S(11)/2)/(32*b**2*(a + b*x**2)**4) - 11*d**5*(d*x)**(sympy.S(7)/2)/(128*b**3*(a + b*x**2)**3) - 77*d**7*(d*x)**(sympy.S(3)/2)/(1024*b**4*(a + b*x**2)**2) + 231*d**7*(d*x)**(sympy.S(3)/2)/(4096*a*b**4*(a + b*x**2)) + 231*sqrt(2)*d**(sympy.S(17)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(5)/4)*b**(sympy.S(19)/4)) - 231*sqrt(2)*d**(sympy.S(17)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(5)/4)*b**(sympy.S(19)/4)) - 231*sqrt(2)*d**(sympy.S(17)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(5)/4)*b**(sympy.S(19)/4)) + 231*sqrt(2)*d**(sympy.S(17)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(5)/4)*b**(sympy.S(19)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_711():
    f = (d*x)**(sympy.S(15)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(13)/2)/(10*b*(a + b*x**2)**5) - 13*d**3*(d*x)**(sympy.S(9)/2)/(160*b**2*(a + b*x**2)**4) - 39*d**5*(d*x)**(sympy.S(5)/2)/(640*b**3*(a + b*x**2)**3) - 39*d**7*sqrt(d*x)/(1024*b**4*(a + b*x**2)**2) + 39*d**7*sqrt(d*x)/(4096*a*b**4*(a + b*x**2)) - 117*sqrt(2)*d**(sympy.S(15)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) + 117*sqrt(2)*d**(sympy.S(15)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) - 117*sqrt(2)*d**(sympy.S(15)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) + 117*sqrt(2)*d**(sympy.S(15)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(7)/4)*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_712():
    f = (d*x)**(sympy.S(13)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(11)/2)/(10*b*(a + b*x**2)**5) - 11*d**3*(d*x)**(sympy.S(7)/2)/(160*b**2*(a + b*x**2)**4) - 77*d**5*(d*x)**(sympy.S(3)/2)/(1920*b**3*(a + b*x**2)**3) + 77*d**5*(d*x)**(sympy.S(3)/2)/(5120*a*b**3*(a + b*x**2)**2) + 77*d**5*(d*x)**(sympy.S(3)/2)/(4096*a**2*b**3*(a + b*x**2)) + 77*sqrt(2)*d**(sympy.S(13)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(9)/4)*b**(sympy.S(15)/4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(9)/4)*b**(sympy.S(15)/4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(9)/4)*b**(sympy.S(15)/4)) + 77*sqrt(2)*d**(sympy.S(13)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(9)/4)*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_713():
    f = (d*x)**(sympy.S(11)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(9)/2)/(10*b*(a + b*x**2)**5) - 9*d**3*(d*x)**(sympy.S(5)/2)/(160*b**2*(a + b*x**2)**4) - 3*d**5*sqrt(d*x)/(128*b**3*(a + b*x**2)**3) + 3*d**5*sqrt(d*x)/(1024*a*b**3*(a + b*x**2)**2) + 21*d**5*sqrt(d*x)/(4096*a**2*b**3*(a + b*x**2)) - 63*sqrt(2)*d**(sympy.S(11)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(11)/4)*b**(sympy.S(13)/4)) + 63*sqrt(2)*d**(sympy.S(11)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(11)/4)*b**(sympy.S(13)/4)) - 63*sqrt(2)*d**(sympy.S(11)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(11)/4)*b**(sympy.S(13)/4)) + 63*sqrt(2)*d**(sympy.S(11)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(11)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_714():
    f = (d*x)**(sympy.S(9)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(7)/2)/(10*b*(a + b*x**2)**5) - 7*d**3*(d*x)**(sympy.S(3)/2)/(160*b**2*(a + b*x**2)**4) + 7*d**3*(d*x)**(sympy.S(3)/2)/(640*a*b**2*(a + b*x**2)**3) + 63*d**3*(d*x)**(sympy.S(3)/2)/(5120*a**2*b**2*(a + b*x**2)**2) + 63*d**3*(d*x)**(sympy.S(3)/2)/(4096*a**3*b**2*(a + b*x**2)) + 63*sqrt(2)*d**(sympy.S(9)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(13)/4)*b**(sympy.S(11)/4)) - 63*sqrt(2)*d**(sympy.S(9)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(13)/4)*b**(sympy.S(11)/4)) - 63*sqrt(2)*d**(sympy.S(9)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(13)/4)*b**(sympy.S(11)/4)) + 63*sqrt(2)*d**(sympy.S(9)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(13)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_715():
    f = (d*x)**(sympy.S(7)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(5)/2)/(10*b*(a + b*x**2)**5) - d**3*sqrt(d*x)/(32*b**2*(a + b*x**2)**4) + d**3*sqrt(d*x)/(384*a*b**2*(a + b*x**2)**3) + 11*d**3*sqrt(d*x)/(3072*a**2*b**2*(a + b*x**2)**2) + 77*d**3*sqrt(d*x)/(12288*a**3*b**2*(a + b*x**2)) - 77*sqrt(2)*d**(sympy.S(7)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(15)/4)*b**(sympy.S(9)/4)) + 77*sqrt(2)*d**(sympy.S(7)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(15)/4)*b**(sympy.S(9)/4)) - 77*sqrt(2)*d**(sympy.S(7)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(15)/4)*b**(sympy.S(9)/4)) + 77*sqrt(2)*d**(sympy.S(7)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(15)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_716():
    f = (d*x)**(sympy.S(5)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*(d*x)**(sympy.S(3)/2)/(10*b*(a + b*x**2)**5) + 3*d*(d*x)**(sympy.S(3)/2)/(160*a*b*(a + b*x**2)**4) + 13*d*(d*x)**(sympy.S(3)/2)/(640*a**2*b*(a + b*x**2)**3) + 117*d*(d*x)**(sympy.S(3)/2)/(5120*a**3*b*(a + b*x**2)**2) + 117*d*(d*x)**(sympy.S(3)/2)/(4096*a**4*b*(a + b*x**2)) + 117*sqrt(2)*d**(sympy.S(5)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(17)/4)*b**(sympy.S(7)/4)) - 117*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(17)/4)*b**(sympy.S(7)/4)) - 117*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(17)/4)*b**(sympy.S(7)/4)) + 117*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(17)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_717():
    f = (d*x)**(sympy.S(3)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = -d*sqrt(d*x)/(10*b*(a + b*x**2)**5) + d*sqrt(d*x)/(160*a*b*(a + b*x**2)**4) + d*sqrt(d*x)/(128*a**2*b*(a + b*x**2)**3) + 11*d*sqrt(d*x)/(1024*a**3*b*(a + b*x**2)**2) + 77*d*sqrt(d*x)/(4096*a**4*b*(a + b*x**2)) - 231*sqrt(2)*d**(sympy.S(3)/2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(19)/4)*b**(sympy.S(5)/4)) + 231*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(19)/4)*b**(sympy.S(5)/4)) - 231*sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(19)/4)*b**(sympy.S(5)/4)) + 231*sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(19)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_718():
    f = sqrt(d*x)/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = (d*x)**(sympy.S(3)/2)/(10*a*d*(a + b*x**2)**5) + 17*(d*x)**(sympy.S(3)/2)/(160*a**2*d*(a + b*x**2)**4) + 221*(d*x)**(sympy.S(3)/2)/(1920*a**3*d*(a + b*x**2)**3) + 663*(d*x)**(sympy.S(3)/2)/(5120*a**4*d*(a + b*x**2)**2) + 663*(d*x)**(sympy.S(3)/2)/(4096*a**5*d*(a + b*x**2)) + 663*sqrt(2)*sqrt(d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(21)/4)*b**(sympy.S(3)/4)) - 663*sqrt(2)*sqrt(d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(21)/4)*b**(sympy.S(3)/4)) - 663*sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(21)/4)*b**(sympy.S(3)/4)) + 663*sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(21)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_719():
    f = 1/(sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = sqrt(d*x)/(10*a*d*(a + b*x**2)**5) + 19*sqrt(d*x)/(160*a**2*d*(a + b*x**2)**4) + 19*sqrt(d*x)/(128*a**3*d*(a + b*x**2)**3) + 209*sqrt(d*x)/(1024*a**4*d*(a + b*x**2)**2) + 1463*sqrt(d*x)/(4096*a**5*d*(a + b*x**2)) - 4389*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(23)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 4389*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(23)/4)*b**(sympy.S(1)/4)*sqrt(d)) - 4389*sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(23)/4)*b**(sympy.S(1)/4)*sqrt(d)) + 4389*sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(23)/4)*b**(sympy.S(1)/4)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_720():
    f = 1/((d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*d*sqrt(d*x)*(a + b*x**2)**5) + 21/(160*a**2*d*sqrt(d*x)*(a + b*x**2)**4) + 119/(640*a**3*d*sqrt(d*x)*(a + b*x**2)**3) + 1547/(5120*a**4*d*sqrt(d*x)*(a + b*x**2)**2) + 13923/(20480*a**5*d*sqrt(d*x)*(a + b*x**2)) - 13923/(4096*a**6*d*sqrt(d*x)) - 13923*sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(25)/4)*d**(sympy.S(3)/2)) + 13923*sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(25)/4)*d**(sympy.S(3)/2)) + 13923*sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(25)/4)*d**(sympy.S(3)/2)) - 13923*sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(25)/4)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_721():
    f = 1/((d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**5) + 23/(160*a**2*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**4) + 437/(1920*a**3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**3) + 437/(1024*a**4*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**2) + 4807/(4096*a**5*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) - 33649/(12288*a**6*d*(d*x)**(sympy.S(3)/2)) + 33649*sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(27)/4)*d**(sympy.S(5)/2)) - 33649*sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(27)/4)*d**(sympy.S(5)/2)) + 33649*sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(27)/4)*d**(sympy.S(5)/2)) - 33649*sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(27)/4)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_722():
    f = 1/((d*x)**(sympy.S(7)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**3)
    F = 1/(10*a*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**5) + 5/(32*a**2*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**4) + 35/(128*a**3*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**3) + 595/(1024*a**4*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**2) + 7735/(4096*a**5*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 13923/(4096*a**6*d*(d*x)**(sympy.S(5)/2)) + 69615*b/(4096*a**7*d**3*sqrt(d*x)) + 69615*sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(29)/4)*d**(sympy.S(7)/2)) - 69615*sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(32768*a**(sympy.S(29)/4)*d**(sympy.S(7)/2)) - 69615*sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(29)/4)*d**(sympy.S(7)/2)) + 69615*sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(16384*a**(sympy.S(29)/4)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_723():
    f = (d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d*(a + b*x**2)) + 2*b*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_724():
    f = (d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(a + b*x**2)) + 2*b*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_725():
    f = sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = 2*a*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(a + b*x**2)) + 2*b*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_726():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/sqrt(d*x)
    F = 2*a*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)) + 2*b*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_727():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*sqrt(d*x)*(a + b*x**2)) + 2*b*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_728():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(5)/2)
    F = -2*a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) + 2*b*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_729():
    f = sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*x)**(sympy.S(7)/2)
    F = -2*a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*sqrt(d*x)*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_730():
    f = (d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = 2*a**3*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d*(a + b*x**2)) + 6*a**2*b*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**3*(a + b*x**2)) + 2*a*b**2*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(19)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_731():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = 2*a**3*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(a + b*x**2)) + 2*a**2*b*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**3*(a + b*x**2)) + 6*a*b**2*(d*x)**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(17)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_732():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = 2*a**3*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(a + b*x**2)) + 6*a**2*b*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**3*(a + b*x**2)) + 6*a*b**2*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_733():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/sqrt(d*x)
    F = 2*a**3*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)) + 6*a**2*b*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d**3*(a + b*x**2)) + 2*a*b**2*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_734():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(d*x)**(sympy.S(3)/2)
    F = -2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*sqrt(d*x)*(a + b*x**2)) + 2*a**2*b*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)) + 6*a*b**2*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_735():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(d*x)**(sympy.S(5)/2)
    F = -2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) + 6*a**2*b*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)) + 6*a*b**2*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_736():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)/(d*x)**(sympy.S(7)/2)
    F = -2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 6*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*sqrt(d*x)*(a + b*x**2)) + 2*a*b**2*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**5*(a + b*x**2)) + 2*b**3*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**7*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_737():
    f = (d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = 2*a**5*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d*(a + b*x**2)) + 10*a**4*b*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**3*(a + b*x**2)) + 4*a**3*b**2*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(19)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(23)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(23*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(27)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(27*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_738():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = 2*a**5*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(a + b*x**2)) + 10*a**4*b*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*d**3*(a + b*x**2)) + 20*a**3*b**2*(d*x)**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(17)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(21)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(25)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(25*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_739():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = 2*a**5*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(a + b*x**2)) + 10*a**4*b*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**3*(a + b*x**2)) + 20*a**3*b**2*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**5*(a + b*x**2)) + 4*a**2*b**3*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(19)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(23)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(23*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_740():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/sqrt(d*x)
    F = 2*a**5*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)) + 2*a**4*b*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)) + 20*a**3*b**2*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(17)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(21)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(21*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_741():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(d*x)**(sympy.S(3)/2)
    F = -2*a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*sqrt(d*x)*(a + b*x**2)) + 10*a**4*b*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**3*(a + b*x**2)) + 20*a**3*b**2*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**7*(a + b*x**2)) + 2*a*b**4*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(19)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(19*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_742():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(d*x)**(sympy.S(5)/2)
    F = -2*a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)) + 10*a**4*b*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)) + 4*a**3*b**2*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(9)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(9*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(13)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(13*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(17)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(17*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_743():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)/(d*x)**(sympy.S(7)/2)
    F = -2*a**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)) - 10*a**4*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*sqrt(d*x)*(a + b*x**2)) + 20*a**3*b**2*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**5*(a + b*x**2)) + 20*a**2*b**3*(d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(7*d**7*(a + b*x**2)) + 10*a*b**4*(d*x)**(sympy.S(11)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(11*d**9*(a + b*x**2)) + 2*b**5*(d*x)**(sympy.S(15)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(15*d**11*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_744():
    f = (d*x)**(sympy.S(7)/2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(7)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(7)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 2*a*d**3*sqrt(d*x)*(a + b*x**2)/(b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)/(5*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_745():
    f = (d*x)**(sympy.S(5)/2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(5)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(5)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)/(3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_746():
    f = (d*x)**(sympy.S(3)/2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(3)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(3)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*d*sqrt(d*x)*(a + b*x**2)/(b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_747():
    f = sqrt(d*x)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = sqrt(2)*sqrt(d)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*sqrt(d)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_748():
    f = 1/(sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -sqrt(2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_749():
    f = 1/((d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (-2*a - 2*b*x**2)/(a*d*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(5)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(5)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(5)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(5)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_750():
    f = 1/((d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (-2*a - 2*b*x**2)/(3*a*d*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(7)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(7)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(7)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(7)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_751():
    f = 1/((d*x)**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = (-2*a - 2*b*x**2)/(5*a*d*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 2*b*(a + b*x**2)/(a**2*d**3*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(9)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(4*a**(sympy.S(9)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(9)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(2*a**(sympy.S(9)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_752():
    f = (d*x)**(sympy.S(15)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -117*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(15)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(15)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 117*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(15)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(15)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 117*a*d**7*sqrt(d*x)*(a + b*x**2)/(16*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(13)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13*d**3*(d*x)**(sympy.S(9)/2)/(16*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*d**5*(d*x)**(sympy.S(5)/2)*(a + b*x**2)/(80*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_753():
    f = (d*x)**(sympy.S(13)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -77*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(13)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(13)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(13)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(13)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(11)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 11*d**3*(d*x)**(sympy.S(7)/2)/(16*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*d**5*(d*x)**(sympy.S(3)/2)*(a + b*x**2)/(48*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_754():
    f = (d*x)**(sympy.S(11)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = 45*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(11)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(11)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(11)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(11)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(9)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 9*d**3*(d*x)**(sympy.S(5)/2)/(16*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*d**5*sqrt(d*x)*(a + b*x**2)/(16*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_755():
    f = (d*x)**(sympy.S(9)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -d*(d*x)**(sympy.S(7)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7*d**3*(d*x)**(sympy.S(3)/2)/(16*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 21*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 21*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 21*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 21*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_756():
    f = (d*x)**(sympy.S(7)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -d*(d*x)**(sympy.S(5)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*d**3*sqrt(d*x)/(16*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_757():
    f = (d*x)**(sympy.S(5)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -d*(d*x)**(sympy.S(3)/2)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*d*(d*x)**(sympy.S(3)/2)/(16*a*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_758():
    f = (d*x)**(sympy.S(3)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -d*sqrt(d*x)/(4*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + d*sqrt(d*x)/(16*a*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_759():
    f = sqrt(d*x)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = (d*x)**(sympy.S(3)/2)/(4*a*d*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*(d*x)**(sympy.S(3)/2)/(16*a**2*d*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*sqrt(2)*sqrt(d)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*sqrt(2)*sqrt(d)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_760():
    f = 1/(sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = sqrt(d*x)/(4*a*d*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7*sqrt(d*x)/(16*a**2*d*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*(21*a + 21*b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(21*a + 21*b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*(21*a + 21*b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(21*a + 21*b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_761():
    f = 1/((d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*d*sqrt(d*x)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 9/(16*a**2*d*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (45*a + 45*b*x**2)/(16*a**3*d*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(13)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(13)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(13)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(13)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_762():
    f = 1/((d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 11/(16*a**2*d*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (77*a + 77*b*x**2)/(48*a**3*d*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(15)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(15)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(15)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(15)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_763():
    f = 1/((d*x)**(sympy.S(7)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = 1/(4*a*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13/(16*a**2*d*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (117*a + 117*b*x**2)/(80*a**3*d*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*b*(a + b*x**2)/(16*a**4*d**3*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(17)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 117*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(128*a**(sympy.S(17)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 117*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(17)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 117*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(64*a**(sympy.S(17)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_764():
    f = (d*x)**(sympy.S(23)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -13923*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(23)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(25)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(23)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(25)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13923*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(23)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(25)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*sqrt(2)*a**(sympy.S(5)/4)*d**(sympy.S(23)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(25)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13923*a*d**11*sqrt(d*x)*(a + b*x**2)/(1024*b**6*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(21)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7*d**3*(d*x)**(sympy.S(17)/2)/(32*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 119*d**5*(d*x)**(sympy.S(13)/2)/(256*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 1547*d**7*(d*x)**(sympy.S(9)/2)/(1024*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*d**9*(d*x)**(sympy.S(5)/2)*(a + b*x**2)/(5120*b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_765():
    f = (d*x)**(sympy.S(21)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -7315*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(21)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(23)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7315*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(21)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(23)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7315*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(21)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(23)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7315*sqrt(2)*a**(sympy.S(3)/4)*d**(sympy.S(21)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(23)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(19)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 19*d**3*(d*x)**(sympy.S(15)/2)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 95*d**5*(d*x)**(sympy.S(11)/2)/(256*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 1045*d**7*(d*x)**(sympy.S(7)/2)/(1024*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7315*d**9*(d*x)**(sympy.S(3)/2)*(a + b*x**2)/(3072*b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_766():
    f = (d*x)**(sympy.S(19)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = 3315*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(19)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(21)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3315*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(19)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*b**(sympy.S(21)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3315*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(19)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(21)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3315*sqrt(2)*a**(sympy.S(1)/4)*d**(sympy.S(19)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*b**(sympy.S(21)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(d*x)**(sympy.S(17)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 17*d**3*(d*x)**(sympy.S(13)/2)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 221*d**5*(d*x)**(sympy.S(9)/2)/(768*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 663*d**7*(d*x)**(sympy.S(5)/2)/(1024*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3315*d**9*sqrt(d*x)*(a + b*x**2)/(1024*b**5*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_767():
    f = (d*x)**(sympy.S(17)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(15)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*d**3*(d*x)**(sympy.S(11)/2)/(32*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 55*d**5*(d*x)**(sympy.S(7)/2)/(256*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 385*d**7*(d*x)**(sympy.S(3)/2)/(1024*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1155*sqrt(2)*d**(sympy.S(17)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 1155*sqrt(2)*d**(sympy.S(17)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 1155*sqrt(2)*d**(sympy.S(17)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1155*sqrt(2)*d**(sympy.S(17)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_768():
    f = (d*x)**(sympy.S(15)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(13)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13*d**3*(d*x)**(sympy.S(9)/2)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 39*d**5*(d*x)**(sympy.S(5)/2)/(256*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 195*d**7*sqrt(d*x)/(1024*b**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 195*sqrt(2)*d**(sympy.S(15)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 195*sqrt(2)*d**(sympy.S(15)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 195*sqrt(2)*d**(sympy.S(15)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 195*sqrt(2)*d**(sympy.S(15)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_769():
    f = (d*x)**(sympy.S(13)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(11)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 11*d**3*(d*x)**(sympy.S(7)/2)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*d**5*(d*x)**(sympy.S(3)/2)/(768*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*d**5*(d*x)**(sympy.S(3)/2)/(1024*a*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*d**(sympy.S(13)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*d**(sympy.S(13)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*d**(sympy.S(13)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_770():
    f = (d*x)**(sympy.S(11)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(9)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3*d**3*(d*x)**(sympy.S(5)/2)/(32*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 15*d**5*sqrt(d*x)/(256*b**3*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 15*d**5*sqrt(d*x)/(1024*a*b**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*d**(sympy.S(11)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*d**(sympy.S(11)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*d**(sympy.S(11)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*d**(sympy.S(11)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_771():
    f = (d*x)**(sympy.S(9)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(7)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7*d**3*(d*x)**(sympy.S(3)/2)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7*d**3*(d*x)**(sympy.S(3)/2)/(256*a*b**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*d**3*(d*x)**(sympy.S(3)/2)/(1024*a**2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 35*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 35*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*sqrt(2)*d**(sympy.S(9)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_772():
    f = (d*x)**(sympy.S(7)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(5)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 5*d**3*sqrt(d*x)/(96*b**2*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*d**3*sqrt(d*x)/(768*a*b**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*d**3*sqrt(d*x)/(3072*a**2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 35*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 35*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 35*sqrt(2)*d**(sympy.S(7)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_773():
    f = (d*x)**(sympy.S(5)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*(d*x)**(sympy.S(3)/2)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + d*(d*x)**(sympy.S(3)/2)/(32*a*b*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 9*d*(d*x)**(sympy.S(3)/2)/(256*a**2*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*d*(d*x)**(sympy.S(3)/2)/(1024*a**3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 45*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 45*sqrt(2)*d**(sympy.S(5)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_774():
    f = (d*x)**(sympy.S(3)/2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = -d*sqrt(d*x)/(8*b*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + d*sqrt(d*x)/(96*a*b*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 11*d*sqrt(d*x)/(768*a**2*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*d*sqrt(d*x)/(3072*a**3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 77*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 77*sqrt(2)*d**(sympy.S(3)/2)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_775():
    f = sqrt(d*x)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = (d*x)**(sympy.S(3)/2)/(8*a*d*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13*(d*x)**(sympy.S(3)/2)/(96*a**2*d*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 39*(d*x)**(sympy.S(3)/2)/(256*a**3*d*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 195*(d*x)**(sympy.S(3)/2)/(1024*a**4*d*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 195*sqrt(2)*sqrt(d)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(17)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 195*sqrt(2)*sqrt(d)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(17)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 195*sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(17)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 195*sqrt(2)*sqrt(d)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(17)/4)*b**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_776():
    f = 1/(sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = sqrt(d*x)/(8*a*d*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 5*sqrt(d*x)/(32*a**2*d*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 55*sqrt(d*x)/(256*a**3*d*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 385*sqrt(d*x)/(1024*a**4*d*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*(1155*a + 1155*b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(19)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(1155*a + 1155*b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(19)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - sqrt(2)*(1155*a + 1155*b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(19)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + sqrt(2)*(1155*a + 1155*b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(19)/4)*b**(sympy.S(1)/4)*sqrt(d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_777():
    f = 1/((d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*d*sqrt(d*x)*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 17/(96*a**2*d*sqrt(d*x)*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 221/(768*a**3*d*sqrt(d*x)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 663/(1024*a**4*d*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (3315*a + 3315*b*x**2)/(1024*a**5*d*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3315*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(21)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3315*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(21)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 3315*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(21)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 3315*sqrt(2)*b**(sympy.S(1)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(21)/4)*d**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_778():
    f = 1/((d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 19/(96*a**2*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 95/(256*a**3*d*(d*x)**(sympy.S(3)/2)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1045/(1024*a**4*d*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (7315*a + 7315*b*x**2)/(3072*a**5*d*(d*x)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7315*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(23)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7315*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(23)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7315*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(23)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 7315*sqrt(2)*b**(sympy.S(3)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(23)/4)*d**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_779():
    f = 1/((d*x)**(sympy.S(7)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2))
    F = 1/(8*a*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 7/(32*a**2*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 119/(256*a**3*d*(d*x)**(sympy.S(5)/2)*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 1547/(1024*a**4*d*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (13923*a + 13923*b*x**2)/(5120*a**5*d*(d*x)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*b*(a + b*x**2)/(1024*a**6*d**3*sqrt(d*x)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(25)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13923*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(d*x) + sqrt(a)*sqrt(d) + sqrt(b)*sqrt(d)*x)/(8192*a**(sympy.S(25)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - 13923*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(25)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + 13923*sqrt(2)*b**(sympy.S(5)/4)*(a + b*x**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(d*x)/(a**(sympy.S(1)/4)*sqrt(d)))/(4096*a**(sympy.S(25)/4)*d**(sympy.S(7)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_780():
    f = (d*x)**m*(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = a**6*(d*x)**(m + 1)/(d*(m + 1)) + 6*a**5*b*(d*x)**(m + 3)/(d**3*(m + 3)) + 15*a**4*b**2*(d*x)**(m + 5)/(d**5*(m + 5)) + 20*a**3*b**3*(d*x)**(m + 7)/(d**7*(m + 7)) + 15*a**2*b**4*(d*x)**(m + 9)/(d**9*(m + 9)) + 6*a*b**5*(d*x)**(m + 11)/(d**11*(m + 11)) + b**6*(d*x)**(m + 13)/(d**13*(m + 13))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_781():
    f = (d*x)**m*(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = a**4*(d*x)**(m + 1)/(d*(m + 1)) + 4*a**3*b*(d*x)**(m + 3)/(d**3*(m + 3)) + 6*a**2*b**2*(d*x)**(m + 5)/(d**5*(m + 5)) + 4*a*b**3*(d*x)**(m + 7)/(d**7*(m + 7)) + b**4*(d*x)**(m + 9)/(d**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_782():
    f = (d*x)**m*(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a**2*(d*x)**(m + 1)/(d*(m + 1)) + 2*a*b*(d*x)**(m + 3)/(d**3*(m + 3)) + b**2*(d*x)**(m + 5)/(d**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_783():
    f = (d*x)**m/(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (d*x)**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**2*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_784():
    f = (d*x)**m/(a**2 + 2*a*b*x**2 + b**2*x**4)**2
    F = (d*x)**(m + 1)*hyper((4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**4*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_785():
    f = (d*x)**m/(a**2 + 2*a*b*x**2 + b**2*x**4)**3
    F = (d*x)**(m + 1)*hyper((6, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**6*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_786():
    f = (d*x)**m*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)*(m + 1)) + 5*a**4*b*(d*x)**(m + 3)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)*(m + 3)) + 10*a**3*b**2*(d*x)**(m + 5)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**5*(a + b*x**2)*(m + 5)) + 10*a**2*b**3*(d*x)**(m + 7)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**7*(a + b*x**2)*(m + 7)) + 5*a*b**4*(d*x)**(m + 9)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**9*(a + b*x**2)*(m + 9)) + b**5*(d*x)**(m + 11)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**11*(a + b*x**2)*(m + 11))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_787():
    f = (d*x)**m*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)*(m + 1)) + 3*a**2*b*(d*x)**(m + 3)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)*(m + 3)) + 3*a*b**2*(d*x)**(m + 5)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**5*(a + b*x**2)*(m + 5)) + b**3*(d*x)**(m + 7)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**7*(a + b*x**2)*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_788():
    f = (d*x)**m*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d*(a + b*x**2)*(m + 1)) + b*(d*x)**(m + 3)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(d**3*(a + b*x**2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_789():
    f = (d*x)**m/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = (d*x)**(m + 1)*(a + b*x**2)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*d*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_790():
    f = (d*x)**m/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*(a + b*x**2)*hyper((3, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**3*d*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_791():
    f = (d*x)**m/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = (d*x)**(m + 1)*(a + b*x**2)*hyper((5, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a**5*d*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_792():
    f = x**7*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = -a**3*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**4*(2*p + 1)) + 3*a**2*(a + b*x**2)**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**4*(p + 1)) - 3*a*(a + b*x**2)**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**4*(2*p + 3)) + (a + b*x**2)**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**4*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_793():
    f = x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = a**2*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**3*(2*p + 1)) - a*(a + b*x**2)**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**3*(p + 1)) + (a + b*x**2)**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**3*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_794():
    f = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = -a*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**2*(2*p + 1)) + (a + b*x**2)**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_795():
    f = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = (a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_796():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/x
    F = -(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((1, 2*p + 1), (2*p + 2,), 1 + b*x**2/a)/(2*a*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_797():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/x**3
    F = b*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((2, 2*p + 1), (2*p + 2,), 1 + b*x**2/a)/(2*a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_798():
    f = x**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = x**5*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(5)/2, -2*p), (sympy.S(7)/2,), -b*x**2/a)/(5*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_799():
    f = x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(3)/2, -2*p), (sympy.S(5)/2,), -b*x**2/a)/(3*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_800():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = x*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S.Half, -2*p), (sympy.S(3)/2,), -b*x**2/a)/(1 + b*x**2/a)**(2*p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_801():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/x**2
    F = -(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(-1)/2, -2*p), (sympy.S.Half,), -b*x**2/a)/(x*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_802():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/x**4
    F = -(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(-3)/2, -2*p), (sympy.S(-1)/2,), -b*x**2/a)/(3*x**3*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_803():
    f = (d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = 2*(d*x)**(sympy.S(5)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(5)/4, -2*p), (sympy.S(9)/4,), -b*x**2/a)/(5*d*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_804():
    f = sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = 2*(d*x)**(sympy.S(3)/2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(3)/4, -2*p), (sympy.S(7)/4,), -b*x**2/a)/(3*d*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_805():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/sqrt(d*x)
    F = 2*sqrt(d*x)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(1)/4, -2*p), (sympy.S(5)/4,), -b*x**2/a)/(d*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_806():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/(d*x)**(sympy.S(3)/2)
    F = -2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(-1)/4, -2*p), (sympy.S(3)/4,), -b*x**2/a)/(d*sqrt(d*x)*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_807():
    f = (a**2 + 2*a*b*x**2 + b**2*x**4)**p/(d*x)**(sympy.S(5)/2)
    F = -2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p*hyper((sympy.S(-3)/4, -2*p), (sympy.S(1)/4,), -b*x**2/a)/(3*d*(d*x)**(sympy.S(3)/2)*(1 + b*x**2/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_808():
    f = x**2*(a + b*x**2 + c*x**4)
    F = a*x**3/3 + b*x**5/5 + c*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_809():
    f = x*(a + b*x**2 + c*x**4)
    F = a*x**2/2 + b*x**4/4 + c*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_810():
    f = a + b*x**2 + c*x**4
    F = a*x + b*x**3/3 + c*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_811():
    f = (a + b*x**2 + c*x**4)/x
    F = a*log(x) + b*x**2/2 + c*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_812():
    f = (a + b*x**2 + c*x**4)/x**2
    F = -a/x + b*x + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_813():
    f = (a + b*x**2 + c*x**4)/x**3
    F = -a/(2*x**2) + b*log(x) + c*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_814():
    f = (a + b*x**2 + c*x**4)/x**4
    F = -a/(3*x**3) - b/x + c*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_815():
    f = (a + b*x**2 + c*x**4)/x**5
    F = -a/(4*x**4) - b/(2*x**2) + c*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_816():
    f = (a + b*x**2 + c*x**4)/x**6
    F = -a/(5*x**5) - b/(3*x**3) - c/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_817():
    f = (a + b*x**2 + c*x**4)/x**7
    F = -a/(6*x**6) - b/(4*x**4) - c/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_818():
    f = (a + b*x**2 + c*x**4)/x**8
    F = -a/(7*x**7) - b/(5*x**5) - c/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_819():
    f = x**2*(a + b*x**2 + c*x**4)**2
    F = a**2*x**3/3 + 2*a*b*x**5/5 + 2*b*c*x**9/9 + c**2*x**11/11 + x**7*(2*a*c/7 + b**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_820():
    f = x*(a + b*x**2 + c*x**4)**2
    F = a**2*x**2/2 + a*b*x**4/2 + b*c*x**8/4 + c**2*x**10/10 + x**6*(a*c/3 + b**2/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_821():
    f = (a + b*x**2 + c*x**4)**2
    F = a**2*x + 2*a*b*x**3/3 + 2*b*c*x**7/7 + c**2*x**9/9 + x**5*(2*a*c/5 + b**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_822():
    f = (a + b*x**2 + c*x**4)**2/x
    F = a**2*log(x) + a*b*x**2 + b*c*x**6/3 + c**2*x**8/8 + x**4*(a*c/2 + b**2/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_823():
    f = (a + b*x**2 + c*x**4)**2/x**2
    F = -a**2/x + 2*a*b*x + 2*b*c*x**5/5 + c**2*x**7/7 + x**3*(2*a*c/3 + b**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_824():
    f = (a + b*x**2 + c*x**4)**2/x**3
    F = -a**2/(2*x**2) + 2*a*b*log(x) + b*c*x**4/2 + c**2*x**6/6 + x**2*(a*c + b**2/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_825():
    f = (a + b*x**2 + c*x**4)**2/x**4
    F = -a**2/(3*x**3) - 2*a*b/x + 2*b*c*x**3/3 + c**2*x**5/5 + x*(2*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_826():
    f = (a + b*x**2 + c*x**4)**2/x**5
    F = -a**2/(4*x**4) - a*b/x**2 + b*c*x**2 + c**2*x**4/4 + (2*a*c + b**2)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_827():
    f = (a + b*x**2 + c*x**4)**2/x**6
    F = -a**2/(5*x**5) - 2*a*b/(3*x**3) + 2*b*c*x + c**2*x**3/3 - (2*a*c + b**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_828():
    f = (a + b*x**2 + c*x**4)**2/x**7
    F = -a**2/(6*x**6) - a*b/(2*x**4) + 2*b*c*log(x) + c**2*x**2/2 - (2*a*c + b**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_829():
    f = (a + b*x**2 + c*x**4)**2/x**8
    F = -a**2/(7*x**7) - 2*a*b/(5*x**5) - 2*b*c/x + c**2*x - (2*a*c + b**2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_830():
    f = (a + b*x**2 + c*x**4)**2/x**9
    F = -a**2/(8*x**8) - a*b/(3*x**6) - b*c/x**2 + c**2*log(x) - (2*a*c + b**2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_831():
    f = (a + b*x**2 + c*x**4)**2/x**10
    F = -a**2/(9*x**9) - 2*a*b/(7*x**7) - 2*b*c/(3*x**3) - c**2/x - (2*a*c + b**2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_832():
    f = (a + b*x**2 + c*x**4)**2/x**11
    F = -a**2/(10*x**10) - a*b/(4*x**8) - b*c/(2*x**4) - c**2/(2*x**2) - (2*a*c + b**2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_833():
    f = (a + b*x**2 + c*x**4)**2/x**12
    F = -a**2/(11*x**11) - 2*a*b/(9*x**9) - 2*b*c/(5*x**5) - c**2/(3*x**3) - (2*a*c + b**2)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_834():
    f = (a + b*x**2 + c*x**4)**2/x**13
    F = -a**2/(12*x**12) - a*b/(5*x**10) - b*c/(3*x**6) - c**2/(4*x**4) - (2*a*c + b**2)/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_835():
    f = x**2*(a + b*x**2 + c*x**4)**3
    F = a**3*x**3/3 + 3*a**2*b*x**5/5 + 3*a*x**7*(a*c + b**2)/7 + 3*b*c**2*x**13/13 + b*x**9*(6*a*c + b**2)/9 + c**3*x**15/15 + 3*c*x**11*(a*c + b**2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_836():
    f = x*(a + b*x**2 + c*x**4)**3
    F = a**3*x**2/2 + 3*a**2*b*x**4/4 + a*x**6*(a*c + b**2)/2 + b*c**2*x**12/4 + b*x**8*(6*a*c + b**2)/8 + c**3*x**14/14 + 3*c*x**10*(a*c + b**2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_837():
    f = (a + b*x**2 + c*x**4)**3
    F = a**3*x + a**2*b*x**3 + 3*a*x**5*(a*c + b**2)/5 + 3*b*c**2*x**11/11 + b*x**7*(6*a*c + b**2)/7 + c**3*x**13/13 + c*x**9*(a*c + b**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_838():
    f = (a + b*x**2 + c*x**4)**3/x
    F = a**3*log(x) + 3*a**2*b*x**2/2 + 3*a*x**4*(a*c + b**2)/4 + 3*b*c**2*x**10/10 + b*x**6*(6*a*c + b**2)/6 + c**3*x**12/12 + 3*c*x**8*(a*c + b**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_839():
    f = (a + b*x**2 + c*x**4)**3/x**2
    F = -a**3/x + 3*a**2*b*x + a*x**3*(a*c + b**2) + b*c**2*x**9/3 + b*x**5*(6*a*c + b**2)/5 + c**3*x**11/11 + 3*c*x**7*(a*c + b**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_840():
    f = (a + b*x**2 + c*x**4)**3/x**3
    F = -a**3/(2*x**2) + 3*a**2*b*log(x) + 3*a*x**2*(a*c + b**2)/2 + 3*b*c**2*x**8/8 + b*x**4*(6*a*c + b**2)/4 + c**3*x**10/10 + c*x**6*(a*c + b**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_841():
    f = (a + b*x**2 + c*x**4)**3/x**4
    F = -a**3/(3*x**3) - 3*a**2*b/x + 3*a*x*(a*c + b**2) + 3*b*c**2*x**7/7 + b*x**3*(6*a*c + b**2)/3 + c**3*x**9/9 + 3*c*x**5*(a*c + b**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_842():
    f = x**7/(a + b*x**2 + c*x**4)
    F = -b*x**2/(2*c**2) + b*(-3*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2)) + x**4/(4*c) + (-a*c + b**2)*log(a + b*x**2 + c*x**4)/(4*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_843():
    f = x**5/(a + b*x**2 + c*x**4)
    F = -b*log(a + b*x**2 + c*x**4)/(4*c**2) + x**2/(2*c) - (-2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_844():
    f = x**3/(a + b*x**2 + c*x**4)
    F = b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + log(a + b*x**2 + c*x**4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_845():
    f = x/(a + b*x**2 + c*x**4)
    F = -atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_846():
    f = 1/(x*(a + b*x**2 + c*x**4))
    F = b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x**2 + c*x**4)/(4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_847():
    f = 1/(x**3*(a + b*x**2 + c*x**4))
    F = -1/(2*a*x**2) - b*log(x)/a**2 + b*log(a + b*x**2 + c*x**4)/(4*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_848():
    f = 1/(x**5*(a + b*x**2 + c*x**4))
    F = -1/(4*a*x**4) + b/(2*a**2*x**2) + b*(-3*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*sqrt(-4*a*c + b**2)) + (-a*c + b**2)*log(x)/a**3 - (-a*c + b**2)*log(a + b*x**2 + c*x**4)/(4*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_849():
    f = x**6/(a + b*x**2 + c*x**4)
    F = -b*x/c**2 + x**3/(3*c) + sqrt(2)*(-a*c + b**2 + b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-a*c + b**2 - b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_850():
    f = x**4/(a + b*x**2 + c*x**4)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_851():
    f = x**2/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_852():
    f = 1/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_853():
    f = 1/(x**2*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_854():
    f = 1/(x**4*(a + b*x**2 + c*x**4))
    F = -1/(3*a*x**3) + b/(a**2*x) + sqrt(2)*sqrt(c)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*sqrt(c)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_855():
    f = x**7/(a + b*x**2 + c*x**4)**2
    F = -b*x**2/(2*c*(-4*a*c + b**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + x**4*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + log(a + b*x**2 + c*x**4)/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_856():
    f = x**5/(a + b*x**2 + c*x**4)**2
    F = 2*a*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + x**2*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_857():
    f = x**3/(a + b*x**2 + c*x**4)**2
    F = -b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_858():
    f = x/(a + b*x**2 + c*x**4)**2
    F = 2*c*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_859():
    f = 1/(x*(a + b*x**2 + c*x**4)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(x)/a**2 - log(a + b*x**2 + c*x**4)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_860():
    f = 1/(x**3*(a + b*x**2 + c*x**4)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-3*a*c + b**2)/(a**2*x**2*(-4*a*c + b**2)) - 2*b*log(x)/a**3 + b*log(a + b*x**2 + c*x**4)/(2*a**3) - (6*a**2*c**2 - 6*a*b**2*c + b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_861():
    f = x**8/(a + b*x**2 + c*x**4)**2
    F = -b*x**3/(2*c*(-4*a*c + b**2)) + x**5*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + x*(-10*a*c + 3*b**2)/(2*c**2*(-4*a*c + b**2)) - sqrt(2)*(-13*a*b*c + 3*b**3 + (20*a**2*c**2 - 19*a*b**2*c + 3*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(-13*a*b*c + 3*b**3 - (20*a**2*c**2 - 19*a*b**2*c + 3*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_862():
    f = x**6/(a + b*x**2 + c*x**4)**2
    F = -b*x/(2*c*(-4*a*c + b**2)) + x**3*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(-6*a*c + b**2 + b*(-8*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(-6*a*c + b**2 - b*(-8*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_863():
    f = x**4/(a + b*x**2 + c*x**4)**2
    F = x*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(b - (4*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_864():
    f = x**2/(a + b*x**2 + c*x**4)**2
    F = -sqrt(2)*sqrt(c)*(2*b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(2*b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - x*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_865():
    f = (a + b*x**2 + c*x**4)**(-2)
    F = -sqrt(2)*sqrt(c)*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + x*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_866():
    f = 1/(x**2*(a + b*x**2 + c*x**4)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 - (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 + (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*a*c + 3*b**2)/(2*a**2*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_867():
    f = x**11/(a + b*x**2 + c*x**4)**3
    F = -b*x**2*(-7*a*c + b**2)/(2*c**2*(-4*a*c + b**2)**2) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*(-4*a*c + b**2)**(sympy.S(5)/2)) + x**8*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x**4*(a*(-16*a*c + b**2) + b*x**2*(-10*a*c + b**2))/(4*c*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + log(a + b*x**2 + c*x**4)/(4*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_868():
    f = x**9/(a + b*x**2 + c*x**4)**3
    F = -6*a**2*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - 3*a*x**2*(2*a + b*x**2)/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + x**6*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_869():
    f = x**7/(a + b*x**2 + c*x**4)**3
    F = 3*a*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*b*x**2*(2*a + b*x**2)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - x**6*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_870():
    f = x**5/(a + b*x**2 + c*x**4)**3
    F = x**2*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + (3*a*b + x**2*(2*a*c + b**2))/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_871():
    f = x**3/(a + b*x**2 + c*x**4)**3
    F = 3*b*c*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - 3*b*(b + 2*c*x**2)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + (2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_872():
    f = x/(a + b*x**2 + c*x**4)**3
    F = -6*c**2*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*c*(b + 2*c*x**2)/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_873():
    f = 1/(x*(a + b*x**2 + c*x**4)**3)
    F = (-2*a*c + b**2 + b*c*x**2)/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + (16*a**2*c**2 - 15*a*b**2*c + 2*b**4 + 2*b*c*x**2*(-7*a*c + b**2))/(4*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*(-4*a*c + b**2)**(sympy.S(5)/2)) + log(x)/a**3 - log(a + b*x**2 + c*x**4)/(4*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_874():
    f = 1/(x**3*(a + b*x**2 + c*x**4)**3)
    F = (-2*a*c + b**2 + b*c*x**2)/(4*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + (20*a**2*c**2 - 20*a*b**2*c + 3*b**4 + 3*b*c*x**2*(-6*a*c + b**2))/(4*a**2*x**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-15*a*c + 3*b**2)*(-2*a*c + b**2)/(2*a**3*x**2*(-4*a*c + b**2)**2) - 3*b*log(x)/a**4 + 3*b*log(a + b*x**2 + c*x**4)/(4*a**4) - (-60*a**3*c**3 + 90*a**2*b**2*c**2 - 30*a*b**4*c + 3*b**6)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**4*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_875():
    f = x**10/(a + b*x**2 + c*x**4)**3
    F = -3*b*x*(-8*a*c + b**2)/(8*c**2*(-4*a*c + b**2)**2) + x**7*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x**5*(12*a*b - x**2*(-28*a*c + b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + x**3*(-28*a*c + b**2)/(8*c*(-4*a*c + b**2)**2) + sqrt(2)*(84*a**2*c**2 - 27*a*b**2*c + 3*b**4 + 3*(44*a**2*b*c**2 - 11*a*b**3*c + b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(84*a**2*c**2 - 27*a*b**2*c + 3*b**4 - 3*(44*a**2*b*c**2 - 11*a*b**3*c + b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_876():
    f = x**8/(a + b*x**2 + c*x**4)**3
    F = x**5*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x**3*(12*a*b + x**2*(20*a*c + b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - x*(20*a*c + b**2)/(8*c*(-4*a*c + b**2)**2) + sqrt(2)*(-16*a*b*c + b**3 + (-40*a**2*c**2 - 18*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(-16*a*b*c + b**3 - (-40*a**2*c**2 - 18*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_877():
    f = x**6/(a + b*x**2 + c*x**4)**3
    F = x**3*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + 3*x*(4*a*b + x**2*(4*a*c + b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + sqrt(2)*(12*a*c + 3*b**2 + 3*b*(12*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(12*a*c + 3*b**2 - 3*b*(12*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_878():
    f = x**4/(a + b*x**2 + c*x**4)**3
    F = -3*sqrt(2)*sqrt(c)*(4*a*c + 3*b**2 + 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(8*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*sqrt(2)*sqrt(c)*(4*a*c + 3*b**2 - 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(8*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + x*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - x*(-4*a*c + 7*b**2 + 12*b*c*x**2)/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_879():
    f = x**2/(a + b*x**2 + c*x**4)**3
    F = -x*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a*c + b**2 - b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(20*a*c + b**2 + b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*(b*(8*a*c + b**2) + c*x**2*(20*a*c + b**2))/(8*a*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_880():
    f = (a + b*x**2 + c*x**4)**(-3)
    F = x*(-2*a*c + b**2 + b*c*x**2)/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + 3*sqrt(2)*sqrt(c)*(-8*a*b*c + b**3 - (56*a**2*c**2 - 10*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + 3*sqrt(2)*sqrt(c)*(56*a**2*c**2 - 10*a*b**2*c + b**4 + b*(-8*a*c + b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + x*(3*b*c*x**2*(-8*a*c + b**2) + (-7*a*c + b**2)*(-4*a*c + 3*b**2))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_881():
    f = 1/(x**2*(a + b*x**2 + c*x**4)**3)
    F = (-2*a*c + b**2 + b*c*x**2)/(4*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + (36*a**2*c**2 - 35*a*b**2*c + 5*b**4 + b*c*x**2*(-32*a*c + 5*b**2))/(8*a**2*x*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - 3*sqrt(2)*sqrt(c)*((-12*a*c + 5*b**2)*(-5*a*c + b**2) - (124*a**2*b*c**2 - 47*a*b**3*c + 5*b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**3*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - 3*sqrt(2)*sqrt(c)*(b*(124*a**2*c**2 - 47*a*b**2*c + 5*b**4)/sqrt(-4*a*c + b**2) + (-12*a*c + 5*b**2)*(-5*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**3*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - (-36*a*c + 15*b**2)*(-5*a*c + b**2)/(8*a**3*x*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_882():
    f = x**5/(a - b*x**2 + c*x**4)
    F = b*log(a - b*x**2 + c*x**4)/(4*c**2) + x**2/(2*c) + (-2*a*c + b**2)*atanh((b - 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_883():
    f = x**3/(a - b*x**2 + c*x**4)
    F = b*atanh((b - 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + log(a - b*x**2 + c*x**4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_884():
    f = x/(a - b*x**2 + c*x**4)
    F = atanh((b - 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_885():
    f = 1/(x*(a - b*x**2 + c*x**4))
    F = b*atanh((b - 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a - b*x**2 + c*x**4)/(4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_886():
    f = 1/(x**3*(a - b*x**2 + c*x**4))
    F = -1/(2*a*x**2) + b*log(x)/a**2 - b*log(a - b*x**2 + c*x**4)/(4*a**2) + (-2*a*c + b**2)*atanh((b - 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_887():
    f = x**4/(a - b*x**2 + c*x**4)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_888():
    f = x**2/(a - b*x**2 + c*x**4)
    F = sqrt(2)*sqrt(b - sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) - sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_889():
    f = 1/(a - b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_890():
    f = 1/(x**2*(a - b*x**2 + c*x**4))
    F = sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atanh(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atanh(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_891():
    f = x**5/(a*x**4 + 2*a*x**2 + a - b)
    F = x**2/(2*a) - log(a*x**4 + 2*a*x**2 + a - b)/(2*a) - (a + b)*atanh(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_892():
    f = x**3/(a*x**4 + 2*a*x**2 + a - b)
    F = log(a*x**4 + 2*a*x**2 + a - b)/(4*a) + atanh(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_893():
    f = x/(a*x**4 + 2*a*x**2 + a - b)
    F = -atanh(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_894():
    f = 1/(x*(a*x**4 + 2*a*x**2 + a - b))
    F = sqrt(a)*atanh(sqrt(a)*(x**2 + 1)/sqrt(b))/(sqrt(b)*(2*a - 2*b)) - log(a*x**4 + 2*a*x**2 + a - b)/(4*a - 4*b) + log(x)/(a - b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_895():
    f = 1/(x**3*(a*x**4 + 2*a*x**2 + a - b))
    F = -sqrt(a)*(a + b)*atanh(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(b)*(a - b)**2) - 2*a*log(x)/(a - b)**2 + a*log(a*x**4 + 2*a*x**2 + a - b)/(2*(a - b)**2) - 1/(x**2*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_896():
    f = x**4/(a*x**4 + 2*a*x**2 + a - b)
    F = x/a + (sqrt(a) - sqrt(b))**(sympy.S(3)/2)*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(5)/4)*sqrt(b)) - (sqrt(a) + sqrt(b))**(sympy.S(3)/2)*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(5)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_897():
    f = x**2/(a*x**4 + 2*a*x**2 + a - b)
    F = -sqrt(sqrt(a) - sqrt(b))*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(3)/4)*sqrt(b)) + sqrt(sqrt(a) + sqrt(b))*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(3)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_898():
    f = 1/(a*x**4 + 2*a*x**2 + a - b)
    F = -atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(1)/4)*sqrt(b)*sqrt(sqrt(a) + sqrt(b))) + atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(1)/4)*sqrt(b)*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_899():
    f = 1/(x**2*(a*x**4 + 2*a*x**2 + a - b))
    F = a**(sympy.S(1)/4)*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) + sqrt(b)))/(2*sqrt(b)*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - a**(sympy.S(1)/4)*atan(a**(sympy.S(1)/4)*x/sqrt(sqrt(a) - sqrt(b)))/(2*sqrt(b)*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) - 1/(x*(a - b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_900():
    f = x**5/(a*x**4 + 2*a*x**2 + a + b)
    F = x**2/(2*a) - log(a*x**4 + 2*a*x**2 + a + b)/(2*a) + (a - b)*atan(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_901():
    f = x**3/(a*x**4 + 2*a*x**2 + a + b)
    F = log(a*x**4 + 2*a*x**2 + a + b)/(4*a) - atan(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_902():
    f = x/(a*x**4 + 2*a*x**2 + a + b)
    F = atan(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_903():
    f = 1/(x*(a*x**4 + 2*a*x**2 + a + b))
    F = -sqrt(a)*atan(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(b)*(a + b)) - log(a*x**4 + 2*a*x**2 + a + b)/(4*a + 4*b) + log(x)/(a + b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_904():
    f = 1/(x**3*(a*x**4 + 2*a*x**2 + a + b))
    F = sqrt(a)*(a - b)*atan(sqrt(a)*(x**2 + 1)/sqrt(b))/(2*sqrt(b)*(a + b)**2) - 2*a*log(x)/(a + b)**2 + a*log(a*x**4 + 2*a*x**2 + a + b)/(2*(a + b)**2) - 1/(x**2*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_905():
    f = x**4/(a*x**4 + 2*a*x**2 + a + b)
    F = x/a + sqrt(2)*(2*sqrt(a)*sqrt(a + b) + a + b)*atan((-sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(5)/4)*sqrt(sqrt(a) + sqrt(a + b))*sqrt(a + b)) - sqrt(2)*(2*sqrt(a)*sqrt(a + b) + a + b)*atan((sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(5)/4)*sqrt(sqrt(a) + sqrt(a + b))*sqrt(a + b)) + sqrt(2)*(-2*sqrt(a)*sqrt(a + b) + a + b)*log(-sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(5)/4)*sqrt(-sqrt(a) + sqrt(a + b))*sqrt(a + b)) - sqrt(2)*(-2*sqrt(a)*sqrt(a + b) + a + b)*log(sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(5)/4)*sqrt(-sqrt(a) + sqrt(a + b))*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_906():
    f = x**2/(a*x**4 + 2*a*x**2 + a + b)
    F = -sqrt(2)*atan((-sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(3)/4)*sqrt(sqrt(a) + sqrt(a + b))) + sqrt(2)*atan((sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(3)/4)*sqrt(sqrt(a) + sqrt(a + b))) + sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(3)/4)*sqrt(-sqrt(a) + sqrt(a + b))) - sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(3)/4)*sqrt(-sqrt(a) + sqrt(a + b)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_907():
    f = 1/(a*x**4 + 2*a*x**2 + a + b)
    F = -sqrt(2)*atan((-sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(1)/4)*sqrt(sqrt(a) + sqrt(a + b))*sqrt(a + b)) + sqrt(2)*atan((sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*a**(sympy.S(1)/4)*sqrt(sqrt(a) + sqrt(a + b))*sqrt(a + b)) - sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(1)/4)*sqrt(-sqrt(a) + sqrt(a + b))*sqrt(a + b)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*a**(sympy.S(1)/4)*sqrt(-sqrt(a) + sqrt(a + b))*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_908():
    f = 1/(x**2*(a*x**4 + 2*a*x**2 + a + b))
    F = sqrt(2)*a**(sympy.S(1)/4)*(2*sqrt(a) + sqrt(a + b))*atan((-sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*sqrt(sqrt(a) + sqrt(a + b))*(a + b)**(sympy.S(3)/2)) - sqrt(2)*a**(sympy.S(1)/4)*(2*sqrt(a) + sqrt(a + b))*atan((sqrt(2)*a**(sympy.S(1)/4)*x + sqrt(-sqrt(a) + sqrt(a + b)))/sqrt(sqrt(a) + sqrt(a + b)))/(4*sqrt(sqrt(a) + sqrt(a + b))*(a + b)**(sympy.S(3)/2)) + sqrt(2)*a**(sympy.S(1)/4)*(2*sqrt(a) - sqrt(a + b))*log(-sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*sqrt(-sqrt(a) + sqrt(a + b))*(a + b)**(sympy.S(3)/2)) - sqrt(2)*a**(sympy.S(1)/4)*(2*sqrt(a) - sqrt(a + b))*log(sqrt(2)*a**(sympy.S(1)/4)*x*sqrt(-sqrt(a) + sqrt(a + b)) + sqrt(a)*x**2 + sqrt(a + b))/(8*sqrt(-sqrt(a) + sqrt(a + b))*(a + b)**(sympy.S(3)/2)) - 1/(x*(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_909():
    f = x/(x**4 + x**2 + 1)
    F = sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_910():
    f = x/(x**4 + 2*x**2 + 10)
    F = atan(x**2/3 + sympy.S(1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_911():
    f = x**2/(x**4 + 9*x**2 + 20)
    F = -2*atan(x/2) + sqrt(5)*atan(sqrt(5)*x/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_912():
    f = x**2/(x**4 - x**2 + 1)
    F = sqrt(3)*log(x**2 - sqrt(3)*x + 1)/12 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/12 + atan(2*x - sqrt(3))/2 + atan(2*x + sqrt(3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_913():
    f = x**2/(x**4 - 2*x**2 + 2)
    F = log(x**2 - x*sqrt(2 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) - log(x**2 + x*sqrt(2 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*x + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2 + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*x + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_914():
    f = x**7*sqrt(a + b*x**2 + c*x**4)
    F = -b*(b + 2*c*x**2)*(-12*a*c + 7*b**2)*sqrt(a + b*x**2 + c*x**4)/(256*c**4) + b*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(512*c**(sympy.S(9)/2)) + x**4*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(10*c) + (a + b*x**2 + c*x**4)**(sympy.S(3)/2)*(-32*a*c + 35*b**2 - 42*b*c*x**2)/(480*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_915():
    f = x**5*sqrt(a + b*x**2 + c*x**4)
    F = -5*b*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(48*c**2) + x**2*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(8*c) + (b + 2*c*x**2)*(-4*a*c + 5*b**2)*sqrt(a + b*x**2 + c*x**4)/(128*c**3) - (-4*a*c + b**2)*(-4*a*c + 5*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(256*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_916():
    f = x**3*sqrt(a + b*x**2 + c*x**4)
    F = -b*(b + 2*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(16*c**2) + b*(-4*a*c + b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(5)/2)) + (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_917():
    f = x*sqrt(a + b*x**2 + c*x**4)
    F = (b + 2*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(8*c) - (-4*a*c + b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_918():
    f = sqrt(a + b*x**2 + c*x**4)/x
    F = -sqrt(a)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/2 + b*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)) + sqrt(a + b*x**2 + c*x**4)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_919():
    f = sqrt(a + b*x**2 + c*x**4)/x**3
    F = sqrt(c)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/2 - sqrt(a + b*x**2 + c*x**4)/(2*x**2) - b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_920():
    f = sqrt(a + b*x**2 + c*x**4)/x**5
    F = -(2*a + b*x**2)*sqrt(a + b*x**2 + c*x**4)/(8*a*x**4) + (-4*a*c + b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(16*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_921():
    f = sqrt(a + b*x**2 + c*x**4)/x**7
    F = -(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*a*x**6) + b*(2*a + b*x**2)*sqrt(a + b*x**2 + c*x**4)/(16*a**2*x**4) - b*(-4*a*c + b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_922():
    f = sqrt(a + b*x**2 + c*x**4)/x**9
    F = -(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(8*a*x**8) + 5*b*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(48*a**2*x**6) - (2*a + b*x**2)*(-4*a*c + 5*b**2)*sqrt(a + b*x**2 + c*x**4)/(128*a**3*x**4) + (-4*a*c + b**2)*(-4*a*c + 5*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(256*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_923():
    f = sqrt(a + b*x**2 + c*x**4)/x**11
    F = -(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(10*a*x**10) + 7*b*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(80*a**2*x**8) - (-32*a*c + 35*b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(480*a**3*x**6) + b*(2*a + b*x**2)*(-12*a*c + 7*b**2)*sqrt(a + b*x**2 + c*x**4)/(256*a**4*x**4) - b*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(512*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_924():
    f = x**4*sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-29*a*c + 8*b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(105*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c)*(-5*a*c + 2*b**2) - 29*a*b*c + 8*b**3)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(210*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + b*x*(-29*a*c + 8*b**2)*sqrt(a + b*x**2 + c*x**4)/(105*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2)) + x**3*(b + 5*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(35*c) - x*(-10*a*c + 4*b**2)*sqrt(a + b*x**2 + c*x**4)/(105*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_925():
    f = x**2*sqrt(a + b*x**2 + c*x**4)
    F = 2*a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) + x*(b + 3*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(15*c) - x*(-6*a*c + 2*b**2)*sqrt(a + b*x**2 + c*x**4)/(15*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_926():
    f = sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + b*x*sqrt(a + b*x**2 + c*x**4)/(3*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + x*sqrt(a + b*x**2 + c*x**4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_927():
    f = sqrt(a + b*x**2 + c*x**4)/x**2
    F = -2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/sqrt(a + b*x**2 + c*x**4) + 2*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2) - sqrt(a + b*x**2 + c*x**4)/x + sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_928():
    f = sqrt(a + b*x**2 + c*x**4)/x**4
    F = -sqrt(a + b*x**2 + c*x**4)/(3*x**3) + b*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(3*a*(sqrt(a) + sqrt(c)*x**2)) - b*sqrt(a + b*x**2 + c*x**4)/(3*a*x) - b*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_929():
    f = sqrt(a + b*x**2 + c*x**4)/x**6
    F = -sqrt(a + b*x**2 + c*x**4)/(5*x**5) - b*sqrt(a + b*x**2 + c*x**4)/(15*a*x**3) - 2*sqrt(c)*x*(-3*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(15*a**2*(sqrt(a) + sqrt(c)*x**2)) + (-6*a*c + 2*b**2)*sqrt(a + b*x**2 + c*x**4)/(15*a**2*x) + 2*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_930():
    f = x**7*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -b*(b + 2*c*x**2)*(-4*a*c + 3*b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(256*c**4) + 3*b*(b + 2*c*x**2)*(-4*a*c + b**2)*(-4*a*c + 3*b**2)*sqrt(a + b*x**2 + c*x**4)/(2048*c**5) - 3*b*(-4*a*c + b**2)**2*(-4*a*c + 3*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4096*c**(sympy.S(11)/2)) + x**4*(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(14*c) + (a + b*x**2 + c*x**4)**(sympy.S(5)/2)*(-16*a*c + 21*b**2 - 30*b*c*x**2)/(560*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_931():
    f = x**5*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -7*b*(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(120*c**2) + x**2*(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(12*c) + (b + 2*c*x**2)*(-4*a*c + 7*b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(384*c**3) - (b + 2*c*x**2)*(-4*a*c + b**2)*(-4*a*c + 7*b**2)*sqrt(a + b*x**2 + c*x**4)/(1024*c**4) + (-4*a*c + b**2)**2*(-4*a*c + 7*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2048*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_932():
    f = x**3*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -b*(b + 2*c*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(32*c**2) + 3*b*(b + 2*c*x**2)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(256*c**3) - 3*b*(-4*a*c + b**2)**2*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(512*c**(sympy.S(7)/2)) + (a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_933():
    f = x*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = (b + 2*c*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(16*c) - (b + 2*c*x**2)*(-12*a*c + 3*b**2)*sqrt(a + b*x**2 + c*x**4)/(128*c**2) + 3*(-4*a*c + b**2)**2*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(256*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_934():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x
    F = -a**(sympy.S(3)/2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/2 - b*(-12*a*c + b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(3)/2)) + (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/6 + sqrt(a + b*x**2 + c*x**4)*(8*a*c + b**2 + 2*b*c*x**2)/(16*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_935():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**3
    F = -3*sqrt(a)*b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/4 + (9*b/8 + 3*c*x**2/4)*sqrt(a + b*x**2 + c*x**4) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(2*x**2) + (12*a*c + 3*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_936():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**5
    F = 3*b*sqrt(c)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/4 - (3*b - 6*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(8*x**2) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*x**4) - (12*a*c + 3*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(16*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_937():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**7
    F = c**(sympy.S(3)/2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/2 - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*x**6) - (2*a*b + x**2*(8*a*c + b**2))*sqrt(a + b*x**2 + c*x**4)/(16*a*x**4) + b*(-12*a*c + b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(32*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_938():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**9
    F = -(2*a + b*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(16*a*x**8) + (2*a + b*x**2)*(-12*a*c + 3*b**2)*sqrt(a + b*x**2 + c*x**4)/(128*a**2*x**4) - 3*(-4*a*c + b**2)**2*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(256*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_939():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**11
    F = -(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*a*x**10) + b*(2*a + b*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(32*a**2*x**8) - 3*b*(2*a + b*x**2)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(256*a**3*x**4) + 3*b*(-4*a*c + b**2)**2*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(512*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_940():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**13
    F = -(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(12*a*x**12) + 7*b*(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(120*a**2*x**10) - (2*a + b*x**2)*(-4*a*c + 7*b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(384*a**3*x**8) + (2*a + b*x**2)*(-4*a*c + b**2)*(-4*a*c + 7*b**2)*sqrt(a + b*x**2 + c*x**4)/(1024*a**4*x**4) - (-4*a*c + b**2)**2*(-4*a*c + 7*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2048*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_941():
    f = x**4*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 8*a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-9*a*c + 2*b**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(1155*c**(sympy.S(15)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c)*(60*a**2*c**2 - 51*a*b**2*c + 8*b**4) + 8*b*(-9*a*c + 2*b**2)*(-3*a*c + b**2))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2310*c**(sympy.S(15)/4)*sqrt(a + b*x**2 + c*x**4)) - 8*b*x*(-9*a*c + 2*b**2)*(-3*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(1155*c**(sympy.S(7)/2)*(sqrt(a) + sqrt(c)*x**2)) + x**3*(b + 3*c*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(33*c) - x**3*(b*(a*c + 2*b**2) + 10*c*x**2*(-3*a*c + b**2))*sqrt(a + b*x**2 + c*x**4)/(385*c**2) + x*sqrt(a + b*x**2 + c*x**4)*(60*a**2*c**2 - 51*a*b**2*c + 8*b**4)/(1155*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_942():
    f = x**2*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(84*a**2*c**2 - 57*a*b**2*c + 8*b**4)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(315*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(4*sqrt(a)*b*sqrt(c)*(-6*a*c + b**2) + 84*a**2*c**2 - 57*a*b**2*c + 8*b**4)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(630*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + x*(3*b + 7*c*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(63*c) - x*(b*(-9*a*c + 4*b**2) + 6*c*x**2*(-7*a*c + 2*b**2))*sqrt(a + b*x**2 + c*x**4)/(315*c**2) + x*sqrt(a + b*x**2 + c*x**4)*(84*a**2*c**2 - 57*a*b**2*c + 8*b**4)/(315*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_943():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-8*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(35*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c)*(-20*a*c + b**2) + 2*b*(-8*a*c + b**2))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(70*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - 2*b*x*(-8*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(35*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)) + x*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/7 + x*sqrt(a + b*x**2 + c*x**4)*(10*a*c + b**2 + 3*b*c*x**2)/(35*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_944():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**2
    F = -a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(12*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(5*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(8*sqrt(a)*b*sqrt(c) + 12*a*c + b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(10*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + x*(7*b + 6*c*x**2)*sqrt(a + b*x**2 + c*x**4)/5 - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x + x*(12*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(5*sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_945():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**4
    F = -8*a**(sympy.S(1)/4)*b*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*sqrt(a + b*x**2 + c*x**4)) + 8*b*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(3*sqrt(a) + 3*sqrt(c)*x**2) - (3*b - 2*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(3*x) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**3) + sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(8*sqrt(a)*b*sqrt(c) + 4*a*c + 3*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_946():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**6
    F = -(b - 6*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(5*x**3) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*x**5) + sqrt(c)*x*(12*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(5*a*(sqrt(a) + sqrt(c)*x**2)) - (12*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(5*a*x) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(12*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(5*a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(8*sqrt(a)*b*sqrt(c) + 12*a*c + b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(10*a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_947():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/x**8
    F = -(3*b + 30*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(35*x**5) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*x**7) - (-20*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(35*a*x**3) - 2*b*sqrt(c)*x*(-8*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(35*a**2*(sqrt(a) + sqrt(c)*x**2)) + 2*b*(-8*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(35*a**2*x) + 2*b*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-8*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(35*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c)*(-20*a*c + b**2) + 2*b*(-8*a*c + b**2))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(70*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_948():
    f = sqrt(-x**4 - 2*x**2 + 3)
    F = x*sqrt(-x**4 - 2*x**2 + 3)/3 - 2*sqrt(3)*elliptic_e(asin(x), sympy.S(-1)/3)/3 + 4*sqrt(3)*elliptic_f(asin(x), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_949():
    f = x**7/sqrt(a + b*x**2 + c*x**4)
    F = -b*(-12*a*c + 5*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(7)/2)) + x**4*sqrt(a + b*x**2 + c*x**4)/(6*c) + sqrt(a + b*x**2 + c*x**4)*(-16*a*c + 15*b**2 - 10*b*c*x**2)/(48*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_950():
    f = x**5/sqrt(a + b*x**2 + c*x**4)
    F = -3*b*sqrt(a + b*x**2 + c*x**4)/(8*c**2) + x**2*sqrt(a + b*x**2 + c*x**4)/(4*c) + (-4*a*c + 3*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_951():
    f = x**3/sqrt(a + b*x**2 + c*x**4)
    F = -b*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*c**(sympy.S(3)/2)) + sqrt(a + b*x**2 + c*x**4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_952():
    f = x/sqrt(a + b*x**2 + c*x**4)
    F = atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_953():
    f = 1/(x*sqrt(a + b*x**2 + c*x**4))
    F = -atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_954():
    f = 1/(x**3*sqrt(a + b*x**2 + c*x**4))
    F = -sqrt(a + b*x**2 + c*x**4)/(2*a*x**2) + b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_955():
    f = 1/(x**5*sqrt(a + b*x**2 + c*x**4))
    F = -sqrt(a + b*x**2 + c*x**4)/(4*a*x**4) + 3*b*sqrt(a + b*x**2 + c*x**4)/(8*a**2*x**2) - (-4*a*c + 3*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_956():
    f = 1/(x**7*sqrt(a + b*x**2 + c*x**4))
    F = -sqrt(a + b*x**2 + c*x**4)/(6*a*x**6) + 5*b*sqrt(a + b*x**2 + c*x**4)/(24*a**2*x**4) - (-16*a*c + 15*b**2)*sqrt(a + b*x**2 + c*x**4)/(48*a**3*x**2) + b*(-12*a*c + 5*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(32*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_957():
    f = x**4/sqrt(a + b*x**2 + c*x**4)
    F = 2*a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c) + 2*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - 2*b*x*sqrt(a + b*x**2 + c*x**4)/(3*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)) + x*sqrt(a + b*x**2 + c*x**4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_958():
    f = x**2/sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + x*sqrt(a + b*x**2 + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_959():
    f = 1/sqrt(a + b*x**2 + c*x**4)
    F = sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_960():
    f = 1/(x**2*sqrt(a + b*x**2 + c*x**4))
    F = sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)) - sqrt(a + b*x**2 + c*x**4)/(a*x) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_961():
    f = 1/(x**4*sqrt(a + b*x**2 + c*x**4))
    F = -sqrt(a + b*x**2 + c*x**4)/(3*a*x**3) - 2*b*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(3*a**2*(sqrt(a) + sqrt(c)*x**2)) + 2*b*sqrt(a + b*x**2 + c*x**4)/(3*a**2*x) + 2*b*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c) + 2*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_962():
    f = x**7/sqrt(a + b*x**2 - c*x**4)
    F = -b*(12*a*c + 5*b**2)*atan((b - 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 - c*x**4)))/(32*c**(sympy.S(7)/2)) - x**4*sqrt(a + b*x**2 - c*x**4)/(6*c) - sqrt(a + b*x**2 - c*x**4)*(16*a*c + 15*b**2 + 10*b*c*x**2)/(48*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_963():
    f = x**5/sqrt(a + b*x**2 - c*x**4)
    F = -3*b*sqrt(a + b*x**2 - c*x**4)/(8*c**2) - x**2*sqrt(a + b*x**2 - c*x**4)/(4*c) - (4*a*c + 3*b**2)*atan((b - 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 - c*x**4)))/(16*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_964():
    f = x**3/sqrt(a + b*x**2 - c*x**4)
    F = -b*atan((b - 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 - c*x**4)))/(4*c**(sympy.S(3)/2)) - sqrt(a + b*x**2 - c*x**4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_965():
    f = x/sqrt(a + b*x**2 - c*x**4)
    F = -atan((b - 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 - c*x**4)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_966():
    f = 1/(x*sqrt(-a + b*x**2 + c*x**4))
    F = -atan((2*a - b*x**2)/(2*sqrt(a)*sqrt(-a + b*x**2 + c*x**4)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_967():
    f = 1/(x**3*sqrt(-a + b*x**2 + c*x**4))
    F = sqrt(-a + b*x**2 + c*x**4)/(2*a*x**2) - b*atan((2*a - b*x**2)/(2*sqrt(a)*sqrt(-a + b*x**2 + c*x**4)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_968():
    f = 1/(x**5*sqrt(-a + b*x**2 + c*x**4))
    F = sqrt(-a + b*x**2 + c*x**4)/(4*a*x**4) + 3*b*sqrt(-a + b*x**2 + c*x**4)/(8*a**2*x**2) - (4*a*c + 3*b**2)*atan((2*a - b*x**2)/(2*sqrt(a)*sqrt(-a + b*x**2 + c*x**4)))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_969():
    f = 1/(x**7*sqrt(-a + b*x**2 + c*x**4))
    F = sqrt(-a + b*x**2 + c*x**4)/(6*a*x**6) + 5*b*sqrt(-a + b*x**2 + c*x**4)/(24*a**2*x**4) + (16*a*c + 15*b**2)*sqrt(-a + b*x**2 + c*x**4)/(48*a**3*x**2) - b*(12*a*c + 5*b**2)*atan((2*a - b*x**2)/(2*sqrt(a)*sqrt(-a + b*x**2 + c*x**4)))/(32*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_970():
    f = x**4/sqrt(a + b*x**2 - c*x**4)
    F = -sqrt(2)*b*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*c**(sympy.S(5)/2)*sqrt(a + b*x**2 - c*x**4)) - x*sqrt(a + b*x**2 - c*x**4)/(3*c) + sqrt(2)*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*(a*c + b**2 - b*sqrt(4*a*c + b**2))*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*c**(sympy.S(5)/2)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_971():
    f = x**2/sqrt(a + b*x**2 - c*x**4)
    F = -sqrt(2)*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(a + b*x**2 - c*x**4)) + sqrt(2)*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_972():
    f = 1/sqrt(a + b*x**2 - c*x**4)
    F = sqrt(2)*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(2*sqrt(c)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_973():
    f = 1/(x**2*sqrt(a + b*x**2 - c*x**4))
    F = -sqrt(a + b*x**2 - c*x**4)/(a*x) + sqrt(2)*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(a + b*x**2 - c*x**4)) - sqrt(2)*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_974():
    f = 1/(x**4*sqrt(a + b*x**2 - c*x**4))
    F = -sqrt(a + b*x**2 - c*x**4)/(3*a*x**3) + 2*b*sqrt(a + b*x**2 - c*x**4)/(3*a**2*x) - sqrt(2)*b*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*a**2*sqrt(c)*sqrt(a + b*x**2 - c*x**4)) + sqrt(2)*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*(a*c + b**2 - b*sqrt(4*a*c + b**2))*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*a**2*sqrt(c)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_975():
    f = x**9/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -b*x**4*sqrt(a + b*x**2 + c*x**4)/(c*(-4*a*c + b**2)) + x**6*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - (b*(-52*a*c + 15*b**2) - 2*c*x**2*(-12*a*c + 5*b**2))*sqrt(a + b*x**2 + c*x**4)/(8*c**3*(-4*a*c + b**2)) + (-12*a*c + 15*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_976():
    f = x**7/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -3*b*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*c**(sympy.S(5)/2)) + x**4*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + sqrt(a + b*x**2 + c*x**4)*(-8*a*c + 3*b**2 - 2*b*c*x**2)/(2*c**2*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_977():
    f = x**5/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -b*sqrt(a + b*x**2 + c*x**4)/(c*(-4*a*c + b**2)) + x**2*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_978():
    f = x**3/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = (2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_979():
    f = x/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(b + 2*c*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_980():
    f = 1/(x*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_981():
    f = 1/(x**3*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*x**2*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - (-8*a*c + 3*b**2)*sqrt(a + b*x**2 + c*x**4)/(2*a**2*x**2*(-4*a*c + b**2)) + 3*b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_982():
    f = 1/(x**5*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*x**4*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - (-12*a*c + 5*b**2)*sqrt(a + b*x**2 + c*x**4)/(4*a**2*x**4*(-4*a*c + b**2)) + b*(-52*a*c + 15*b**2)*sqrt(a + b*x**2 + c*x**4)/(8*a**3*x**2*(-4*a*c + b**2)) - (-12*a*c + 15*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(16*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_983():
    f = x**6/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - b*x*sqrt(a + b*x**2 + c*x**4)/(c*(-4*a*c + b**2)) + x**3*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + x*(-6*a*c + 2*b**2)*sqrt(a + b*x**2 + c*x**4)/(c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_984():
    f = x**4/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = a**(sympy.S(1)/4)*b*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*(-4*sqrt(a)*sqrt(c) + 2*b)*sqrt(a + b*x**2 + c*x**4)) - b*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)) + x*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_985():
    f = x**2/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + 2*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/((sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)) - x*(b + 2*c*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-2*sqrt(a)*sqrt(c) + b)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_986():
    f = (a + b*x**2 + c*x**4)**(sympy.S(-3)/2)
    F = -b*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)) + x*(-2*a*c + b**2 + b*c*x**2)/(a*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + b*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c) + b)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_987():
    f = 1/(x**2*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*x*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + 2*sqrt(c)*x*(-3*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)/(a**2*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)) - (-6*a*c + 2*b**2)*sqrt(a + b*x**2 + c*x**4)/(a**2*x*(-4*a*c + b**2)) - 2*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_988():
    f = x**4/sqrt(b*x**2 + c*x**4)
    F = -2*b*sqrt(b*x**2 + c*x**4)/(3*c**2*x) + x*sqrt(b*x**2 + c*x**4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_989():
    f = x**3/sqrt(b*x**2 + c*x**4)
    F = -b*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*c**(sympy.S(3)/2)) + sqrt(b*x**2 + c*x**4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_990():
    f = x**2/sqrt(b*x**2 + c*x**4)
    F = sqrt(b*x**2 + c*x**4)/(c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_991():
    f = x/sqrt(b*x**2 + c*x**4)
    F = atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_992():
    f = 1/sqrt(b*x**2 + c*x**4)
    F = -atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_993():
    f = 1/(x*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_994():
    f = 1/(x**2*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(2*b*x**3) + c*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_995():
    f = 1/(x**3*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(3*b*x**4) + 2*c*sqrt(b*x**2 + c*x**4)/(3*b**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_996():
    f = 1/(x**4*sqrt(b*x**2 + c*x**4))
    F = -sqrt(b*x**2 + c*x**4)/(4*b*x**5) + 3*c*sqrt(b*x**2 + c*x**4)/(8*b**2*x**3) - 3*c**2*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_997():
    f = x**4/sqrt(a + c*x**4)
    F = -a**(sympy.S(3)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*c**(sympy.S(5)/4)*sqrt(a + c*x**4)) + x*sqrt(a + c*x**4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_998():
    f = x**3/sqrt(a + c*x**4)
    F = sqrt(a + c*x**4)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_999():
    f = x**2/sqrt(a + c*x**4)
    F = -a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + x*sqrt(a + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1000():
    f = x/sqrt(a + c*x**4)
    F = atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1001():
    f = 1/sqrt(a + c*x**4)
    F = sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1002():
    f = 1/(x*sqrt(a + c*x**4))
    F = -atanh(sqrt(a + c*x**4)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1003():
    f = 1/(x**2*sqrt(a + c*x**4))
    F = sqrt(c)*x*sqrt(a + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)) - sqrt(a + c*x**4)/(a*x) - c**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1004():
    f = 1/(x**3*sqrt(a + c*x**4))
    F = -sqrt(a + c*x**4)/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1005():
    f = 1/(x**4*sqrt(a + c*x**4))
    F = -sqrt(a + c*x**4)/(3*a*x**3) - c**(sympy.S(3)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*a**(sympy.S(5)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1006():
    f = x**4/sqrt(a + b*x**2)
    F = 3*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2)) - 3*a*x*sqrt(a + b*x**2)/(8*b**2) + x**3*sqrt(a + b*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1007():
    f = x**3/sqrt(a + b*x**2)
    F = -a*sqrt(a + b*x**2)/b**2 + (a + b*x**2)**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1008():
    f = x**2/sqrt(a + b*x**2)
    F = -a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2)) + x*sqrt(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1009():
    f = x/sqrt(a + b*x**2)
    F = sqrt(a + b*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1010():
    f = 1/sqrt(a + b*x**2)
    F = atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1011():
    f = 1/(x*sqrt(a + b*x**2))
    F = -atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1012():
    f = 1/(x**2*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1013():
    f = 1/(x**3*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(2*a*x**2) + b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1014():
    f = 1/(x**4*sqrt(a + b*x**2))
    F = -sqrt(a + b*x**2)/(3*a*x**3) + 2*b*sqrt(a + b*x**2)/(3*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1015():
    f = x**4/sqrt(c*x**4)
    F = x**5/(3*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1016():
    f = x**3/sqrt(c*x**4)
    F = x**4/(2*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1017():
    f = x**2/sqrt(c*x**4)
    F = x**3/sqrt(c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1018():
    f = x/sqrt(c*x**4)
    F = x**2*log(x)/sqrt(c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1019():
    f = 1/sqrt(c*x**4)
    F = -x/sqrt(c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1020():
    f = 1/(x*sqrt(c*x**4))
    F = -1/(2*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1021():
    f = 1/(x**2*sqrt(c*x**4))
    F = -1/(3*x*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1022():
    f = 1/(x**3*sqrt(c*x**4))
    F = -1/(4*x**2*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1023():
    f = 1/(x**4*sqrt(c*x**4))
    F = -1/(5*x**3*sqrt(c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1024():
    f = x**4/sqrt(a)
    F = x**5/(5*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1025():
    f = x**3/sqrt(a)
    F = x**4/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1026():
    f = x**2/sqrt(a)
    F = x**3/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1027():
    f = x/sqrt(a)
    F = x**2/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1028():
    f = 1/sqrt(a)
    F = x/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1029():
    f = 1/(sqrt(a)*x)
    F = log(x)/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1030():
    f = 1/(sqrt(a)*x**2)
    F = -1/(sqrt(a)*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1031():
    f = 1/(sqrt(a)*x**3)
    F = -1/(2*sqrt(a)*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1032():
    f = 1/(sqrt(a)*x**4)
    F = -1/(3*sqrt(a)*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1033():
    f = 1/sqrt(-x**4 - 2*x**2 + 3)
    F = sqrt(3)*elliptic_f(asin(x), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1034():
    f = 1/sqrt(-x**4 + 5*x**2 - 1)
    F = -21**(sympy.S(3)/4)*elliptic_f(acos(sqrt(2)*x/sqrt(sqrt(21) + 5)), sympy.S.Half + 5*sqrt(21)/42)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1035():
    f = x**(sympy.S(5)/2)*(a + b*x**2 + c*x**4)
    F = 2*a*x**(sympy.S(7)/2)/7 + 2*b*x**(sympy.S(11)/2)/11 + 2*c*x**(sympy.S(15)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1036():
    f = x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)
    F = 2*a*x**(sympy.S(5)/2)/5 + 2*b*x**(sympy.S(9)/2)/9 + 2*c*x**(sympy.S(13)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1037():
    f = sqrt(x)*(a + b*x**2 + c*x**4)
    F = 2*a*x**(sympy.S(3)/2)/3 + 2*b*x**(sympy.S(7)/2)/7 + 2*c*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1038():
    f = (a + b*x**2 + c*x**4)/sqrt(x)
    F = 2*a*sqrt(x) + 2*b*x**(sympy.S(5)/2)/5 + 2*c*x**(sympy.S(9)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1039():
    f = (a + b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    F = -2*a/sqrt(x) + 2*b*x**(sympy.S(3)/2)/3 + 2*c*x**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1040():
    f = (a + b*x**2 + c*x**4)/x**(sympy.S(5)/2)
    F = -2*a/(3*x**(sympy.S(3)/2)) + 2*b*sqrt(x) + 2*c*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1041():
    f = (a + b*x**2 + c*x**4)/x**(sympy.S(7)/2)
    F = -2*a/(5*x**(sympy.S(5)/2)) - 2*b/sqrt(x) + 2*c*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1042():
    f = x**(sympy.S(5)/2)*(a + b*x**2 + c*x**4)**2
    F = 2*a**2*x**(sympy.S(7)/2)/7 + 4*a*b*x**(sympy.S(11)/2)/11 + 4*b*c*x**(sympy.S(19)/2)/19 + 2*c**2*x**(sympy.S(23)/2)/23 + x**(sympy.S(15)/2)*(4*a*c/15 + 2*b**2/15)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1043():
    f = x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**2
    F = 2*a**2*x**(sympy.S(5)/2)/5 + 4*a*b*x**(sympy.S(9)/2)/9 + 4*b*c*x**(sympy.S(17)/2)/17 + 2*c**2*x**(sympy.S(21)/2)/21 + x**(sympy.S(13)/2)*(4*a*c/13 + 2*b**2/13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1044():
    f = sqrt(x)*(a + b*x**2 + c*x**4)**2
    F = 2*a**2*x**(sympy.S(3)/2)/3 + 4*a*b*x**(sympy.S(7)/2)/7 + 4*b*c*x**(sympy.S(15)/2)/15 + 2*c**2*x**(sympy.S(19)/2)/19 + x**(sympy.S(11)/2)*(4*a*c/11 + 2*b**2/11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1045():
    f = (a + b*x**2 + c*x**4)**2/sqrt(x)
    F = 2*a**2*sqrt(x) + 4*a*b*x**(sympy.S(5)/2)/5 + 4*b*c*x**(sympy.S(13)/2)/13 + 2*c**2*x**(sympy.S(17)/2)/17 + x**(sympy.S(9)/2)*(4*a*c/9 + 2*b**2/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1046():
    f = (a + b*x**2 + c*x**4)**2/x**(sympy.S(3)/2)
    F = -2*a**2/sqrt(x) + 4*a*b*x**(sympy.S(3)/2)/3 + 4*b*c*x**(sympy.S(11)/2)/11 + 2*c**2*x**(sympy.S(15)/2)/15 + x**(sympy.S(7)/2)*(4*a*c/7 + 2*b**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1047():
    f = (a + b*x**2 + c*x**4)**2/x**(sympy.S(5)/2)
    F = -2*a**2/(3*x**(sympy.S(3)/2)) + 4*a*b*sqrt(x) + 4*b*c*x**(sympy.S(9)/2)/9 + 2*c**2*x**(sympy.S(13)/2)/13 + x**(sympy.S(5)/2)*(4*a*c/5 + 2*b**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1048():
    f = (a + b*x**2 + c*x**4)**2/x**(sympy.S(7)/2)
    F = -2*a**2/(5*x**(sympy.S(5)/2)) - 4*a*b/sqrt(x) + 4*b*c*x**(sympy.S(7)/2)/7 + 2*c**2*x**(sympy.S(11)/2)/11 + x**(sympy.S(3)/2)*(4*a*c/3 + 2*b**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1049():
    f = x**(sympy.S(5)/2)*(a + b*x**2 + c*x**4)**3
    F = 2*a**3*x**(sympy.S(7)/2)/7 + 6*a**2*b*x**(sympy.S(11)/2)/11 + 2*a*x**(sympy.S(15)/2)*(a*c + b**2)/5 + 2*b*c**2*x**(sympy.S(27)/2)/9 + 2*b*x**(sympy.S(19)/2)*(6*a*c + b**2)/19 + 2*c**3*x**(sympy.S(31)/2)/31 + 6*c*x**(sympy.S(23)/2)*(a*c + b**2)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1050():
    f = x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**3
    F = 2*a**3*x**(sympy.S(5)/2)/5 + 2*a**2*b*x**(sympy.S(9)/2)/3 + 6*a*x**(sympy.S(13)/2)*(a*c + b**2)/13 + 6*b*c**2*x**(sympy.S(25)/2)/25 + 2*b*x**(sympy.S(17)/2)*(6*a*c + b**2)/17 + 2*c**3*x**(sympy.S(29)/2)/29 + 2*c*x**(sympy.S(21)/2)*(a*c + b**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1051():
    f = sqrt(x)*(a + b*x**2 + c*x**4)**3
    F = 2*a**3*x**(sympy.S(3)/2)/3 + 6*a**2*b*x**(sympy.S(7)/2)/7 + 6*a*x**(sympy.S(11)/2)*(a*c + b**2)/11 + 6*b*c**2*x**(sympy.S(23)/2)/23 + 2*b*x**(sympy.S(15)/2)*(6*a*c + b**2)/15 + 2*c**3*x**(sympy.S(27)/2)/27 + 6*c*x**(sympy.S(19)/2)*(a*c + b**2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1052():
    f = (a + b*x**2 + c*x**4)**3/sqrt(x)
    F = 2*a**3*sqrt(x) + 6*a**2*b*x**(sympy.S(5)/2)/5 + 2*a*x**(sympy.S(9)/2)*(a*c + b**2)/3 + 2*b*c**2*x**(sympy.S(21)/2)/7 + 2*b*x**(sympy.S(13)/2)*(6*a*c + b**2)/13 + 2*c**3*x**(sympy.S(25)/2)/25 + 6*c*x**(sympy.S(17)/2)*(a*c + b**2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1053():
    f = (a + b*x**2 + c*x**4)**3/x**(sympy.S(3)/2)
    F = -2*a**3/sqrt(x) + 2*a**2*b*x**(sympy.S(3)/2) + 6*a*x**(sympy.S(7)/2)*(a*c + b**2)/7 + 6*b*c**2*x**(sympy.S(19)/2)/19 + 2*b*x**(sympy.S(11)/2)*(6*a*c + b**2)/11 + 2*c**3*x**(sympy.S(23)/2)/23 + 2*c*x**(sympy.S(15)/2)*(a*c + b**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1054():
    f = (a + b*x**2 + c*x**4)**3/x**(sympy.S(5)/2)
    F = -2*a**3/(3*x**(sympy.S(3)/2)) + 6*a**2*b*sqrt(x) + 6*a*x**(sympy.S(5)/2)*(a*c + b**2)/5 + 6*b*c**2*x**(sympy.S(17)/2)/17 + 2*b*x**(sympy.S(9)/2)*(6*a*c + b**2)/9 + 2*c**3*x**(sympy.S(21)/2)/21 + 6*c*x**(sympy.S(13)/2)*(a*c + b**2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1055():
    f = (a + b*x**2 + c*x**4)**3/x**(sympy.S(7)/2)
    F = -2*a**3/(5*x**(sympy.S(5)/2)) - 6*a**2*b/sqrt(x) + 2*a*x**(sympy.S(3)/2)*(a*c + b**2) + 2*b*c**2*x**(sympy.S(15)/2)/5 + 2*b*x**(sympy.S(7)/2)*(6*a*c + b**2)/7 + 2*c**3*x**(sympy.S(19)/2)/19 + 6*c*x**(sympy.S(11)/2)*(a*c + b**2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1056():
    f = x**(sympy.S(9)/2)/(a + b*x**2 + c*x**4)
    F = 2*x**(sympy.S(3)/2)/(3*c) - 2**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1057():
    f = x**(sympy.S(7)/2)/(a + b*x**2 + c*x**4)
    F = 2*sqrt(x)/c + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1058():
    f = x**(sympy.S(5)/2)/(a + b*x**2 + c*x**4)
    F = -2**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1059():
    f = x**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = 2**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1060():
    f = sqrt(x)/(a + b*x**2 + c*x**4)
    F = 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1061():
    f = 1/(sqrt(x)*(a + b*x**2 + c*x**4))
    F = -2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/((-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1062():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4))
    F = -2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2/(a*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1063():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2 + c*x**4))
    F = 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2/(3*a*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1064():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2 + c*x**4))
    F = -2/(5*a*x**(sympy.S(5)/2)) + 2*b/(a**2*sqrt(x)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1065():
    f = x**(sympy.S(13)/2)/(a + b*x**2 + c*x**4)**2
    F = -b*x**(sympy.S(3)/2)/(2*c*(-4*a*c + b**2)) + x**(sympy.S(7)/2)*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - 2**(sympy.S(1)/4)*(-20*a*b*c + 3*b**3 - (-14*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(1)/4)*(-20*a*b*c + 3*b**3 - (-14*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(1)/4)*(-20*a*b*c + 3*b**3 + (-14*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(1)/4)*(-20*a*b*c + 3*b**3 + (-14*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1066():
    f = x**(sympy.S(11)/2)/(a + b*x**2 + c*x**4)**2
    F = -b*sqrt(x)/(2*c*(-4*a*c + b**2)) + x**(sympy.S(5)/2)*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - 2**(sympy.S(3)/4)*(-10*a*c + b**2 - b*(-12*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-10*a*c + b**2 - b*(-12*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-10*a*c + b**2 + b*(-12*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-10*a*c + b**2 + b*(-12*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1067():
    f = x**(sympy.S(9)/2)/(a + b*x**2 + c*x**4)**2
    F = x**(sympy.S(3)/2)*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + 2**(sympy.S(1)/4)*(b - (12*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*(b - (12*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*(12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(1)/4)*(12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1068():
    f = x**(sympy.S(7)/2)/(a + b*x**2 + c*x**4)**2
    F = sqrt(x)*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + 2**(sympy.S(3)/4)*(4*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(3)/4)*(4*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(3)/4)*(4*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(3)/4)*(4*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1069():
    f = x**(sympy.S(5)/2)/(a + b*x**2 + c*x**4)**2
    F = 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(4*b - sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(4*b - sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(4*b + sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(4*b + sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - x**(sympy.S(3)/2)*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1070():
    f = x**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)**2
    F = 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-4*b/sqrt(-4*a*c + b**2) + 3)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-4*b/sqrt(-4*a*c + b**2) + 3)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(4*b/sqrt(-4*a*c + b**2) + 3)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(4*b/sqrt(-4*a*c + b**2) + 3)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)) - sqrt(x)*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1071():
    f = sqrt(x)/(a + b*x**2 + c*x**4)**2
    F = 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b + (-20*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b + (-20*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b - (-20*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b - (-20*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)) + x**(sympy.S(3)/2)*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1072():
    f = 1/(sqrt(x)*(a + b*x**2 + c*x**4)**2)
    F = -2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-28*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-28*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-28*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-28*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(x)*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1073():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*sqrt(x)*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-28*a*b*c + 5*b**3 + (-18*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-28*a*b*c + 5*b**3 + (-18*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-28*a*b*c + 5*b**3 - (-18*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-28*a*b*c + 5*b**3 - (-18*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(8*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-18*a*c + 5*b**2)/(2*a**2*sqrt(x)*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1074():
    f = x**(sympy.S(15)/2)/(a + b*x**2 + c*x**4)**3
    F = x**(sympy.S(9)/2)*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + 3*x**(sympy.S(5)/2)*(8*a*b + x**2*(12*a*c + b**2))/(16*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - sqrt(x)*(36*a*c + 3*b**2)/(16*c*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(-84*a*b*c + 3*b**3 - 3*(-24*a**2*c**2 - 30*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(-84*a*b*c + 3*b**3 - 3*(-24*a**2*c**2 - 30*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(-84*a*b*c + 3*b**3 + 3*(-24*a**2*c**2 - 30*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(-84*a*b*c + 3*b**3 + 3*(-24*a**2*c**2 - 30*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1075():
    f = x**(sympy.S(13)/2)/(a + b*x**2 + c*x**4)**3
    F = x**(sympy.S(7)/2)*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x**(sympy.S(3)/2)*(24*a*b + x**2*(28*a*c + 5*b**2))/(16*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + 2**(sympy.S(1)/4)*(28*a*c + 5*b**2 - (172*a*b*c + 5*b**3)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(1)/4)*(28*a*c + 5*b**2 - (172*a*b*c + 5*b**3)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**2) + 2**(sympy.S(1)/4)*(172*a*b*c + 5*b**3 + sqrt(-4*a*c + b**2)*(28*a*c + 5*b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 2**(sympy.S(1)/4)*(172*a*b*c + 5*b**3 + sqrt(-4*a*c + b**2)*(28*a*c + 5*b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1076():
    f = x**(sympy.S(11)/2)/(a + b*x**2 + c*x**4)**3
    F = x**(sympy.S(5)/2)*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(x)*(24*a*b + x**2*(20*a*c + 7*b**2))/(16*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - 2**(sympy.S(3)/4)*(60*a*c + 21*b**2 - 3*(36*a*b*c + 7*b**3)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(60*a*c + 21*b**2 - 3*(36*a*b*c + 7*b**3)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 2**(sympy.S(3)/4)*(108*a*b*c + 21*b**3 + 3*sqrt(-4*a*c + b**2)*(20*a*c + 7*b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 2**(sympy.S(3)/4)*(108*a*b*c + 21*b**3 + 3*sqrt(-4*a*c + b**2)*(20*a*c + 7*b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1077():
    f = x**(sympy.S(9)/2)/(a + b*x**2 + c*x**4)**3
    F = 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(20*a*c + 11*b**2 - 4*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(20*a*c + 11*b**2 - 4*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(20*a*c + 11*b**2 + 4*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(20*a*c + 11*b**2 + 4*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + x**(sympy.S(3)/2)*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - 3*x**(sympy.S(3)/2)*(-4*a*c + 5*b**2 + 8*b*c*x**2)/(16*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1078():
    f = x**(sympy.S(7)/2)/(a + b*x**2 + c*x**4)**3
    F = -2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(28*a*c + 41*b**2 - 36*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(28*a*c + 41*b**2 - 36*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(28*a*c + 41*b**2 + 36*b*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(28*a*c + 41*b**2 + 36*b*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(32*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + sqrt(x)*(2*a + b*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - sqrt(x)*(-4*a*c + 13*b**2 + 24*b*c*x**2)/(16*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1079():
    f = x**(sympy.S(5)/2)/(a + b*x**2 + c*x**4)**3
    F = -x**(sympy.S(3)/2)*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-68*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(12*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-68*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(12*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(68*a*b*c/sqrt(-4*a*c + b**2) + 12*a*c - b**3/sqrt(-4*a*c + b**2) + b**2)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**2) - 3*2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(68*a*b*c/sqrt(-4*a*c + b**2) + 12*a*c - b**3/sqrt(-4*a*c + b**2) + b**2)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**2) + 3*x**(sympy.S(3)/2)*(b*(4*a*c + b**2) + c*x**2*(12*a*c + b**2))/(16*a*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1080():
    f = x**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)**3
    F = -sqrt(x)*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-68*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(44*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-68*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(44*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(68*a*b*c/sqrt(-4*a*c + b**2) + 44*a*c - b**3/sqrt(-4*a*c + b**2) + b**2)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(68*a*b*c/sqrt(-4*a*c + b**2) + 44*a*c - b**3/sqrt(-4*a*c + b**2) + b**2)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**2) + sqrt(x)*(b*(20*a*c + b**2) + c*x**2*(44*a*c + b**2))/(16*a*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1081():
    f = sqrt(x)/(a + b*x**2 + c*x**4)**3
    F = x**(sympy.S(3)/2)*(-2*a*c + b**2 + b*c*x**2)/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(520*a**2*c**2 - 54*a*b**2*c + 5*b**4 + b*(-44*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(520*a**2*c**2 - 54*a*b**2*c + 5*b**4 + b*(-44*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(520*a**2*c**2 - 54*a*b**2*c + 5*b**4 - b*(-44*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(520*a**2*c**2 - 54*a*b**2*c + 5*b**4 - b*(-44*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + x**(sympy.S(3)/2)*(52*a**2*c**2 - 45*a*b**2*c + 5*b**4 + b*c*x**2*(-44*a*c + 5*b**2))/(16*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1082():
    f = 1/(sqrt(x)*(a + b*x**2 + c*x**4)**3)
    F = sqrt(x)*(-2*a*c + b**2 + b*c*x**2)/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(280*a**2*c**2 - 66*a*b**2*c + 7*b**4 + b*(-52*a*c + 7*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(280*a**2*c**2 - 66*a*b**2*c + 7*b**4 + b*(-52*a*c + 7*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(280*a**2*c**2 - 66*a*b**2*c + 7*b**4 - b*(-52*a*c + 7*b**2)*sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(280*a**2*c**2 - 66*a*b**2*c + 7*b**4 - b*(-52*a*c + 7*b**2)*sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(64*a**2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*(-4*a*c + b**2)**(sympy.S(5)/2)) + sqrt(x)*(60*a**2*c**2 - 55*a*b**2*c + 7*b**4 + b*c*x**2*(-52*a*c + 7*b**2))/(16*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1083():
    f = (d*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)
    F = 2*(d*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1084():
    f = sqrt(d*x)*sqrt(a + b*x**2 + c*x**4)
    F = 2*(d*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1085():
    f = sqrt(a + b*x**2 + c*x**4)/sqrt(d*x)
    F = 2*sqrt(d*x)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(1)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1086():
    f = sqrt(a + b*x**2 + c*x**4)/(d*x)**(sympy.S(3)/2)
    F = -2*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(-1)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(d*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1087():
    f = (d*x)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*a*(d*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1088():
    f = sqrt(d*x)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*a*(d*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1089():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/sqrt(d*x)
    F = 2*a*sqrt(d*x)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(1)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1090():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(-1)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(d*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1091():
    f = (d*x)**(sympy.S(3)/2)/sqrt(a + b*x**2 + c*x**4)
    F = 2*(d*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S.Half, sympy.S.Half, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1092():
    f = sqrt(d*x)/sqrt(a + b*x**2 + c*x**4)
    F = 2*(d*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S.Half, sympy.S.Half, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1093():
    f = 1/(sqrt(d*x)*sqrt(a + b*x**2 + c*x**4))
    F = 2*sqrt(d*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/4, sympy.S.Half, sympy.S.Half, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1094():
    f = 1/((d*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4))
    F = -2*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/4, sympy.S.Half, sympy.S.Half, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*sqrt(d*x)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1095():
    f = (d*x)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*(d*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*a*d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1096():
    f = sqrt(d*x)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*(d*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*a*d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1097():
    f = 1/(sqrt(d*x)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 2*sqrt(d*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*d*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1098():
    f = 1/((d*x)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -2*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*d*sqrt(d*x)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1099():
    f = (d*x)**m*(a + b*x**2 + c*x**4)**3
    F = a**3*(d*x)**(m + 1)/(d*(m + 1)) + 3*a**2*b*(d*x)**(m + 3)/(d**3*(m + 3)) + 3*a*(d*x)**(m + 5)*(a*c + b**2)/(d**5*(m + 5)) + 3*b*c**2*(d*x)**(m + 11)/(d**11*(m + 11)) + b*(d*x)**(m + 7)*(6*a*c + b**2)/(d**7*(m + 7)) + c**3*(d*x)**(m + 13)/(d**13*(m + 13)) + 3*c*(d*x)**(m + 9)*(a*c + b**2)/(d**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1100():
    f = (d*x)**m*(a + b*x**2 + c*x**4)**2
    F = a**2*(d*x)**(m + 1)/(d*(m + 1)) + 2*a*b*(d*x)**(m + 3)/(d**3*(m + 3)) + 2*b*c*(d*x)**(m + 7)/(d**7*(m + 7)) + c**2*(d*x)**(m + 9)/(d**9*(m + 9)) + (d*x)**(m + 5)*(2*a*c + b**2)/(d**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1101():
    f = (d*x)**m*(a + b*x**2 + c*x**4)
    F = a*(d*x)**(m + 1)/(d*(m + 1)) + b*(d*x)**(m + 3)/(d**3*(m + 3)) + c*(d*x)**(m + 5)/(d**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1102():
    f = (d*x)**m/(a + b*x**2 + c*x**4)
    F = -2*c*(d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*(d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(d*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1103():
    f = (d*x)**m/(a + b*x**2 + c*x**4)**2
    F = -c*(d*x)**(m + 1)*(-4*a*c*(3 - m) + b**2*(1 - m) - b*(1 - m)*sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(2*a*d*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + c*(d*x)**(m + 1)*(-4*a*c*(3 - m) + b**2*(1 - m) + b*(1 - m)*sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(2*a*d*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + (d*x)**(m + 1)*(-2*a*c + b**2 + b*c*x**2)/(2*a*d*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1104():
    f = (d*x)**m*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = a*(d*x)**(m + 1)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S.Half, sympy.S(-3)/2, sympy.S(-3)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1105():
    f = (d*x)**m*sqrt(a + b*x**2 + c*x**4)
    F = (d*x)**(m + 1)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S.Half, sympy.S(-1)/2, sympy.S(-1)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1106():
    f = (d*x)**m/sqrt(a + b*x**2 + c*x**4)
    F = (d*x)**(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S.Half, sympy.S.Half, sympy.S.Half, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1107():
    f = (d*x)**m/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S.Half, sympy.S(3)/2, sympy.S(3)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*d*(m + 1)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1108():
    f = (d*x)**m*(a + b*x**2 + c*x**4)**p
    F = (d*x)**(m + 1)*(a + b*x**2 + c*x**4)**p*appellf1(m/2 + sympy.S.Half, -p, -p, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1109():
    f = x**7*(a + b*x**2 + c*x**4)**p
    F = -2**(p - 2)*b*(-(b + 2*c*x**2 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(6*a*c - b**2*(p + 3))*(a + b*x**2 + c*x**4)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(c**3*(p + 1)*(2*p + 3)*sqrt(-4*a*c + b**2)) + x**4*(a + b*x**2 + c*x**4)**(p + 1)/(4*c*(p + 2)) + (a + b*x**2 + c*x**4)**(p + 1)*(-2*a*c*(2*p + 3) + b**2*(p + 2)*(p + 3) - 2*b*c*x**2*(p + 1)*(p + 3))/(8*c**3*(p + 1)*(p + 2)*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1110():
    f = x**5*(a + b*x**2 + c*x**4)**p
    F = 2**(p - 1)*(-(b + 2*c*x**2 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(2*a*c - b**2*(p + 2))*(a + b*x**2 + c*x**4)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(c**2*(p + 1)*(2*p + 3)*sqrt(-4*a*c + b**2)) - b*(p + 2)*(a + b*x**2 + c*x**4)**(p + 1)/(4*c**2*(p + 1)*(2*p + 3)) + x**2*(a + b*x**2 + c*x**4)**(p + 1)/(2*c*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1111():
    f = x**3*(a + b*x**2 + c*x**4)**p
    F = 2**(p - 1)*b*(-(b + 2*c*x**2 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(a + b*x**2 + c*x**4)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(c*(p + 1)*sqrt(-4*a*c + b**2)) + (a + b*x**2 + c*x**4)**(p + 1)/(4*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1112():
    f = x*(a + b*x**2 + c*x**4)**p
    F = -2**p*(-(b + 2*c*x**2 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(a + b*x**2 + c*x**4)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/((p + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1113():
    f = (a + b*x**2 + c*x**4)**p/x
    F = 4**(p - 1)*(a + b*x**2 + c*x**4)**p*appellf1(-2*p, -p, -p, 1 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**2), -(b + sqrt(-4*a*c + b**2))/(2*c*x**2))/(p*((b + 2*c*x**2 - sqrt(-4*a*c + b**2))/(c*x**2))**p*((b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(c*x**2))**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1114():
    f = (a + b*x**2 + c*x**4)**p/x**3
    F = -2**(2*p - 1)*(a + b*x**2 + c*x**4)**p*appellf1(1 - 2*p, -p, -p, 2 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**2), -(b + sqrt(-4*a*c + b**2))/(2*c*x**2))/(x**2*((b + 2*c*x**2 - sqrt(-4*a*c + b**2))/(c*x**2))**p*((b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(c*x**2))**p*(1 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1115():
    f = (a + b*x**2 + c*x**4)**p/x**5
    F = -4**(p - 1)*(a + b*x**2 + c*x**4)**p*appellf1(2 - 2*p, -p, -p, 3 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**2), -(b + sqrt(-4*a*c + b**2))/(2*c*x**2))/(x**4*((b + 2*c*x**2 - sqrt(-4*a*c + b**2))/(c*x**2))**p*((b + 2*c*x**2 + sqrt(-4*a*c + b**2))/(c*x**2))**p*(1 - p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1116():
    f = x**4*(a + b*x**2 + c*x**4)**p
    F = x**5*(a + b*x**2 + c*x**4)**p*appellf1(sympy.S(5)/2, -p, -p, sympy.S(7)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1117():
    f = x**2*(a + b*x**2 + c*x**4)**p
    F = x**3*(a + b*x**2 + c*x**4)**p*appellf1(sympy.S(3)/2, -p, -p, sympy.S(5)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1118():
    f = (a + b*x**2 + c*x**4)**p
    F = x*(a + b*x**2 + c*x**4)**p*appellf1(sympy.S.Half, -p, -p, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/((2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1119():
    f = (a + b*x**2 + c*x**4)**p/x**2
    F = -(a + b*x**2 + c*x**4)**p*appellf1(sympy.S(-1)/2, -p, -p, sympy.S.Half, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(x*(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_2_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1120():
    f = (a + b*x**2 + c*x**4)**p/x**4
    F = -(a + b*x**2 + c*x**4)**p*appellf1(sympy.S(-3)/2, -p, -p, sympy.S(-1)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*x**3*(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F

