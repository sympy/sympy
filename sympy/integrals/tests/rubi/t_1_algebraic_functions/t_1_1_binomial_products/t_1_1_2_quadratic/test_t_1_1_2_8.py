"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.2 Quadratic/1.1.2.8 P(x) (c x)^m (a+b x^2)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, D, F, a, b, c, d, e, f, m = symbols('A B C D F a b c d e f m')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_1():
    f = x**3*(A + B*x)*sqrt(a + b*x**2)
    F = A*x**2*(a + b*x**2)**(sympy.S(3)/2)/(5*b) + B*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(5)/2)) + B*a**2*x*sqrt(a + b*x**2)/(16*b**2) + B*x**3*(a + b*x**2)**(sympy.S(3)/2)/(6*b) - a*(16*A + 15*B*x)*(a + b*x**2)**(sympy.S(3)/2)/(120*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_2():
    f = x**2*(A + B*x)*sqrt(a + b*x**2)
    F = -A*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(3)/2)) - A*a*x*sqrt(a + b*x**2)/(8*b) + B*x**2*(a + b*x**2)**(sympy.S(3)/2)/(5*b) - (a + b*x**2)**(sympy.S(3)/2)*(-15*A*b*x + 8*B*a)/(60*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_3():
    f = x*(A + B*x)*sqrt(a + b*x**2)
    F = -B*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(3)/2)) - B*a*x*sqrt(a + b*x**2)/(8*b) + (4*A + 3*B*x)*(a + b*x**2)**(sympy.S(3)/2)/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_4():
    f = (A + B*x)*sqrt(a + b*x**2)
    F = A*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)) + A*x*sqrt(a + b*x**2)/2 + B*(a + b*x**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_5():
    f = (A + B*x)*sqrt(a + b*x**2)/x
    F = -A*sqrt(a)*atanh(sqrt(a + b*x**2)/sqrt(a)) + B*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)) + (2*A + B*x)*sqrt(a + b*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_6():
    f = (A + B*x)*sqrt(a + b*x**2)/x**2
    F = A*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - B*sqrt(a)*atanh(sqrt(a + b*x**2)/sqrt(a)) - (A - B*x)*sqrt(a + b*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_7():
    f = (A + B*x)*sqrt(a + b*x**2)/x**3
    F = -A*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*sqrt(a)) + B*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - (A + 2*B*x)*sqrt(a + b*x**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_8():
    f = x**3*(A + B*x)*(a + b*x**2)**(sympy.S(3)/2)
    F = A*x**2*(a + b*x**2)**(sympy.S(5)/2)/(7*b) + 3*B*a**4*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(5)/2)) + 3*B*a**3*x*sqrt(a + b*x**2)/(128*b**2) + B*a**2*x*(a + b*x**2)**(sympy.S(3)/2)/(64*b**2) + B*x**3*(a + b*x**2)**(sympy.S(5)/2)/(8*b) - a*(32*A + 35*B*x)*(a + b*x**2)**(sympy.S(5)/2)/(560*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_9():
    f = x**2*(A + B*x)*(a + b*x**2)**(sympy.S(3)/2)
    F = -A*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(3)/2)) - A*a**2*x*sqrt(a + b*x**2)/(16*b) - A*a*x*(a + b*x**2)**(sympy.S(3)/2)/(24*b) + B*x**2*(a + b*x**2)**(sympy.S(5)/2)/(7*b) - (a + b*x**2)**(sympy.S(5)/2)*(-35*A*b*x + 12*B*a)/(210*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_10():
    f = x*(A + B*x)*(a + b*x**2)**(sympy.S(3)/2)
    F = -B*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(3)/2)) - B*a**2*x*sqrt(a + b*x**2)/(16*b) - B*a*x*(a + b*x**2)**(sympy.S(3)/2)/(24*b) + (6*A + 5*B*x)*(a + b*x**2)**(sympy.S(5)/2)/(30*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_11():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(3)/2)
    F = 3*A*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)) + 3*A*a*x*sqrt(a + b*x**2)/8 + A*x*(a + b*x**2)**(sympy.S(3)/2)/4 + B*(a + b*x**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_12():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(3)/2)/x
    F = -A*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + 3*B*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)) + a*(8*A + 3*B*x)*sqrt(a + b*x**2)/8 + (4*A + 3*B*x)*(a + b*x**2)**(sympy.S(3)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_13():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(3)/2)/x**2
    F = 3*A*a*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 - B*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + sqrt(a + b*x**2)*(3*A*b*x + 2*B*a)/2 - (3*A - B*x)*(a + b*x**2)**(sympy.S(3)/2)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_14():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(3)/2)/x**3
    F = -3*A*sqrt(a)*b*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + 3*B*a*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + sqrt(a + b*x**2)*(3*A*b*x - 3*B*a)/(2*x) - (A - B*x)*(a + b*x**2)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_15():
    f = x**3*(A + B*x)*(a + b*x**2)**(sympy.S(5)/2)
    F = A*x**2*(a + b*x**2)**(sympy.S(7)/2)/(9*b) + 3*B*a**5*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(5)/2)) + 3*B*a**4*x*sqrt(a + b*x**2)/(256*b**2) + B*a**3*x*(a + b*x**2)**(sympy.S(3)/2)/(128*b**2) + B*a**2*x*(a + b*x**2)**(sympy.S(5)/2)/(160*b**2) + B*x**3*(a + b*x**2)**(sympy.S(7)/2)/(10*b) - a*(160*A + 189*B*x)*(a + b*x**2)**(sympy.S(7)/2)/(5040*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_16():
    f = x**2*(A + B*x)*(a + b*x**2)**(sympy.S(5)/2)
    F = -5*A*a**4*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(3)/2)) - 5*A*a**3*x*sqrt(a + b*x**2)/(128*b) - 5*A*a**2*x*(a + b*x**2)**(sympy.S(3)/2)/(192*b) - A*a*x*(a + b*x**2)**(sympy.S(5)/2)/(48*b) + B*x**2*(a + b*x**2)**(sympy.S(7)/2)/(9*b) - (a + b*x**2)**(sympy.S(7)/2)*(-63*A*b*x + 16*B*a)/(504*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_17():
    f = x*(A + B*x)*(a + b*x**2)**(sympy.S(5)/2)
    F = -5*B*a**4*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(3)/2)) - 5*B*a**3*x*sqrt(a + b*x**2)/(128*b) - 5*B*a**2*x*(a + b*x**2)**(sympy.S(3)/2)/(192*b) - B*a*x*(a + b*x**2)**(sympy.S(5)/2)/(48*b) + (8*A + 7*B*x)*(a + b*x**2)**(sympy.S(7)/2)/(56*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_18():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(5)/2)
    F = 5*A*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*sqrt(b)) + 5*A*a**2*x*sqrt(a + b*x**2)/16 + 5*A*a*x*(a + b*x**2)**(sympy.S(3)/2)/24 + A*x*(a + b*x**2)**(sympy.S(5)/2)/6 + B*(a + b*x**2)**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_19():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(5)/2)/x
    F = -A*a**(sympy.S(5)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + 5*B*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*sqrt(b)) + a**2*(16*A + 5*B*x)*sqrt(a + b*x**2)/16 + a*(8*A + 5*B*x)*(a + b*x**2)**(sympy.S(3)/2)/24 + (6*A + 5*B*x)*(a + b*x**2)**(sympy.S(5)/2)/30
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_20():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(5)/2)/x**2
    F = 15*A*a**2*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/8 - B*a**(sympy.S(5)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + a*sqrt(a + b*x**2)*(15*A*b*x + 8*B*a)/8 + (a + b*x**2)**(sympy.S(3)/2)*(15*A*b*x + 4*B*a)/12 - (5*A - B*x)*(a + b*x**2)**(sympy.S(5)/2)/(5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_21():
    f = (A + B*x)*(a + b*x**2)**(sympy.S(5)/2)/x**3
    F = -5*A*a**(sympy.S(3)/2)*b*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + 15*B*a**2*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/8 + 5*a*b*(4*A + 3*B*x)*sqrt(a + b*x**2)/8 - (a + b*x**2)**(sympy.S(3)/2)*(-10*A*b*x + 15*B*a)/(12*x) - (2*A - B*x)*(a + b*x**2)**(sympy.S(5)/2)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_22():
    f = x**3*(A + B*x)/sqrt(a + b*x**2)
    F = A*x**2*sqrt(a + b*x**2)/(3*b) + 3*B*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2)) + B*x**3*sqrt(a + b*x**2)/(4*b) - a*(16*A + 9*B*x)*sqrt(a + b*x**2)/(24*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_23():
    f = x**2*(A + B*x)/sqrt(a + b*x**2)
    F = -A*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2)) + B*x**2*sqrt(a + b*x**2)/(3*b) - sqrt(a + b*x**2)*(-3*A*b*x + 4*B*a)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_24():
    f = x*(A + B*x)/sqrt(a + b*x**2)
    F = -B*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2)) + (2*A + B*x)*sqrt(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_25():
    f = (A + B*x)/sqrt(a + b*x**2)
    F = A*atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b) + B*sqrt(a + b*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_26():
    f = (A + B*x)/(x*sqrt(a + b*x**2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a) + B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_27():
    f = (A + B*x)/(x**2*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(a*x) - B*atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_28():
    f = (A + B*x)/(x**3*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(2*a*x**2) + A*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2)) - B*sqrt(a + b*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_29():
    f = x**3*(A + B*x)/(a + b*x**2)**(sympy.S(3)/2)
    F = -3*B*a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(5)/2)) - x**2*(A + B*x)/(b*sqrt(a + b*x**2)) + (4*A + 3*B*x)*sqrt(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_30():
    f = x**2*(A + B*x)/(a + b*x**2)**(sympy.S(3)/2)
    F = A*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2) + 2*B*sqrt(a + b*x**2)/b**2 - x*(A + B*x)/(b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_31():
    f = x*(A + B*x)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2) - (A + B*x)/(b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_32():
    f = (A + B*x)/(a + b*x**2)**(sympy.S(3)/2)
    F = -(-A*b*x + B*a)/(a*b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_33():
    f = (A + B*x)/(x*(a + b*x**2)**(sympy.S(3)/2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(3)/2) + (A + B*x)/(a*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_34():
    f = (A + B*x)/(x**2*(a + b*x**2)**(sympy.S(3)/2))
    F = -2*A*sqrt(a + b*x**2)/(a**2*x) - B*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(3)/2) + (A + B*x)/(a*x*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_35():
    f = (A + B*x)/(x**3*(a + b*x**2)**(sympy.S(3)/2))
    F = -3*A*sqrt(a + b*x**2)/(2*a**2*x**2) + 3*A*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(5)/2)) - 2*B*sqrt(a + b*x**2)/(a**2*x) + (A + B*x)/(a*x**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_36():
    f = x**3*(A + B*x)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(5)/2) - x**2*(A + B*x)/(3*b*(a + b*x**2)**(sympy.S(3)/2)) - (2*A + 3*B*x)/(3*b**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_37():
    f = x**2*(A + B*x)/(a + b*x**2)**(sympy.S(5)/2)
    F = -2*B/(3*b**2*sqrt(a + b*x**2)) - x**2*(-A*b*x + B*a)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_38():
    f = x*(A + B*x)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*x/(3*a*b*sqrt(a + b*x**2)) + (-A - B*x)/(3*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_39():
    f = (A + B*x)/(a + b*x**2)**(sympy.S(5)/2)
    F = 2*A*x/(3*a**2*sqrt(a + b*x**2)) + (A*b*x - B*a)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_40():
    f = (A + B*x)/(x*(a + b*x**2)**(sympy.S(5)/2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(5)/2) + (A + B*x)/(3*a*(a + b*x**2)**(sympy.S(3)/2)) + (3*A + 2*B*x)/(3*a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_41():
    f = (A + B*x)/(x**2*(a + b*x**2)**(sympy.S(5)/2))
    F = -8*A*sqrt(a + b*x**2)/(3*a**3*x) - B*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(5)/2) + (A + B*x)/(3*a*x*(a + b*x**2)**(sympy.S(3)/2)) + (4*A + 3*B*x)/(3*a**2*x*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_42():
    f = (A + B*x)/(x**3*(a + b*x**2)**(sympy.S(5)/2))
    F = -5*A*sqrt(a + b*x**2)/(2*a**3*x**2) + 5*A*b*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(7)/2)) - 8*B*sqrt(a + b*x**2)/(3*a**3*x) + (A + B*x)/(3*a*x**2*(a + b*x**2)**(sympy.S(3)/2)) + (5*A + 4*B*x)/(3*a**2*x**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_43():
    f = x*(1 - x)/sqrt(1 - x**2)
    F = sqrt(1 - x**2)*(x/2 - 1) - asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_44():
    f = (-x**2 + x)/sqrt(1 - x**2)
    F = sqrt(1 - x**2)*(x/2 - 1) - asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_45():
    f = (x**2 + 3)/(x**2 - 3)
    F = x - 2*sqrt(3)*atanh(sqrt(3)*x/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_46():
    f = (x**2 - 1)/(x**2 + 1)
    F = x - 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_47():
    f = x**7*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(9)/2) - x**7*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - x**5*(7*B*a - x*(A*b - 8*C*a))/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2)) - x**3*(35*B*a - x*(6*A*b - 48*C*a))/(105*a*b**3*(a + b*x**2)**(sympy.S(3)/2)) - x*(35*B*a - x*(8*A*b - 64*C*a))/(35*a*b**4*sqrt(a + b*x**2)) - sqrt(a + b*x**2)*(16*A*b - 128*C*a)/(35*a*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_48():
    f = x**6*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = C*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(9)/2) - x**4*(6*B + 7*C*x)/(35*b**2*(a + b*x**2)**(sympy.S(5)/2)) - x**2*(24*B + 35*C*x)/(105*b**3*(a + b*x**2)**(sympy.S(3)/2)) - (16*B + 35*C*x)/(35*b**4*sqrt(a + b*x**2)) - x**6*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_49():
    f = x**5*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = (4*A*b + 24*C*a)/(105*b**4*(a + b*x**2)**(sympy.S(3)/2)) - x**5*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - x**4*(A*b - 5*B*b*x + 6*C*a)/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2)) - (4*A*b + 24*C*a)/(35*a*b**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_50():
    f = x**4*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = -x**4*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - x**2*(4*B*a + x*(2*A*b + 5*C*a))/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2)) - (8*B*a + x*(6*A*b + 15*C*a))/(105*a*b**3*(a + b*x**2)**(sympy.S(3)/2)) + x*(2*A*b + 5*C*a)/(35*a**2*b**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_51():
    f = x**3*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = 2*B*x/(35*a**2*b**2*sqrt(a + b*x**2)) - x**3*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - x*(3*B*a + x*(3*A*b + 4*C*a))/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2)) - (6*A*b - 3*B*b*x + 8*C*a)/(105*a*b**3*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_52():
    f = x**2*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = -x**2*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - (2*B*a + x*(4*A*b + 3*C*a))/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2)) + x*(4*A*b + 3*C*a)/(105*a**2*b**2*(a + b*x**2)**(sympy.S(3)/2)) + x*(8*A*b + 6*C*a)/(105*a**3*b**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_53():
    f = x*(A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = 4*B*x/(105*a**2*b*(a + b*x**2)**(sympy.S(3)/2)) + 8*B*x/(105*a**3*b*sqrt(a + b*x**2)) - x*(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) - (5*A*b - B*b*x + 2*C*a)/(35*a*b**2*(a + b*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_54():
    f = (A + B*x + C*x**2)/(a + b*x**2)**(sympy.S(9)/2)
    F = -(B*a - x*(A*b - C*a))/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) + x*(6*A*b + C*a)/(35*a**2*b*(a + b*x**2)**(sympy.S(5)/2)) + x*(24*A*b + 4*C*a)/(105*a**3*b*(a + b*x**2)**(sympy.S(3)/2)) + x*(48*A*b + 8*C*a)/(105*a**4*b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_55():
    f = (A + B*x + C*x**2)/(x*(a + b*x**2)**(sympy.S(9)/2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(9)/2) + (A*b + B*b*x - C*a)/(7*a*b*(a + b*x**2)**(sympy.S(7)/2)) + (7*A + 6*B*x)/(35*a**2*(a + b*x**2)**(sympy.S(5)/2)) + (35*A + 24*B*x)/(105*a**3*(a + b*x**2)**(sympy.S(3)/2)) + (35*A + 16*B*x)/(35*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_56():
    f = (A + B*x + C*x**2)/(x**2*(a + b*x**2)**(sympy.S(9)/2))
    F = -A*sqrt(a + b*x**2)/(a**5*x) - B*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(9)/2) + (B - x*(A*b/a - C))/(7*a*(a + b*x**2)**(sympy.S(7)/2)) + (7*B - x*(13*A*b/a - 6*C))/(35*a**2*(a + b*x**2)**(sympy.S(5)/2)) + (35*B - x*(87*A*b/a - 24*C))/(105*a**3*(a + b*x**2)**(sympy.S(3)/2)) + (35*B - x*(93*A*b/a - 16*C))/(35*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_57():
    f = (A + B*x + C*x**2)/(x**3*(a + b*x**2)**(sympy.S(9)/2))
    F = -A*sqrt(a + b*x**2)/(2*a**5*x**2) - B*sqrt(a + b*x**2)/(a**5*x) - (B*b*x + a*(A*b/a - C))/(7*a**2*(a + b*x**2)**(sympy.S(7)/2)) - (14*A*b + 13*B*b*x - 7*C*a)/(35*a**3*(a + b*x**2)**(sympy.S(5)/2)) - (105*A*b + 87*B*b*x - 35*C*a)/(105*a**4*(a + b*x**2)**(sympy.S(3)/2)) - (140*A*b + 93*B*b*x - 35*C*a)/(35*a**5*sqrt(a + b*x**2)) + (9*A*b - 2*C*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_58():
    f = A*(c*x)**m/(a + b*x**2)
    F = A*(c*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_59():
    f = (c*x)**m*(A + B*x)/(a + b*x**2)
    F = A*(c*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*c*(m + 1)) + B*(c*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -b*x**2/a)/(a*c**2*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_60():
    f = (c*x)**m*(A + C*x**2)/(a + b*x**2)
    F = C*(c*x)**(m + 1)/(b*c*(m + 1)) + (c*x)**(m + 1)*(A*b - C*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_61():
    f = (c*x)**m*(A + B*x + C*x**2)/(a + b*x**2)
    F = B*(c*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -b*x**2/a)/(a*c**2*(m + 2)) + C*(c*x)**(m + 1)/(b*c*(m + 1)) + (c*x)**(m + 1)*(A*b - C*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_62():
    f = x**3*(a + b*x**2)*(A + B*x + C*x**2 + D*x**3)
    F = A*a*x**4/4 + B*a*x**5/5 + C*b*x**8/8 + D*b*x**9/9 + x**7*(B*b + D*a)/7 + x**6*(A*b + C*a)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_63():
    f = x**2*(a + b*x**2)*(A + B*x + C*x**2 + D*x**3)
    F = A*a*x**3/3 + B*a*x**4/4 + C*b*x**7/7 + D*b*x**8/8 + x**6*(B*b + D*a)/6 + x**5*(A*b + C*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_64():
    f = x*(a + b*x**2)*(A + B*x + C*x**2 + D*x**3)
    F = A*a*x**2/2 + B*a*x**3/3 + C*b*x**6/6 + D*b*x**7/7 + x**5*(B*b + D*a)/5 + x**4*(A*b + C*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_65():
    f = (a + b*x**2)*(A + B*x + C*x**2 + D*x**3)
    F = A*a*x + B*a*x**2/2 + C*b*x**5/5 + D*b*x**6/6 + x**4*(B*b + D*a)/4 + x**3*(A*b + C*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_66():
    f = (a + b*x**2)*(A + B*x + C*x**2 + D*x**3)/x
    F = A*a*log(x) + B*a*x + C*b*x**4/4 + D*b*x**5/5 + x**3*(B*b + D*a)/3 + x**2*(A*b + C*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_67():
    f = (a + b*x**2)*(A + B*x + C*x**2 + D*x**3)/x**2
    F = -A*a/x + B*a*log(x) + C*b*x**3/3 + D*b*x**4/4 + x**2*(B*b + D*a)/2 + x*(A*b + C*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_68():
    f = (a + b*x**2)*(A + B*x + C*x**2 + D*x**3)/x**3
    F = -A*a/(2*x**2) - B*a/x + C*b*x**2/2 + D*b*x**3/3 + x*(B*b + D*a) + (A*b + C*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_69():
    f = (a + b*x**2)*(A + B*x + C*x**2 + D*x**3)/x**4
    F = -A*a/(3*x**3) - B*a/(2*x**2) + C*b*x + D*b*x**2/2 + (B*b + D*a)*log(x) - (A*b + C*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_70():
    f = x**3*(a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)
    F = A*a**2*x**4/4 + B*a**2*x**5/5 + C*b**2*x**10/10 + D*b**2*x**11/11 + a*x**7*(2*B*b + D*a)/7 + a*x**6*(2*A*b + C*a)/6 + b*x**9*(B*b + 2*D*a)/9 + b*x**8*(A*b + 2*C*a)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_71():
    f = x**2*(a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)
    F = A*a**2*x**3/3 + B*a**2*x**4/4 + C*b**2*x**9/9 + D*b**2*x**10/10 + a*x**6*(2*B*b + D*a)/6 + a*x**5*(2*A*b + C*a)/5 + b*x**8*(B*b + 2*D*a)/8 + b*x**7*(A*b + 2*C*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_72():
    f = x*(a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)
    F = A*(a + b*x**2)**3/(6*b) + B*a**2*x**3/3 + C*a**2*x**4/4 + C*a*b*x**6/3 + C*b**2*x**8/8 + D*b**2*x**9/9 + a*x**5*(2*B*b + D*a)/5 + b*x**7*(B*b + 2*D*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_73():
    f = (a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)
    F = A*a**2*x + B*(a + b*x**2)**3/(6*b) + C*b**2*x**7/7 + D*a**2*x**4/4 + D*a*b*x**6/3 + D*b**2*x**8/8 + a*x**3*(2*A*b + C*a)/3 + b*x**5*(A*b + 2*C*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_74():
    f = (a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)/x
    F = A*a**2*log(x) + A*a*b*x**2 + A*b**2*x**4/4 + B*a**2*x + C*(a + b*x**2)**3/(6*b) + D*b**2*x**7/7 + a*x**3*(2*B*b + D*a)/3 + b*x**5*(B*b + 2*D*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_75():
    f = (a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)/x**2
    F = -A*a**2/x + B*a**2*log(x) + B*a*b*x**2 + B*b**2*x**4/4 + C*b**2*x**5/5 + D*(a + b*x**2)**3/(6*b) + a*x*(2*A*b + C*a) + b*x**3*(A*b + 2*C*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_76():
    f = (a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)/x**3
    F = -A*a**2/(2*x**2) - B*a**2/x + C*b**2*x**4/4 + D*b**2*x**5/5 + a*x*(2*B*b + D*a) + a*(2*A*b + C*a)*log(x) + b*x**3*(B*b + 2*D*a)/3 + b*x**2*(A*b + 2*C*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_77():
    f = (a + b*x**2)**2*(A + B*x + C*x**2 + D*x**3)/x**4
    F = -A*a**2/(3*x**3) - B*a**2/(2*x**2) + C*b**2*x**3/3 + D*b**2*x**4/4 + a*(2*B*b + D*a)*log(x) - a*(2*A*b + C*a)/x + b*x**2*(B*b + 2*D*a)/2 + b*x*(A*b + 2*C*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_78():
    f = x**3*(a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)
    F = A*a**3*x**4/4 + B*a**3*x**5/5 + C*b**3*x**12/12 + D*b**3*x**13/13 + a**2*x**7*(3*B*b + D*a)/7 + a**2*x**6*(3*A*b + C*a)/6 + a*b*x**9*(B*b + D*a)/3 + 3*a*b*x**8*(A*b + C*a)/8 + b**2*x**11*(B*b + 3*D*a)/11 + b**2*x**10*(A*b + 3*C*a)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_79():
    f = x**2*(a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)
    F = A*a**3*x**3/3 + B*a**3*x**4/4 + C*b**3*x**11/11 + D*b**3*x**12/12 + a**2*x**6*(3*B*b + D*a)/6 + a**2*x**5*(3*A*b + C*a)/5 + 3*a*b*x**8*(B*b + D*a)/8 + 3*a*b*x**7*(A*b + C*a)/7 + b**2*x**10*(B*b + 3*D*a)/10 + b**2*x**9*(A*b + 3*C*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_80():
    f = x*(a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)
    F = A*(a + b*x**2)**4/(8*b) + B*a**3*x**3/3 + C*a**3*x**4/4 + C*a**2*b*x**6/2 + 3*C*a*b**2*x**8/8 + C*b**3*x**10/10 + D*b**3*x**11/11 + a**2*x**5*(3*B*b + D*a)/5 + 3*a*b*x**7*(B*b + D*a)/7 + b**2*x**9*(B*b + 3*D*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_81():
    f = (a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)
    F = A*a**3*x + B*(a + b*x**2)**4/(8*b) + C*b**3*x**9/9 + D*a**3*x**4/4 + D*a**2*b*x**6/2 + 3*D*a*b**2*x**8/8 + D*b**3*x**10/10 + a**2*x**3*(3*A*b + C*a)/3 + 3*a*b*x**5*(A*b + C*a)/5 + b**2*x**7*(A*b + 3*C*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_82():
    f = (a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)/x
    F = A*a**3*log(x) + 3*A*a**2*b*x**2/2 + 3*A*a*b**2*x**4/4 + A*b**3*x**6/6 + B*a**3*x + C*(a + b*x**2)**4/(8*b) + D*b**3*x**9/9 + a**2*x**3*(3*B*b + D*a)/3 + 3*a*b*x**5*(B*b + D*a)/5 + b**2*x**7*(B*b + 3*D*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_83():
    f = (a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)/x**2
    F = -A*a**3/x + B*a**3*log(x) + 3*B*a**2*b*x**2/2 + 3*B*a*b**2*x**4/4 + B*b**3*x**6/6 + C*b**3*x**7/7 + D*(a + b*x**2)**4/(8*b) + a**2*x*(3*A*b + C*a) + a*b*x**3*(A*b + C*a) + b**2*x**5*(A*b + 3*C*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_84():
    f = (a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)/x**3
    F = -A*a**3/(2*x**2) - B*a**3/x + C*b**3*x**6/6 + D*b**3*x**7/7 + a**2*x*(3*B*b + D*a) + a**2*(3*A*b + C*a)*log(x) + a*b*x**3*(B*b + D*a) + 3*a*b*x**2*(A*b + C*a)/2 + b**2*x**5*(B*b + 3*D*a)/5 + b**2*x**4*(A*b + 3*C*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_85():
    f = (a + b*x**2)**3*(A + B*x + C*x**2 + D*x**3)/x**4
    F = -A*a**3/(3*x**3) - B*a**3/(2*x**2) + C*b**3*x**5/5 + D*b**3*x**6/6 + a**2*(3*B*b + D*a)*log(x) - a**2*(3*A*b + C*a)/x + 3*a*b*x**2*(B*b + D*a)/2 + 3*a*b*x*(A*b + C*a) + b**2*x**4*(B*b + 3*D*a)/4 + b**2*x**3*(A*b + 3*C*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_86():
    f = x**4*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)
    F = C*x**5/(5*b) + D*x**6/(6*b) + a**(sympy.S(3)/2)*(A*b - C*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) + a**2*(B*b - D*a)*log(a + b*x**2)/(2*b**4) - a*x**2*(B*b - D*a)/(2*b**3) - a*x*(A*b - C*a)/b**3 + x**4*(B*b - D*a)/(4*b**2) + x**3*(A*b - C*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_87():
    f = x**3*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)
    F = C*x**4/(4*b) + D*x**5/(5*b) + a**(sympy.S(3)/2)*(B*b - D*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) - a*x*(B*b - D*a)/b**3 - a*(A*b - C*a)*log(a + b*x**2)/(2*b**3) + x**3*(B*b - D*a)/(3*b**2) + x**2*(A*b - C*a)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_88():
    f = x**2*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)
    F = C*x**3/(3*b) + D*x**4/(4*b) - sqrt(a)*(A*b - C*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(5)/2) - a*(B*b - D*a)*log(a + b*x**2)/(2*b**3) + x**2*(B*b - D*a)/(2*b**2) + x*(A*b - C*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_89():
    f = x*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)
    F = C*x**2/(2*b) + D*x**3/(3*b) - sqrt(a)*(B*b - D*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(5)/2) + x*(B*b - D*a)/b**2 + (A*b - C*a)*log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_90():
    f = (A + B*x + C*x**2 + D*x**3)/(a + b*x**2)
    F = C*x/b + D*x**2/(2*b) + (B*b - D*a)*log(a + b*x**2)/(2*b**2) + (A*b - C*a)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_91():
    f = (A + B*x + C*x**2 + D*x**3)/(x*(a + b*x**2))
    F = A*log(x)/a + D*x/b - (A*b - C*a)*log(a + b*x**2)/(2*a*b) + (B*b - D*a)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_92():
    f = (A + B*x + C*x**2 + D*x**3)/(x**2*(a + b*x**2))
    F = -A/(a*x) + B*log(x)/a - (B*b - D*a)*log(a + b*x**2)/(2*a*b) - (A*b - C*a)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_93():
    f = (A + B*x + C*x**2 + D*x**3)/(x**3*(a + b*x**2))
    F = -A/(2*a*x**2) - B/(a*x) - (A*b - C*a)*log(x)/a**2 + (A*b - C*a)*log(a + b*x**2)/(2*a**2) - (B*b - D*a)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_94():
    f = x**4*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**2
    F = D*x**4/(4*b**2) - sqrt(a)*(3*A*b - 5*C*a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) - a*(2*B*b - 3*D*a)*log(a + b*x**2)/(2*b**4) + x**2*(2*B*b - 3*D*a)/(2*b**3) + x*(3*A*b - 5*C*a)/(2*b**3) - x**4*(a*(B - D*a/b) - x*(A*b - C*a))/(2*a*b*(a + b*x**2)) - x**3*(3*A*b - 5*C*a)/(6*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_95():
    f = x**3*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**2
    F = D*x**3/(3*b**2) - sqrt(a)*(3*B*b - 5*D*a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) + x*(3*B*b - 5*D*a)/(2*b**3) + (A*b - 2*C*a)*log(a + b*x**2)/(2*b**3) - x**3*(a*(B - D*a/b) - x*(A*b - C*a))/(2*a*b*(a + b*x**2)) - x**2*(A*b - 2*C*a)/(2*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_96():
    f = x**2*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**2
    F = D*x**2/(2*b**2) + (B*b - 2*D*a)*log(a + b*x**2)/(2*b**3) - x**2*(a*(B - D*a/b) - x*(A*b - C*a))/(2*a*b*(a + b*x**2)) - x*(A*b - 3*C*a)/(2*a*b**2) + (A*b - 3*C*a)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_97():
    f = x*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**2
    F = C*log(a + b*x**2)/(2*b**2) + D*x/b**2 - x*(a*(B - D*a/b) - x*(A*b - C*a))/(2*a*b*(a + b*x**2)) + (B*b - 3*D*a)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_98():
    f = (A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**2
    F = D*log(a + b*x**2)/(2*b**2) + (-a*(B - D*a/b) + x*(A*b - C*a))/(2*a*b*(a + b*x**2)) + (A*b + C*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_99():
    f = (A + B*x + C*x**2 + D*x**3)/(x*(a + b*x**2)**2)
    F = A*log(x)/a**2 - A*log(a + b*x**2)/(2*a**2) + (A*b - C*a + x*(B*b - D*a))/(2*a*b*(a + b*x**2)) + (B*b + D*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_100():
    f = (A + B*x + C*x**2 + D*x**3)/(x**2*(a + b*x**2)**2)
    F = -A/(a**2*x) + B*log(x)/a**2 - B*log(a + b*x**2)/(2*a**2) + (B*b - D*a - b*x*(A*b/a - C))/(2*a*b*(a + b*x**2)) - (3*A*b - C*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_101():
    f = (A + B*x + C*x**2 + D*x**3)/(x**3*(a + b*x**2)**2)
    F = -A/(2*a**2*x**2) - B/(a**2*x) - (A*b/a - C + x*(B*b/a - D))/(2*a*(a + b*x**2)) - (2*A*b - C*a)*log(x)/a**3 + (2*A*b - C*a)*log(a + b*x**2)/(2*a**3) - (3*B*b - D*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_102():
    f = x**4*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**3
    F = (B*b - 3*D*a)*log(a + b*x**2)/(2*b**4) - x**4*(a*(B - D*a/b) - x*(A*b - C*a))/(4*a*b*(a + b*x**2)**2) + x**3*(A*b - 5*C*a + x*(4*B*b - 8*D*a))/(8*a*b**2*(a + b*x**2)) - x**2*(B*b - 3*D*a)/(2*a*b**3) + x*(-3*A*b + 15*C*a)/(8*a*b**3) + (3*A*b - 15*C*a)*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_103():
    f = x**3*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**3
    F = C*log(a + b*x**2)/(2*b**3) - x**3*(a*(B - D*a/b) - x*(A*b - C*a))/(4*a*b*(a + b*x**2)**2) - x**2*(4*C*a - x*(3*B*b - 7*D*a))/(8*a*b**2*(a + b*x**2)) + x*(-3*B*b + 15*D*a)/(8*a*b**3) + (3*B*b - 15*D*a)*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_104():
    f = x**2*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**3
    F = D*log(a + b*x**2)/(2*b**3) - x**2*(a*(B - D*a/b) - x*(A*b - C*a))/(4*a*b*(a + b*x**2)**2) - x*(A*b + 3*C*a - x*(2*B*b - 6*D*a))/(8*a*b**2*(a + b*x**2)) + (A*b + 3*C*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_105():
    f = x*(A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**3
    F = -x*(a*(B - D*a/b) - x*(A*b - C*a))/(4*a*b*(a + b*x**2)**2) - (2*A*b + 2*C*a - x*(B*b - 5*D*a))/(8*a*b**2*(a + b*x**2)) + (B*b + 3*D*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_106():
    f = (A + B*x + C*x**2 + D*x**3)/(a + b*x**2)**3
    F = (-a*(B - D*a/b) + x*(A*b - C*a))/(4*a*b*(a + b*x**2)**2) - (4*D*a**2 - b*x*(3*A*b + C*a))/(8*a**2*b**2*(a + b*x**2)) + (3*A*b + C*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_107():
    f = (A + B*x + C*x**2 + D*x**3)/(x*(a + b*x**2)**3)
    F = A*log(x)/a**3 - A*log(a + b*x**2)/(2*a**3) + (A*b - C*a + x*(B*b - D*a))/(4*a*b*(a + b*x**2)**2) + (4*A*b + x*(3*B*b + D*a))/(8*a**2*b*(a + b*x**2)) + (3*B*b + D*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_108():
    f = (A + B*x + C*x**2 + D*x**3)/(x**2*(a + b*x**2)**3)
    F = -A/(a**3*x) + B*log(x)/a**3 - B*log(a + b*x**2)/(2*a**3) + (B*b - D*a - b*x*(A*b/a - C))/(4*a*b*(a + b*x**2)**2) + (4*B - x*(7*A*b/a - 3*C))/(8*a**2*(a + b*x**2)) - (15*A*b - 3*C*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_109():
    f = (A + B*x + C*x**2 + D*x**3)/(x**3*(a + b*x**2)**3)
    F = -A/(2*a**3*x**2) - B/(a**3*x) - (A*b/a - C + x*(B*b/a - D))/(4*a*(a + b*x**2)**2) - (8*A*b - 4*C*a + x*(7*B*b - 3*D*a))/(8*a**3*(a + b*x**2)) - (3*A*b - C*a)*log(x)/a**4 + (3*A*b - C*a)*log(a + b*x**2)/(2*a**4) - (15*B*b - 3*D*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_110():
    f = (4*x**3 - x)/(x**2 + 5)**2
    F = 2*log(x**2 + 5) + 21/(2*x**2 + 10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_111():
    f = (x**3 - x)/sqrt(x**2 - 2)
    F = (x**2 - 2)**(sympy.S(3)/2)/3 + sqrt(x**2 - 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_112():
    f = (2*x**4 - x**2)/(2*x**2 + 1)
    F = x**3/3 - x + sqrt(2)*atan(sqrt(2)*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_113():
    f = (x**4 + x**3)/(x**2 + 1)
    F = x**3/3 + x**2/2 - x - log(x**2 + 1)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_114():
    f = x**6*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)
    F = -a**(sympy.S(5)/2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(13)/2) + a**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**6 - a*x**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**5) + f*x**11/(11*b) + x**9*(-a*f + b*e)/(9*b**2) + x**7*(a**2*f - a*b*e + b**2*d)/(7*b**3) + x**5*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_115():
    f = x**4*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(11)/2) - a*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**5 + f*x**9/(9*b) + x**7*(-a*f + b*e)/(7*b**2) + x**5*(a**2*f - a*b*e + b**2*d)/(5*b**3) + x**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_116():
    f = x**2*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)
    F = -sqrt(a)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(9)/2) + f*x**7/(7*b) + x**5*(-a*f + b*e)/(5*b**2) + x**3*(a**2*f - a*b*e + b**2*d)/(3*b**3) + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_117():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)
    F = f*x**5/(5*b) + x**3*(-a*f + b*e)/(3*b**2) + x*(a**2*f - a*b*e + b**2*d)/b**3 + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_118():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**2*(a + b*x**2))
    F = f*x**3/(3*b) + x*(-a*f + b*e)/b**2 - c/(a*x) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_119():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**4*(a + b*x**2))
    F = f*x/b - c/(3*a*x**3) + (-a*d + b*c)/(a**2*x) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_120():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**6*(a + b*x**2))
    F = -c/(5*a*x**5) + (-a*d + b*c)/(3*a**2*x**3) - (a**2*e - a*b*d + b**2*c)/(a**3*x) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_121():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**8*(a + b*x**2))
    F = -c/(7*a*x**7) + (-a*d + b*c)/(5*a**2*x**5) - (a**2*e - a*b*d + b**2*c)/(3*a**3*x**3) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**4*x) + sqrt(b)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_122():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**10*(a + b*x**2))
    F = -c/(9*a*x**9) + (-a*d + b*c)/(7*a**2*x**7) - (a**2*e - a*b*d + b**2*c)/(5*a**3*x**5) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**4*x**3) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**5*x) - b**(sympy.S(3)/2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(11)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_123():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**12*(a + b*x**2))
    F = -c/(11*a*x**11) + (-a*d + b*c)/(9*a**2*x**9) - (a**2*e - a*b*d + b**2*c)/(7*a**3*x**7) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(5*a**4*x**5) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**5*x**3) + b**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**6*x) + b**(sympy.S(5)/2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(13)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_124():
    f = x**6*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**2
    F = a**(sympy.S(3)/2)*(-11*a**3*f + 9*a**2*b*e - 7*a*b**2*d + 5*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(13)/2)) - a*x*(-11*a**3*f + 9*a**2*b*e - 7*a*b**2*d + 5*b**3*c)/(2*b**6) + f*x**9/(9*b**2) + x**7*(-2*a*f + b*e)/(7*b**3) + x**3*(-11*a**3*f + 9*a**2*b*e - 7*a*b**2*d + 5*b**3*c)/(6*b**5) + x**7*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(2*a*(a + b*x**2)) - x**5*(-11*a**3*f + 9*a**2*b*e - 7*a*b**2*d + 5*b**3*c)/(10*a*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_125():
    f = x**4*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**2
    F = -sqrt(a)*(-9*a**3*f + 7*a**2*b*e - 5*a*b**2*d + 3*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(11)/2)) + f*x**7/(7*b**2) + x**5*(-2*a*f + b*e)/(5*b**3) + x*(-9*a**3*f + 7*a**2*b*e - 5*a*b**2*d + 3*b**3*c)/(2*b**5) + x**5*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(2*a*(a + b*x**2)) - x**3*(-9*a**3*f + 7*a**2*b*e - 5*a*b**2*d + 3*b**3*c)/(6*a*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_126():
    f = x**2*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**2
    F = f*x**5/(5*b**2) + x**3*(-2*a*f + b*e)/(3*b**3) + x**3*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(2*a*(a + b*x**2)) - x*(-7*a**3*f + 5*a**2*b*e - 3*a*b**2*d + b**3*c)/(2*a*b**4) + (-7*a**3*f + 5*a**2*b*e - 3*a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_127():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**2
    F = f*x**3/(3*b**2) + x*(-2*a*f + b*e)/b**3 + x*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(2*a*(a + b*x**2)) + (5*a**3*f - 3*a**2*b*e + a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_128():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**2*(a + b*x**2)**2)
    F = f*x/b**2 - x*(-a**2*f/b**2 + a*e/b - d + b*c/a)/(2*a*(a + b*x**2)) - c/(a**2*x) - (3*a**3*f - a**2*b*e - a*b**2*d + 3*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_129():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**4*(a + b*x**2)**2)
    F = x*(-a*f/b + e - b*d/a + b**2*c/a**2)/(2*a*(a + b*x**2)) - c/(3*a**2*x**3) + (-a*d + 2*b*c)/(a**3*x) + (a**3*f + a**2*b*e - 3*a*b**2*d + 5*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_130():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**6*(a + b*x**2)**2)
    F = -c/(5*a**2*x**5) + (-a*d + 2*b*c)/(3*a**3*x**3) - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*a**4*(a + b*x**2)) - (a**2*e - 2*a*b*d + 3*b**2*c)/(a**4*x) - (-a**3*f + 3*a**2*b*e - 5*a*b**2*d + 7*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(9)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_131():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**8*(a + b*x**2)**2)
    F = -c/(7*a**2*x**7) + (-a*d + 2*b*c)/(5*a**3*x**5) - (a**2*e - 2*a*b*d + 3*b**2*c)/(3*a**4*x**3) + b*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*a**5*(a + b*x**2)) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(a**5*x) + sqrt(b)*(-3*a**3*f + 5*a**2*b*e - 7*a*b**2*d + 9*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_132():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**10*(a + b*x**2)**2)
    F = -c/(9*a**2*x**9) + (-a*d + 2*b*c)/(7*a**3*x**7) - (a**2*e - 2*a*b*d + 3*b**2*c)/(5*a**4*x**5) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(3*a**5*x**3) - b**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*a**6*(a + b*x**2)) - b*(-2*a**3*f + 3*a**2*b*e - 4*a*b**2*d + 5*b**3*c)/(a**6*x) - b**(sympy.S(3)/2)*(-5*a**3*f + 7*a**2*b*e - 9*a*b**2*d + 11*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_133():
    f = x**8*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**3
    F = a**(sympy.S(3)/2)*(-143*a**3*f + 99*a**2*b*e - 63*a*b**2*d + 35*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(15)/2)) - a**2*x*(-17*a**3*f + 13*a**2*b*e - 9*a*b**2*d + 5*b**3*c)/(8*b**7*(a + b*x**2)) - a*x*(-63*a**3*f + 43*a**2*b*e - 27*a*b**2*d + 15*b**3*c)/(4*b**7) + f*x**9/(9*b**3) + x**7*(-3*a*f + b*e)/(7*b**4) + x**3*(-23*a**3*f + 15*a**2*b*e - 9*a*b**2*d + 5*b**3*c)/(6*b**6) + x**9*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(4*a*(a + b*x**2)**2) - x**5*(-29*a**3*f + 17*a**2*b*e - 9*a*b**2*d + 5*b**3*c)/(20*a*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_134():
    f = x**6*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**3
    F = -sqrt(a)*(-99*a**3*f + 63*a**2*b*e - 35*a*b**2*d + 15*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(13)/2)) + a*x*(-15*a**3*f + 11*a**2*b*e - 7*a*b**2*d + 3*b**3*c)/(8*b**6*(a + b*x**2)) + f*x**7/(7*b**3) + x**5*(-3*a*f + b*e)/(5*b**4) + x*(-21*a**3*f + 13*a**2*b*e - 7*a*b**2*d + 3*b**3*c)/(2*b**6) + x**7*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(4*a*(a + b*x**2)**2) - x**3*(-27*a**3*f + 15*a**2*b*e - 7*a*b**2*d + 3*b**3*c)/(12*a*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_135():
    f = x**4*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**3
    F = f*x**5/(5*b**3) + x**3*(-3*a*f + b*e)/(3*b**4) - x*(-13*a**3*f + 9*a**2*b*e - 5*a*b**2*d + b**3*c)/(8*b**5*(a + b*x**2)) + x**5*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(4*a*(a + b*x**2)**2) - x*(-25*a**3*f + 13*a**2*b*e - 5*a*b**2*d + b**3*c)/(4*a*b**5) + (-63*a**3*f + 35*a**2*b*e - 15*a*b**2*d + 3*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_136():
    f = x**2*(c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**3
    F = f*x**3/(3*b**3) + x*(-3*a*f + b*e)/b**4 + x**3*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(4*a*(a + b*x**2)**2) - x*(11*a**3*f - 7*a**2*b*e + 3*a*b**2*d + b**3*c)/(8*a*b**4*(a + b*x**2)) + (35*a**3*f - 15*a**2*b*e + 3*a*b**2*d + b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_137():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(a + b*x**2)**3
    F = f*x/b**3 + x*(-a*(a**2*f - a*b*e + b**2*d)/b**3 + c)/(4*a*(a + b*x**2)**2) + x*(9*a**3*f - 5*a**2*b*e + a*b**2*d + 3*b**3*c)/(8*a**2*b**3*(a + b*x**2)) + (-15*a**3*f + 3*a**2*b*e + a*b**2*d + 3*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_138():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**2*(a + b*x**2)**3)
    F = -x*(-a**2*f/b**2 + a*e/b - d + b*c/a)/(4*a*(a + b*x**2)**2) - c/(a**3*x) - x*(5*a**3*f - a**2*b*e - 3*a*b**2*d + 7*b**3*c)/(8*a**3*b**2*(a + b*x**2)) - (-3*a**3*f - a**2*b*e - 3*a*b**2*d + 15*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_139():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**4*(a + b*x**2)**3)
    F = x*(-a*f/b + e - b*d/a + b**2*c/a**2)/(4*a*(a + b*x**2)**2) - c/(3*a**3*x**3) + (-a*d + 3*b*c)/(a**4*x) + x*(a**3*f + 3*a**2*b*e - 7*a*b**2*d + 11*b**3*c)/(8*a**4*b*(a + b*x**2)) + (a**3*f + 3*a**2*b*e - 15*a*b**2*d + 35*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(9)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_140():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**6*(a + b*x**2)**3)
    F = -c/(5*a**3*x**5) - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*a**4*(a + b*x**2)**2) + (-a*d + 3*b*c)/(3*a**4*x**3) - x*(-3*a**3*f + 7*a**2*b*e - 11*a*b**2*d + 15*b**3*c)/(8*a**5*(a + b*x**2)) - (a**2*e - 3*a*b*d + 6*b**2*c)/(a**5*x) - (-3*a**3*f + 15*a**2*b*e - 35*a*b**2*d + 63*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(11)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_141():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**8*(a + b*x**2)**3)
    F = -c/(7*a**3*x**7) + (-a*d + 3*b*c)/(5*a**4*x**5) + b*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*a**5*(a + b*x**2)**2) - (a**2*e - 3*a*b*d + 6*b**2*c)/(3*a**5*x**3) + b*x*(-7*a**3*f + 11*a**2*b*e - 15*a*b**2*d + 19*b**3*c)/(8*a**6*(a + b*x**2)) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(a**6*x) + sqrt(b)*(-15*a**3*f + 35*a**2*b*e - 63*a*b**2*d + 99*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_142():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**10*(a + b*x**2)**3)
    F = -c/(9*a**3*x**9) + (-a*d + 3*b*c)/(7*a**4*x**7) - (a**2*e - 3*a*b*d + 6*b**2*c)/(5*a**5*x**5) - b**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*a**6*(a + b*x**2)**2) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(3*a**6*x**3) - b**2*x*(-11*a**3*f + 15*a**2*b*e - 19*a*b**2*d + 23*b**3*c)/(8*a**7*(a + b*x**2)) - b*(-3*a**3*f + 6*a**2*b*e - 10*a*b**2*d + 15*b**3*c)/(a**7*x) - b**(sympy.S(3)/2)*(-35*a**3*f + 63*a**2*b*e - 99*a*b**2*d + 143*b**3*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_143():
    f = x**5*(c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = a**2*sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**6 - a*(a + b*x**2)**(sympy.S(3)/2)*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)/(3*b**6) + f*(a + b*x**2)**(sympy.S(11)/2)/(11*b**6) + (a + b*x**2)**(sympy.S(9)/2)*(-5*a*f + b*e)/(9*b**6) + (a + b*x**2)**(sympy.S(7)/2)*(10*a**2*f - 4*a*b*e + b**2*d)/(7*b**6) + (a + b*x**2)**(sympy.S(5)/2)*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(5*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_144():
    f = x**3*(c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = -a*sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**5 + f*(a + b*x**2)**(sympy.S(9)/2)/(9*b**5) + (a + b*x**2)**(sympy.S(7)/2)*(-4*a*f + b*e)/(7*b**5) + (a + b*x**2)**(sympy.S(5)/2)*(6*a**2*f - 3*a*b*e + b**2*d)/(5*b**5) + (a + b*x**2)**(sympy.S(3)/2)*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_145():
    f = x*(c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = f*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) + (a + b*x**2)**(sympy.S(5)/2)*(-3*a*f + b*e)/(5*b**4) + (a + b*x**2)**(sympy.S(3)/2)*(3*a**2*f - 2*a*b*e + b**2*d)/(3*b**4) + sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_146():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x*sqrt(a + b*x**2))
    F = f*(a + b*x**2)**(sympy.S(5)/2)/(5*b**3) + (a + b*x**2)**(sympy.S(3)/2)*(-2*a*f + b*e)/(3*b**3) + sqrt(a + b*x**2)*(a**2*f - a*b*e + b**2*d)/b**3 - c*atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_147():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**3*sqrt(a + b*x**2))
    F = f*(a + b*x**2)**(sympy.S(3)/2)/(3*b**2) + sqrt(a + b*x**2)*(-a*f + b*e)/b**2 - c*sqrt(a + b*x**2)/(2*a*x**2) + (-2*a*d + b*c)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_148():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**5*sqrt(a + b*x**2))
    F = f*sqrt(a + b*x**2)/b - c*sqrt(a + b*x**2)/(4*a*x**4) + sqrt(a + b*x**2)*(-4*a*d + 3*b*c)/(8*a**2*x**2) - (8*a**2*e - 4*a*b*d + 3*b**2*c)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_149():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**7*sqrt(a + b*x**2))
    F = -c*sqrt(a + b*x**2)/(6*a*x**6) + sqrt(a + b*x**2)*(-6*a*d + 5*b*c)/(24*a**2*x**4) - sqrt(a + b*x**2)*(8*a**2*e - 6*a*b*d + 5*b**2*c)/(16*a**3*x**2) + (-16*a**3*f + 8*a**2*b*e - 6*a*b**2*d + 5*b**3*c)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_150():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**9*sqrt(a + b*x**2))
    F = -c*sqrt(a + b*x**2)/(8*a*x**8) + sqrt(a + b*x**2)*(-8*a*d + 7*b*c)/(48*a**2*x**6) - sqrt(a + b*x**2)*(48*a**2*e - 40*a*b*d + 35*b**2*c)/(192*a**3*x**4) + sqrt(a + b*x**2)*(-64*a**3*f + 48*a**2*b*e - 40*a*b**2*d + 35*b**3*c)/(128*a**4*x**2) - b*(-64*a**3*f + 48*a**2*b*e - 40*a*b**2*d + 35*b**3*c)*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_151():
    f = x**4*(c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = a**2*(-63*a**3*f + 70*a**2*b*e - 80*a*b**2*d + 96*b**3*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(11)/2)) - a*x*sqrt(a + b*x**2)*(-63*a**3*f + 70*a**2*b*e - 80*a*b**2*d + 96*b**3*c)/(256*b**5) + f*x**9*sqrt(a + b*x**2)/(10*b) + x**7*sqrt(a + b*x**2)*(-9*a*f + 10*b*e)/(80*b**2) + x**5*sqrt(a + b*x**2)*(63*a**2*f - 70*a*b*e + 80*b**2*d)/(480*b**3) + x**3*sqrt(a + b*x**2)*(-63*a**3*f + 70*a**2*b*e - 80*a*b**2*d + 96*b**3*c)/(384*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_152():
    f = x**2*(c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = -a*(-35*a**3*f + 40*a**2*b*e - 48*a*b**2*d + 64*b**3*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(9)/2)) + f*x**7*sqrt(a + b*x**2)/(8*b) + x**5*sqrt(a + b*x**2)*(-7*a*f + 8*b*e)/(48*b**2) + x**3*sqrt(a + b*x**2)*(35*a**2*f - 40*a*b*e + 48*b**2*d)/(192*b**3) + x*sqrt(a + b*x**2)*(-35*a**3*f + 40*a**2*b*e - 48*a*b**2*d + 64*b**3*c)/(128*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_153():
    f = (c + d*x**2 + e*x**4 + f*x**6)/sqrt(a + b*x**2)
    F = f*x**5*sqrt(a + b*x**2)/(6*b) + x**3*sqrt(a + b*x**2)*(-5*a*f + 6*b*e)/(24*b**2) + x*sqrt(a + b*x**2)*(5*a**2*f - 6*a*b*e + 8*b**2*d)/(16*b**3) + (-5*a**3*f + 6*a**2*b*e - 8*a*b**2*d + 16*b**3*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_154():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**2*sqrt(a + b*x**2))
    F = f*x**3*sqrt(a + b*x**2)/(4*b) + x*sqrt(a + b*x**2)*(-3*a*f + 4*b*e)/(8*b**2) + (3*a**2*f - 4*a*b*e + 8*b**2*d)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2)) - c*sqrt(a + b*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_155():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**4*sqrt(a + b*x**2))
    F = f*x*sqrt(a + b*x**2)/(2*b) + (-a*f + 2*b*e)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2)) - c*sqrt(a + b*x**2)/(3*a*x**3) + sqrt(a + b*x**2)*(-3*a*d + 2*b*c)/(3*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_156():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**6*sqrt(a + b*x**2))
    F = f*atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b) - c*sqrt(a + b*x**2)/(5*a*x**5) + sqrt(a + b*x**2)*(-5*a*d + 4*b*c)/(15*a**2*x**3) - sqrt(a + b*x**2)*(15*a**2*e - 10*a*b*d + 8*b**2*c)/(15*a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_157():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**8*sqrt(a + b*x**2))
    F = -c*sqrt(a + b*x**2)/(7*a*x**7) + sqrt(a + b*x**2)*(-7*a*d + 6*b*c)/(35*a**2*x**5) - sqrt(a + b*x**2)*(35*a**2*e - 28*a*b*d + 24*b**2*c)/(105*a**3*x**3) + sqrt(a + b*x**2)*(-105*a**3*f + 70*a**2*b*e - 56*a*b**2*d + 48*b**3*c)/(105*a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_158():
    f = (c + d*x**2 + e*x**4 + f*x**6)/(x**10*sqrt(a + b*x**2))
    F = -c*sqrt(a + b*x**2)/(9*a*x**9) + sqrt(a + b*x**2)*(-9*a*d + 8*b*c)/(63*a**2*x**7) - sqrt(a + b*x**2)*(21*a**2*e - 18*a*b*d + 16*b**2*c)/(105*a**3*x**5) + sqrt(a + b*x**2)*(-105*a**3*f + 84*a**2*b*e - 72*a*b**2*d + 64*b**3*c)/(315*a**4*x**3) - 2*b*sqrt(a + b*x**2)*(-105*a**3*f + 84*a**2*b*e - 72*a*b**2*d + 64*b**3*c)/(315*a**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_159():
    f = x**8*(A + B*x**2 + C*x**4 + D*x**6)/(a + b*x**2)**(sympy.S(9)/2)
    F = D*x**9/(6*b**3*(a + b*x**2)**(sympy.S(3)/2)) + (16*A*b**3 - 72*B*a*b**2 + 198*C*a**2*b - 429*D*a**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(15)/2)) + x**9*(A - a*(B*b**2 - C*a*b + D*a**2)/b**3)/(7*a*(a + b*x**2)**(sympy.S(7)/2)) - x*sqrt(a + b*x**2)*(16*A*b**3 - 3*a*(24*B*b**2 - 66*C*a*b + 143*D*a**2))/(16*a*b**7) - x**9*(2*A*b**3 - a*(9*B*b**2 - 16*C*a*b + 23*D*a**2))/(35*a**2*b**3*(a + b*x**2)**(sympy.S(5)/2)) - x**7*(16*A*b**3 - 3*a*(24*B*b**2 - 66*C*a*b + 143*D*a**2))/(210*a**2*b**4*(a + b*x**2)**(sympy.S(3)/2)) - x**5*(16*A*b**3 - 3*a*(24*B*b**2 - 66*C*a*b + 143*D*a**2))/(30*a**2*b**5*sqrt(a + b*x**2)) + x**3*sqrt(a + b*x**2)*(16*A*b**3 - 3*a*(24*B*b**2 - 66*C*a*b + 143*D*a**2))/(24*a**2*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_160():
    f = x**6*(A + B*x**2 + C*x**4 + D*x**6)/(a + b*x**2)**(sympy.S(9)/2)
    F = D*x**7/(4*b**3*(a + b*x**2)**(sympy.S(3)/2)) + (8*B*b**2 - 36*C*a*b + 99*D*a**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(13)/2)) + x**7*(A - a*(B*b**2 - C*a*b + D*a**2)/b**3)/(7*a*(a + b*x**2)**(sympy.S(7)/2)) + x**7*(B*b**2 - 2*C*a*b + 3*D*a**2)/(5*a*b**3*(a + b*x**2)**(sympy.S(5)/2)) + x**5*(8*B*b**2 - 36*C*a*b + 99*D*a**2)/(60*a*b**4*(a + b*x**2)**(sympy.S(3)/2)) + x**3*(8*B*b**2 - 36*C*a*b + 99*D*a**2)/(12*a*b**5*sqrt(a + b*x**2)) - x*sqrt(a + b*x**2)*(8*B*b**2 - 36*C*a*b + 99*D*a**2)/(8*a*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_161():
    f = x**4*(A + B*x**2 + C*x**4 + D*x**6)/(a + b*x**2)**(sympy.S(9)/2)
    F = D*x*sqrt(a + b*x**2)/(2*b**5) + a*x*(C*b - 3*D*a)/(3*b**5*(a + b*x**2)**(sympy.S(3)/2)) - x*(4*C*b - 15*D*a)/(3*b**5*sqrt(a + b*x**2)) + (2*C*b - 9*D*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(11)/2)) + x**5*(A - a*(B*b**2 - C*a*b + D*a**2)/b**3)/(7*a*(a + b*x**2)**(sympy.S(7)/2)) + x**5*(2*A*b**3 + a*(5*B*b**2 - 12*C*a*b + 19*D*a**2))/(35*a**2*b**3*(a + b*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_162():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(a + b*x**2)**(sympy.S(9)/2)
    F = A*x/(a*(a + b*x**2)**(sympy.S(7)/2)) + x**3*(6*A*b + B*a)/(3*a**2*(a + b*x**2)**(sympy.S(7)/2)) + x**5*(24*A*b**2 + a*(4*B*b + 3*C*a))/(15*a**3*(a + b*x**2)**(sympy.S(7)/2)) + x**7*(48*A*b**3 + a*(8*B*b**2 + 6*C*a*b + 15*D*a**2))/(105*a**4*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_163():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(x**2*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(a*x*(a + b*x**2)**(sympy.S(7)/2)) - x*(8*A*b - B*a)/(a**2*(a + b*x**2)**(sympy.S(7)/2)) - x**3*(48*A*b**2 - a*(6*B*b + C*a))/(3*a**3*(a + b*x**2)**(sympy.S(7)/2)) - x**5*(-3*D*a**3 + 4*b*(48*A*b**2 - a*(6*B*b + C*a)))/(15*a**4*(a + b*x**2)**(sympy.S(7)/2)) - 2*b*x**7*(-3*D*a**3 + 4*b*(48*A*b**2 - a*(6*B*b + C*a)))/(105*a**5*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_164():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(x**4*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(3*a*x**3*(a + b*x**2)**(sympy.S(7)/2)) + (10*A*b - 3*B*a)/(3*a**2*x*(a + b*x**2)**(sympy.S(7)/2)) + x*(80*A*b**2 - 3*a*(8*B*b - C*a))/(3*a**3*(a + b*x**2)**(sympy.S(7)/2)) + x**3*(160*A*b**3 - a*(48*B*b**2 - 6*C*a*b - D*a**2))/(3*a**4*(a + b*x**2)**(sympy.S(7)/2)) + 4*b*x**5*(160*A*b**3 - a*(48*B*b**2 - 6*C*a*b - D*a**2))/(15*a**5*(a + b*x**2)**(sympy.S(7)/2)) + 8*b**2*x**7*(160*A*b**3 - a*(48*B*b**2 - 6*C*a*b - D*a**2))/(105*a**6*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_165():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(x**6*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(5*a*x**5*(a + b*x**2)**(sympy.S(7)/2)) + (12*A*b - 5*B*a)/(15*a**2*x**3*(a + b*x**2)**(sympy.S(7)/2)) - (24*A*b**2 - a*(10*B*b - 3*C*a))/(3*a**3*x*(a + b*x**2)**(sympy.S(7)/2)) - x*(192*A*b**3 - a*(80*B*b**2 - 24*C*a*b + 3*D*a**2))/(21*a**4*(a + b*x**2)**(sympy.S(7)/2)) - x*(384*A*b**3 - 2*a*(80*B*b**2 - 24*C*a*b + 3*D*a**2))/(35*a**5*(a + b*x**2)**(sympy.S(5)/2)) - x*(1536*A*b**3 - 8*a*(80*B*b**2 - 24*C*a*b + 3*D*a**2))/(105*a**6*(a + b*x**2)**(sympy.S(3)/2)) - x*(3072*A*b**3 - 16*a*(80*B*b**2 - 24*C*a*b + 3*D*a**2))/(105*a**7*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_166():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(x**8*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(7*a*x**7*(a + b*x**2)**(sympy.S(7)/2)) + (2*A*b - B*a)/(5*a**2*x**5*(a + b*x**2)**(sympy.S(7)/2)) - (24*A*b**2 - a*(12*B*b - 5*C*a))/(15*a**3*x**3*(a + b*x**2)**(sympy.S(7)/2)) + (48*A*b**3 - a*(24*B*b**2 - 10*C*a*b + 3*D*a**2))/(3*a**4*x*(a + b*x**2)**(sympy.S(7)/2)) + 8*b*x*(48*A*b**3 - a*(24*B*b**2 - 10*C*a*b + 3*D*a**2))/(21*a**5*(a + b*x**2)**(sympy.S(7)/2)) + 16*b*x*(48*A*b**3 - a*(24*B*b**2 - 10*C*a*b + 3*D*a**2))/(35*a**6*(a + b*x**2)**(sympy.S(5)/2)) + 64*b*x*(48*A*b**3 - a*(24*B*b**2 - 10*C*a*b + 3*D*a**2))/(105*a**7*(a + b*x**2)**(sympy.S(3)/2)) + 128*b*x*(48*A*b**3 - a*(24*B*b**2 - 10*C*a*b + 3*D*a**2))/(105*a**8*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_167():
    f = (A + B*x**2 + C*x**4 + D*x**6)/(x**10*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(9*a*x**9*(a + b*x**2)**(sympy.S(7)/2)) + (16*A*b - 9*B*a)/(63*a**2*x**7*(a + b*x**2)**(sympy.S(7)/2)) - (32*A*b**2 - 9*a*(2*B*b - C*a))/(45*a**3*x**5*(a + b*x**2)**(sympy.S(7)/2)) + (128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(45*a**4*x**3*(a + b*x**2)**(sympy.S(7)/2)) - 2*b*(128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(9*a**5*x*(a + b*x**2)**(sympy.S(7)/2)) - 16*b**2*x*(128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(63*a**6*(a + b*x**2)**(sympy.S(7)/2)) - 32*b**2*x*(128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(105*a**7*(a + b*x**2)**(sympy.S(5)/2)) - 128*b**2*x*(128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(315*a**8*(a + b*x**2)**(sympy.S(3)/2)) - 256*b**2*x*(128*A*b**3 - 3*a*(24*B*b**2 - 12*C*a*b + 5*D*a**2))/(315*a**9*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_168():
    f = (c*x**5 + d*x**7 + e*x**9 + f*x**11)/sqrt(a + b*x**2)
    F = a**2*sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**6 - a*(a + b*x**2)**(sympy.S(3)/2)*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)/(3*b**6) + f*(a + b*x**2)**(sympy.S(11)/2)/(11*b**6) + (a + b*x**2)**(sympy.S(9)/2)*(-5*a*f + b*e)/(9*b**6) + (a + b*x**2)**(sympy.S(7)/2)*(10*a**2*f - 4*a*b*e + b**2*d)/(7*b**6) + (a + b*x**2)**(sympy.S(5)/2)*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(5*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_169():
    f = (c*x**3 + d*x**5 + e*x**7 + f*x**9)/sqrt(a + b*x**2)
    F = -a*sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**5 + f*(a + b*x**2)**(sympy.S(9)/2)/(9*b**5) + (a + b*x**2)**(sympy.S(7)/2)*(-4*a*f + b*e)/(7*b**5) + (a + b*x**2)**(sympy.S(5)/2)*(6*a**2*f - 3*a*b*e + b**2*d)/(5*b**5) + (a + b*x**2)**(sympy.S(3)/2)*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_170():
    f = (c*x + d*x**3 + e*x**5 + f*x**7)/sqrt(a + b*x**2)
    F = f*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) + (a + b*x**2)**(sympy.S(5)/2)*(-3*a*f + b*e)/(5*b**4) + (a + b*x**2)**(sympy.S(3)/2)*(3*a**2*f - 2*a*b*e + b**2*d)/(3*b**4) + sqrt(a + b*x**2)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_171():
    f = x**2*(A + B*x**2 + C*x**4 + D*x**6 + F*x**8)/(a + b*x**2)**(sympy.S(9)/2)
    F = F*x*sqrt(a + b*x**2)/(2*b**5) - x*(D*b - 4*F*a)/(b**5*sqrt(a + b*x**2)) + (2*D*b - 9*F*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(11)/2)) + x**3*(A*b**4 - a*(B*b**3 - C*a*b**2 + D*a**2*b - F*a**3))/(7*a*b**4*(a + b*x**2)**(sympy.S(7)/2)) + x**3*(4*A*b**4 + a*(3*B*b**3 - 10*C*a*b**2 + 17*D*a**2*b - 24*F*a**3))/(35*a**2*b**4*(a + b*x**2)**(sympy.S(5)/2)) + x**3*(8*A*b**4 + a*(6*B*b**3 + 15*C*a*b**2 - 71*D*a**2*b + 162*F*a**3))/(105*a**3*b**4*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_8_P_x_c_x_pow_m_a_plus_b_x_pow_2_pow_p_172():
    f = (A + B*x**2 + C*x**4 + D*x**6 + F*x**8)/(x**2*(a + b*x**2)**(sympy.S(9)/2))
    F = -A/(a*x*(a + b*x**2)**(sympy.S(7)/2)) - x*(8*A*b - B*a)/(a**2*(a + b*x**2)**(sympy.S(7)/2)) - x**3*(48*A*b**2 - a*(6*B*b + C*a))/(3*a**3*(a + b*x**2)**(sympy.S(7)/2)) - x**5*(192*A*b**3 - a*(24*B*b**2 + 4*C*a*b + 3*D*a**2))/(15*a**4*(a + b*x**2)**(sympy.S(7)/2)) - x**7*(384*A*b**4 - a*(48*B*b**3 + 8*C*a*b**2 + 6*D*a**2*b + 15*F*a**3))/(105*a**5*(a + b*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F

