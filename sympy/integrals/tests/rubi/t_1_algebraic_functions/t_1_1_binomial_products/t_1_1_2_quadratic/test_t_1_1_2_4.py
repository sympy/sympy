"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.2 Quadratic/1.1.2.4 (e x)^m (a+b x^2)^p (c+d x^2)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, m, p, q = symbols('A B a b c d e m p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1():
    f = x**2*(A + B*x**2)*(a + b*x**2)
    F = A*a*x**3/3 + B*b*x**7/7 + x**5*(A*b + B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_2():
    f = x*(A + B*x**2)*(a + b*x**2)
    F = A*a*x**2/2 + B*b*x**6/6 + x**4*(A*b + B*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_3():
    f = (A + B*x**2)*(a + b*x**2)
    F = A*a*x + B*b*x**5/5 + x**3*(A*b/3 + B*a/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_4():
    f = (A + B*x**2)*(a + b*x**2)/x
    F = A*a*log(x) + B*b*x**4/4 + x**2*(A*b + B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_5():
    f = (A + B*x**2)*(a + b*x**2)/x**2
    F = -A*a/x + B*b*x**3/3 + x*(A*b + B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_6():
    f = (A + B*x**2)*(a + b*x**2)/x**3
    F = -A*a/(2*x**2) + B*b*x**2/2 + (A*b + B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_7():
    f = (A + B*x**2)*(a + b*x**2)/x**4
    F = -A*a/(3*x**3) + B*b*x - (A*b + B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_8():
    f = (A + B*x**2)*(a + b*x**2)/x**5
    F = -A*a/(4*x**4) + B*b*log(x) - (A*b + B*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_9():
    f = (A + B*x**2)*(a + b*x**2)/x**6
    F = -A*a/(5*x**5) - B*b/x - (A*b + B*a)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_10():
    f = (A + B*x**2)*(a + b*x**2)/x**7
    F = -A*a/(6*x**6) - B*b/(2*x**2) - (A*b + B*a)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_11():
    f = x**2*(A + B*x**2)*(a + b*x**2)**2
    F = A*a**2*x**3/3 + B*b**2*x**9/9 + a*x**5*(2*A*b + B*a)/5 + b*x**7*(A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_12():
    f = x*(A + B*x**2)*(a + b*x**2)**2
    F = B*(a + b*x**2)**4/(8*b**2) + (a + b*x**2)**3*(A*b - B*a)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_13():
    f = (A + B*x**2)*(a + b*x**2)**2
    F = A*a**2*x + B*b**2*x**7/7 + a*x**3*(2*A*b + B*a)/3 + b*x**5*(A*b + 2*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_14():
    f = (A + B*x**2)*(a + b*x**2)**2/x
    F = A*a**2*log(x) + A*a*b*x**2 + A*b**2*x**4/4 + B*(a + b*x**2)**3/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_15():
    f = (A + B*x**2)*(a + b*x**2)**2/x**2
    F = -A*a**2/x + B*b**2*x**5/5 + a*x*(2*A*b + B*a) + b*x**3*(A*b + 2*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_16():
    f = (A + B*x**2)*(a + b*x**2)**2/x**3
    F = -A*a**2/(2*x**2) + B*b**2*x**4/4 + a*(2*A*b + B*a)*log(x) + b*x**2*(A*b + 2*B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_17():
    f = (A + B*x**2)*(a + b*x**2)**2/x**4
    F = -A*a**2/(3*x**3) + B*b**2*x**3/3 - a*(2*A*b + B*a)/x + b*x*(A*b + 2*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_18():
    f = (A + B*x**2)*(a + b*x**2)**2/x**5
    F = -A*a**2/(4*x**4) + B*b**2*x**2/2 - a*(2*A*b + B*a)/(2*x**2) + b*(A*b + 2*B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_19():
    f = (A + B*x**2)*(a + b*x**2)**2/x**6
    F = -A*a**2/(5*x**5) + B*b**2*x - a*(2*A*b + B*a)/(3*x**3) - b*(A*b + 2*B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_20():
    f = (A + B*x**2)*(a + b*x**2)**2/x**7
    F = -A*a**2/(6*x**6) + B*b**2*log(x) - a*(2*A*b + B*a)/(4*x**4) - b*(A*b + 2*B*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_21():
    f = (A + B*x**2)*(a + b*x**2)**2/x**8
    F = -A*a**2/(7*x**7) - B*b**2/x - a*(2*A*b + B*a)/(5*x**5) - b*(A*b + 2*B*a)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_22():
    f = (A + B*x**2)*(a + b*x**2)**2/x**9
    F = -A*(a + b*x**2)**3/(8*a*x**8) + (a + b*x**2)**3*(A*b - 4*B*a)/(24*a**2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_23():
    f = x**9*(A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x**10/10 + B*b**5*x**22/22 + a**4*x**12*(5*A*b + B*a)/12 + 5*a**3*b*x**14*(2*A*b + B*a)/14 + 5*a**2*b**2*x**16*(A*b + B*a)/8 + 5*a*b**3*x**18*(A*b + 2*B*a)/18 + b**4*x**20*(A*b + 5*B*a)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_24():
    f = x**8*(A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x**9/9 + B*b**5*x**21/21 + a**4*x**11*(5*A*b + B*a)/11 + 5*a**3*b*x**13*(2*A*b + B*a)/13 + 2*a**2*b**2*x**15*(A*b + B*a)/3 + 5*a*b**3*x**17*(A*b + 2*B*a)/17 + b**4*x**19*(A*b + 5*B*a)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_25():
    f = x**7*(A + B*x**2)*(a + b*x**2)**5
    F = B*(a + b*x**2)**10/(20*b**5) - a**3*(a + b*x**2)**6*(A*b - B*a)/(12*b**5) + a**2*(a + b*x**2)**7*(3*A*b - 4*B*a)/(14*b**5) - 3*a*(a + b*x**2)**8*(A*b - 2*B*a)/(16*b**5) + (a + b*x**2)**9*(A*b - 4*B*a)/(18*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_26():
    f = x**6*(A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x**7/7 + B*b**5*x**19/19 + a**4*x**9*(5*A*b + B*a)/9 + 5*a**3*b*x**11*(2*A*b + B*a)/11 + 10*a**2*b**2*x**13*(A*b + B*a)/13 + a*b**3*x**15*(A*b + 2*B*a)/3 + b**4*x**17*(A*b + 5*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_27():
    f = x**5*(A + B*x**2)*(a + b*x**2)**5
    F = B*(a + b*x**2)**9/(18*b**4) + a**2*(a + b*x**2)**6*(A*b - B*a)/(12*b**4) - a*(a + b*x**2)**7*(2*A*b - 3*B*a)/(14*b**4) + (a + b*x**2)**8*(A*b - 3*B*a)/(16*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_28():
    f = x**4*(A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x**5/5 + B*b**5*x**17/17 + a**4*x**7*(5*A*b + B*a)/7 + 5*a**3*b*x**9*(2*A*b + B*a)/9 + 10*a**2*b**2*x**11*(A*b + B*a)/11 + 5*a*b**3*x**13*(A*b + 2*B*a)/13 + b**4*x**15*(A*b + 5*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_29():
    f = x**3*(A + B*x**2)*(a + b*x**2)**5
    F = B*(a + b*x**2)**8/(16*b**3) - a*(a + b*x**2)**6*(A*b - B*a)/(12*b**3) + (a + b*x**2)**7*(A*b - 2*B*a)/(14*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_30():
    f = x**2*(A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x**3/3 + B*b**5*x**15/15 + a**4*x**5*(5*A*b + B*a)/5 + 5*a**3*b*x**7*(2*A*b + B*a)/7 + 10*a**2*b**2*x**9*(A*b + B*a)/9 + 5*a*b**3*x**11*(A*b + 2*B*a)/11 + b**4*x**13*(A*b + 5*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_31():
    f = x*(A + B*x**2)*(a + b*x**2)**5
    F = B*(a + b*x**2)**7/(14*b**2) + (a + b*x**2)**6*(A*b - B*a)/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_32():
    f = (A + B*x**2)*(a + b*x**2)**5
    F = A*a**5*x + B*b**5*x**13/13 + a**4*x**3*(5*A*b + B*a)/3 + a**3*b*x**5*(2*A*b + B*a) + 10*a**2*b**2*x**7*(A*b + B*a)/7 + 5*a*b**3*x**9*(A*b + 2*B*a)/9 + b**4*x**11*(A*b + 5*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_33():
    f = (A + B*x**2)*(a + b*x**2)**5/x
    F = A*a**5*log(x) + 5*A*a**4*b*x**2/2 + 5*A*a**3*b**2*x**4/2 + 5*A*a**2*b**3*x**6/3 + 5*A*a*b**4*x**8/8 + A*b**5*x**10/10 + B*(a + b*x**2)**6/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_34():
    f = (A + B*x**2)*(a + b*x**2)**5/x**2
    F = -A*a**5/x + B*b**5*x**11/11 + a**4*x*(5*A*b + B*a) + 5*a**3*b*x**3*(2*A*b + B*a)/3 + 2*a**2*b**2*x**5*(A*b + B*a) + 5*a*b**3*x**7*(A*b + 2*B*a)/7 + b**4*x**9*(A*b + 5*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_35():
    f = (A + B*x**2)*(a + b*x**2)**5/x**3
    F = -A*a**5/(2*x**2) + B*b**5*x**10/10 + a**4*(5*A*b + B*a)*log(x) + 5*a**3*b*x**2*(2*A*b + B*a)/2 + 5*a**2*b**2*x**4*(A*b + B*a)/2 + 5*a*b**3*x**6*(A*b + 2*B*a)/6 + b**4*x**8*(A*b + 5*B*a)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_36():
    f = (A + B*x**2)*(a + b*x**2)**5/x**4
    F = -A*a**5/(3*x**3) + B*b**5*x**9/9 - a**4*(5*A*b + B*a)/x + 5*a**3*b*x*(2*A*b + B*a) + 10*a**2*b**2*x**3*(A*b + B*a)/3 + a*b**3*x**5*(A*b + 2*B*a) + b**4*x**7*(A*b + 5*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_37():
    f = (A + B*x**2)*(a + b*x**2)**5/x**5
    F = -A*a**5/(4*x**4) + B*b**5*x**8/8 - a**4*(5*A*b + B*a)/(2*x**2) + 5*a**3*b*(2*A*b + B*a)*log(x) + 5*a**2*b**2*x**2*(A*b + B*a) + 5*a*b**3*x**4*(A*b + 2*B*a)/4 + b**4*x**6*(A*b + 5*B*a)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_38():
    f = (A + B*x**2)*(a + b*x**2)**5/x**6
    F = -A*a**5/(5*x**5) + B*b**5*x**7/7 - a**4*(5*A*b + B*a)/(3*x**3) - 5*a**3*b*(2*A*b + B*a)/x + 10*a**2*b**2*x*(A*b + B*a) + 5*a*b**3*x**3*(A*b + 2*B*a)/3 + b**4*x**5*(A*b + 5*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_39():
    f = (A + B*x**2)*(a + b*x**2)**5/x**7
    F = -A*a**5/(6*x**6) + B*b**5*x**6/6 - a**4*(5*A*b + B*a)/(4*x**4) - 5*a**3*b*(2*A*b + B*a)/(2*x**2) + 10*a**2*b**2*(A*b + B*a)*log(x) + 5*a*b**3*x**2*(A*b + 2*B*a)/2 + b**4*x**4*(A*b + 5*B*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_40():
    f = (A + B*x**2)*(a + b*x**2)**5/x**8
    F = -A*a**5/(7*x**7) + B*b**5*x**5/5 - a**4*(5*A*b + B*a)/(5*x**5) - 5*a**3*b*(2*A*b + B*a)/(3*x**3) - 10*a**2*b**2*(A*b + B*a)/x + 5*a*b**3*x*(A*b + 2*B*a) + b**4*x**3*(A*b + 5*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_41():
    f = (A + B*x**2)*(a + b*x**2)**5/x**9
    F = -A*a**5/(8*x**8) + B*b**5*x**4/4 - a**4*(5*A*b + B*a)/(6*x**6) - 5*a**3*b*(2*A*b + B*a)/(4*x**4) - 5*a**2*b**2*(A*b + B*a)/x**2 + 5*a*b**3*(A*b + 2*B*a)*log(x) + b**4*x**2*(A*b + 5*B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_42():
    f = (A + B*x**2)*(a + b*x**2)**5/x**10
    F = -A*a**5/(9*x**9) + B*b**5*x**3/3 - a**4*(5*A*b + B*a)/(7*x**7) - a**3*b*(2*A*b + B*a)/x**5 - 10*a**2*b**2*(A*b + B*a)/(3*x**3) - 5*a*b**3*(A*b + 2*B*a)/x + b**4*x*(A*b + 5*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_43():
    f = (A + B*x**2)*(a + b*x**2)**5/x**11
    F = -A*a**5/(10*x**10) + B*b**5*x**2/2 - a**4*(5*A*b + B*a)/(8*x**8) - 5*a**3*b*(2*A*b + B*a)/(6*x**6) - 5*a**2*b**2*(A*b + B*a)/(2*x**4) - 5*a*b**3*(A*b + 2*B*a)/(2*x**2) + b**4*(A*b + 5*B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_44():
    f = (A + B*x**2)*(a + b*x**2)**5/x**12
    F = -A*a**5/(11*x**11) + B*b**5*x - a**4*(5*A*b + B*a)/(9*x**9) - 5*a**3*b*(2*A*b + B*a)/(7*x**7) - 2*a**2*b**2*(A*b + B*a)/x**5 - 5*a*b**3*(A*b + 2*B*a)/(3*x**3) - b**4*(A*b + 5*B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_45():
    f = (A + B*x**2)*(a + b*x**2)**5/x**13
    F = -A*(a + b*x**2)**6/(12*a*x**12) - B*a**5/(10*x**10) - 5*B*a**4*b/(8*x**8) - 5*B*a**3*b**2/(3*x**6) - 5*B*a**2*b**3/(2*x**4) - 5*B*a*b**4/(2*x**2) + B*b**5*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_46():
    f = (A + B*x**2)*(a + b*x**2)**5/x**14
    F = -A*a**5/(13*x**13) - B*b**5/x - a**4*(5*A*b + B*a)/(11*x**11) - 5*a**3*b*(2*A*b + B*a)/(9*x**9) - 10*a**2*b**2*(A*b + B*a)/(7*x**7) - a*b**3*(A*b + 2*B*a)/x**5 - b**4*(A*b + 5*B*a)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_47():
    f = (A + B*x**2)*(a + b*x**2)**5/x**15
    F = -A*(a + b*x**2)**6/(14*a*x**14) + (a + b*x**2)**6*(A*b - 7*B*a)/(84*a**2*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_48():
    f = (A + B*x**2)*(a + b*x**2)**5/x**16
    F = -A*a**5/(15*x**15) - B*b**5/(3*x**3) - a**4*(5*A*b + B*a)/(13*x**13) - 5*a**3*b*(2*A*b + B*a)/(11*x**11) - 10*a**2*b**2*(A*b + B*a)/(9*x**9) - 5*a*b**3*(A*b + 2*B*a)/(7*x**7) - b**4*(A*b + 5*B*a)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_49():
    f = (A + B*x**2)*(a + b*x**2)**5/x**17
    F = -A*(a + b*x**2)**6/(16*a*x**16) + (a + b*x**2)**6*(A*b - 4*B*a)/(56*a**2*x**14) - b*(a + b*x**2)**6*(A*b - 4*B*a)/(336*a**3*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_50():
    f = (A + B*x**2)*(a + b*x**2)**5/x**18
    F = -A*a**5/(17*x**17) - B*b**5/(5*x**5) - a**4*(5*A*b + B*a)/(15*x**15) - 5*a**3*b*(2*A*b + B*a)/(13*x**13) - 10*a**2*b**2*(A*b + B*a)/(11*x**11) - 5*a*b**3*(A*b + 2*B*a)/(9*x**9) - b**4*(A*b + 5*B*a)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_51():
    f = (A + B*x**2)*(a + b*x**2)**5/x**19
    F = -A*a**5/(18*x**18) - B*b**5/(6*x**6) - a**4*(5*A*b + B*a)/(16*x**16) - 5*a**3*b*(2*A*b + B*a)/(14*x**14) - 5*a**2*b**2*(A*b + B*a)/(6*x**12) - a*b**3*(A*b + 2*B*a)/(2*x**10) - b**4*(A*b + 5*B*a)/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_52():
    f = (A + B*x**2)*(a + b*x**2)**5/x**20
    F = -A*a**5/(19*x**19) - B*b**5/(7*x**7) - a**4*(5*A*b + B*a)/(17*x**17) - a**3*b*(2*A*b + B*a)/(3*x**15) - 10*a**2*b**2*(A*b + B*a)/(13*x**13) - 5*a*b**3*(A*b + 2*B*a)/(11*x**11) - b**4*(A*b + 5*B*a)/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_53():
    f = (A + B*x**2)*(a + b*x**2)**5/x**21
    F = -A*a**5/(20*x**20) - B*b**5/(8*x**8) - a**4*(5*A*b + B*a)/(18*x**18) - 5*a**3*b*(2*A*b + B*a)/(16*x**16) - 5*a**2*b**2*(A*b + B*a)/(7*x**14) - 5*a*b**3*(A*b + 2*B*a)/(12*x**12) - b**4*(A*b + 5*B*a)/(10*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_54():
    f = (A + B*x**2)*(a + b*x**2)**5/x**22
    F = -A*a**5/(21*x**21) - B*b**5/(9*x**9) - a**4*(5*A*b + B*a)/(19*x**19) - 5*a**3*b*(2*A*b + B*a)/(17*x**17) - 2*a**2*b**2*(A*b + B*a)/(3*x**15) - 5*a*b**3*(A*b + 2*B*a)/(13*x**13) - b**4*(A*b + 5*B*a)/(11*x**11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_55():
    f = (A + B*x**2)*(a + b*x**2)**5/x**23
    F = -A*a**5/(22*x**22) - B*b**5/(10*x**10) - a**4*(5*A*b + B*a)/(20*x**20) - 5*a**3*b*(2*A*b + B*a)/(18*x**18) - 5*a**2*b**2*(A*b + B*a)/(8*x**16) - 5*a*b**3*(A*b + 2*B*a)/(14*x**14) - b**4*(A*b + 5*B*a)/(12*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_56():
    f = x**6*(A + B*x**2)/(a + b*x**2)
    F = B*x**7/(7*b) - a**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(9)/2) + a**2*x*(A*b - B*a)/b**4 - a*x**3*(A*b - B*a)/(3*b**3) + x**5*(A*b - B*a)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_57():
    f = x**5*(A + B*x**2)/(a + b*x**2)
    F = B*x**6/(6*b) + a**2*(A*b - B*a)*log(a + b*x**2)/(2*b**4) - a*x**2*(A*b - B*a)/(2*b**3) + x**4*(A*b - B*a)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_58():
    f = x**4*(A + B*x**2)/(a + b*x**2)
    F = B*x**5/(5*b) + a**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) - a*x*(A*b - B*a)/b**3 + x**3*(A*b - B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_59():
    f = x**3*(A + B*x**2)/(a + b*x**2)
    F = B*x**4/(4*b) - a*(A*b - B*a)*log(a + b*x**2)/(2*b**3) + x**2*(A*b - B*a)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_60():
    f = x**2*(A + B*x**2)/(a + b*x**2)
    F = B*x**3/(3*b) - sqrt(a)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(5)/2) + x*(A*b - B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_61():
    f = x*(A + B*x**2)/(a + b*x**2)
    F = B*x**2/(2*b) + (A*b - B*a)*log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_62():
    f = (A + B*x**2)/(a + b*x**2)
    F = B*x/b + (A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_63():
    f = (A + B*x**2)/(x*(a + b*x**2))
    F = A*log(x)/a - (A*b - B*a)*log(a + b*x**2)/(2*a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_64():
    f = (A + B*x**2)/(x**2*(a + b*x**2))
    F = -A/(a*x) - (A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_65():
    f = (A + B*x**2)/(x**3*(a + b*x**2))
    F = -A/(2*a*x**2) - (A*b - B*a)*log(x)/a**2 + (A*b - B*a)*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_66():
    f = (A + B*x**2)/(x**4*(a + b*x**2))
    F = -A/(3*a*x**3) + (A*b - B*a)/(a**2*x) + sqrt(b)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_67():
    f = (A + B*x**2)/(x**5*(a + b*x**2))
    F = -A/(4*a*x**4) + (A*b - B*a)/(2*a**2*x**2) + b*(A*b - B*a)*log(x)/a**3 - b*(A*b - B*a)*log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_68():
    f = (A + B*x**2)/(x**6*(a + b*x**2))
    F = -A/(5*a*x**5) + (A*b - B*a)/(3*a**2*x**3) - b*(A*b - B*a)/(a**3*x) - b**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_69():
    f = (A + B*x**2)/(x**7*(a + b*x**2))
    F = -A/(6*a*x**6) + (A*b - B*a)/(4*a**2*x**4) - b*(A*b - B*a)/(2*a**3*x**2) - b**2*(A*b - B*a)*log(x)/a**4 + b**2*(A*b - B*a)*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_70():
    f = (A + B*x**2)/(x**8*(a + b*x**2))
    F = -A/(7*a*x**7) + (A*b - B*a)/(5*a**2*x**5) - b*(A*b - B*a)/(3*a**3*x**3) + b**2*(A*b - B*a)/(a**4*x) + b**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_71():
    f = x**9*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**8/(8*b**2) - a**4*(A*b - B*a)/(2*b**6*(a + b*x**2)) - a**3*(4*A*b - 5*B*a)*log(a + b*x**2)/(2*b**6) + a**2*x**2*(3*A*b - 4*B*a)/(2*b**5) - a*x**4*(2*A*b - 3*B*a)/(4*b**4) + x**6*(A*b - 2*B*a)/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_72():
    f = x**8*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**7/(7*b**2) - a**(sympy.S(5)/2)*(7*A*b - 9*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(11)/2)) + a**3*x*(A*b - B*a)/(2*b**5*(a + b*x**2)) + a**2*x*(3*A*b - 4*B*a)/b**5 - a*x**3*(2*A*b - 3*B*a)/(3*b**4) + x**5*(A*b - 2*B*a)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_73():
    f = x**7*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**6/(6*b**2) + a**3*(A*b - B*a)/(2*b**5*(a + b*x**2)) + a**2*(3*A*b - 4*B*a)*log(a + b*x**2)/(2*b**5) - a*x**2*(2*A*b - 3*B*a)/(2*b**4) + x**4*(A*b - 2*B*a)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_74():
    f = x**6*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**5/(5*b**2) + a**(sympy.S(3)/2)*(5*A*b - 7*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(9)/2)) - a**2*x*(A*b - B*a)/(2*b**4*(a + b*x**2)) - a*x*(2*A*b - 3*B*a)/b**4 + x**3*(A*b - 2*B*a)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_75():
    f = x**5*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**4/(4*b**2) - a**2*(A*b - B*a)/(2*b**4*(a + b*x**2)) - a*(2*A*b - 3*B*a)*log(a + b*x**2)/(2*b**4) + x**2*(A*b - 2*B*a)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_76():
    f = x**4*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**3/(3*b**2) - sqrt(a)*(3*A*b - 5*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) + a*x*(A*b - B*a)/(2*b**3*(a + b*x**2)) + x*(A*b - 2*B*a)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_77():
    f = x**3*(A + B*x**2)/(a + b*x**2)**2
    F = B*x**2/(2*b**2) + a*(A*b - B*a)/(2*b**3*(a + b*x**2)) + (A*b - 2*B*a)*log(a + b*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_78():
    f = x**2*(A + B*x**2)/(a + b*x**2)**2
    F = B*x/b**2 - x*(A*b - B*a)/(2*b**2*(a + b*x**2)) + (A*b - 3*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_79():
    f = x*(A + B*x**2)/(a + b*x**2)**2
    F = B*log(a + b*x**2)/(2*b**2) - (A*b - B*a)/(2*b**2*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_80():
    f = (A + B*x**2)/(a + b*x**2)**2
    F = x*(A*b - B*a)/(2*a*b*(a + b*x**2)) + (A*b + B*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_81():
    f = (A + B*x**2)/(x*(a + b*x**2)**2)
    F = A*log(x)/a**2 - A*log(a + b*x**2)/(2*a**2) + (A*b - B*a)/(2*a*b*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_82():
    f = (A + B*x**2)/(x**2*(a + b*x**2)**2)
    F = -A/(a**2*x) - x*(A*b - B*a)/(2*a**2*(a + b*x**2)) - (3*A*b - B*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_83():
    f = (A + B*x**2)/(x**3*(a + b*x**2)**2)
    F = -A/(2*a**2*x**2) - (A*b - B*a)/(2*a**2*(a + b*x**2)) - (2*A*b - B*a)*log(x)/a**3 + (2*A*b - B*a)*log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_84():
    f = (A + B*x**2)/(x**4*(a + b*x**2)**2)
    F = -A/(3*a**2*x**3) + b*x*(A*b - B*a)/(2*a**3*(a + b*x**2)) + (2*A*b - B*a)/(a**3*x) + sqrt(b)*(5*A*b - 3*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_85():
    f = (A + B*x**2)/(x**5*(a + b*x**2)**2)
    F = -A/(4*a**2*x**4) + b*(A*b - B*a)/(2*a**3*(a + b*x**2)) + (2*A*b - B*a)/(2*a**3*x**2) + b*(3*A*b - 2*B*a)*log(x)/a**4 - b*(3*A*b - 2*B*a)*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_86():
    f = (A + B*x**2)/(x**6*(a + b*x**2)**2)
    F = -A/(5*a**2*x**5) + (2*A*b - B*a)/(3*a**3*x**3) - b**2*x*(A*b - B*a)/(2*a**4*(a + b*x**2)) - b*(3*A*b - 2*B*a)/(a**4*x) - b**(sympy.S(3)/2)*(7*A*b - 5*B*a)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_87():
    f = (A + B*x**2)/(x**7*(a + b*x**2)**2)
    F = -A/(6*a**2*x**6) + (2*A*b - B*a)/(4*a**3*x**4) - b**2*(A*b - B*a)/(2*a**4*(a + b*x**2)) - b*(3*A*b - 2*B*a)/(2*a**4*x**2) - b**2*(4*A*b - 3*B*a)*log(x)/a**5 + b**2*(4*A*b - 3*B*a)*log(a + b*x**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_88():
    f = x**11*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**8/(8*b**3) + a**5*(A*b - B*a)/(4*b**7*(a + b*x**2)**2) - a**4*(5*A*b - 6*B*a)/(2*b**7*(a + b*x**2)) - 5*a**3*(2*A*b - 3*B*a)*log(a + b*x**2)/(2*b**7) + a**2*x**2*(3*A*b - 5*B*a)/b**6 - 3*a*x**4*(A*b - 2*B*a)/(4*b**5) + x**6*(A*b - 3*B*a)/(6*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_89():
    f = x**9*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**6/(6*b**3) - a**4*(A*b - B*a)/(4*b**6*(a + b*x**2)**2) + a**3*(4*A*b - 5*B*a)/(2*b**6*(a + b*x**2)) + a**2*(3*A*b - 5*B*a)*log(a + b*x**2)/b**6 - 3*a*x**2*(A*b - 2*B*a)/(2*b**5) + x**4*(A*b - 3*B*a)/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_90():
    f = x**7*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**4/(4*b**3) + a**3*(A*b - B*a)/(4*b**5*(a + b*x**2)**2) - a**2*(3*A*b - 4*B*a)/(2*b**5*(a + b*x**2)) - 3*a*(A*b - 2*B*a)*log(a + b*x**2)/(2*b**5) + x**2*(A*b - 3*B*a)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_91():
    f = x**5*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**2/(2*b**3) - a**2*(A*b - B*a)/(4*b**4*(a + b*x**2)**2) + a*(2*A*b - 3*B*a)/(2*b**4*(a + b*x**2)) + (A*b - 3*B*a)*log(a + b*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_92():
    f = x**3*(A + B*x**2)/(a + b*x**2)**3
    F = B*log(a + b*x**2)/(2*b**3) + a*(A*b - B*a)/(4*b**3*(a + b*x**2)**2) - (A*b - 2*B*a)/(2*b**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_93():
    f = x*(A + B*x**2)/(a + b*x**2)**3
    F = -(A + B*x**2)**2/((a + b*x**2)**2*(4*A*b - 4*B*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_94():
    f = (A + B*x**2)/(x*(a + b*x**2)**3)
    F = A/(2*a**2*(a + b*x**2)) + A*log(x)/a**3 - A*log(a + b*x**2)/(2*a**3) + (A*b - B*a)/(4*a*b*(a + b*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_95():
    f = (A + B*x**2)/(x**3*(a + b*x**2)**3)
    F = -A/(2*a**3*x**2) - (A*b - B*a)/(4*a**2*(a + b*x**2)**2) - (2*A*b - B*a)/(2*a**3*(a + b*x**2)) - (3*A*b - B*a)*log(x)/a**4 + (3*A*b - B*a)*log(a + b*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_96():
    f = (A + B*x**2)/(x**5*(a + b*x**2)**3)
    F = -A/(4*a**3*x**4) + b*(A*b - B*a)/(4*a**3*(a + b*x**2)**2) + b*(3*A*b - 2*B*a)/(2*a**4*(a + b*x**2)) + (3*A*b - B*a)/(2*a**4*x**2) + 3*b*(2*A*b - B*a)*log(x)/a**5 - 3*b*(2*A*b - B*a)*log(a + b*x**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_97():
    f = (A + B*x**2)/(x**7*(a + b*x**2)**3)
    F = -A/(6*a**3*x**6) - b**2*(A*b - B*a)/(4*a**4*(a + b*x**2)**2) + (3*A*b - B*a)/(4*a**4*x**4) - b**2*(4*A*b - 3*B*a)/(2*a**5*(a + b*x**2)) - 3*b*(2*A*b - B*a)/(2*a**5*x**2) - 2*b**2*(5*A*b - 3*B*a)*log(x)/a**6 + b**2*(5*A*b - 3*B*a)*log(a + b*x**2)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_98():
    f = x**10*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**7/(7*b**3) - 9*a**(sympy.S(5)/2)*(7*A*b - 11*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(13)/2)) - a**4*x*(A*b - B*a)/(4*b**6*(a + b*x**2)**2) + a**3*x*(17*A*b - 21*B*a)/(8*b**6*(a + b*x**2)) + 2*a**2*x*(3*A*b - 5*B*a)/b**6 - a*x**3*(A*b - 2*B*a)/b**5 + x**5*(A*b - 3*B*a)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_99():
    f = x**8*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**5/(5*b**3) + 7*a**(sympy.S(3)/2)*(5*A*b - 9*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(11)/2)) + a**3*x*(A*b - B*a)/(4*b**5*(a + b*x**2)**2) - a**2*x*(13*A*b - 17*B*a)/(8*b**5*(a + b*x**2)) - 3*a*x*(A*b - 2*B*a)/b**5 + x**3*(A*b - 3*B*a)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_100():
    f = x**6*(A + B*x**2)/(a + b*x**2)**3
    F = B*x**3/(3*b**3) - 5*sqrt(a)*(3*A*b - 7*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*b**(sympy.S(9)/2)) - a**2*x*(A*b - B*a)/(4*b**4*(a + b*x**2)**2) + a*x*(9*A*b - 13*B*a)/(8*b**4*(a + b*x**2)) + x*(A*b - 3*B*a)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_101():
    f = x**4*(A + B*x**2)/(a + b*x**2)**3
    F = B*x/b**3 + a*x*(A*b - B*a)/(4*b**3*(a + b*x**2)**2) - x*(5*A*b - 9*B*a)/(8*b**3*(a + b*x**2)) + (3*A*b - 15*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_102():
    f = x**2*(A + B*x**2)/(a + b*x**2)**3
    F = -x*(A*b - B*a)/(4*b**2*(a + b*x**2)**2) + x*(A*b - 5*B*a)/(8*a*b**2*(a + b*x**2)) + (A*b + 3*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_103():
    f = (A + B*x**2)/(a + b*x**2)**3
    F = x*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) + x*(3*A*b + B*a)/(8*a**2*b*(a + b*x**2)) + (3*A*b + B*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_104():
    f = (A + B*x**2)/(x**2*(a + b*x**2)**3)
    F = -A/(a**3*x) - x*(A*b - B*a)/(4*a**2*(a + b*x**2)**2) - x*(7*A*b - 3*B*a)/(8*a**3*(a + b*x**2)) - (15*A*b - 3*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_105():
    f = (A + B*x**2)/(x**4*(a + b*x**2)**3)
    F = -A/(3*a**3*x**3) + b*x*(A*b - B*a)/(4*a**3*(a + b*x**2)**2) + b*x*(11*A*b - 7*B*a)/(8*a**4*(a + b*x**2)) + (3*A*b - B*a)/(a**4*x) + 5*sqrt(b)*(7*A*b - 3*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_106():
    f = (A + B*x**2)/(x**6*(a + b*x**2)**3)
    F = -A/(5*a**3*x**5) - b**2*x*(A*b - B*a)/(4*a**4*(a + b*x**2)**2) + (3*A*b - B*a)/(3*a**4*x**3) - b**2*x*(15*A*b - 11*B*a)/(8*a**5*(a + b*x**2)) - 3*b*(2*A*b - B*a)/(a**5*x) - 7*b**(sympy.S(3)/2)*(9*A*b - 5*B*a)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_107():
    f = (a + b*x**2)/(x**2 + 1)
    F = b*x + (a - b)*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_108():
    f = (a + b*x**2)/(1 - x**2)
    F = -b*x + (a + b)*atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_109():
    f = (x**2 + 1)/(x**2 - 1)**2
    F = x/(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_110():
    f = (1 - x**2)/(x**2 + 1)**2
    F = x/(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_111():
    f = (2*x**2 + 3)/(x**2 + 1)**2
    F = x/(2*x**2 + 2) + 5*atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_112():
    f = (x**2 - 2)/(x**2 + 1)**2
    F = -3*x/(2*x**2 + 2) - atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_113():
    f = (x**2 + 3)/(x**2 + 1)**2
    F = x/(x**2 + 1) + 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_114():
    f = (a + b*x**2)/(-a + b*x**2)**2
    F = x/(a - b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_115():
    f = (a + b*x**2)/(a - b*x**2)**2
    F = x/(a - b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_116():
    f = (A + B*x**2)/(a - b*x**2)
    F = -B*x/b + (A*b + B*a)*atanh(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_117():
    f = (x**2 + 1)/(x**2 + 16)**3
    F = 19*x/(2048*x**2 + 32768) - 15*x/(64*(x**2 + 16)**2) + 19*atan(x/4)/8192
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_118():
    f = (2*x**2 + 1)/(x**5*(x**2 + 1)**3)
    F = -1/(4*x**4*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_119():
    f = (1 - x**2)**2/(x**2 - 1)**2
    F = x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_120():
    f = x**3*(a*c + b*c*x**2)/(a + b*x**2)
    F = c*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_121():
    f = x**2*(a*c + b*c*x**2)/(a + b*x**2)
    F = c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_122():
    f = x*(a*c + b*c*x**2)/(a + b*x**2)
    F = c*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_123():
    f = (a*c + b*c*x**2)/(a + b*x**2)
    F = c*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_124():
    f = (a*c + b*c*x**2)/(x*(a + b*x**2))
    F = c*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_125():
    f = (a*c + b*c*x**2)/(x**2*(a + b*x**2))
    F = -c/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_126():
    f = (a*c + b*c*x**2)/(x**3*(a + b*x**2))
    F = -c/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_127():
    f = x**3*(a*c + b*c*x**2)/(a + b*x**2)**2
    F = -a*c*log(a + b*x**2)/(2*b**2) + c*x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_128():
    f = x**2*(a*c + b*c*x**2)/(a + b*x**2)**2
    F = -sqrt(a)*c*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(3)/2) + c*x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_129():
    f = x*(a*c + b*c*x**2)/(a + b*x**2)**2
    F = c*log(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_130():
    f = (a*c + b*c*x**2)/(a + b*x**2)**2
    F = c*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_131():
    f = (a*c + b*c*x**2)/(x*(a + b*x**2)**2)
    F = c*log(x)/a - c*log(a + b*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_132():
    f = (a*c + b*c*x**2)/(x**2*(a + b*x**2)**2)
    F = -c/(a*x) - sqrt(b)*c*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_133():
    f = (a*c + b*c*x**2)/(x**3*(a + b*x**2)**2)
    F = -c/(2*a*x**2) - b*c*log(x)/a**2 + b*c*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_134():
    f = x**3*(a*c + b*c*x**2)/(a + b*x**2)**3
    F = a*c/(2*b**2*(a + b*x**2)) + c*log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_135():
    f = x**2*(a*c + b*c*x**2)/(a + b*x**2)**3
    F = -c*x/(2*b*(a + b*x**2)) + c*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_136():
    f = x*(a*c + b*c*x**2)/(a + b*x**2)**3
    F = -c/(2*b*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_137():
    f = (a*c + b*c*x**2)/(a + b*x**2)**3
    F = c*x/(2*a*(a + b*x**2)) + c*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_138():
    f = (a*c + b*c*x**2)/(x*(a + b*x**2)**3)
    F = c/(2*a*(a + b*x**2)) + c*log(x)/a**2 - c*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_139():
    f = (a*c + b*c*x**2)/(x**2*(a + b*x**2)**3)
    F = c/(2*a*x*(a + b*x**2)) - 3*c/(2*a**2*x) - 3*sqrt(b)*c*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_140():
    f = (a*c + b*c*x**2)/(x**3*(a + b*x**2)**3)
    F = -b*c/(2*a**2*(a + b*x**2)) - c/(2*a**2*x**2) - 2*b*c*log(x)/a**3 + b*c*log(a + b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_141():
    f = x**4*(a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x**5/5 + a*x**7*(a*d + 2*b*c)/7 + b**2*d*x**11/11 + b*x**9*(2*a*d + b*c)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_142():
    f = x**3*(a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x**4/4 + a*x**6*(a*d + 2*b*c)/6 + b**2*d*x**10/10 + b*x**8*(2*a*d + b*c)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_143():
    f = x**2*(a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x**3/3 + a*x**5*(a*d + 2*b*c)/5 + b**2*d*x**9/9 + b*x**7*(2*a*d + b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_144():
    f = x*(a + b*x**2)**2*(c + d*x**2)
    F = d*(a + b*x**2)**4/(8*b**2) + (a + b*x**2)**3*(-a*d + b*c)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_145():
    f = (a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x + a*x**3*(a*d + 2*b*c)/3 + b**2*d*x**7/7 + b*x**5*(2*a*d + b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_146():
    f = (a + b*x**2)**2*(c + d*x**2)/x
    F = a**2*c*log(x) + a*b*c*x**2 + b**2*c*x**4/4 + d*(a + b*x**2)**3/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_147():
    f = (a + b*x**2)**2*(c + d*x**2)/x**2
    F = -a**2*c/x + a*x*(a*d + 2*b*c) + b**2*d*x**5/5 + b*x**3*(2*a*d + b*c)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_148():
    f = (a + b*x**2)**2*(c + d*x**2)/x**3
    F = -a**2*c/(2*x**2) + a*(a*d + 2*b*c)*log(x) + b**2*d*x**4/4 + b*x**2*(2*a*d + b*c)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_149():
    f = (a + b*x**2)**2*(c + d*x**2)/x**4
    F = -a**2*c/(3*x**3) - a*(a*d + 2*b*c)/x + b**2*d*x**3/3 + b*x*(2*a*d + b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_150():
    f = x**4*(a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x**5/5 + 2*a*c*x**7*(a*d + b*c)/7 + b**2*d**2*x**13/13 + 2*b*d*x**11*(a*d + b*c)/11 + x**9*(a**2*d**2/9 + 4*a*b*c*d/9 + b**2*c**2/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_151():
    f = x**3*(a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x**4/4 + a*c*x**6*(a*d + b*c)/3 + b**2*d**2*x**12/12 + b*d*x**10*(a*d + b*c)/5 + x**8*(a**2*d**2/8 + a*b*c*d/2 + b**2*c**2/8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_152():
    f = x**2*(a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x**3/3 + 2*a*c*x**5*(a*d + b*c)/5 + b**2*d**2*x**11/11 + 2*b*d*x**9*(a*d + b*c)/9 + x**7*(a**2*d**2/7 + 4*a*b*c*d/7 + b**2*c**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_153():
    f = x*(a + b*x**2)**2*(c + d*x**2)**2
    F = d**2*(a + b*x**2)**5/(10*b**3) + d*(a + b*x**2)**4*(-a*d + b*c)/(4*b**3) + (a + b*x**2)**3*(-a*d + b*c)**2/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_154():
    f = (a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x + 2*a*c*x**3*(a*d + b*c)/3 + b**2*d**2*x**9/9 + 2*b*d*x**7*(a*d + b*c)/7 + x**5*(a**2*d**2 + 4*a*b*c*d + b**2*c**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_155():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x
    F = a**2*c**2*log(x) + a*c*x**2*(a*d + b*c) + b**2*d**2*x**8/8 + b*d*x**6*(a*d + b*c)/3 + x**4*(a**2*d**2/4 + a*b*c*d + b**2*c**2/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_156():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**2
    F = -a**2*c**2/x + 2*a*c*x*(a*d + b*c) + b**2*d**2*x**7/7 + 2*b*d*x**5*(a*d + b*c)/5 + x**3*(a**2*d**2/3 + 4*a*b*c*d/3 + b**2*c**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_157():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**3
    F = -a**2*c**2/(2*x**2) + 2*a*c*(a*d + b*c)*log(x) + b**2*d**2*x**6/6 + b*d*x**4*(a*d + b*c)/2 + x**2*(a**2*d**2/2 + 2*a*b*c*d + b**2*c**2/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_158():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**4
    F = -a**2*c**2/(3*x**3) - 2*a*c*(a*d + b*c)/x + b**2*d**2*x**5/5 + 2*b*d*x**3*(a*d + b*c)/3 + x*(a**2*d**2 + 4*a*b*c*d + b**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_159():
    f = x**4*(a + b*x**2)**2*(c + d*x**2)**3
    F = a**2*c**3*x**5/5 + a*c**2*x**7*(3*a*d + 2*b*c)/7 + b**2*d**3*x**15/15 + b*d**2*x**13*(2*a*d + 3*b*c)/13 + c*x**9*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/9 + d*x**11*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_160():
    f = x**3*(a + b*x**2)**2*(c + d*x**2)**3
    F = b**2*(c + d*x**2)**7/(14*d**4) - b*(c + d*x**2)**6*(-2*a*d + 3*b*c)/(12*d**4) - c*(c + d*x**2)**4*(-a*d + b*c)**2/(8*d**4) + (c + d*x**2)**5*(-a*d + b*c)*(-a*d + 3*b*c)/(10*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_161():
    f = x**2*(a + b*x**2)**2*(c + d*x**2)**3
    F = a**2*c**3*x**3/3 + a*c**2*x**5*(3*a*d + 2*b*c)/5 + b**2*d**3*x**13/13 + b*d**2*x**11*(2*a*d + 3*b*c)/11 + c*x**7*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/7 + d*x**9*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_162():
    f = x*(a + b*x**2)**2*(c + d*x**2)**3
    F = b**2*(c + d*x**2)**6/(12*d**3) - b*(c + d*x**2)**5*(-a*d + b*c)/(5*d**3) + (c + d*x**2)**4*(-a*d + b*c)**2/(8*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_163():
    f = (a + b*x**2)**2*(c + d*x**2)**3
    F = a**2*c**3*x + a*c**2*x**3*(3*a*d + 2*b*c)/3 + b**2*d**3*x**11/11 + b*d**2*x**9*(2*a*d + 3*b*c)/9 + c*x**5*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/5 + d*x**7*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_164():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x
    F = a**2*c**3*log(x) + a*c**2*x**2*(3*a*d + 2*b*c)/2 + b**2*d**3*x**10/10 + b*d**2*x**8*(2*a*d + 3*b*c)/8 + c*x**4*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/4 + d*x**6*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_165():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**2
    F = -a**2*c**3/x + a*c**2*x*(3*a*d + 2*b*c) + b**2*d**3*x**9/9 + b*d**2*x**7*(2*a*d + 3*b*c)/7 + c*x**3*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/3 + d*x**5*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_166():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**3
    F = -a**2*c**3/(2*x**2) + a*c**2*(3*a*d + 2*b*c)*log(x) + b**2*d**3*x**8/8 + b*d**2*x**6*(2*a*d + 3*b*c)/6 + c*x**2*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/2 + d*x**4*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_167():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**4
    F = -a**2*c**3/(3*x**3) - a*c**2*(3*a*d + 2*b*c)/x + b**2*d**3*x**7/7 + b*d**2*x**5*(2*a*d + 3*b*c)/5 + c*x*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2) + d*x**3*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_168():
    f = x**4*(a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**7/(7*d) - b*x**5*(-2*a*d + b*c)/(5*d**2) + c**(sympy.S(3)/2)*(-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/d**(sympy.S(9)/2) - c*x*(-a*d + b*c)**2/d**4 + x**3*(-a*d + b*c)**2/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_169():
    f = x**3*(a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**6/(6*d) - b*x**4*(-2*a*d + b*c)/(4*d**2) - c*(-a*d + b*c)**2*log(c + d*x**2)/(2*d**4) + x**2*(-a*d + b*c)**2/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_170():
    f = x**2*(a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**5/(5*d) - b*x**3*(-2*a*d + b*c)/(3*d**2) - sqrt(c)*(-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/d**(sympy.S(7)/2) + x*(-a*d + b*c)**2/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_171():
    f = x*(a + b*x**2)**2/(c + d*x**2)
    F = -b*x**2*(-a*d + b*c)/(2*d**2) + (a + b*x**2)**2/(4*d) + (-a*d + b*c)**2*log(c + d*x**2)/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_172():
    f = (a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**3/(3*d) - b*x*(-2*a*d + b*c)/d**2 + (-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_173():
    f = (a + b*x**2)**2/(x*(c + d*x**2))
    F = a**2*log(x)/c + b**2*x**2/(2*d) - (-a*d + b*c)**2*log(c + d*x**2)/(2*c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_174():
    f = (a + b*x**2)**2/(x**2*(c + d*x**2))
    F = -a**2/(c*x) + b**2*x/d - (-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_175():
    f = (a + b*x**2)**2/(x**3*(c + d*x**2))
    F = -a**2/(2*c*x**2) + a*(-a*d + 2*b*c)*log(x)/c**2 + (-a*d + b*c)**2*log(c + d*x**2)/(2*c**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_176():
    f = (a + b*x**2)**2/(x**4*(c + d*x**2))
    F = -a**2/(3*c*x**3) - a*(-a*d + 2*b*c)/(c**2*x) + (-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(5)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_177():
    f = (a + b*x**2)**2/(x**5*(c + d*x**2))
    F = -a**2/(4*c*x**4) - a*(-a*d + 2*b*c)/(2*c**2*x**2) + (-a*d + b*c)**2*log(x)/c**3 - (-a*d + b*c)**2*log(c + d*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_178():
    f = (a + b*x**2)**2/(x**6*(c + d*x**2))
    F = -a**2/(5*c*x**5) - a*(-a*d + 2*b*c)/(3*c**2*x**3) - (-a*d + b*c)**2/(c**3*x) - sqrt(d)*(-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/c**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_179():
    f = (a + b*x**2)**2/(x**7*(c + d*x**2))
    F = -a**2/(6*c*x**6) - a*(-a*d + 2*b*c)/(4*c**2*x**4) - (-a*d + b*c)**2/(2*c**3*x**2) - d*(-a*d + b*c)**2*log(x)/c**4 + d*(-a*d + b*c)**2*log(c + d*x**2)/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_180():
    f = x**4*(a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x**5/(5*d**2) - sqrt(c)*(-3*a*d + 7*b*c)*(-a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*d**(sympy.S(9)/2)) + x*(-3*a*d + 7*b*c)*(-a*d + b*c)/(2*d**4) + x**5*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - x**3*(-3*a*d + 7*b*c)*(-a*d + b*c)/(6*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_181():
    f = x**3*(a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x**4/(4*d**2) - b*x**2*(-a*d + b*c)/d**3 + c*(-a*d + b*c)**2/(2*d**4*(c + d*x**2)) + (-a*d + b*c)*(-a*d + 3*b*c)*log(c + d*x**2)/(2*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_182():
    f = x**2*(a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x**3/(3*d**2) + x**3*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - x*(-a*d + b*c)*(-a*d + 5*b*c)/(2*c*d**3) + (-a*d + b*c)*(-a*d + 5*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*sqrt(c)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_183():
    f = x*(a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x**2/(2*d**2) - b*(-a*d + b*c)*log(c + d*x**2)/d**3 - (-a*d + b*c)**2/(2*d**3*(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_184():
    f = (a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x/d**2 + x*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - (-a*d + b*c)*(a*d + 3*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_185():
    f = (a + b*x**2)**2/(x*(c + d*x**2)**2)
    F = a**2*log(x)/c**2 - (a**2/(2*c**2) - b**2/(2*d**2))*log(c + d*x**2) + (-a*d + b*c)**2/(2*c*d**2*(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_186():
    f = (a + b*x**2)**2/(x**2*(c + d*x**2)**2)
    F = -a**2/(c*x*(c + d*x**2)) - x*(3*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(2*c**2*d*(c + d*x**2)) + (-a*d + b*c)*(3*a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(5)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_187():
    f = (a + b*x**2)**2/(x**3*(c + d*x**2)**2)
    F = -a**2/(2*c**2*x**2) + 2*a*(-a*d + b*c)*log(x)/c**3 - a*(-a*d + b*c)*log(c + d*x**2)/c**3 - (-a*d + b*c)**2/(2*c**2*d*(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_188():
    f = (a + b*x**2)**2/(x**4*(c + d*x**2)**2)
    F = -a**2/(3*c*x**3*(c + d*x**2)) - a*(-5*a*d + 6*b*c)/(3*c**3*x) + x*(5*a**2*d**2 - 6*a*b*c*d + 3*b**2*c**2)/(6*c**3*(c + d*x**2)) + (-5*a*d + b*c)*(-a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(7)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_189():
    f = x**4*(a + b*x**2)**2/(c + d*x**2)**3
    F = b**2*x**3/(3*d**3) - x*(-a*d + b*c)*(-a*d + 9*b*c)/(8*d**4*(c + d*x**2)) + x**5*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - x*(a**2*d**2 - 10*a*b*c*d + 13*b**2*c**2)/(4*c*d**4) + (3*a**2*d**2 - 30*a*b*c*d + 35*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*sqrt(c)*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_190():
    f = x**3*(a + b*x**2)**2/(c + d*x**2)**3
    F = b**2*x**2/(2*d**3) - b*(-2*a*d + 3*b*c)*log(c + d*x**2)/(2*d**4) + c*(-a*d + b*c)**2/(4*d**4*(c + d*x**2)**2) - (-a*d + b*c)*(-a*d + 3*b*c)/(2*d**4*(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_191():
    f = x**2*(a + b*x**2)**2/(c + d*x**2)**3
    F = b**2*x/d**3 + x**3*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) + x*(-a*d + b*c)*(a*d + 7*b*c)/(8*c*d**3*(c + d*x**2)) - (-a**2*d**2 - 6*a*b*c*d + 15*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(3)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_192():
    f = x*(a + b*x**2)**2/(c + d*x**2)**3
    F = b**2*log(c + d*x**2)/(2*d**3) + b*(-a*d + b*c)/(d**3*(c + d*x**2)) - (-a*d + b*c)**2/(4*d**3*(c + d*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_193():
    f = (a + b*x**2)**2/(c + d*x**2)**3
    F = x*(3*a**2/c**2 - 3*b**2/d**2)/(8*c + 8*d*x**2) - x*(a + b*x**2)*(-a*d + b*c)/(4*c*d*(c + d*x**2)**2) + (3*a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_194():
    f = (a + b*x**2)**2/(x*(c + d*x**2)**3)
    F = a**2*log(x)/c**3 - a**2*log(c + d*x**2)/(2*c**3) + (a**2/c**2 - b**2/d**2)/(2*c + 2*d*x**2) + (-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_195():
    f = (a + b*x**2)**2/(x**2*(c + d*x**2)**3)
    F = -a**2/(c*x*(c + d*x**2)**2) - x*(5*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(4*c**2*d*(c + d*x**2)**2) + x*(3*a*d*(-5*a*d + 2*b*c) + b**2*c**2)/(8*c**3*d*(c + d*x**2)) + (3*a*d*(-5*a*d + 2*b*c) + b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(7)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_196():
    f = (a + b*x**2)**2/(x**3*(c + d*x**2)**3)
    F = -a**2/(2*c**3*x**2) + a*(-a*d + b*c)/(c**3*(c + d*x**2)) + a*(-3*a*d + 2*b*c)*log(x)/c**4 - a*(-3*a*d + 2*b*c)*log(c + d*x**2)/(2*c**4) - (-a*d + b*c)**2/(4*c**2*d*(c + d*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_197():
    f = (a + b*x**2)**2/(x**4*(c + d*x**2)**3)
    F = -a**2/(3*c*x**3*(c + d*x**2)**2) - a*(-7*a*d + 6*b*c)/(3*c**4*x) + x*(7*a**2*d**2 - 6*a*b*c*d + 3*b**2*c**2)/(12*c**3*(c + d*x**2)**2) + x*(-7*a*d + 3*b*c)**2/(24*c**4*(c + d*x**2)) + (35*a**2*d**2 - 30*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(9)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_198():
    f = x**5*(c + d*x**2)/(a + b*x**2)
    F = a**2*(-a*d + b*c)*log(a + b*x**2)/(2*b**4) - a*x**2*(-a*d + b*c)/(2*b**3) + d*x**6/(6*b) + x**4*(-a*d + b*c)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_199():
    f = x**4*(c + d*x**2)/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) - a*x*(-a*d + b*c)/b**3 + d*x**5/(5*b) + x**3*(-a*d + b*c)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_200():
    f = x**3*(c + d*x**2)/(a + b*x**2)
    F = -a*(-a*d + b*c)*log(a + b*x**2)/(2*b**3) + d*x**4/(4*b) + x**2*(-a*d + b*c)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_201():
    f = x**2*(c + d*x**2)/(a + b*x**2)
    F = -sqrt(a)*(-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(5)/2) + d*x**3/(3*b) + x*(-a*d + b*c)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_202():
    f = x*(c + d*x**2)/(a + b*x**2)
    F = d*x**2/(2*b) + (-a*d + b*c)*log(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_203():
    f = (c + d*x**2)/(a + b*x**2)
    F = d*x/b + (-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_204():
    f = (c + d*x**2)/(x*(a + b*x**2))
    F = c*log(x)/a - (-a*d + b*c)*log(a + b*x**2)/(2*a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_205():
    f = (c + d*x**2)/(x**2*(a + b*x**2))
    F = -c/(a*x) - (-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_206():
    f = (c + d*x**2)/(x**3*(a + b*x**2))
    F = -c/(2*a*x**2) - (-a*d + b*c)*log(x)/a**2 + (-a*d + b*c)*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_207():
    f = (c + d*x**2)/(x**4*(a + b*x**2))
    F = -c/(3*a*x**3) + (-a*d + b*c)/(a**2*x) + sqrt(b)*(-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_208():
    f = x**5*(c + d*x**2)**2/(a + b*x**2)
    F = a**2*(-a*d + b*c)**2*log(a + b*x**2)/(2*b**5) - a*x**2*(-a*d + b*c)**2/(2*b**4) + d**2*x**8/(8*b) + d*x**6*(-a*d + 2*b*c)/(6*b**2) + x**4*(-a*d + b*c)**2/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_209():
    f = x**4*(c + d*x**2)**2/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(9)/2) - a*x*(-a*d + b*c)**2/b**4 + d**2*x**7/(7*b) + d*x**5*(-a*d + 2*b*c)/(5*b**2) + x**3*(-a*d + b*c)**2/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_210():
    f = x**3*(c + d*x**2)**2/(a + b*x**2)
    F = -a*(-a*d + b*c)**2*log(a + b*x**2)/(2*b**4) + d**2*x**6/(6*b) + d*x**4*(-a*d + 2*b*c)/(4*b**2) + x**2*(-a*d + b*c)**2/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_211():
    f = x**2*(c + d*x**2)**2/(a + b*x**2)
    F = -sqrt(a)*(-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(7)/2) + d**2*x**5/(5*b) + d*x**3*(-a*d + 2*b*c)/(3*b**2) + x*(-a*d + b*c)**2/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_212():
    f = x*(c + d*x**2)**2/(a + b*x**2)
    F = (c + d*x**2)**2/(4*b) + d*x**2*(-a*d + b*c)/(2*b**2) + (-a*d + b*c)**2*log(a + b*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_213():
    f = (c + d*x**2)**2/(a + b*x**2)
    F = d**2*x**3/(3*b) + d*x*(-a*d + 2*b*c)/b**2 + (-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_214():
    f = (c + d*x**2)**2/(x*(a + b*x**2))
    F = d**2*x**2/(2*b) + c**2*log(x)/a - (-a*d + b*c)**2*log(a + b*x**2)/(2*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_215():
    f = (c + d*x**2)**2/(x**2*(a + b*x**2))
    F = d**2*x/b - c**2/(a*x) - (-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_216():
    f = (c + d*x**2)**2/(x**3*(a + b*x**2))
    F = -c**2/(2*a*x**2) - c*(-2*a*d + b*c)*log(x)/a**2 + (-a*d + b*c)**2*log(a + b*x**2)/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_217():
    f = (c + d*x**2)**2/(x**4*(a + b*x**2))
    F = -c**2/(3*a*x**3) + c*(-2*a*d + b*c)/(a**2*x) + (-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_218():
    f = x**5*(c + d*x**2)**3/(a + b*x**2)
    F = a**2*(-a*d + b*c)**3*log(a + b*x**2)/(2*b**6) - a*x**2*(-a*d + b*c)**3/(2*b**5) + d**3*x**10/(10*b) + d**2*x**8*(-a*d + 3*b*c)/(8*b**2) + d*x**6*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(6*b**3) + x**4*(-a*d + b*c)**3/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_219():
    f = x**4*(c + d*x**2)**3/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(11)/2) - a*x*(-a*d + b*c)**3/b**5 + d**3*x**9/(9*b) + d**2*x**7*(-a*d + 3*b*c)/(7*b**2) + d*x**5*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(5*b**3) + x**3*(-a*d + b*c)**3/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_220():
    f = x**3*(c + d*x**2)**3/(a + b*x**2)
    F = -a*(-a*d + b*c)**3*log(a + b*x**2)/(2*b**5) + d**3*x**8/(8*b) + d**2*x**6*(-a*d + 3*b*c)/(6*b**2) + d*x**4*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(4*b**3) + x**2*(-a*d + b*c)**3/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_221():
    f = x**2*(c + d*x**2)**3/(a + b*x**2)
    F = -sqrt(a)*(-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(9)/2) + d**3*x**7/(7*b) + d**2*x**5*(-a*d + 3*b*c)/(5*b**2) + d*x**3*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(3*b**3) + x*(-a*d + b*c)**3/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_222():
    f = x*(c + d*x**2)**3/(a + b*x**2)
    F = (c + d*x**2)**3/(6*b) + (c + d*x**2)**2*(-a*d + b*c)/(4*b**2) + d*x**2*(-a*d + b*c)**2/(2*b**3) + (-a*d + b*c)**3*log(a + b*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_223():
    f = (c + d*x**2)**3/(a + b*x**2)
    F = d**3*x**5/(5*b) + d**2*x**3*(-a*d + 3*b*c)/(3*b**2) + d*x*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/b**3 + (-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_224():
    f = (c + d*x**2)**3/(x*(a + b*x**2))
    F = d**3*x**4/(4*b) + d**2*x**2*(-a*d + 3*b*c)/(2*b**2) + c**3*log(x)/a - (-a*d + b*c)**3*log(a + b*x**2)/(2*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_225():
    f = (c + d*x**2)**3/(x**2*(a + b*x**2))
    F = d**3*x**3/(3*b) + d**2*x*(-a*d + 3*b*c)/b**2 - c**3/(a*x) - (-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_226():
    f = (c + d*x**2)**3/(x**3*(a + b*x**2))
    F = d**3*x**2/(2*b) - c**3/(2*a*x**2) - c**2*(-3*a*d + b*c)*log(x)/a**2 + (-a*d + b*c)**3*log(a + b*x**2)/(2*a**2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_227():
    f = (c + d*x**2)**3/(x**4*(a + b*x**2))
    F = d**3*x/b - c**3/(3*a*x**3) + c**2*(-3*a*d + b*c)/(a**2*x) + (-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_228():
    f = x**5/((a + b*x**2)*(c + d*x**2))
    F = a**2*log(a + b*x**2)/(2*b**2*(-a*d + b*c)) - c**2*log(c + d*x**2)/(2*d**2*(-a*d + b*c)) + x**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_229():
    f = x**4/((a + b*x**2)*(c + d*x**2))
    F = a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(3)/2)*(-a*d + b*c)) - c**(sympy.S(3)/2)*atan(sqrt(d)*x/sqrt(c))/(d**(sympy.S(3)/2)*(-a*d + b*c)) + x/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_230():
    f = x**3/((a + b*x**2)*(c + d*x**2))
    F = -a*log(a + b*x**2)/(2*b*(-a*d + b*c)) + c*log(c + d*x**2)/(2*d*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_231():
    f = x**2/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(a)*atan(sqrt(b)*x/sqrt(a))/(sqrt(b)*(-a*d + b*c)) + sqrt(c)*atan(sqrt(d)*x/sqrt(c))/(sqrt(d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_232():
    f = x/((a + b*x**2)*(c + d*x**2))
    F = log(a + b*x**2)/(-2*a*d + 2*b*c) - log(c + d*x**2)/(-2*a*d + 2*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_233():
    f = 1/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(d)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*(-a*d + b*c)) + sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_234():
    f = 1/(x*(a + b*x**2)*(c + d*x**2))
    F = d*log(c + d*x**2)/(2*c*(-a*d + b*c)) - b*log(a + b*x**2)/(2*a*(-a*d + b*c)) + log(x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_235():
    f = 1/(x**2*(a + b*x**2)*(c + d*x**2))
    F = d**(sympy.S(3)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(3)/2)*(-a*d + b*c)) - 1/(a*c*x) - b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_236():
    f = 1/(x**3*(a + b*x**2)*(c + d*x**2))
    F = -d**2*log(c + d*x**2)/(2*c**2*(-a*d + b*c)) - 1/(2*a*c*x**2) + b**2*log(a + b*x**2)/(2*a**2*(-a*d + b*c)) - (a*d + b*c)*log(x)/(a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_237():
    f = 1/(x**4*(a + b*x**2)*(c + d*x**2))
    F = -d**(sympy.S(5)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(5)/2)*(-a*d + b*c)) - 1/(3*a*c*x**3) + (a*d + b*c)/(a**2*c**2*x) + b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_238():
    f = 1/(x**5*(a + b*x**2)*(c + d*x**2))
    F = d**3*log(c + d*x**2)/(2*c**3*(-a*d + b*c)) - 1/(4*a*c*x**4) + (a*d + b*c)/(2*a**2*c**2*x**2) - b**3*log(a + b*x**2)/(2*a**3*(-a*d + b*c)) + (a**2*d**2 + a*b*c*d + b**2*c**2)*log(x)/(a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_239():
    f = 1/(x**6*(a + b*x**2)*(c + d*x**2))
    F = d**(sympy.S(7)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(7)/2)*(-a*d + b*c)) - 1/(5*a*c*x**5) + (a*d + b*c)/(3*a**2*c**2*x**3) - (a**2*d**2 + a*b*c*d + b**2*c**2)/(a**3*c**3*x) - b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(7)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_240():
    f = 1/(x**7*(a + b*x**2)*(c + d*x**2))
    F = -d**4*log(c + d*x**2)/(2*c**4*(-a*d + b*c)) - 1/(6*a*c*x**6) + (a*d + b*c)/(4*a**2*c**2*x**4) - (a**2*d**2 + a*b*c*d + b**2*c**2)/(2*a**3*c**3*x**2) + b**4*log(a + b*x**2)/(2*a**4*(-a*d + b*c)) - (a*d + b*c)*(a**2*d**2 + b**2*c**2)*log(x)/(a**4*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_241():
    f = x**5/((a + b*x**2)**2*(c + d*x**2))
    F = -a**2/(2*b**2*(a + b*x**2)*(-a*d + b*c)) - a*(-a*d + 2*b*c)*log(a + b*x**2)/(2*b**2*(-a*d + b*c)**2) + c**2*log(c + d*x**2)/(2*d*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_242():
    f = x**4/((a + b*x**2)*(c + d*x**2)**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(b)*(-a*d + b*c)**2) + sqrt(c)*(-3*a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*d**(sympy.S(3)/2)*(-a*d + b*c)**2) - c*x/(2*d*(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_243():
    f = x**3/((a + b*x**2)*(c + d*x**2)**2)
    F = -a*log(a + b*x**2)/(2*(-a*d + b*c)**2) + a*log(c + d*x**2)/(2*(-a*d + b*c)**2) - c/(2*d*(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_244():
    f = x**2/((a + b*x**2)*(c + d*x**2)**2)
    F = -sqrt(a)*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(-a*d + b*c)**2 + x/((c + d*x**2)*(-2*a*d + 2*b*c)) + (a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*sqrt(c)*sqrt(d)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_245():
    f = x/((a + b*x**2)*(c + d*x**2)**2)
    F = b*log(a + b*x**2)/(2*(-a*d + b*c)**2) - b*log(c + d*x**2)/(2*(-a*d + b*c)**2) + 1/((c + d*x**2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_246():
    f = 1/((a + b*x**2)*(c + d*x**2)**2)
    F = -d*x/(2*c*(c + d*x**2)*(-a*d + b*c)) - sqrt(d)*(-a*d + 3*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_247():
    f = 1/(x*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*(c + d*x**2)*(-a*d + b*c)) + d*(-a*d + 2*b*c)*log(c + d*x**2)/(2*c**2*(-a*d + b*c)**2) - b**2*log(a + b*x**2)/(2*a*(-a*d + b*c)**2) + log(x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_248():
    f = 1/(x**2*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*x*(c + d*x**2)*(-a*d + b*c)) + d**(sympy.S(3)/2)*(-3*a*d + 5*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(5)/2)*(-a*d + b*c)**2) - (-3*a*d + 2*b*c)/(2*a*c**2*x*(-a*d + b*c)) - b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_249():
    f = 1/(x**3*(a + b*x**2)*(c + d*x**2)**2)
    F = d**2/(2*c**2*(c + d*x**2)*(-a*d + b*c)) - d**2*(-2*a*d + 3*b*c)*log(c + d*x**2)/(2*c**3*(-a*d + b*c)**2) - 1/(2*a*c**2*x**2) + b**3*log(a + b*x**2)/(2*a**2*(-a*d + b*c)**2) - (2*a*d + b*c)*log(x)/(a**2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_250():
    f = 1/(x**4*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*x**3*(c + d*x**2)*(-a*d + b*c)) - d**(sympy.S(5)/2)*(-5*a*d + 7*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(7)/2)*(-a*d + b*c)**2) - (-5*a*d + 2*b*c)/(6*a*c**2*x**3*(-a*d + b*c)) + (-5*a**2*d**2 + 2*a*b*c*d + 2*b**2*c**2)/(2*a**2*c**3*x*(-a*d + b*c)) + b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_251():
    f = x**5/((a + b*x**2)**3*(c + d*x**2))
    F = -a**2/(4*b**2*(a + b*x**2)**2*(-a*d + b*c)) + a*(-a*d + 2*b*c)/(2*b**2*(a + b*x**2)*(-a*d + b*c)**2) + c**2*log(a + b*x**2)/(2*(-a*d + b*c)**3) - c**2*log(c + d*x**2)/(2*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_252():
    f = x**4/((a + b*x**2)*(c + d*x**2)**3)
    F = a**(sympy.S(3)/2)*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(-a*d + b*c)**3 - c*x/(4*d*(c + d*x**2)**2*(-a*d + b*c)) + x*(-5*a*d + b*c)/(8*d*(c + d*x**2)*(-a*d + b*c)**2) + (-3*a**2*d**2 - 6*a*b*c*d + b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*sqrt(c)*d**(sympy.S(3)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_253():
    f = x**3/((a + b*x**2)*(c + d*x**2)**3)
    F = -a*b*log(a + b*x**2)/(2*(-a*d + b*c)**3) + a*b*log(c + d*x**2)/(2*(-a*d + b*c)**3) - a/(2*(c + d*x**2)*(-a*d + b*c)**2) - c/(4*d*(c + d*x**2)**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_254():
    f = x**2/((a + b*x**2)*(c + d*x**2)**3)
    F = -sqrt(a)*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(-a*d + b*c)**3 + x/((c + d*x**2)**2*(-4*a*d + 4*b*c)) + x*(a*d + 3*b*c)/(8*c*(c + d*x**2)*(-a*d + b*c)**2) + (-a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(3)/2)*sqrt(d)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_255():
    f = x/((a + b*x**2)*(c + d*x**2)**3)
    F = b**2*log(a + b*x**2)/(2*(-a*d + b*c)**3) - b**2*log(c + d*x**2)/(2*(-a*d + b*c)**3) + b/(2*(c + d*x**2)*(-a*d + b*c)**2) + 1/((c + d*x**2)**2*(-4*a*d + 4*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_256():
    f = 1/((a + b*x**2)*(c + d*x**2)**3)
    F = -d*x/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - d*x*(-3*a*d + 7*b*c)/(8*c**2*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(d)*(3*a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**3) + b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_257():
    f = 1/(x*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - d*(-a*d + 2*b*c)/(2*c**2*(c + d*x**2)*(-a*d + b*c)**2) + d*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)*log(c + d*x**2)/(2*c**3*(-a*d + b*c)**3) - b**3*log(a + b*x**2)/(2*a*(-a*d + b*c)**3) + log(x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_258():
    f = 1/(x**2*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*x*(c + d*x**2)**2*(-a*d + b*c)) - d*(-5*a*d + 9*b*c)/(8*c**2*x*(c + d*x**2)*(-a*d + b*c)**2) + d**(sympy.S(3)/2)*(15*a**2*d**2 - 42*a*b*c*d + 35*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(7)/2)*(-a*d + b*c)**3) - (15*a**2*d**2 - 27*a*b*c*d + 8*b**2*c**2)/(8*a*c**3*x*(-a*d + b*c)**2) - b**(sympy.S(7)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_259():
    f = 1/(x**3*(a + b*x**2)*(c + d*x**2)**3)
    F = d**2/(4*c**2*(c + d*x**2)**2*(-a*d + b*c)) + d**2*(-2*a*d + 3*b*c)/(2*c**3*(c + d*x**2)*(-a*d + b*c)**2) - d**2*(3*a**2*d**2 - 8*a*b*c*d + 6*b**2*c**2)*log(c + d*x**2)/(2*c**4*(-a*d + b*c)**3) - 1/(2*a*c**3*x**2) + b**4*log(a + b*x**2)/(2*a**2*(-a*d + b*c)**3) - (3*a*d + b*c)*log(x)/(a**2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_260():
    f = 1/(x**4*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*x**3*(c + d*x**2)**2*(-a*d + b*c)) - d*(-7*a*d + 11*b*c)/(8*c**2*x**3*(c + d*x**2)*(-a*d + b*c)**2) - d**(sympy.S(5)/2)*(35*a**2*d**2 - 90*a*b*c*d + 63*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(9)/2)*(-a*d + b*c)**3) - (35*a**2*d**2 - 55*a*b*c*d + 8*b**2*c**2)/(24*a*c**3*x**3*(-a*d + b*c)**2) + (35*a**3*d**3 - 55*a**2*b*c*d**2 + 8*a*b**2*c**2*d + 8*b**3*c**3)/(8*a**2*c**4*x*(-a*d + b*c)**2) + b**(sympy.S(9)/2)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_261():
    f = x/((x**2 + 1)*(x**2 + 4))
    F = log(x**2 + 1)/6 - log(x**2 + 4)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_262():
    f = x**4*(c + d*x**2)/(a + b*x**2)**2
    F = -sqrt(a)*(-5*a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(7)/2)) + a*x*(-a*d + b*c)/(2*b**3*(a + b*x**2)) + d*x**3/(3*b**2) + x*(-2*a*d + b*c)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_263():
    f = x**3*(c + d*x**2)/(a + b*x**2)**2
    F = a*(-a*d + b*c)/(2*b**3*(a + b*x**2)) + d*x**2/(2*b**2) + (-2*a*d + b*c)*log(a + b*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_264():
    f = x**2*(c + d*x**2)/(a + b*x**2)**2
    F = d*x/b**2 - x*(-a*d + b*c)/(2*b**2*(a + b*x**2)) + (-3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_265():
    f = x*(c + d*x**2)/(a + b*x**2)**2
    F = d*log(a + b*x**2)/(2*b**2) + (a*d - b*c)/(2*b**2*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_266():
    f = (c + d*x**2)/(a + b*x**2)**2
    F = x*(-a*d + b*c)/(2*a*b*(a + b*x**2)) + (a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_267():
    f = (c + d*x**2)/(x*(a + b*x**2)**2)
    F = (-a*d + b*c)/(2*a*b*(a + b*x**2)) + c*log(x)/a**2 - c*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_268():
    f = (c + d*x**2)/(x**2*(a + b*x**2)**2)
    F = -c/(a**2*x) - x*(-a*d + b*c)/(2*a**2*(a + b*x**2)) - (-a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_269():
    f = (c + d*x**2)/(x**3*(a + b*x**2)**2)
    F = -c/(2*a**2*x**2) - (-a*d + b*c)/(2*a**2*(a + b*x**2)) - (-a*d + 2*b*c)*log(x)/a**3 + (-a*d + 2*b*c)*log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_270():
    f = (c + d*x**2)/(x**4*(a + b*x**2)**2)
    F = -c/(3*a**2*x**3) + b*x*(-a*d + b*c)/(2*a**3*(a + b*x**2)) + (-a*d + 2*b*c)/(a**3*x) + sqrt(b)*(-3*a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_271():
    f = x**4*(c + d*x**2)**2/(a + b*x**2)**2
    F = -sqrt(a)*(-7*a*d + 3*b*c)*(-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(9)/2)) + d**2*x**5/(5*b**2) + x*(-7*a*d + 3*b*c)*(-a*d + b*c)/(2*b**4) + x**5*(-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) - x**3*(-7*a*d + 3*b*c)*(-a*d + b*c)/(6*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_272():
    f = x**3*(c + d*x**2)**2/(a + b*x**2)**2
    F = a*(-a*d + b*c)**2/(2*b**4*(a + b*x**2)) + d**2*x**4/(4*b**2) + d*x**2*(-a*d + b*c)/b**3 + (-3*a*d + b*c)*(-a*d + b*c)*log(a + b*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_273():
    f = x**2*(c + d*x**2)**2/(a + b*x**2)**2
    F = d**2*x**3/(3*b**2) + x**3*(-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) - x*(-5*a*d + b*c)*(-a*d + b*c)/(2*a*b**3) + (-5*a*d + b*c)*(-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_274():
    f = x*(c + d*x**2)**2/(a + b*x**2)**2
    F = d**2*x**2/(2*b**2) + d*(-a*d + b*c)*log(a + b*x**2)/b**3 - (-a*d + b*c)**2/(2*b**3*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_275():
    f = (c + d*x**2)**2/(a + b*x**2)**2
    F = d**2*x/b**2 + x*(-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) + (-a*d + b*c)*(3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_276():
    f = (c + d*x**2)**2/(x*(a + b*x**2)**2)
    F = -(-d**2/(2*b**2) + c**2/(2*a**2))*log(a + b*x**2) + (-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) + c**2*log(x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_277():
    f = (c + d*x**2)**2/(x**2*(a + b*x**2)**2)
    F = -c**2/(a*x*(a + b*x**2)) - x*(a*d**2/b - 2*c*d + 3*b*c**2/a)/(2*a*(a + b*x**2)) - (-a*d + b*c)*(a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_278():
    f = (c + d*x**2)**2/(x**3*(a + b*x**2)**2)
    F = -c**2/(2*a**2*x**2) - (-a*d + b*c)**2/(2*a**2*b*(a + b*x**2)) - 2*c*(-a*d + b*c)*log(x)/a**3 + c*(-a*d + b*c)*log(a + b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_279():
    f = (c + d*x**2)**2/(x**4*(a + b*x**2)**2)
    F = -c**2/(3*a*x**3*(a + b*x**2)) + c*(-6*a*d + 5*b*c)/(3*a**3*x) + x*(3*a**2*d**2 - 6*a*b*c*d + 5*b**2*c**2)/(6*a**3*(a + b*x**2)) + (-a*d + b*c)*(-a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_280():
    f = x**4*(c + d*x**2)**3/(a + b*x**2)**2
    F = -3*sqrt(a)*(-3*a*d + b*c)*(-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(11)/2)) - x**3*(c + d*x**2)**3/(2*b*(a + b*x**2)) + 9*d**3*x**7/(14*b**2) + 3*d**2*x**5*(-3*a*d + 7*b*c)/(10*b**3) + d*x**3*(3*a**2*d**2 - 7*a*b*c*d + 5*b**2*c**2)/(2*b**4) + x*(-9*a*d + 3*b*c)*(-a*d + b*c)**2/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_281():
    f = x**3*(c + d*x**2)**3/(a + b*x**2)**2
    F = a*(-a*d + b*c)**3/(2*b**5*(a + b*x**2)) + d**3*x**6/(6*b**2) + d**2*x**4*(-2*a*d + 3*b*c)/(4*b**3) + 3*d*x**2*(-a*d + b*c)**2/(2*b**4) + (-4*a*d + b*c)*(-a*d + b*c)**2*log(a + b*x**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_282():
    f = x**2*(c + d*x**2)**3/(a + b*x**2)**2
    F = -x*(c + d*x**2)**3/(2*b*(a + b*x**2)) + 7*d*x*(c + d*x**2)**2/(10*b**2) + d*x*(c + d*x**2)*(-35*a*d + 33*b*c)/(30*b**3) + d*x*(105*a**2*d**2 - 190*a*b*c*d + 81*b**2*c**2)/(30*b**4) + (-7*a*d + b*c)*(-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_283():
    f = x*(c + d*x**2)**3/(a + b*x**2)**2
    F = d**3*x**4/(4*b**2) + d**2*x**2*(-2*a*d + 3*b*c)/(2*b**3) + 3*d*(-a*d + b*c)**2*log(a + b*x**2)/(2*b**4) - (-a*d + b*c)**3/(2*b**4*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_284():
    f = (c + d*x**2)**3/(a + b*x**2)**2
    F = d**3*x**3/(3*b**2) + d**2*x*(-2*a*d + 3*b*c)/b**3 + x*(-a*d + b*c)**3/(2*a*b**3*(a + b*x**2)) + (-a*d + b*c)**2*(5*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_285():
    f = (c + d*x**2)**3/(x*(a + b*x**2)**2)
    F = d**3*x**2/(2*b**2) + (-a*d + b*c)**3/(2*a*b**3*(a + b*x**2)) + c**3*log(x)/a**2 - (-a*d + b*c)**2*(2*a*d + b*c)*log(a + b*x**2)/(2*a**2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_286():
    f = (c + d*x**2)**3/(x**2*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*x*(a + b*x**2)) - d**2*x*(-3*a*d + b*c)/(2*a*b**2) - c**2*(-a*d + 3*b*c)/(2*a**2*b*x) - 3*(-a*d + b*c)**2*(a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_287():
    f = (c + d*x**2)**3/(x**3*(a + b*x**2)**2)
    F = -c**3/(2*a**2*x**2) - (-a*d + b*c)**3/(2*a**2*b**2*(a + b*x**2)) - c**2*(-3*a*d + 2*b*c)*log(x)/a**3 + (-a*d + b*c)**2*(a*d + 2*b*c)*log(a + b*x**2)/(2*a**3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_288():
    f = (c + d*x**2)**3/(x**4*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*x**3*(a + b*x**2)) - c**2*(-3*a*d + 5*b*c)/(6*a**2*b*x**3) + c*(2*a**2*d**2 - 9*a*b*c*d + 5*b**2*c**2)/(2*a**3*b*x) + (-a*d + b*c)**2*(a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_289():
    f = x**4/((a + b*x**2)**2*(c + d*x**2))
    F = -sqrt(a)*(-a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*b**(sympy.S(3)/2)*(-a*d + b*c)**2) + a*x/(2*b*(a + b*x**2)*(-a*d + b*c)) + c**(sympy.S(3)/2)*atan(sqrt(d)*x/sqrt(c))/(sqrt(d)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_290():
    f = x**3/((a + b*x**2)**2*(c + d*x**2))
    F = a/(2*b*(a + b*x**2)*(-a*d + b*c)) + c*log(a + b*x**2)/(2*(-a*d + b*c)**2) - c*log(c + d*x**2)/(2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_291():
    f = x**2/((a + b*x**2)**2*(c + d*x**2))
    F = -sqrt(c)*sqrt(d)*atan(sqrt(d)*x/sqrt(c))/(-a*d + b*c)**2 - x/((a + b*x**2)*(-2*a*d + 2*b*c)) + (a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*sqrt(b)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_292():
    f = x/((a + b*x**2)**2*(c + d*x**2))
    F = -d*log(a + b*x**2)/(2*(-a*d + b*c)**2) + d*log(c + d*x**2)/(2*(-a*d + b*c)**2) - 1/((a + b*x**2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_293():
    f = 1/((a + b*x**2)**2*(c + d*x**2))
    F = d**(sympy.S(3)/2)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*(-a*d + b*c)**2) + b*x/(2*a*(a + b*x**2)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_294():
    f = 1/(x*(a + b*x**2)**2*(c + d*x**2))
    F = -d**2*log(c + d*x**2)/(2*c*(-a*d + b*c)**2) + b/(2*a*(a + b*x**2)*(-a*d + b*c)) - b*(-2*a*d + b*c)*log(a + b*x**2)/(2*a**2*(-a*d + b*c)**2) + log(x)/(a**2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_295():
    f = 1/(x**2*(a + b*x**2)**2*(c + d*x**2))
    F = -d**(sympy.S(5)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(3)/2)*(-a*d + b*c)**2) + b/(2*a*x*(a + b*x**2)*(-a*d + b*c)) - (-2*a*d + 3*b*c)/(2*a**2*c*x*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_296():
    f = 1/(x**3*(a + b*x**2)**2*(c + d*x**2))
    F = d**3*log(c + d*x**2)/(2*c**2*(-a*d + b*c)**2) - b**2/(2*a**2*(a + b*x**2)*(-a*d + b*c)) - 1/(2*a**2*c*x**2) + b**2*(-3*a*d + 2*b*c)*log(a + b*x**2)/(2*a**3*(-a*d + b*c)**2) - (a*d + 2*b*c)*log(x)/(a**3*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_297():
    f = 1/(x**4*(a + b*x**2)**2*(c + d*x**2))
    F = d**(sympy.S(7)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(5)/2)*(-a*d + b*c)**2) + b/(2*a*x**3*(a + b*x**2)*(-a*d + b*c)) - (-2*a*d + 5*b*c)/(6*a**2*c*x**3*(-a*d + b*c)) + (-2*a**2*d**2 - 2*a*b*c*d + 5*b**2*c**2)/(2*a**3*c**2*x*(-a*d + b*c)) + b**(sympy.S(5)/2)*(-7*a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_298():
    f = 1/(x**5*(a + b*x**2)**2*(c + d*x**2))
    F = -d**4*log(c + d*x**2)/(2*c**3*(-a*d + b*c)**2) - 1/(4*a**2*c*x**4) + b**3/(2*a**3*(a + b*x**2)*(-a*d + b*c)) + (a*d + 2*b*c)/(2*a**3*c**2*x**2) - b**3*(-4*a*d + 3*b*c)*log(a + b*x**2)/(2*a**4*(-a*d + b*c)**2) + (a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)*log(x)/(a**4*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_299():
    f = 1/(x**6*(a + b*x**2)**2*(c + d*x**2))
    F = -d**(sympy.S(9)/2)*atan(sqrt(d)*x/sqrt(c))/(c**(sympy.S(7)/2)*(-a*d + b*c)**2) + b/(2*a*x**5*(a + b*x**2)*(-a*d + b*c)) - (-2*a*d + 7*b*c)/(10*a**2*c*x**5*(-a*d + b*c)) + (-2*a**2*d**2 - 2*a*b*c*d + 7*b**2*c**2)/(6*a**3*c**2*x**3*(-a*d + b*c)) - (-2*a**3*d**3 - 2*a**2*b*c*d**2 - 2*a*b**2*c**2*d + 7*b**3*c**3)/(2*a**4*c**3*x*(-a*d + b*c)) - b**(sympy.S(7)/2)*(-9*a*d + 7*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(9)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_300():
    f = 1/(x**7*(a + b*x**2)**2*(c + d*x**2))
    F = d**5*log(c + d*x**2)/(2*c**4*(-a*d + b*c)**2) - 1/(6*a**2*c*x**6) + (a*d + 2*b*c)/(4*a**3*c**2*x**4) - b**4/(2*a**4*(a + b*x**2)*(-a*d + b*c)) - (a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)/(2*a**4*c**3*x**2) + b**4*(-5*a*d + 4*b*c)*log(a + b*x**2)/(2*a**5*(-a*d + b*c)**2) - (a**3*d**3 + 2*a**2*b*c*d**2 + 3*a*b**2*c**2*d + 4*b**3*c**3)*log(x)/(a**5*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_301():
    f = x**4/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -sqrt(a)*(a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(b)*(-a*d + b*c)**3) + a*x/(2*b*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + sqrt(c)*(3*a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*sqrt(d)*(-a*d + b*c)**3) + x*(a*d + b*c)/(2*b*(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_302():
    f = x**3/((a + b*x**2)**2*(c + d*x**2)**2)
    F = a/(2*(a + b*x**2)*(-a*d + b*c)**2) + c/(2*(c + d*x**2)*(-a*d + b*c)**2) + (a*d + b*c)*log(a + b*x**2)/(2*(-a*d + b*c)**3) - (a*d + b*c)*log(c + d*x**2)/(2*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_303():
    f = x**2/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -d*x/((c + d*x**2)*(-a*d + b*c)**2) - x/((a + b*x**2)*(c + d*x**2)*(-2*a*d + 2*b*c)) - sqrt(d)*(a*d + 3*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*sqrt(c)*(-a*d + b*c)**3) + sqrt(b)*(3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_304():
    f = x/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -b*d*log(a + b*x**2)/(-a*d + b*c)**3 + b*d*log(c + d*x**2)/(-a*d + b*c)**3 - b/(2*(a + b*x**2)*(-a*d + b*c)**2) - d/(2*(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_305():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**2)
    F = d**(sympy.S(3)/2)*(-a*d + 5*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**3) + b*x/(2*a*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*x*(a*d + b*c)/(2*a*c*(c + d*x**2)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*(-5*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_306():
    f = 1/(x*(a + b*x**2)**2*(c + d*x**2)**2)
    F = d**2/(2*c*(c + d*x**2)*(-a*d + b*c)**2) - d**2*(-a*d + 3*b*c)*log(c + d*x**2)/(2*c**2*(-a*d + b*c)**3) + b**2/(2*a*(a + b*x**2)*(-a*d + b*c)**2) - b**2*(-3*a*d + b*c)*log(a + b*x**2)/(2*a**2*(-a*d + b*c)**3) + log(x)/(a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_307():
    f = 1/(x**2*(a + b*x**2)**2*(c + d*x**2)**2)
    F = -d**(sympy.S(5)/2)*(-3*a*d + 7*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(5)/2)*(-a*d + b*c)**3) + b/(2*a*x*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(a*d + b*c)/(2*a*c*x*(c + d*x**2)*(-a*d + b*c)**2) - (3*a**2*d**2 - 4*a*b*c*d + 3*b**2*c**2)/(2*a**2*c**2*x*(-a*d + b*c)**2) - b**(sympy.S(5)/2)*(-7*a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_308():
    f = 1/(x**3*(a + b*x**2)**2*(c + d*x**2)**2)
    F = -d**3/(2*c**2*(c + d*x**2)*(-a*d + b*c)**2) + d**3*(-a*d + 2*b*c)*log(c + d*x**2)/(c**3*(-a*d + b*c)**3) - b**3/(2*a**2*(a + b*x**2)*(-a*d + b*c)**2) - 1/(2*a**2*c**2*x**2) + b**3*(-2*a*d + b*c)*log(a + b*x**2)/(a**3*(-a*d + b*c)**3) - (2*a*d + 2*b*c)*log(x)/(a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_309():
    f = 1/(x**4*(a + b*x**2)**2*(c + d*x**2)**2)
    F = d**(sympy.S(7)/2)*(-5*a*d + 9*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(7)/2)*(-a*d + b*c)**3) + b/(2*a*x**3*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(a*d + b*c)/(2*a*c*x**3*(c + d*x**2)*(-a*d + b*c)**2) - (5*a**2*d**2 - 4*a*b*c*d + 5*b**2*c**2)/(6*a**2*c**2*x**3*(-a*d + b*c)**2) + (a*d + b*c)*(5*a**2*d**2 - 9*a*b*c*d + 5*b**2*c**2)/(2*a**3*c**3*x*(-a*d + b*c)**2) + b**(sympy.S(7)/2)*(-9*a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_310():
    f = x**4/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*sqrt(a)*sqrt(b)*(a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*(-a*d + b*c)**4) + a*x/(2*b*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + x*(9*a*d + 3*b*c)/(8*(c + d*x**2)*(-a*d + b*c)**3) + (3*a**2*d**2 + 18*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*sqrt(c)*sqrt(d)*(-a*d + b*c)**4) + x*(2*a*d + b*c)/(4*b*(c + d*x**2)**2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_311():
    f = x**3/((a + b*x**2)**2*(c + d*x**2)**3)
    F = a*b/(2*(a + b*x**2)*(-a*d + b*c)**3) + b*(2*a*d + b*c)*log(a + b*x**2)/(2*(-a*d + b*c)**4) - b*(2*a*d + b*c)*log(c + d*x**2)/(2*(-a*d + b*c)**4) + c/(4*(c + d*x**2)**2*(-a*d + b*c)**2) + (a*d + b*c)/(2*(c + d*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_312():
    f = x**2/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*d*x/(4*(c + d*x**2)**2*(-a*d + b*c)**2) - x/((a + b*x**2)*(c + d*x**2)**2*(-2*a*d + 2*b*c)) - d*x*(a*d + 11*b*c)/(8*c*(c + d*x**2)*(-a*d + b*c)**3) - sqrt(d)*(-a**2*d**2 + 10*a*b*c*d + 15*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(3)/2)*(-a*d + b*c)**4) + b**(sympy.S(3)/2)*(5*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*sqrt(a)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_313():
    f = x/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*b**2*d*log(a + b*x**2)/(2*(-a*d + b*c)**4) + 3*b**2*d*log(c + d*x**2)/(2*(-a*d + b*c)**4) - b**2/(2*(a + b*x**2)*(-a*d + b*c)**3) - b*d/((c + d*x**2)*(-a*d + b*c)**3) - d/(4*(c + d*x**2)**2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_314():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**3)
    F = d**(sympy.S(3)/2)*(3*a**2*d**2 - 14*a*b*c*d + 35*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**4) + b*x/(2*a*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*x*(a*d + 2*b*c)/(4*a*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d*x*(-a*d + 4*b*c)*(3*a*d + b*c)/(8*a*c**2*(c + d*x**2)*(-a*d + b*c)**3) + b**(sympy.S(5)/2)*(-7*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_315():
    f = 1/(x*(a + b*x**2)**2*(c + d*x**2)**3)
    F = d**2/(4*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d**2*(-a*d + 3*b*c)/(2*c**2*(c + d*x**2)*(-a*d + b*c)**3) - d**2*(a**2*d**2 - 4*a*b*c*d + 6*b**2*c**2)*log(c + d*x**2)/(2*c**3*(-a*d + b*c)**4) + b**3/(2*a*(a + b*x**2)*(-a*d + b*c)**3) - b**3*(-4*a*d + b*c)*log(a + b*x**2)/(2*a**2*(-a*d + b*c)**4) + log(x)/(a**2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_316():
    f = 1/(x**2*(a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*d**(sympy.S(5)/2)*(5*a**2*d**2 - 18*a*b*c*d + 21*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(7)/2)*(-a*d + b*c)**4) + b/(2*a*x*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(a*d + 2*b*c)/(4*a*c*x*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(-5*a**2*d**2 + 13*a*b*c*d + 4*b**2*c**2)/(8*a*c**2*x*(c + d*x**2)*(-a*d + b*c)**3) - (-3*a*d + 6*b*c)*(5*a**2*d**2 - 3*a*b*c*d + 2*b**2*c**2)/(8*a**2*c**3*x*(-a*d + b*c)**3) - 3*b**(sympy.S(7)/2)*(-3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_317():
    f = 1/(x**3*(a + b*x**2)**2*(c + d*x**2)**3)
    F = -d**3/(4*c**2*(c + d*x**2)**2*(-a*d + b*c)**2) - d**3*(-a*d + 2*b*c)/(c**3*(c + d*x**2)*(-a*d + b*c)**3) + d**3*(3*a**2*d**2 - 10*a*b*c*d + 10*b**2*c**2)*log(c + d*x**2)/(2*c**4*(-a*d + b*c)**4) - b**4/(2*a**2*(a + b*x**2)*(-a*d + b*c)**3) - 1/(2*a**2*c**3*x**2) + b**4*(-5*a*d + 2*b*c)*log(a + b*x**2)/(2*a**3*(-a*d + b*c)**4) - (3*a*d + 2*b*c)*log(x)/(a**3*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_318():
    f = 1/(x**4*(a + b*x**2)**2*(c + d*x**2)**3)
    F = d**(sympy.S(7)/2)*(35*a**2*d**2 - 110*a*b*c*d + 99*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(9)/2)*(-a*d + b*c)**4) + b/(2*a*x**3*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(a*d + 2*b*c)/(4*a*c*x**3*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(-7*a**2*d**2 + 15*a*b*c*d + 4*b**2*c**2)/(8*a*c**2*x**3*(c + d*x**2)*(-a*d + b*c)**3) - (-35*a**3*d**3 + 75*a**2*b*c*d**2 - 24*a*b**2*c**2*d + 20*b**3*c**3)/(24*a**2*c**3*x**3*(-a*d + b*c)**3) + (-35*a**4*d**4 + 75*a**3*b*c*d**3 - 24*a**2*b**2*c**2*d**2 - 24*a*b**3*c**3*d + 20*b**4*c**4)/(8*a**3*c**4*x*(-a*d + b*c)**3) + b**(sympy.S(9)/2)*(-11*a*d + 5*b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_319():
    f = x**m*(A + B*x**2)*(a + b*x**2)**3
    F = A*a**3*x**(m + 1)/(m + 1) + B*b**3*x**(m + 9)/(m + 9) + a**2*x**(m + 3)*(3*A*b + B*a)/(m + 3) + 3*a*b*x**(m + 5)*(A*b + B*a)/(m + 5) + b**2*x**(m + 7)*(A*b + 3*B*a)/(m + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_320():
    f = x**m*(A + B*x**2)*(a + b*x**2)**2
    F = A*a**2*x**(m + 1)/(m + 1) + B*b**2*x**(m + 7)/(m + 7) + a*x**(m + 3)*(2*A*b + B*a)/(m + 3) + b*x**(m + 5)*(A*b + 2*B*a)/(m + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_321():
    f = x**m*(A + B*x**2)*(a + b*x**2)
    F = A*a*x**(m + 1)/(m + 1) + B*b*x**(m + 5)/(m + 5) + x**(m + 3)*(A*b + B*a)/(m + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_322():
    f = x**m*(A + B*x**2)/(a + b*x**2)
    F = B*x**(m + 1)/(b*(m + 1)) + x**(m + 1)*(A*b - B*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_323():
    f = x**m*(A + B*x**2)/(a + b*x**2)**2
    F = x**(m + 1)*(A*b - B*a)/(2*a*b*(a + b*x**2)) + x**(m + 1)*(A*(-b*m + b) + B*a*(m + 1))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_324():
    f = x**m*(A + B*x**2)/(a + b*x**2)**3
    F = x**(m + 1)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) + x**(m + 1)*(A*b*(3 - m) + B*a*(m + 1))*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(4*a**3*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_325():
    f = x**m*(a + b*x**2)**2*(c + d*x**2)**3
    F = a**2*c**3*x**(m + 1)/(m + 1) + a*c**2*x**(m + 3)*(3*a*d + 2*b*c)/(m + 3) + b**2*d**3*x**(m + 11)/(m + 11) + b*d**2*x**(m + 9)*(2*a*d + 3*b*c)/(m + 9) + c*x**(m + 5)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/(m + 5) + d*x**(m + 7)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/(m + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_326():
    f = x**m*(a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x**(m + 1)/(m + 1) + 2*a*c*x**(m + 3)*(a*d + b*c)/(m + 3) + b**2*d**2*x**(m + 9)/(m + 9) + 2*b*d*x**(m + 7)*(a*d + b*c)/(m + 7) + x**(m + 5)*(a**2*d**2 + 4*a*b*c*d + b**2*c**2)/(m + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_327():
    f = x**m*(a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x**(m + 1)/(m + 1) + a*x**(m + 3)*(a*d + 2*b*c)/(m + 3) + b**2*d*x**(m + 7)/(m + 7) + b*x**(m + 5)*(2*a*d + b*c)/(m + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_328():
    f = x**m*(a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**(m + 3)/(d*(m + 3)) - b*x**(m + 1)*(-2*a*d + b*c)/(d**2*(m + 1)) + x**(m + 1)*(-a*d + b*c)**2*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_329():
    f = x**m*(a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x**(m + 1)/(d**2*(m + 1)) + x**(m + 1)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - x**(m + 1)*(-a*d + b*c)*(a*d*(1 - m) + b*c*(m + 3))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*d**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_330():
    f = x**m*(a + b*x**2)**2/(c + d*x**2)**3
    F = x**(m + 1)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - x**(m + 1)*(-a*d + b*c)*(a*d*(3 - m) + b*c*(m + 5))/(8*c**2*d**2*(c + d*x**2)) + x**(m + 1)*(a**2*d**2*(m**2 - 4*m + 3) + 2*a*b*c*d*(1 - m**2) + b**2*c**2*(m**2 + 4*m + 3))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*d**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_331():
    f = x**m*(c + d*x**2)**3/(a + b*x**2)
    F = d**3*x**(m + 5)/(b*(m + 5)) + d**2*x**(m + 3)*(-a*d + 3*b*c)/(b**2*(m + 3)) + d*x**(m + 1)*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(b**3*(m + 1)) + x**(m + 1)*(-a*d + b*c)**3*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_332():
    f = x**m*(c + d*x**2)**2/(a + b*x**2)
    F = d**2*x**(m + 3)/(b*(m + 3)) + d*x**(m + 1)*(-a*d + 2*b*c)/(b**2*(m + 1)) + x**(m + 1)*(-a*d + b*c)**2*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_333():
    f = x**m*(c + d*x**2)/(a + b*x**2)
    F = d*x**(m + 1)/(b*(m + 1)) + x**(m + 1)*(-a*d + b*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_334():
    f = x**m/((a + b*x**2)*(c + d*x**2))
    F = -d*x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*(m + 1)*(-a*d + b*c)) + b*x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_335():
    f = x**m/((a + b*x**2)**2*(c + d*x**2))
    F = d**2*x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*(m + 1)*(-a*d + b*c)**2) + b*x**(m + 1)/(2*a*(a + b*x**2)*(-a*d + b*c)) + b*x**(m + 1)*(-a*d*(3 - m) + b*c*(1 - m))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_336():
    f = x**m/((a + b*x**2)**3*(c + d*x**2))
    F = -d**3*x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*(m + 1)*(-a*d + b*c)**3) + b*x**(m + 1)/(4*a*(a + b*x**2)**2*(-a*d + b*c)) + b*x**(m + 1)*(-a*d*(7 - m) + b*c*(3 - m))/(8*a**2*(a + b*x**2)*(-a*d + b*c)**2) + b*x**(m + 1)*(a**2*d**2*(m**2 - 8*m + 15) - 2*a*b*c*d*(m**2 - 6*m + 5) + b**2*c**2*(m**2 - 4*m + 3))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_337():
    f = x**m*(c + d*x**2)**3/(a + b*x**2)**2
    F = x**(m + 1)*(c + d*x**2)**2*(-a*d + b*c)/(2*a*b*(a + b*x**2)) - d**2*x**(m + 3)*(-a*d*(m + 5) + b*c*(m + 3))/(2*a*b**2*(m + 3)) - d*x**(m + 1)*(a**2*d**2*(m + 5) - 3*a*b*c*d*(m + 3) + 2*b**2*c**2*(m + 1))/(2*a*b**3*(m + 1)) + x**(m + 1)*(-a*d + b*c)**2*(a*d*(m + 5) + b*(-c*m + c))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_338():
    f = x**m*(c + d*x**2)**2/(a + b*x**2)**2
    F = d**2*x**(m + 1)/(b**2*(m + 1)) + x**(m + 1)*(-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) + x**(m + 1)*(-a*d + b*c)*(a*d*(m + 3) + b*(-c*m + c))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_339():
    f = x**m*(c + d*x**2)/(a + b*x**2)**2
    F = x**(m + 1)*(-a*d + b*c)/(2*a*b*(a + b*x**2)) + x**(m + 1)*(a*d*(m + 1) + b*(-c*m + c))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_340():
    f = x**m/((a + b*x**2)**2*(c + d*x**2))
    F = d**2*x**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*(m + 1)*(-a*d + b*c)**2) + b*x**(m + 1)/(2*a*(a + b*x**2)*(-a*d + b*c)) + b*x**(m + 1)*(-a*d*(3 - m) + b*c*(1 - m))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_341():
    f = x**m/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -d**2*x**(m + 1)*(a*d*(1 - m) - b*c*(5 - m))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*(m + 1)*(-a*d + b*c)**3) + b*x**(m + 1)/(2*a*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*x**(m + 1)*(a*d + b*c)/(2*a*c*(c + d*x**2)*(-a*d + b*c)**2) - b**2*x**(m + 1)*(a*d*(5 - m) - b*(-c*m + c))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_342():
    f = x**m/((a + b*x**2)**2*(c + d*x**2)**3)
    F = d**2*x**(m + 1)*(a**2*d**2*(m**2 - 4*m + 3) - 2*a*b*c*d*(m**2 - 8*m + 7) + b**2*c**2*(m**2 - 12*m + 35))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*(m + 1)*(-a*d + b*c)**4) + b*x**(m + 1)/(2*a*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*x**(m + 1)*(a*d + 2*b*c)/(4*a*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d*x**(m + 1)*(-a**2*d**2*(3 - m) + a*b*c*d*(11 - m) + 4*b**2*c**2)/(8*a*c**2*(c + d*x**2)*(-a*d + b*c)**3) - b**3*x**(m + 1)*(a*d*(7 - m) - b*(-c*m + c))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*(m + 1)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_343():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(a + b*x**2)
    F = 2*A*a*x**(sympy.S(9)/2)/9 + 2*B*b*x**(sympy.S(17)/2)/17 + x**(sympy.S(13)/2)*(2*A*b + 2*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_344():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(a + b*x**2)
    F = 2*A*a*x**(sympy.S(7)/2)/7 + 2*B*b*x**(sympy.S(15)/2)/15 + x**(sympy.S(11)/2)*(2*A*b + 2*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_345():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(a + b*x**2)
    F = 2*A*a*x**(sympy.S(5)/2)/5 + 2*B*b*x**(sympy.S(13)/2)/13 + x**(sympy.S(9)/2)*(2*A*b + 2*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_346():
    f = sqrt(x)*(A + B*x**2)*(a + b*x**2)
    F = 2*A*a*x**(sympy.S(3)/2)/3 + 2*B*b*x**(sympy.S(11)/2)/11 + x**(sympy.S(7)/2)*(2*A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_347():
    f = (A + B*x**2)*(a + b*x**2)/sqrt(x)
    F = 2*A*a*sqrt(x) + 2*B*b*x**(sympy.S(9)/2)/9 + x**(sympy.S(5)/2)*(2*A*b + 2*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_348():
    f = (A + B*x**2)*(a + b*x**2)/x**(sympy.S(3)/2)
    F = -2*A*a/sqrt(x) + 2*B*b*x**(sympy.S(7)/2)/7 + x**(sympy.S(3)/2)*(2*A*b + 2*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_349():
    f = (A + B*x**2)*(a + b*x**2)/x**(sympy.S(5)/2)
    F = -2*A*a/(3*x**(sympy.S(3)/2)) + 2*B*b*x**(sympy.S(5)/2)/5 + sqrt(x)*(2*A*b + 2*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_350():
    f = (A + B*x**2)*(a + b*x**2)/x**(sympy.S(7)/2)
    F = -2*A*a/(5*x**(sympy.S(5)/2)) + 2*B*b*x**(sympy.S(3)/2)/3 - (2*A*b + 2*B*a)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_351():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(a + b*x**2)**2
    F = 2*A*a**2*x**(sympy.S(9)/2)/9 + 2*B*b**2*x**(sympy.S(21)/2)/21 + 2*a*x**(sympy.S(13)/2)*(2*A*b + B*a)/13 + 2*b*x**(sympy.S(17)/2)*(A*b + 2*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_352():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(a + b*x**2)**2
    F = 2*A*a**2*x**(sympy.S(7)/2)/7 + 2*B*b**2*x**(sympy.S(19)/2)/19 + 2*a*x**(sympy.S(11)/2)*(2*A*b + B*a)/11 + 2*b*x**(sympy.S(15)/2)*(A*b + 2*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_353():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(a + b*x**2)**2
    F = 2*A*a**2*x**(sympy.S(5)/2)/5 + 2*B*b**2*x**(sympy.S(17)/2)/17 + 2*a*x**(sympy.S(9)/2)*(2*A*b + B*a)/9 + 2*b*x**(sympy.S(13)/2)*(A*b + 2*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_354():
    f = sqrt(x)*(A + B*x**2)*(a + b*x**2)**2
    F = 2*A*a**2*x**(sympy.S(3)/2)/3 + 2*B*b**2*x**(sympy.S(15)/2)/15 + 2*a*x**(sympy.S(7)/2)*(2*A*b + B*a)/7 + 2*b*x**(sympy.S(11)/2)*(A*b + 2*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_355():
    f = (A + B*x**2)*(a + b*x**2)**2/sqrt(x)
    F = 2*A*a**2*sqrt(x) + 2*B*b**2*x**(sympy.S(13)/2)/13 + 2*a*x**(sympy.S(5)/2)*(2*A*b + B*a)/5 + 2*b*x**(sympy.S(9)/2)*(A*b + 2*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_356():
    f = (A + B*x**2)*(a + b*x**2)**2/x**(sympy.S(3)/2)
    F = -2*A*a**2/sqrt(x) + 2*B*b**2*x**(sympy.S(11)/2)/11 + 2*a*x**(sympy.S(3)/2)*(2*A*b + B*a)/3 + 2*b*x**(sympy.S(7)/2)*(A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_357():
    f = (A + B*x**2)*(a + b*x**2)**2/x**(sympy.S(5)/2)
    F = -2*A*a**2/(3*x**(sympy.S(3)/2)) + 2*B*b**2*x**(sympy.S(9)/2)/9 + 2*a*sqrt(x)*(2*A*b + B*a) + 2*b*x**(sympy.S(5)/2)*(A*b + 2*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_358():
    f = (A + B*x**2)*(a + b*x**2)**2/x**(sympy.S(7)/2)
    F = -2*A*a**2/(5*x**(sympy.S(5)/2)) + 2*B*b**2*x**(sympy.S(7)/2)/7 - 2*a*(2*A*b + B*a)/sqrt(x) + 2*b*x**(sympy.S(3)/2)*(A*b + 2*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_359():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(a + b*x**2)**3
    F = 2*A*a**3*x**(sympy.S(9)/2)/9 + 2*B*b**3*x**(sympy.S(25)/2)/25 + 2*a**2*x**(sympy.S(13)/2)*(3*A*b + B*a)/13 + 6*a*b*x**(sympy.S(17)/2)*(A*b + B*a)/17 + 2*b**2*x**(sympy.S(21)/2)*(A*b + 3*B*a)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_360():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(a + b*x**2)**3
    F = 2*A*a**3*x**(sympy.S(7)/2)/7 + 2*B*b**3*x**(sympy.S(23)/2)/23 + 2*a**2*x**(sympy.S(11)/2)*(3*A*b + B*a)/11 + 2*a*b*x**(sympy.S(15)/2)*(A*b + B*a)/5 + 2*b**2*x**(sympy.S(19)/2)*(A*b + 3*B*a)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_361():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(a + b*x**2)**3
    F = 2*A*a**3*x**(sympy.S(5)/2)/5 + 2*B*b**3*x**(sympy.S(21)/2)/21 + 2*a**2*x**(sympy.S(9)/2)*(3*A*b + B*a)/9 + 6*a*b*x**(sympy.S(13)/2)*(A*b + B*a)/13 + 2*b**2*x**(sympy.S(17)/2)*(A*b + 3*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_362():
    f = sqrt(x)*(A + B*x**2)*(a + b*x**2)**3
    F = 2*A*a**3*x**(sympy.S(3)/2)/3 + 2*B*b**3*x**(sympy.S(19)/2)/19 + 2*a**2*x**(sympy.S(7)/2)*(3*A*b + B*a)/7 + 6*a*b*x**(sympy.S(11)/2)*(A*b + B*a)/11 + 2*b**2*x**(sympy.S(15)/2)*(A*b + 3*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_363():
    f = (A + B*x**2)*(a + b*x**2)**3/sqrt(x)
    F = 2*A*a**3*sqrt(x) + 2*B*b**3*x**(sympy.S(17)/2)/17 + 2*a**2*x**(sympy.S(5)/2)*(3*A*b + B*a)/5 + 2*a*b*x**(sympy.S(9)/2)*(A*b + B*a)/3 + 2*b**2*x**(sympy.S(13)/2)*(A*b + 3*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_364():
    f = (A + B*x**2)*(a + b*x**2)**3/x**(sympy.S(3)/2)
    F = -2*A*a**3/sqrt(x) + 2*B*b**3*x**(sympy.S(15)/2)/15 + 2*a**2*x**(sympy.S(3)/2)*(3*A*b + B*a)/3 + 6*a*b*x**(sympy.S(7)/2)*(A*b + B*a)/7 + 2*b**2*x**(sympy.S(11)/2)*(A*b + 3*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_365():
    f = (A + B*x**2)*(a + b*x**2)**3/x**(sympy.S(5)/2)
    F = -2*A*a**3/(3*x**(sympy.S(3)/2)) + 2*B*b**3*x**(sympy.S(13)/2)/13 + 2*a**2*sqrt(x)*(3*A*b + B*a) + 6*a*b*x**(sympy.S(5)/2)*(A*b + B*a)/5 + 2*b**2*x**(sympy.S(9)/2)*(A*b + 3*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_366():
    f = (A + B*x**2)*(a + b*x**2)**3/x**(sympy.S(7)/2)
    F = -2*A*a**3/(5*x**(sympy.S(5)/2)) + 2*B*b**3*x**(sympy.S(11)/2)/11 - 2*a**2*(3*A*b + B*a)/sqrt(x) + 2*a*b*x**(sympy.S(3)/2)*(A*b + B*a) + 2*b**2*x**(sympy.S(7)/2)*(A*b + 3*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_367():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(a + b*x**2)
    F = 2*B*x**(sympy.S(9)/2)/(9*b) - sqrt(2)*a**(sympy.S(5)/4)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(13)/4)) + sqrt(2)*a**(sympy.S(5)/4)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(13)/4)) - sqrt(2)*a**(sympy.S(5)/4)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4)) + sqrt(2)*a**(sympy.S(5)/4)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4)) - 2*a*sqrt(x)*(A*b - B*a)/b**3 + x**(sympy.S(5)/2)*(2*A*b - 2*B*a)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_368():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(a + b*x**2)
    F = 2*B*x**(sympy.S(7)/2)/(7*b) - sqrt(2)*a**(sympy.S(3)/4)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(11)/4)) + sqrt(2)*a**(sympy.S(3)/4)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(11)/4)) + sqrt(2)*a**(sympy.S(3)/4)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4)) - sqrt(2)*a**(sympy.S(3)/4)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4)) + x**(sympy.S(3)/2)*(2*A*b - 2*B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_369():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(a + b*x**2)
    F = 2*B*x**(sympy.S(5)/2)/(5*b) + sqrt(2)*a**(sympy.S(1)/4)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(9)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(9)/4)) + sqrt(2)*a**(sympy.S(1)/4)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) + sqrt(x)*(2*A*b - 2*B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_370():
    f = sqrt(x)*(A + B*x**2)/(a + b*x**2)
    F = 2*B*x**(sympy.S(3)/2)/(3*b) + sqrt(2)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_371():
    f = (A + B*x**2)/(sqrt(x)*(a + b*x**2))
    F = 2*B*sqrt(x)/b - sqrt(2)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_372():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(a + b*x**2))
    F = -2*A/(a*sqrt(x)) - sqrt(2)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_373():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*(a + b*x**2))
    F = -2*A/(3*a*x**(sympy.S(3)/2)) + sqrt(2)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + sqrt(2)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_374():
    f = (A + B*x**2)/(x**(sympy.S(7)/2)*(a + b*x**2))
    F = -2*A/(5*a*x**(sympy.S(5)/2)) + (2*A*b - 2*B*a)/(a**2*sqrt(x)) + sqrt(2)*b**(sympy.S(1)/4)*(A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_375():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(a + b*x**2)**2
    F = sqrt(2)*a**(sympy.S(1)/4)*(5*A*b - 9*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(13)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(5*A*b - 9*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(13)/4)) + sqrt(2)*a**(sympy.S(1)/4)*(5*A*b - 9*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(5*A*b - 9*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4)) + sqrt(x)*(5*A*b - 9*B*a)/(2*b**3) + x**(sympy.S(9)/2)*(A*b - B*a)/(2*a*b*(a + b*x**2)) - x**(sympy.S(5)/2)*(5*A*b - 9*B*a)/(10*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_376():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(a + b*x**2)**2
    F = x**(sympy.S(7)/2)*(A*b - B*a)/(2*a*b*(a + b*x**2)) - x**(sympy.S(3)/2)*(3*A*b - 7*B*a)/(6*a*b**2) + sqrt(2)*(3*A*b - 7*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(3*A*b - 7*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(3*A*b - 7*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(3*A*b - 7*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_377():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(a + b*x**2)**2
    F = x**(sympy.S(5)/2)*(A*b - B*a)/(2*a*b*(a + b*x**2)) - sqrt(x)*(A*b - 5*B*a)/(2*a*b**2) - sqrt(2)*(A*b - 5*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(A*b - 5*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(A*b - 5*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(A*b - 5*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_378():
    f = sqrt(x)*(A + B*x**2)/(a + b*x**2)**2
    F = x**(sympy.S(3)/2)*(A*b - B*a)/(2*a*b*(a + b*x**2)) + sqrt(2)*(A*b + 3*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(A*b + 3*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(A*b + 3*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(A*b + 3*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_379():
    f = (A + B*x**2)/(sqrt(x)*(a + b*x**2)**2)
    F = sqrt(x)*(A*b - B*a)/(2*a*b*(a + b*x**2)) - sqrt(2)*(3*A*b + B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(3*A*b + B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(3*A*b + B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(3*A*b + B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_380():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(a + b*x**2)**2)
    F = (A*b - B*a)/(2*a*b*sqrt(x)*(a + b*x**2)) - (5*A*b - B*a)/(2*a**2*b*sqrt(x)) - sqrt(2)*(5*A*b - B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(5*A*b - B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(5*A*b - B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(5*A*b - B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_381():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*(a + b*x**2)**2)
    F = (A*b - B*a)/(2*a*b*x**(sympy.S(3)/2)*(a + b*x**2)) - (7*A*b - 3*B*a)/(6*a**2*b*x**(sympy.S(3)/2)) + sqrt(2)*(7*A*b - 3*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(7*A*b - 3*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + sqrt(2)*(7*A*b - 3*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(7*A*b - 3*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_382():
    f = (A + B*x**2)/(x**(sympy.S(7)/2)*(a + b*x**2)**2)
    F = (A*b - B*a)/(2*a*b*x**(sympy.S(5)/2)*(a + b*x**2)) - (9*A*b - 5*B*a)/(10*a**2*b*x**(sympy.S(5)/2)) + (9*A*b - 5*B*a)/(2*a**3*sqrt(x)) + sqrt(2)*b**(sympy.S(1)/4)*(9*A*b - 5*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(9*A*b - 5*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(9*A*b - 5*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(9*A*b - 5*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_383():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(a + b*x**2)**3
    F = x**(sympy.S(9)/2)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) + x**(sympy.S(5)/2)*(A*b - 9*B*a)/(16*a*b**2*(a + b*x**2)) - sqrt(x)*(5*A*b - 45*B*a)/(16*a*b**3) - sqrt(2)*(5*A*b - 45*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(5*A*b - 45*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) - sqrt(2)*(5*A*b - 45*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(5*A*b - 45*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_384():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(a + b*x**2)**3
    F = x**(sympy.S(7)/2)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) - x**(sympy.S(3)/2)*(A*b + 7*B*a)/(16*a*b**2*(a + b*x**2)) + sqrt(2)*(3*A*b + 21*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(3*A*b + 21*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(3*A*b + 21*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(3*A*b + 21*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_385():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(a + b*x**2)**3
    F = x**(sympy.S(5)/2)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) - sqrt(x)*(3*A*b + 5*B*a)/(16*a*b**2*(a + b*x**2)) - sqrt(2)*(3*A*b + 5*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(3*A*b + 5*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(3*A*b + 5*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(3*A*b + 5*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_386():
    f = sqrt(x)*(A + B*x**2)/(a + b*x**2)**3
    F = x**(sympy.S(3)/2)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) + x**(sympy.S(3)/2)*(5*A*b + 3*B*a)/(16*a**2*b*(a + b*x**2)) + sqrt(2)*(5*A*b + 3*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(5*A*b + 3*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(5*A*b + 3*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(5*A*b + 3*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_387():
    f = (A + B*x**2)/(sqrt(x)*(a + b*x**2)**3)
    F = sqrt(x)*(A*b - B*a)/(4*a*b*(a + b*x**2)**2) + sqrt(x)*(7*A*b + B*a)/(16*a**2*b*(a + b*x**2)) - sqrt(2)*(21*A*b + 3*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(21*A*b + 3*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(21*A*b + 3*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(21*A*b + 3*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_388():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(a + b*x**2)**3)
    F = (A*b - B*a)/(4*a*b*sqrt(x)*(a + b*x**2)**2) + (9*A*b - B*a)/(16*a**2*b*sqrt(x)*(a + b*x**2)) - (45*A*b - 5*B*a)/(16*a**3*b*sqrt(x)) - sqrt(2)*(45*A*b - 5*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(45*A*b - 5*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(45*A*b - 5*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(45*A*b - 5*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_389():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*(a + b*x**2)**3)
    F = (A*b - B*a)/(4*a*b*x**(sympy.S(3)/2)*(a + b*x**2)**2) + (11*A*b - 3*B*a)/(16*a**2*b*x**(sympy.S(3)/2)*(a + b*x**2)) - (77*A*b - 21*B*a)/(48*a**3*b*x**(sympy.S(3)/2)) + sqrt(2)*(77*A*b - 21*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(77*A*b - 21*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) + sqrt(2)*(77*A*b - 21*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(77*A*b - 21*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(15)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_390():
    f = (A + B*x**2)/(x**(sympy.S(7)/2)*(a + b*x**2)**3)
    F = (A*b - B*a)/(4*a*b*x**(sympy.S(5)/2)*(a + b*x**2)**2) + (13*A*b - 5*B*a)/(16*a**2*b*x**(sympy.S(5)/2)*(a + b*x**2)) - (117*A*b - 45*B*a)/(80*a**3*b*x**(sympy.S(5)/2)) + (117*A*b - 45*B*a)/(16*a**4*sqrt(x)) + 9*sqrt(2)*b**(sympy.S(1)/4)*(13*A*b - 5*B*a)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(17)/4)) - 9*sqrt(2)*b**(sympy.S(1)/4)*(13*A*b - 5*B*a)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(128*a**(sympy.S(17)/4)) - 9*sqrt(2)*b**(sympy.S(1)/4)*(13*A*b - 5*B*a)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(17)/4)) + 9*sqrt(2)*b**(sympy.S(1)/4)*(13*A*b - 5*B*a)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_391():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2*(c + d*x**2)
    F = 2*a**2*c*x**(sympy.S(9)/2)/9 + 2*a*x**(sympy.S(13)/2)*(a*d + 2*b*c)/13 + 2*b**2*d*x**(sympy.S(21)/2)/21 + 2*b*x**(sympy.S(17)/2)*(2*a*d + b*c)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_392():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)
    F = 2*a**2*c*x**(sympy.S(7)/2)/7 + 2*a*x**(sympy.S(11)/2)*(a*d + 2*b*c)/11 + 2*b**2*d*x**(sympy.S(19)/2)/19 + 2*b*x**(sympy.S(15)/2)*(2*a*d + b*c)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_393():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)
    F = 2*a**2*c*x**(sympy.S(5)/2)/5 + 2*a*x**(sympy.S(9)/2)*(a*d + 2*b*c)/9 + 2*b**2*d*x**(sympy.S(17)/2)/17 + 2*b*x**(sympy.S(13)/2)*(2*a*d + b*c)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_394():
    f = sqrt(x)*(a + b*x**2)**2*(c + d*x**2)
    F = 2*a**2*c*x**(sympy.S(3)/2)/3 + 2*a*x**(sympy.S(7)/2)*(a*d + 2*b*c)/7 + 2*b**2*d*x**(sympy.S(15)/2)/15 + 2*b*x**(sympy.S(11)/2)*(2*a*d + b*c)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_395():
    f = (a + b*x**2)**2*(c + d*x**2)/sqrt(x)
    F = 2*a**2*c*sqrt(x) + 2*a*x**(sympy.S(5)/2)*(a*d + 2*b*c)/5 + 2*b**2*d*x**(sympy.S(13)/2)/13 + 2*b*x**(sympy.S(9)/2)*(2*a*d + b*c)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_396():
    f = (a + b*x**2)**2*(c + d*x**2)/x**(sympy.S(3)/2)
    F = -2*a**2*c/sqrt(x) + 2*a*x**(sympy.S(3)/2)*(a*d + 2*b*c)/3 + 2*b**2*d*x**(sympy.S(11)/2)/11 + 2*b*x**(sympy.S(7)/2)*(2*a*d + b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_397():
    f = (a + b*x**2)**2*(c + d*x**2)/x**(sympy.S(5)/2)
    F = -2*a**2*c/(3*x**(sympy.S(3)/2)) + 2*a*sqrt(x)*(a*d + 2*b*c) + 2*b**2*d*x**(sympy.S(9)/2)/9 + 2*b*x**(sympy.S(5)/2)*(2*a*d + b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_398():
    f = (a + b*x**2)**2*(c + d*x**2)/x**(sympy.S(7)/2)
    F = -2*a**2*c/(5*x**(sympy.S(5)/2)) - 2*a*(a*d + 2*b*c)/sqrt(x) + 2*b**2*d*x**(sympy.S(7)/2)/7 + 2*b*x**(sympy.S(3)/2)*(2*a*d + b*c)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_399():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2*(c + d*x**2)**2
    F = 2*a**2*c**2*x**(sympy.S(9)/2)/9 + 4*a*c*x**(sympy.S(13)/2)*(a*d + b*c)/13 + 2*b**2*d**2*x**(sympy.S(25)/2)/25 + 4*b*d*x**(sympy.S(21)/2)*(a*d + b*c)/21 + x**(sympy.S(17)/2)*(2*a**2*d**2/17 + 8*a*b*c*d/17 + 2*b**2*c**2/17)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_400():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)**2
    F = 2*a**2*c**2*x**(sympy.S(7)/2)/7 + 4*a*c*x**(sympy.S(11)/2)*(a*d + b*c)/11 + 2*b**2*d**2*x**(sympy.S(23)/2)/23 + 4*b*d*x**(sympy.S(19)/2)*(a*d + b*c)/19 + x**(sympy.S(15)/2)*(2*a**2*d**2/15 + 8*a*b*c*d/15 + 2*b**2*c**2/15)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_401():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)**2
    F = 2*a**2*c**2*x**(sympy.S(5)/2)/5 + 4*a*c*x**(sympy.S(9)/2)*(a*d + b*c)/9 + 2*b**2*d**2*x**(sympy.S(21)/2)/21 + 4*b*d*x**(sympy.S(17)/2)*(a*d + b*c)/17 + x**(sympy.S(13)/2)*(2*a**2*d**2/13 + 8*a*b*c*d/13 + 2*b**2*c**2/13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_402():
    f = sqrt(x)*(a + b*x**2)**2*(c + d*x**2)**2
    F = 2*a**2*c**2*x**(sympy.S(3)/2)/3 + 4*a*c*x**(sympy.S(7)/2)*(a*d + b*c)/7 + 2*b**2*d**2*x**(sympy.S(19)/2)/19 + 4*b*d*x**(sympy.S(15)/2)*(a*d + b*c)/15 + x**(sympy.S(11)/2)*(2*a**2*d**2/11 + 8*a*b*c*d/11 + 2*b**2*c**2/11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_403():
    f = (a + b*x**2)**2*(c + d*x**2)**2/sqrt(x)
    F = 2*a**2*c**2*sqrt(x) + 4*a*c*x**(sympy.S(5)/2)*(a*d + b*c)/5 + 2*b**2*d**2*x**(sympy.S(17)/2)/17 + 4*b*d*x**(sympy.S(13)/2)*(a*d + b*c)/13 + x**(sympy.S(9)/2)*(2*a**2*d**2/9 + 8*a*b*c*d/9 + 2*b**2*c**2/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_404():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**(sympy.S(3)/2)
    F = -2*a**2*c**2/sqrt(x) + 4*a*c*x**(sympy.S(3)/2)*(a*d + b*c)/3 + 2*b**2*d**2*x**(sympy.S(15)/2)/15 + 4*b*d*x**(sympy.S(11)/2)*(a*d + b*c)/11 + x**(sympy.S(7)/2)*(2*a**2*d**2/7 + 8*a*b*c*d/7 + 2*b**2*c**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_405():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**(sympy.S(5)/2)
    F = -2*a**2*c**2/(3*x**(sympy.S(3)/2)) + 4*a*c*sqrt(x)*(a*d + b*c) + 2*b**2*d**2*x**(sympy.S(13)/2)/13 + 4*b*d*x**(sympy.S(9)/2)*(a*d + b*c)/9 + x**(sympy.S(5)/2)*(2*a**2*d**2/5 + 8*a*b*c*d/5 + 2*b**2*c**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_406():
    f = (a + b*x**2)**2*(c + d*x**2)**2/x**(sympy.S(7)/2)
    F = -2*a**2*c**2/(5*x**(sympy.S(5)/2)) - 4*a*c*(a*d + b*c)/sqrt(x) + 2*b**2*d**2*x**(sympy.S(11)/2)/11 + 4*b*d*x**(sympy.S(7)/2)*(a*d + b*c)/7 + x**(sympy.S(3)/2)*(2*a**2*d**2/3 + 8*a*b*c*d/3 + 2*b**2*c**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_407():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2*(c + d*x**2)**3
    F = 2*a**2*c**3*x**(sympy.S(9)/2)/9 + 2*a*c**2*x**(sympy.S(13)/2)*(3*a*d + 2*b*c)/13 + 2*b**2*d**3*x**(sympy.S(29)/2)/29 + 2*b*d**2*x**(sympy.S(25)/2)*(2*a*d + 3*b*c)/25 + 2*c*x**(sympy.S(17)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/17 + 2*d*x**(sympy.S(21)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_408():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)**3
    F = 2*a**2*c**3*x**(sympy.S(7)/2)/7 + 2*a*c**2*x**(sympy.S(11)/2)*(3*a*d + 2*b*c)/11 + 2*b**2*d**3*x**(sympy.S(27)/2)/27 + 2*b*d**2*x**(sympy.S(23)/2)*(2*a*d + 3*b*c)/23 + 2*c*x**(sympy.S(15)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/15 + 2*d*x**(sympy.S(19)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_409():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)**3
    F = 2*a**2*c**3*x**(sympy.S(5)/2)/5 + 2*a*c**2*x**(sympy.S(9)/2)*(3*a*d + 2*b*c)/9 + 2*b**2*d**3*x**(sympy.S(25)/2)/25 + 2*b*d**2*x**(sympy.S(21)/2)*(2*a*d + 3*b*c)/21 + 2*c*x**(sympy.S(13)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/13 + 2*d*x**(sympy.S(17)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_410():
    f = sqrt(x)*(a + b*x**2)**2*(c + d*x**2)**3
    F = 2*a**2*c**3*x**(sympy.S(3)/2)/3 + 2*a*c**2*x**(sympy.S(7)/2)*(3*a*d + 2*b*c)/7 + 2*b**2*d**3*x**(sympy.S(23)/2)/23 + 2*b*d**2*x**(sympy.S(19)/2)*(2*a*d + 3*b*c)/19 + 2*c*x**(sympy.S(11)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/11 + 2*d*x**(sympy.S(15)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_411():
    f = (a + b*x**2)**2*(c + d*x**2)**3/sqrt(x)
    F = 2*a**2*c**3*sqrt(x) + 2*a*c**2*x**(sympy.S(5)/2)*(3*a*d + 2*b*c)/5 + 2*b**2*d**3*x**(sympy.S(21)/2)/21 + 2*b*d**2*x**(sympy.S(17)/2)*(2*a*d + 3*b*c)/17 + 2*c*x**(sympy.S(9)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/9 + 2*d*x**(sympy.S(13)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_412():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**(sympy.S(3)/2)
    F = -2*a**2*c**3/sqrt(x) + 2*a*c**2*x**(sympy.S(3)/2)*(3*a*d + 2*b*c)/3 + 2*b**2*d**3*x**(sympy.S(19)/2)/19 + 2*b*d**2*x**(sympy.S(15)/2)*(2*a*d + 3*b*c)/15 + 2*c*x**(sympy.S(7)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/7 + 2*d*x**(sympy.S(11)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_413():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**(sympy.S(5)/2)
    F = -2*a**2*c**3/(3*x**(sympy.S(3)/2)) + 2*a*c**2*sqrt(x)*(3*a*d + 2*b*c) + 2*b**2*d**3*x**(sympy.S(17)/2)/17 + 2*b*d**2*x**(sympy.S(13)/2)*(2*a*d + 3*b*c)/13 + 2*c*x**(sympy.S(5)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/5 + 2*d*x**(sympy.S(9)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_414():
    f = (a + b*x**2)**2*(c + d*x**2)**3/x**(sympy.S(7)/2)
    F = -2*a**2*c**3/(5*x**(sympy.S(5)/2)) - 2*a*c**2*(3*a*d + 2*b*c)/sqrt(x) + 2*b**2*d**3*x**(sympy.S(15)/2)/15 + 2*b*d**2*x**(sympy.S(11)/2)*(2*a*d + 3*b*c)/11 + 2*c*x**(sympy.S(3)/2)*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/3 + 2*d*x**(sympy.S(7)/2)*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_415():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2/(c + d*x**2)
    F = 2*b**2*x**(sympy.S(13)/2)/(13*d) - 2*b*x**(sympy.S(9)/2)*(-2*a*d + b*c)/(9*d**2) - sqrt(2)*c**(sympy.S(5)/4)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(17)/4)) + sqrt(2)*c**(sympy.S(5)/4)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(17)/4)) - sqrt(2)*c**(sympy.S(5)/4)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(17)/4)) + sqrt(2)*c**(sympy.S(5)/4)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(17)/4)) - 2*c*sqrt(x)*(-a*d + b*c)**2/d**4 + 2*x**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_416():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2/(c + d*x**2)
    F = 2*b**2*x**(sympy.S(11)/2)/(11*d) - 2*b*x**(sympy.S(7)/2)*(-2*a*d + b*c)/(7*d**2) - sqrt(2)*c**(sympy.S(3)/4)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(15)/4)) + sqrt(2)*c**(sympy.S(3)/4)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(15)/4)) + sqrt(2)*c**(sympy.S(3)/4)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(15)/4)) - sqrt(2)*c**(sympy.S(3)/4)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(15)/4)) + 2*x**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_417():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2/(c + d*x**2)
    F = 2*b**2*x**(sympy.S(9)/2)/(9*d) - 2*b*x**(sympy.S(5)/2)*(-2*a*d + b*c)/(5*d**2) + sqrt(2)*c**(sympy.S(1)/4)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(13)/4)) + 2*sqrt(x)*(-a*d + b*c)**2/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_418():
    f = sqrt(x)*(a + b*x**2)**2/(c + d*x**2)
    F = 2*b**2*x**(sympy.S(7)/2)/(7*d) - 2*b*x**(sympy.S(3)/2)*(-2*a*d + b*c)/(3*d**2) + sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(1)/4)*d**(sympy.S(11)/4)) - sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(1)/4)*d**(sympy.S(11)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*d**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*d**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_419():
    f = (a + b*x**2)**2/(sqrt(x)*(c + d*x**2))
    F = 2*b**2*x**(sympy.S(5)/2)/(5*d) - 2*b*sqrt(x)*(-2*a*d + b*c)/d**2 - sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_420():
    f = (a + b*x**2)**2/(x**(sympy.S(3)/2)*(c + d*x**2))
    F = -2*a**2/(c*sqrt(x)) + 2*b**2*x**(sympy.S(3)/2)/(3*d) - sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(5)/4)*d**(sympy.S(7)/4)) + sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(5)/4)*d**(sympy.S(7)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*d**(sympy.S(7)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*d**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_421():
    f = (a + b*x**2)**2/(x**(sympy.S(5)/2)*(c + d*x**2))
    F = -2*a**2/(3*c*x**(sympy.S(3)/2)) + 2*b**2*sqrt(x)/d + sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_422():
    f = (a + b*x**2)**2/(x**(sympy.S(7)/2)*(c + d*x**2))
    F = -2*a**2/(5*c*x**(sympy.S(5)/2)) - 2*a*(-a*d + 2*b*c)/(c**2*sqrt(x)) + sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(9)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(9)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)*d**(sympy.S(3)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)*d**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_423():
    f = (a + b*x**2)**2/(x**(sympy.S(9)/2)*(c + d*x**2))
    F = -2*a**2/(7*c*x**(sympy.S(7)/2)) - 2*a*(-a*d + 2*b*c)/(3*c**2*x**(sympy.S(3)/2)) - sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(11)/4)*d**(sympy.S(1)/4)) + sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(11)/4)*d**(sympy.S(1)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)*d**(sympy.S(1)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)*d**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_424():
    f = (c + d*x**2)**2/(x**(sympy.S(11)/2)*(a + b*x**2))
    F = -2*c**2/(9*a*x**(sympy.S(9)/2)) + 2*c*(-2*a*d + b*c)/(5*a**2*x**(sympy.S(5)/2)) - 2*(-a*d + b*c)**2/(a**3*sqrt(x)) - sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**2*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**2*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**2*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**2*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_425():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2/(c + d*x**2)**2
    F = 2*b**2*x**(sympy.S(9)/2)/(9*d**2) + sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + 13*b*c)*(-a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(17)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + 13*b*c)*(-a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(17)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + 13*b*c)*(-a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(17)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + 13*b*c)*(-a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(17)/4)) + sqrt(x)*(-5*a*d + 13*b*c)*(-a*d + b*c)/(2*d**4) + x**(sympy.S(9)/2)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - x**(sympy.S(5)/2)*(-5*a*d + 13*b*c)*(-a*d + b*c)/(10*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_426():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2/(c + d*x**2)**2
    F = 2*b**2*x**(sympy.S(7)/2)/(7*d**2) + x**(sympy.S(7)/2)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - x**(sympy.S(3)/2)*(-3*a*d + 11*b*c)*(-a*d + b*c)/(6*c*d**3) + sqrt(2)*(-3*a*d + 11*b*c)*(-a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*d**(sympy.S(15)/4)) - sqrt(2)*(-3*a*d + 11*b*c)*(-a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*d**(sympy.S(15)/4)) - sqrt(2)*(-3*a*d + 11*b*c)*(-a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*d**(sympy.S(15)/4)) + sqrt(2)*(-3*a*d + 11*b*c)*(-a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*d**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_427():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2/(c + d*x**2)**2
    F = 2*b**2*x**(sympy.S(5)/2)/(5*d**2) + x**(sympy.S(5)/2)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - sqrt(x)*(-a*d + b*c)*(-a*d + 9*b*c)/(2*c*d**3) - sqrt(2)*(-a*d + b*c)*(-a*d + 9*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*d**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)*(-a*d + 9*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*d**(sympy.S(13)/4)) - sqrt(2)*(-a*d + b*c)*(-a*d + 9*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*d**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)*(-a*d + 9*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*d**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_428():
    f = sqrt(x)*(a + b*x**2)**2/(c + d*x**2)**2
    F = 2*b**2*x**(sympy.S(3)/2)/(3*d**2) + x**(sympy.S(3)/2)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - sqrt(2)*(-a*d + b*c)*(a*d + 7*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*d**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)*(a*d + 7*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*d**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)*(a*d + 7*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*d**(sympy.S(11)/4)) - sqrt(2)*(-a*d + b*c)*(a*d + 7*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*d**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_429():
    f = (a + b*x**2)**2/(sqrt(x)*(c + d*x**2)**2)
    F = 2*b**2*sqrt(x)/d**2 + sqrt(x)*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) + sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_430():
    f = (a + b*x**2)**2/(x**(sympy.S(3)/2)*(c + d*x**2)**2)
    F = -2*a**2/(c*sqrt(x)*(c + d*x**2)) - x**(sympy.S(3)/2)*(5*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(2*c**2*d*(c + d*x**2)) + sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*d**(sympy.S(7)/4)) - sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*d**(sympy.S(7)/4)) - sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*d**(sympy.S(7)/4)) + sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*d**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_431():
    f = (a + b*x**2)**2/(x**(sympy.S(5)/2)*(c + d*x**2)**2)
    F = -2*a**2/(3*c*x**(sympy.S(3)/2)*(c + d*x**2)) - sqrt(x)*(7*a**2*d**2 - 6*a*b*c*d + 3*b**2*c**2)/(6*c**2*d*(c + d*x**2)) - sqrt(2)*(-a*d + b*c)*(7*a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)*(7*a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)*(7*a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)*(7*a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_432():
    f = (a + b*x**2)**2/(x**(sympy.S(7)/2)*(c + d*x**2)**2)
    F = -2*a**2/(5*c*x**(sympy.S(5)/2)*(c + d*x**2)) - (9*a**2*d**2 - 10*a*b*c*d + 5*b**2*c**2)/(10*c**2*d*sqrt(x)*(c + d*x**2)) + (-9*a*d + b*c)*(-a*d + b*c)/(2*c**3*d*sqrt(x)) + sqrt(2)*(-9*a*d + b*c)*(-a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-9*a*d + b*c)*(-a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-9*a*d + b*c)*(-a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*d**(sympy.S(3)/4)) + sqrt(2)*(-9*a*d + b*c)*(-a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*d**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_433():
    f = x**(sympy.S(7)/2)*(a + b*x**2)**2/(c + d*x**2)**3
    F = x**(sympy.S(9)/2)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - sqrt(x)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)/(16*c*d**4) - x**(sympy.S(9)/2)*(-a*d + b*c)*(-a*d + 17*b*c)/(16*c**2*d**2*(c + d*x**2)) + x**(sympy.S(5)/2)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)/(80*c**2*d**3) - sqrt(2)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(17)/4)) + sqrt(2)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(17)/4)) - sqrt(2)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(17)/4)) + sqrt(2)*(5*a**2*d**2 - 90*a*b*c*d + 117*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_434():
    f = x**(sympy.S(5)/2)*(a + b*x**2)**2/(c + d*x**2)**3
    F = x**(sympy.S(7)/2)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - x**(sympy.S(3)/2)*(3*a**2*d/c + 42*a*b - 77*b**2*c/d)/(48*c*d**2) - x**(sympy.S(7)/2)*(-a*d + b*c)*(a*d + 15*b*c)/(16*c**2*d**2*(c + d*x**2)) - sqrt(2)*(-3*a**2*d**2 - 42*a*b*c*d + 77*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*d**(sympy.S(15)/4)) + sqrt(2)*(-3*a**2*d**2 - 42*a*b*c*d + 77*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*d**(sympy.S(15)/4)) + sqrt(2)*(-3*a**2*d**2 - 42*a*b*c*d + 77*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*d**(sympy.S(15)/4)) - sqrt(2)*(-3*a**2*d**2 - 42*a*b*c*d + 77*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*d**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_435():
    f = x**(sympy.S(3)/2)*(a + b*x**2)**2/(c + d*x**2)**3
    F = x**(sympy.S(5)/2)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - sqrt(x)*(3*a**2*d/c + 10*a*b - 45*b**2*c/d)/(16*c*d**2) - x**(sympy.S(5)/2)*(-a*d + b*c)*(3*a*d + 13*b*c)/(16*c**2*d**2*(c + d*x**2)) + sqrt(2)*(-3*a**2*d**2 - 10*a*b*c*d + 45*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*d**(sympy.S(13)/4)) - sqrt(2)*(-3*a**2*d**2 - 10*a*b*c*d + 45*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*d**(sympy.S(13)/4)) + sqrt(2)*(-3*a**2*d**2 - 10*a*b*c*d + 45*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*d**(sympy.S(13)/4)) - sqrt(2)*(-3*a**2*d**2 - 10*a*b*c*d + 45*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*d**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_436():
    f = sqrt(x)*(a + b*x**2)**2/(c + d*x**2)**3
    F = x**(sympy.S(3)/2)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - x**(sympy.S(3)/2)*(-a*d + b*c)*(5*a*d + 11*b*c)/(16*c**2*d**2*(c + d*x**2)) + sqrt(2)*(5*a**2*d**2 + 6*a*b*c*d + 21*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*d**(sympy.S(11)/4)) - sqrt(2)*(5*a**2*d**2 + 6*a*b*c*d + 21*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*d**(sympy.S(11)/4)) - sqrt(2)*(5*a**2*d**2 + 6*a*b*c*d + 21*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*d**(sympy.S(11)/4)) + sqrt(2)*(5*a**2*d**2 + 6*a*b*c*d + 21*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*d**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_437():
    f = (a + b*x**2)**2/(sqrt(x)*(c + d*x**2)**3)
    F = sqrt(x)*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**2)**2) - sqrt(x)*(-a*d + b*c)*(7*a*d + 9*b*c)/(16*c**2*d**2*(c + d*x**2)) - sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_438():
    f = (a + b*x**2)**2/(x**(sympy.S(3)/2)*(c + d*x**2)**3)
    F = -2*a**2/(c*sqrt(x)*(c + d*x**2)**2) - x**(sympy.S(3)/2)*(9*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(4*c**2*d*(c + d*x**2)**2) + x**(sympy.S(3)/2)*(5*a*d*(-9*a*d + 2*b*c) + 3*b**2*c**2)/(16*c**3*d*(c + d*x**2)) + sqrt(2)*(5*a*d*(-9*a*d + 2*b*c) + 3*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*d**(sympy.S(7)/4)) - sqrt(2)*(5*a*d*(-9*a*d + 2*b*c) + 3*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*d**(sympy.S(7)/4)) - sqrt(2)*(5*a*d*(-9*a*d + 2*b*c) + 3*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*d**(sympy.S(7)/4)) + sqrt(2)*(5*a*d*(-9*a*d + 2*b*c) + 3*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*d**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_439():
    f = (a + b*x**2)**2/(x**(sympy.S(5)/2)*(c + d*x**2)**3)
    F = -2*a**2/(3*c*x**(sympy.S(3)/2)*(c + d*x**2)**2) - sqrt(x)*(11*a**2*d**2 - 6*a*b*c*d + 3*b**2*c**2)/(12*c**2*d*(c + d*x**2)**2) + sqrt(x)*(7*a*d*(-11*a*d + 6*b*c) + 3*b**2*c**2)/(48*c**3*d*(c + d*x**2)) - sqrt(2)*(7*a*d*(-11*a*d + 6*b*c) + 3*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(7*a*d*(-11*a*d + 6*b*c) + 3*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(7*a*d*(-11*a*d + 6*b*c) + 3*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(7*a*d*(-11*a*d + 6*b*c) + 3*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_440():
    f = (a + b*x**2)**2/(x**(sympy.S(7)/2)*(c + d*x**2)**3)
    F = -2*a**2/(5*c*x**(sympy.S(5)/2)*(c + d*x**2)**2) - (13*a**2*d**2 - 10*a*b*c*d + 5*b**2*c**2)/(20*c**2*d*sqrt(x)*(c + d*x**2)**2) - (-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)/(80*c**3*d*sqrt(x)*(c + d*x**2)) + (-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)/(16*c**4*d*sqrt(x)) + sqrt(2)*(-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*d**(sympy.S(3)/4)) - sqrt(2)*(-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*d**(sympy.S(3)/4)) + sqrt(2)*(-9*a*d*(-13*a*d + 10*b*c) + 5*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*d**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_441():
    f = x**(sympy.S(5)/2)*(c + d*x**2)**3/(a + b*x**2)
    F = -sqrt(2)*a**(sympy.S(3)/4)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(19)/4)) + sqrt(2)*a**(sympy.S(3)/4)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(19)/4)) + sqrt(2)*a**(sympy.S(3)/4)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(19)/4)) - sqrt(2)*a**(sympy.S(3)/4)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(19)/4)) + 2*d**3*x**(sympy.S(15)/2)/(15*b) + 2*d**2*x**(sympy.S(11)/2)*(-a*d + 3*b*c)/(11*b**2) + 2*d*x**(sympy.S(7)/2)*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(7*b**3) + 2*x**(sympy.S(3)/2)*(-a*d + b*c)**3/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_442():
    f = x**(sympy.S(3)/2)*(c + d*x**2)**3/(a + b*x**2)
    F = sqrt(2)*a**(sympy.S(1)/4)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(17)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(17)/4)) + sqrt(2)*a**(sympy.S(1)/4)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(17)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(17)/4)) + 2*d**3*x**(sympy.S(13)/2)/(13*b) + 2*d**2*x**(sympy.S(9)/2)*(-a*d + 3*b*c)/(9*b**2) + 2*d*x**(sympy.S(5)/2)*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(5*b**3) + 2*sqrt(x)*(-a*d + b*c)**3/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_443():
    f = sqrt(x)*(c + d*x**2)**3/(a + b*x**2)
    F = 2*d**3*x**(sympy.S(11)/2)/(11*b) + 2*d**2*x**(sympy.S(7)/2)*(-a*d + 3*b*c)/(7*b**2) + 2*d*x**(sympy.S(3)/2)*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/(3*b**3) + sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) - sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(15)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_444():
    f = (c + d*x**2)**3/(sqrt(x)*(a + b*x**2))
    F = 2*d**3*x**(sympy.S(9)/2)/(9*b) + 2*d**2*x**(sympy.S(5)/2)*(-a*d + 3*b*c)/(5*b**2) + 2*d*sqrt(x)*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/b**3 - sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_445():
    f = (c + d*x**2)**3/(x**(sympy.S(3)/2)*(a + b*x**2))
    F = 2*d**3*x**(sympy.S(7)/2)/(7*b) + 2*d**2*x**(sympy.S(3)/2)*(-a*d + 3*b*c)/(3*b**2) - 2*c**3/(a*sqrt(x)) - sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_446():
    f = (c + d*x**2)**3/(x**(sympy.S(5)/2)*(a + b*x**2))
    F = 2*d**3*x**(sympy.S(5)/2)/(5*b) + 2*d**2*sqrt(x)*(-a*d + 3*b*c)/b**2 - 2*c**3/(3*a*x**(sympy.S(3)/2)) + sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_447():
    f = (c + d*x**2)**3/(x**(sympy.S(7)/2)*(a + b*x**2))
    F = 2*d**3*x**(sympy.S(3)/2)/(3*b) - 2*c**3/(5*a*x**(sympy.S(5)/2)) + 2*c**2*(-3*a*d + b*c)/(a**2*sqrt(x)) + sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_448():
    f = (c + d*x**2)**3/(x**(sympy.S(9)/2)*(a + b*x**2))
    F = 2*d**3*sqrt(x)/b - 2*c**3/(7*a*x**(sympy.S(7)/2)) + 2*c**2*(-3*a*d + b*c)/(3*a**2*x**(sympy.S(3)/2)) - sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_449():
    f = (c + d*x**2)**3/(x**(sympy.S(11)/2)*(a + b*x**2))
    F = -2*c**3/(9*a*x**(sympy.S(9)/2)) + 2*c**2*(-3*a*d + b*c)/(5*a**2*x**(sympy.S(5)/2)) - 2*c*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/(a**3*sqrt(x)) - sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(13)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(13)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_450():
    f = (c + d*x**2)**3/(x**(sympy.S(13)/2)*(a + b*x**2))
    F = -2*c**3/(11*a*x**(sympy.S(11)/2)) + 2*c**2*(-3*a*d + b*c)/(7*a**2*x**(sympy.S(7)/2)) - 2*c*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/(3*a**3*x**(sympy.S(3)/2)) + sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(15)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_451():
    f = (c + d*x**2)**3/(x**(sympy.S(15)/2)*(a + b*x**2))
    F = -2*c**3/(13*a*x**(sympy.S(13)/2)) + 2*c**2*(-3*a*d + b*c)/(9*a**2*x**(sympy.S(9)/2)) - 2*c*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/(5*a**3*x**(sympy.S(5)/2)) + 2*(-a*d + b*c)**3/(a**4*sqrt(x)) + sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(17)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(17)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(17)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_452():
    f = x**(sympy.S(7)/2)*(c + d*x**2)**3/(a + b*x**2)**2
    F = sqrt(2)*a**(sympy.S(1)/4)*(-17*a*d + 5*b*c)*(-a*d + b*c)**2*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(21)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(-17*a*d + 5*b*c)*(-a*d + b*c)**2*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(21)/4)) + sqrt(2)*a**(sympy.S(1)/4)*(-17*a*d + 5*b*c)*(-a*d + b*c)**2*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(21)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(-17*a*d + 5*b*c)*(-a*d + b*c)**2*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(21)/4)) - x**(sympy.S(5)/2)*(c + d*x**2)**3/(2*b*(a + b*x**2)) + 17*d**3*x**(sympy.S(13)/2)/(26*b**2) + d**2*x**(sympy.S(9)/2)*(-17*a*d + 39*b*c)/(18*b**3) + d*x**(sympy.S(5)/2)*(17*a**2*d**2 - 39*a*b*c*d + 27*b**2*c**2)/(10*b**4) + sqrt(x)*(-17*a*d + 5*b*c)*(-a*d + b*c)**2/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_453():
    f = x**(sympy.S(5)/2)*(c + d*x**2)**3/(a + b*x**2)**2
    F = -x**(sympy.S(3)/2)*(c + d*x**2)**3/(2*b*(a + b*x**2)) + 15*d**3*x**(sympy.S(11)/2)/(22*b**2) + 3*d**2*x**(sympy.S(7)/2)*(-5*a*d + 11*b*c)/(14*b**3) + d*x**(sympy.S(3)/2)*(5*a**2*d**2 - 11*a*b*c*d + 7*b**2*c**2)/(2*b**4) + sqrt(2)*(-15*a*d + 3*b*c)*(-a*d + b*c)**2*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)) - sqrt(2)*(-15*a*d + 3*b*c)*(-a*d + b*c)**2*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)) - sqrt(2)*(-15*a*d + 3*b*c)*(-a*d + b*c)**2*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(19)/4)) + sqrt(2)*(-15*a*d + 3*b*c)*(-a*d + b*c)**2*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(19)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_454():
    f = x**(sympy.S(3)/2)*(c + d*x**2)**3/(a + b*x**2)**2
    F = -sqrt(x)*(c + d*x**2)**3/(2*b*(a + b*x**2)) + 13*d*sqrt(x)*(c + d*x**2)**2/(18*b**2) + d*sqrt(x)*(c + d*x**2)*(-117*a*d + 113*b*c)/(90*b**3) + d*sqrt(x)*(585*a**2*d**2 - 1098*a*b*c*d + 497*b**2*c**2)/(90*b**4) - sqrt(2)*(-13*a*d + b*c)*(-a*d + b*c)**2*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-13*a*d + b*c)*(-a*d + b*c)**2*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) - sqrt(2)*(-13*a*d + b*c)*(-a*d + b*c)**2*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-13*a*d + b*c)*(-a*d + b*c)**2*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_455():
    f = sqrt(x)*(c + d*x**2)**3/(a + b*x**2)**2
    F = x**(sympy.S(3)/2)*(c + d*x**2)**2*(-a*d + b*c)/(2*a*b*(a + b*x**2)) - d**2*x**(sympy.S(7)/2)*(-11*a*d + 7*b*c)/(14*a*b**2) - d*x**(sympy.S(3)/2)*(11*a**2*d**2 - 21*a*b*c*d + 6*b**2*c**2)/(6*a*b**3) + sqrt(2)*(-a*d + b*c)**2*(11*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)) - sqrt(2)*(-a*d + b*c)**2*(11*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)) - sqrt(2)*(-a*d + b*c)**2*(11*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(15)/4)) + sqrt(2)*(-a*d + b*c)**2*(11*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_456():
    f = (c + d*x**2)**3/(sqrt(x)*(a + b*x**2)**2)
    F = 2*d**3*x**(sympy.S(5)/2)/(5*b**2) + 2*d**2*sqrt(x)*(-2*a*d + 3*b*c)/b**3 + sqrt(x)*(-a*d + b*c)**3/(2*a*b**3*(a + b*x**2)) - 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) + 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) - 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) + 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_457():
    f = (c + d*x**2)**3/(x**(sympy.S(3)/2)*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*sqrt(x)*(a + b*x**2)) - d**2*x**(sympy.S(3)/2)*(-7*a*d + 3*b*c)/(6*a*b**2) - c**2*(-a*d + 5*b*c)/(2*a**2*b*sqrt(x)) - sqrt(2)*(-a*d + b*c)**2*(7*a*d + 5*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)**2*(7*a*d + 5*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)) + sqrt(2)*(-a*d + b*c)**2*(7*a*d + 5*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*b**(sympy.S(11)/4)) - sqrt(2)*(-a*d + b*c)**2*(7*a*d + 5*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_458():
    f = (c + d*x**2)**3/(x**(sympy.S(5)/2)*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*x**(sympy.S(3)/2)*(a + b*x**2)) - d**2*sqrt(x)*(-5*a*d + b*c)/(2*a*b**2) - c**2*(-3*a*d + 7*b*c)/(6*a**2*b*x**(sympy.S(3)/2)) + sqrt(2)*(-a*d + b*c)**2*(5*a*d + 7*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**2*(5*a*d + 7*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*(5*a*d + 7*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**2*(5*a*d + 7*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_459():
    f = (c + d*x**2)**3/(x**(sympy.S(7)/2)*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*x**(sympy.S(5)/2)*(a + b*x**2)) - c**2*(-5*a*d + 9*b*c)/(10*a**2*b*x**(sympy.S(5)/2)) + c*(2*a**2*d**2 - 15*a*b*c*d + 9*b**2*c**2)/(2*a**3*b*sqrt(x)) + 3*sqrt(2)*(-a*d + b*c)**2*(a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*(-a*d + b*c)**2*(a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)) - 3*sqrt(2)*(-a*d + b*c)**2*(a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)) + 3*sqrt(2)*(-a*d + b*c)**2*(a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_460():
    f = (c + d*x**2)**3/(x**(sympy.S(9)/2)*(a + b*x**2)**2)
    F = (c + d*x**2)**2*(-a*d + b*c)/(2*a*b*x**(sympy.S(7)/2)*(a + b*x**2)) - c**2*(-7*a*d + 11*b*c)/(14*a**2*b*x**(sympy.S(7)/2)) + c*(6*a**2*d**2 - 21*a*b*c*d + 11*b**2*c**2)/(6*a**3*b*x**(sympy.S(3)/2)) - sqrt(2)*(-a*d + b*c)**2*(a*d + 11*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)**2*(a*d + 11*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)**2*(a*d + 11*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)**2*(a*d + 11*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(15)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_461():
    f = x**(sympy.S(9)/2)/((a + b*x**2)*(c + d*x**2))
    F = sqrt(2)*a**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(7)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(7)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(7)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(7)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(7)/4)*(-a*d + b*c)) + 2*x**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_462():
    f = x**(sympy.S(7)/2)/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(2)*a**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(5)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(5)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(5)/4)*(-a*d + b*c)) + 2*sqrt(x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_463():
    f = x**(sympy.S(5)/2)/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(2)*a**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_464():
    f = x**(sympy.S(3)/2)/((a + b*x**2)*(c + d*x**2))
    F = sqrt(2)*a**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(1)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(1)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*d**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(1)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(1)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*d**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_465():
    f = sqrt(x)/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(2)*d**(sympy.S(1)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(1)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(1)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(1)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_466():
    f = 1/(sqrt(x)*(a + b*x**2)*(c + d*x**2))
    F = sqrt(2)*d**(sympy.S(3)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(3)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_467():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)*(c + d*x**2))
    F = sqrt(2)*d**(sympy.S(5)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(5)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(5)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(5)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(5)/4)*(-a*d + b*c)) - 2/(a*c*sqrt(x)) - sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_468():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)*(c + d*x**2))
    F = -sqrt(2)*d**(sympy.S(7)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(7)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(7)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(7)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(7)/4)*(-a*d + b*c)) - 2/(3*a*c*x**(sympy.S(3)/2)) + sqrt(2)*b**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_469():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)*(c + d*x**2))
    F = -sqrt(2)*d**(sympy.S(9)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(9)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(4*c**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(9)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(9)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)*(-a*d + b*c)) - 2/(5*a*c*x**(sympy.S(5)/2)) + (2*a*d + 2*b*c)/(a**2*c**2*sqrt(x)) + sqrt(2)*b**(sympy.S(9)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(9)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(9)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(9)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_470():
    f = x**(sympy.S(11)/2)/((a + b*x**2)*(c + d*x**2)**2)
    F = sqrt(2)*a**(sympy.S(9)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(9)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(9)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(9)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(5)/4)*(-9*a*d + 5*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(5)/4)*(-9*a*d + 5*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(9)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(5)/4)*(-9*a*d + 5*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(5)/4)*(-9*a*d + 5*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(9)/4)*(-a*d + b*c)**2) - c*x**(sympy.S(5)/2)/(2*d*(c + d*x**2)*(-a*d + b*c)) + sqrt(x)*(-4*a*d + 5*b*c)/(2*b*d**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_471():
    f = x**(sympy.S(9)/2)/((a + b*x**2)*(c + d*x**2)**2)
    F = sqrt(2)*a**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(7)/4)*(-a*d + b*c)**2) - c*x**(sympy.S(3)/2)/(2*d*(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_472():
    f = x**(sympy.S(7)/2)/((a + b*x**2)*(c + d*x**2)**2)
    F = -sqrt(2)*a**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(1)/4)*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*b**(sympy.S(1)/4)*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(5)/4)*(-a*d + b*c)**2) - sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*c**(sympy.S(1)/4)*(-5*a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(5)/4)*(-a*d + b*c)**2) - c*sqrt(x)/(2*d*(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_473():
    f = x**(sympy.S(5)/2)/((a + b*x**2)*(c + d*x**2)**2)
    F = -sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**2) + x**(sympy.S(3)/2)/((c + d*x**2)*(-2*a*d + 2*b*c)) + sqrt(2)*(3*a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*(3*a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*(3*a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*(3*a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_474():
    f = x**(sympy.S(3)/2)/((a + b*x**2)*(c + d*x**2)**2)
    F = sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**2) + sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**2) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**2) + sqrt(x)/((c + d*x**2)*(-2*a*d + 2*b*c)) - sqrt(2)*(a*d + 3*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**2) + sqrt(2)*(a*d + 3*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**2) - sqrt(2)*(a*d + 3*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**2) + sqrt(2)*(a*d + 3*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_475():
    f = sqrt(x)/((a + b*x**2)*(c + d*x**2)**2)
    F = -d*x**(sympy.S(3)/2)/(2*c*(c + d*x**2)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(1)/4)*(-a*d + 5*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(1)/4)*(-a*d + 5*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(1)/4)*(-a*d + 5*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(1)/4)*(-a*d + 5*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_476():
    f = 1/(sqrt(x)*(a + b*x**2)*(c + d*x**2)**2)
    F = -d*sqrt(x)/(2*c*(c + d*x**2)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_477():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*sqrt(x)*(c + d*x**2)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(5)/4)*(-5*a*d + 9*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(5)/4)*(-5*a*d + 9*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(5)/4)*(-5*a*d + 9*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(5)/4)*(-5*a*d + 9*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*(-a*d + b*c)**2) - (-5*a*d + 4*b*c)/(2*a*c**2*sqrt(x)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(9)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(9)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(9)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(9)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_478():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*x**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(7)/4)*(-7*a*d + 11*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(7)/4)*(-7*a*d + 11*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(7)/4)*(-7*a*d + 11*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(7)/4)*(-7*a*d + 11*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*(-a*d + b*c)**2) - (-7*a*d + 4*b*c)/(6*a*c**2*x**(sympy.S(3)/2)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(11)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(11)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(11)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(11)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_479():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)*(c + d*x**2)**2)
    F = -d/(2*c*x**(sympy.S(5)/2)*(c + d*x**2)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(9)/4)*(-9*a*d + 13*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(9)/4)*(-9*a*d + 13*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(9)/4)*(-9*a*d + 13*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(9)/4)*(-9*a*d + 13*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*(-a*d + b*c)**2) - (-9*a*d + 4*b*c)/(10*a*c**2*x**(sympy.S(5)/2)*(-a*d + b*c)) + (-9*a**2*d**2 + 4*a*b*c*d + 4*b**2*c**2)/(2*a**2*c**3*sqrt(x)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(13)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(13)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(13)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(13)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_480():
    f = x**(sympy.S(7)/2)/((a + b*x**2)*(c + d*x**2)**3)
    F = -sqrt(2)*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) - c*sqrt(x)/(4*d*(c + d*x**2)**2*(-a*d + b*c)) + sqrt(x)*(-9*a*d + b*c)/(16*d*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*(-5*a**2*d**2 - 30*a*b*c*d + 3*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*(-5*a**2*d**2 - 30*a*b*c*d + 3*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*(-5*a**2*d**2 - 30*a*b*c*d + 3*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*(-5*a**2*d**2 - 30*a*b*c*d + 3*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_481():
    f = x**(sympy.S(5)/2)/((a + b*x**2)*(c + d*x**2)**3)
    F = -sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) + x**(sympy.S(3)/2)/((c + d*x**2)**2*(-4*a*d + 4*b*c)) + x**(sympy.S(3)/2)*(3*a*d + 5*b*c)/(16*c*(c + d*x**2)*(-a*d + b*c)**2) + sqrt(2)*(-3*a**2*d**2 + 30*a*b*c*d + 5*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*(-3*a**2*d**2 + 30*a*b*c*d + 5*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*(-3*a**2*d**2 + 30*a*b*c*d + 5*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*(-3*a**2*d**2 + 30*a*b*c*d + 5*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*d**(sympy.S(3)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_482():
    f = x**(sympy.S(3)/2)/((a + b*x**2)*(c + d*x**2)**3)
    F = sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*(-a*d + b*c)**3) + sqrt(x)/((c + d*x**2)**2*(-4*a*d + 4*b*c)) + sqrt(x)*(a*d + 7*b*c)/(16*c*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*(-3*a**2*d**2 + 14*a*b*c*d + 21*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*(-3*a**2*d**2 + 14*a*b*c*d + 21*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*(-3*a**2*d**2 + 14*a*b*c*d + 21*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*(-3*a**2*d**2 + 14*a*b*c*d + 21*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_483():
    f = sqrt(x)/((a + b*x**2)*(c + d*x**2)**3)
    F = -d*x**(sympy.S(3)/2)/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - d*x**(sympy.S(3)/2)*(-5*a*d + 13*b*c)/(16*c**2*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(1)/4)*(5*a**2*d**2 - 18*a*b*c*d + 45*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(1)/4)*(5*a**2*d**2 - 18*a*b*c*d + 45*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(1)/4)*(5*a**2*d**2 - 18*a*b*c*d + 45*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(1)/4)*(5*a**2*d**2 - 18*a*b*c*d + 45*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(9)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(9)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(9)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(9)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_484():
    f = 1/(sqrt(x)*(a + b*x**2)*(c + d*x**2)**3)
    F = -d*sqrt(x)/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - d*sqrt(x)*(-7*a*d + 15*b*c)/(16*c**2*(c + d*x**2)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(3)/4)*(21*a**2*d**2 - 66*a*b*c*d + 77*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(3)/4)*(21*a**2*d**2 - 66*a*b*c*d + 77*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(3)/4)*(21*a**2*d**2 - 66*a*b*c*d + 77*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(3)/4)*(21*a**2*d**2 - 66*a*b*c*d + 77*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(11)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(11)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(11)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(11)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_485():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*sqrt(x)*(c + d*x**2)**2*(-a*d + b*c)) - d*(-9*a*d + 17*b*c)/(16*c**2*sqrt(x)*(c + d*x**2)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(5)/4)*(45*a**2*d**2 - 130*a*b*c*d + 117*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(5)/4)*(45*a**2*d**2 - 130*a*b*c*d + 117*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(5)/4)*(45*a**2*d**2 - 130*a*b*c*d + 117*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(5)/4)*(45*a**2*d**2 - 130*a*b*c*d + 117*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*(-a*d + b*c)**3) - (45*a**2*d**2 - 85*a*b*c*d + 32*b**2*c**2)/(16*a*c**3*sqrt(x)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(13)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(13)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(13)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(13)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_486():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*x**(sympy.S(3)/2)*(c + d*x**2)**2*(-a*d + b*c)) - d*(-11*a*d + 19*b*c)/(16*c**2*x**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(7)/4)*(77*a**2*d**2 - 210*a*b*c*d + 165*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(77*a**2*d**2 - 210*a*b*c*d + 165*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(7)/4)*(77*a**2*d**2 - 210*a*b*c*d + 165*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(77*a**2*d**2 - 210*a*b*c*d + 165*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*(-a*d + b*c)**3) - (77*a**2*d**2 - 133*a*b*c*d + 32*b**2*c**2)/(48*a*c**3*x**(sympy.S(3)/2)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(15)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(15)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(15)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(15)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_487():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)*(c + d*x**2)**3)
    F = -d/(4*c*x**(sympy.S(5)/2)*(c + d*x**2)**2*(-a*d + b*c)) - d*(-13*a*d + 21*b*c)/(16*c**2*x**(sympy.S(5)/2)*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(9)/4)*(117*a**2*d**2 - 306*a*b*c*d + 221*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(9)/4)*(117*a**2*d**2 - 306*a*b*c*d + 221*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(9)/4)*(117*a**2*d**2 - 306*a*b*c*d + 221*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(9)/4)*(117*a**2*d**2 - 306*a*b*c*d + 221*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*(-a*d + b*c)**3) - (117*a**2*d**2 - 189*a*b*c*d + 32*b**2*c**2)/(80*a*c**3*x**(sympy.S(5)/2)*(-a*d + b*c)**2) + (117*a**3*d**3 - 189*a**2*b*c*d**2 + 32*a*b**2*c**2*d + 32*b**3*c**3)/(16*a**2*c**4*sqrt(x)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(17)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(17)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(4*a**(sympy.S(9)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(17)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(17)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(9)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_488():
    f = x**(sympy.S(7)/2)/((a + b*x**2)**2*(c + d*x**2)**2)
    F = sqrt(2)*a**(sympy.S(1)/4)*(3*a*d + 5*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(1)/4)*(3*a*d + 5*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*b**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*a**(sympy.S(1)/4)*(3*a*d + 5*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*a**(sympy.S(1)/4)*(3*a*d + 5*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*(-a*d + b*c)**3) + a*sqrt(x)/(2*b*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(1)/4)*(5*a*d + 3*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*c**(sympy.S(1)/4)*(5*a*d + 3*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*d**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*c**(sympy.S(1)/4)*(5*a*d + 3*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*c**(sympy.S(1)/4)*(5*a*d + 3*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*d**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(x)*(a*d + b*c)/(2*b*(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_489():
    f = x**(sympy.S(5)/2)/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -d*x**(sympy.S(3)/2)/((c + d*x**2)*(-a*d + b*c)**2) - x**(sympy.S(3)/2)/((a + b*x**2)*(c + d*x**2)*(-2*a*d + 2*b*c)) - sqrt(2)*d**(sympy.S(1)/4)*(3*a*d + 5*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(1)/4)*(3*a*d + 5*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(1)/4)*(3*a*d + 5*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(1)/4)*(3*a*d + 5*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(1)/4)*(5*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(1)/4)*(5*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(1)/4)*(5*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(1)/4)*(5*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_490():
    f = x**(sympy.S(3)/2)/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -d*sqrt(x)/((c + d*x**2)*(-a*d + b*c)**2) - sqrt(x)/((a + b*x**2)*(c + d*x**2)*(-2*a*d + 2*b*c)) + sqrt(2)*d**(sympy.S(3)/4)*(a*d + 7*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(3)/4)*(a*d + 7*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(3)/4)*(a*d + 7*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(3)/4)*(a*d + 7*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(3)/4)*(7*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(3)/4)*(7*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(3)/4)*(7*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(3)/4)*(7*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_491():
    f = sqrt(x)/((a + b*x**2)**2*(c + d*x**2)**2)
    F = sqrt(2)*d**(sympy.S(5)/4)*(-a*d + 9*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(5)/4)*(-a*d + 9*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(5)/4)*(-a*d + 9*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(5)/4)*(-a*d + 9*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(5)/4)*(-a*d + b*c)**3) + b*x**(sympy.S(3)/2)/(2*a*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*x**(sympy.S(3)/2)*(a*d + b*c)/(2*a*c*(c + d*x**2)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(5)/4)*(-9*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(5)/4)*(-9*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(5)/4)*(-9*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(5)/4)*(-9*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_492():
    f = 1/(sqrt(x)*(a + b*x**2)**2*(c + d*x**2)**2)
    F = -sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + b*sqrt(x)/(2*a*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*sqrt(x)*(a*d + b*c)/(2*a*c*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_493():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)**2)
    F = -sqrt(2)*d**(sympy.S(9)/4)*(-5*a*d + 13*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(9)/4)*(-5*a*d + 13*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(9)/4)*(-5*a*d + 13*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(9)/4)*(-5*a*d + 13*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(9)/4)*(-a*d + b*c)**3) + b/(2*a*sqrt(x)*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(a*d + b*c)/(2*a*c*sqrt(x)*(c + d*x**2)*(-a*d + b*c)**2) + (-5*a**2*d**2 + 8*a*b*c*d - 5*b**2*c**2)/(2*a**2*c**2*sqrt(x)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + 5*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + 5*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + 5*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + 5*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_494():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)**2)
    F = sqrt(2)*d**(sympy.S(11)/4)*(-7*a*d + 15*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(11)/4)*(-7*a*d + 15*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(11)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(11)/4)*(-7*a*d + 15*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(11)/4)*(-7*a*d + 15*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(11)/4)*(-a*d + b*c)**3) + b/(2*a*x**(sympy.S(3)/2)*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(a*d + b*c)/(2*a*c*x**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)**2) + (-7*a**2*d**2 + 8*a*b*c*d - 7*b**2*c**2)/(6*a**2*c**2*x**(sympy.S(3)/2)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(11)/4)*(-15*a*d + 7*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(11)/4)*(-15*a*d + 7*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(11)/4)*(-15*a*d + 7*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(11)/4)*(-15*a*d + 7*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_495():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)**2*(c + d*x**2)**2)
    F = sqrt(2)*d**(sympy.S(13)/4)*(-9*a*d + 17*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(13)/4)*(-9*a*d + 17*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(16*c**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(13)/4)*(-9*a*d + 17*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(13)/4)*(-9*a*d + 17*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)*(-a*d + b*c)**3) + b/(2*a*x**(sympy.S(5)/2)*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(a*d + b*c)/(2*a*c*x**(sympy.S(5)/2)*(c + d*x**2)*(-a*d + b*c)**2) + (-9*a**2*d**2 + 8*a*b*c*d - 9*b**2*c**2)/(10*a**2*c**2*x**(sympy.S(5)/2)*(-a*d + b*c)**2) + (a*d + b*c)*(9*a**2*d**2 - 17*a*b*c*d + 9*b**2*c**2)/(2*a**3*c**3*sqrt(x)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 9*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 9*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 9*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 9*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_496():
    f = x**(sympy.S(7)/2)/((a + b*x**2)**2*(c + d*x**2)**3)
    F = sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*(7*a*d + 5*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*(-a*d + b*c)**4) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*(7*a*d + 5*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*(-a*d + b*c)**4) + sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*(7*a*d + 5*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*(-a*d + b*c)**4) - sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*(7*a*d + 5*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*(-a*d + b*c)**4) + a*sqrt(x)/(2*b*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + sqrt(x)*(17*a*d + 7*b*c)/(16*(c + d*x**2)*(-a*d + b*c)**3) - sqrt(2)*(5*a**2*d**2 + 70*a*b*c*d + 21*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**4) + sqrt(2)*(5*a**2*d**2 + 70*a*b*c*d + 21*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**4) - sqrt(2)*(5*a**2*d**2 + 70*a*b*c*d + 21*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**4) + sqrt(2)*(5*a**2*d**2 + 70*a*b*c*d + 21*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(3)/4)*d**(sympy.S(1)/4)*(-a*d + b*c)**4) + sqrt(x)*(2*a*d + b*c)/(4*b*(c + d*x**2)**2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_497():
    f = x**(sympy.S(5)/2)/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*d*x**(sympy.S(3)/2)/(4*(c + d*x**2)**2*(-a*d + b*c)**2) - x**(sympy.S(3)/2)/((a + b*x**2)*(c + d*x**2)**2*(-2*a*d + 2*b*c)) - 3*d*x**(sympy.S(3)/2)*(a*d + 7*b*c)/(16*c*(c + d*x**2)*(-a*d + b*c)**3) - 3*sqrt(2)*d**(sympy.S(1)/4)*(-a**2*d**2 + 18*a*b*c*d + 15*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*d**(sympy.S(1)/4)*(-a**2*d**2 + 18*a*b*c*d + 15*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(5)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*d**(sympy.S(1)/4)*(-a**2*d**2 + 18*a*b*c*d + 15*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*d**(sympy.S(1)/4)*(-a**2*d**2 + 18*a*b*c*d + 15*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(5)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*b**(sympy.S(5)/4)*(3*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*b**(sympy.S(5)/4)*(3*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(1)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*b**(sympy.S(5)/4)*(3*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*b**(sympy.S(5)/4)*(3*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_498():
    f = x**(sympy.S(3)/2)/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*d*sqrt(x)/(4*(c + d*x**2)**2*(-a*d + b*c)**2) - sqrt(x)/((a + b*x**2)*(c + d*x**2)**2*(-2*a*d + 2*b*c)) - d*sqrt(x)*(a*d + 23*b*c)/(16*c*(c + d*x**2)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a**2*d**2 + 22*a*b*c*d + 77*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a**2*d**2 + 22*a*b*c*d + 77*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(7)/4)*(-a*d + b*c)**4) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a**2*d**2 + 22*a*b*c*d + 77*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a**2*d**2 + 22*a*b*c*d + 77*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(7)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(7)/4)*(11*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(7)/4)*(11*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(3)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(7)/4)*(11*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(7)/4)*(11*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_499():
    f = sqrt(x)/((a + b*x**2)**2*(c + d*x**2)**3)
    F = sqrt(2)*d**(sympy.S(5)/4)*(5*a**2*d**2 - 26*a*b*c*d + 117*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(5)/4)*(5*a**2*d**2 - 26*a*b*c*d + 117*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(9)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(5)/4)*(5*a**2*d**2 - 26*a*b*c*d + 117*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*(-a*d + b*c)**4) + sqrt(2)*d**(sympy.S(5)/4)*(5*a**2*d**2 - 26*a*b*c*d + 117*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(9)/4)*(-a*d + b*c)**4) + b*x**(sympy.S(3)/2)/(2*a*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*x**(sympy.S(3)/2)*(a*d + 2*b*c)/(4*a*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d*x**(sympy.S(3)/2)*(-5*a**2*d**2 + 21*a*b*c*d + 8*b**2*c**2)/(16*a*c**2*(c + d*x**2)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(5)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(9)/4)*(-13*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_500():
    f = 1/(sqrt(x)*(a + b*x**2)**2*(c + d*x**2)**3)
    F = -3*sqrt(2)*d**(sympy.S(7)/4)*(7*a**2*d**2 - 30*a*b*c*d + 55*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*d**(sympy.S(7)/4)*(7*a**2*d**2 - 30*a*b*c*d + 55*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(11)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*d**(sympy.S(7)/4)*(7*a**2*d**2 - 30*a*b*c*d + 55*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*d**(sympy.S(7)/4)*(7*a**2*d**2 - 30*a*b*c*d + 55*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(11)/4)*(-a*d + b*c)**4) + b*sqrt(x)/(2*a*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*sqrt(x)*(a*d + 2*b*c)/(4*a*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d*sqrt(x)*(-7*a**2*d**2 + 23*a*b*c*d + 8*b**2*c**2)/(16*a*c**2*(c + d*x**2)*(-a*d + b*c)**3) - 3*sqrt(2)*b**(sympy.S(11)/4)*(-5*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*b**(sympy.S(11)/4)*(-5*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*b**(sympy.S(11)/4)*(-5*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*b**(sympy.S(11)/4)*(-5*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_501():
    f = 1/(x**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)**3)
    F = -sqrt(2)*d**(sympy.S(9)/4)*(45*a**2*d**2 - 170*a*b*c*d + 221*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*(-a*d + b*c)**4) + sqrt(2)*d**(sympy.S(9)/4)*(45*a**2*d**2 - 170*a*b*c*d + 221*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(13)/4)*(-a*d + b*c)**4) + sqrt(2)*d**(sympy.S(9)/4)*(45*a**2*d**2 - 170*a*b*c*d + 221*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(9)/4)*(45*a**2*d**2 - 170*a*b*c*d + 221*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(13)/4)*(-a*d + b*c)**4) + b/(2*a*sqrt(x)*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(a*d + 2*b*c)/(4*a*c*sqrt(x)*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(-9*a**2*d**2 + 25*a*b*c*d + 8*b**2*c**2)/(16*a*c**2*sqrt(x)*(c + d*x**2)*(-a*d + b*c)**3) + (45*a**3*d**3 - 125*a**2*b*c*d**2 + 96*a*b**2*c**2*d - 40*b**3*c**3)/(16*a**2*c**3*sqrt(x)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 5*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 5*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(9)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 5*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(13)/4)*(-17*a*d + 5*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_502():
    f = 1/(x**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)**3)
    F = sqrt(2)*d**(sympy.S(11)/4)*(77*a**2*d**2 - 266*a*b*c*d + 285*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(11)/4)*(77*a**2*d**2 - 266*a*b*c*d + 285*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(15)/4)*(-a*d + b*c)**4) + sqrt(2)*d**(sympy.S(11)/4)*(77*a**2*d**2 - 266*a*b*c*d + 285*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*(-a*d + b*c)**4) - sqrt(2)*d**(sympy.S(11)/4)*(77*a**2*d**2 - 266*a*b*c*d + 285*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(15)/4)*(-a*d + b*c)**4) + b/(2*a*x**(sympy.S(3)/2)*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(a*d + 2*b*c)/(4*a*c*x**(sympy.S(3)/2)*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(-11*a**2*d**2 + 27*a*b*c*d + 8*b**2*c**2)/(16*a*c**2*x**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)**3) + (77*a**3*d**3 - 189*a**2*b*c*d**2 + 96*a*b**2*c**2*d - 56*b**3*c**3)/(48*a**2*c**3*x**(sympy.S(3)/2)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(15)/4)*(-19*a*d + 7*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(15)/4)*(-19*a*d + 7*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(11)/4)*(-a*d + b*c)**4) + sqrt(2)*b**(sympy.S(15)/4)*(-19*a*d + 7*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*(-a*d + b*c)**4) - sqrt(2)*b**(sympy.S(15)/4)*(-19*a*d + 7*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(11)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_503():
    f = 1/(x**(sympy.S(7)/2)*(a + b*x**2)**2*(c + d*x**2)**3)
    F = 3*sqrt(2)*d**(sympy.S(13)/4)*(39*a**2*d**2 - 126*a*b*c*d + 119*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*d**(sympy.S(13)/4)*(39*a**2*d**2 - 126*a*b*c*d + 119*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*sqrt(x) + sqrt(c) + sqrt(d)*x)/(128*c**(sympy.S(17)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*d**(sympy.S(13)/4)*(39*a**2*d**2 - 126*a*b*c*d + 119*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*d**(sympy.S(13)/4)*(39*a**2*d**2 - 126*a*b*c*d + 119*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)*(-a*d + b*c)**4) + b/(2*a*x**(sympy.S(5)/2)*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(a*d + 2*b*c)/(4*a*c*x**(sympy.S(5)/2)*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(-13*a**2*d**2 + 29*a*b*c*d + 8*b**2*c**2)/(16*a*c**2*x**(sympy.S(5)/2)*(c + d*x**2)*(-a*d + b*c)**3) + (117*a**3*d**3 - 261*a**2*b*c*d**2 + 96*a*b**2*c**2*d - 72*b**3*c**3)/(80*a**2*c**3*x**(sympy.S(5)/2)*(-a*d + b*c)**3) + (-117*a**4*d**4 + 261*a**3*b*c*d**3 - 96*a**2*b**2*c**2*d**2 - 96*a*b**3*c**3*d + 72*b**4*c**4)/(16*a**3*c**4*sqrt(x)*(-a*d + b*c)**3) + 3*sqrt(2)*b**(sympy.S(17)/4)*(-7*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*b**(sympy.S(17)/4)*(-7*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x) + sqrt(a) + sqrt(b)*x)/(16*a**(sympy.S(13)/4)*(-a*d + b*c)**4) - 3*sqrt(2)*b**(sympy.S(17)/4)*(-7*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*(-a*d + b*c)**4) + 3*sqrt(2)*b**(sympy.S(17)/4)*(-7*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(13)/4)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_504():
    f = x**5*(A + B*x**2)*sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(9)/2)/(9*b**4) + a**2*(a + b*x**2)**(sympy.S(3)/2)*(A*b - B*a)/(3*b**4) - a*(a + b*x**2)**(sympy.S(5)/2)*(2*A*b - 3*B*a)/(5*b**4) + (a + b*x**2)**(sympy.S(7)/2)*(A*b - 3*B*a)/(7*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_505():
    f = x**4*(A + B*x**2)*sqrt(a + b*x**2)
    F = B*x**5*(a + b*x**2)**(sympy.S(3)/2)/(8*b) + a**3*(8*A*b - 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(7)/2)) - a**2*x*sqrt(a + b*x**2)*(8*A*b - 5*B*a)/(128*b**3) + a*x**3*sqrt(a + b*x**2)*(8*A*b - 5*B*a)/(192*b**2) + x**5*sqrt(a + b*x**2)*(8*A*b - 5*B*a)/(48*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_506():
    f = x**3*(A + B*x**2)*sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(7)/2)/(7*b**3) - a*(a + b*x**2)**(sympy.S(3)/2)*(A*b - B*a)/(3*b**3) + (a + b*x**2)**(sympy.S(5)/2)*(A*b - 2*B*a)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_507():
    f = x**2*(A + B*x**2)*sqrt(a + b*x**2)
    F = B*x**3*(a + b*x**2)**(sympy.S(3)/2)/(6*b) - a**2*(2*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(5)/2)) + a*x*sqrt(a + b*x**2)*(2*A*b - B*a)/(16*b**2) + x**3*sqrt(a + b*x**2)*(2*A*b - B*a)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_508():
    f = x*(A + B*x**2)*sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(5)/2)/(5*b**2) + (a + b*x**2)**(sympy.S(3)/2)*(A*b - B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_509():
    f = (A + B*x**2)*sqrt(a + b*x**2)
    F = B*x*(a + b*x**2)**(sympy.S(3)/2)/(4*b) + a*(4*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(3)/2)) + x*sqrt(a + b*x**2)*(4*A*b - B*a)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_510():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x
    F = -A*sqrt(a)*atanh(sqrt(a + b*x**2)/sqrt(a)) + A*sqrt(a + b*x**2) + B*(a + b*x**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_511():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**2
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(a*x) + (2*A*b + B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)) + x*sqrt(a + b*x**2)*(2*A*b + B*a)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_512():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**3
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(2*a*x**2) + sqrt(a + b*x**2)*(A*b + 2*B*a)/(2*a) - (A*b + 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_513():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**4
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(3*a*x**3) + B*sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - B*sqrt(a + b*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_514():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**5
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(4*a*x**4) + sqrt(a + b*x**2)*(A*b - 4*B*a)/(8*a*x**2) + b*(A*b - 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_515():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**6
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(5*a*x**5) + (a + b*x**2)**(sympy.S(3)/2)*(2*A*b - 5*B*a)/(15*a**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_516():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**7
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(6*a*x**6) + sqrt(a + b*x**2)*(A*b - 2*B*a)/(8*a*x**4) + b*sqrt(a + b*x**2)*(A*b - 2*B*a)/(16*a**2*x**2) - b**2*(A*b - 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_517():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**8
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(7*a*x**7) + (a + b*x**2)**(sympy.S(3)/2)*(4*A*b - 7*B*a)/(35*a**2*x**5) - 2*b*(a + b*x**2)**(sympy.S(3)/2)*(4*A*b - 7*B*a)/(105*a**3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_518():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**9
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(8*a*x**8) + sqrt(a + b*x**2)*(5*A*b - 8*B*a)/(48*a*x**6) + b*sqrt(a + b*x**2)*(5*A*b - 8*B*a)/(192*a**2*x**4) - b**2*sqrt(a + b*x**2)*(5*A*b - 8*B*a)/(128*a**3*x**2) + b**3*(5*A*b - 8*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_519():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**10
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(9*a*x**9) + (a + b*x**2)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(21*a**2*x**7) - 4*b*(a + b*x**2)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(105*a**3*x**5) + 8*b**2*(a + b*x**2)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(315*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_520():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**11
    F = -A*(a + b*x**2)**(sympy.S(3)/2)/(10*a*x**10) + sqrt(a + b*x**2)*(7*A*b - 10*B*a)/(80*a*x**8) + b*sqrt(a + b*x**2)*(7*A*b - 10*B*a)/(480*a**2*x**6) - b**2*sqrt(a + b*x**2)*(7*A*b - 10*B*a)/(384*a**3*x**4) + b**3*sqrt(a + b*x**2)*(7*A*b - 10*B*a)/(256*a**4*x**2) - b**4*(7*A*b - 10*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(256*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_521():
    f = x**5*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*(a + b*x**2)**(sympy.S(11)/2)/(11*b**4) + a**2*(a + b*x**2)**(sympy.S(5)/2)*(A*b - B*a)/(5*b**4) - a*(a + b*x**2)**(sympy.S(7)/2)*(2*A*b - 3*B*a)/(7*b**4) + (a + b*x**2)**(sympy.S(9)/2)*(A*b - 3*B*a)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_522():
    f = x**4*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*x**5*(a + b*x**2)**(sympy.S(5)/2)/(10*b) + 3*a**4*(2*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(7)/2)) - 3*a**3*x*sqrt(a + b*x**2)*(2*A*b - B*a)/(256*b**3) + a**2*x**3*sqrt(a + b*x**2)*(2*A*b - B*a)/(128*b**2) + a*x**5*sqrt(a + b*x**2)*(2*A*b - B*a)/(32*b) + x**5*(a + b*x**2)**(sympy.S(3)/2)*(2*A*b - B*a)/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_523():
    f = x**3*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*(a + b*x**2)**(sympy.S(9)/2)/(9*b**3) - a*(a + b*x**2)**(sympy.S(5)/2)*(A*b - B*a)/(5*b**3) + (a + b*x**2)**(sympy.S(7)/2)*(A*b - 2*B*a)/(7*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_524():
    f = x**2*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*x**3*(a + b*x**2)**(sympy.S(5)/2)/(8*b) - a**3*(8*A*b - 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(5)/2)) + a**2*x*sqrt(a + b*x**2)*(8*A*b - 3*B*a)/(128*b**2) + a*x**3*sqrt(a + b*x**2)*(8*A*b - 3*B*a)/(64*b) + x**3*(a + b*x**2)**(sympy.S(3)/2)*(8*A*b - 3*B*a)/(48*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_525():
    f = x*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*(a + b*x**2)**(sympy.S(7)/2)/(7*b**2) + (a + b*x**2)**(sympy.S(5)/2)*(A*b - B*a)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_526():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = B*x*(a + b*x**2)**(sympy.S(5)/2)/(6*b) + a**2*(6*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(3)/2)) + a*x*sqrt(a + b*x**2)*(6*A*b - B*a)/(16*b) + x*(a + b*x**2)**(sympy.S(3)/2)*(6*A*b - B*a)/(24*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_527():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x
    F = -A*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + A*a*sqrt(a + b*x**2) + A*(a + b*x**2)**(sympy.S(3)/2)/3 + B*(a + b*x**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_528():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**2
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(a*x) + 3*a*(4*A*b + B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)) + x*sqrt(a + b*x**2)*(3*A*b/2 + 3*B*a/8) + x*(a + b*x**2)**(sympy.S(3)/2)*(4*A*b + B*a)/(4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_529():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**3
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(2*a*x**2) - sqrt(a)*(3*A*b + 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + sqrt(a + b*x**2)*(3*A*b/2 + B*a) + (a + b*x**2)**(sympy.S(3)/2)*(3*A*b + 2*B*a)/(6*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_530():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**4
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(3*a*x**3) + sqrt(b)*(2*A*b + 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + b*x*sqrt(a + b*x**2)*(2*A*b + 3*B*a)/(2*a) - (a + b*x**2)**(sympy.S(3)/2)*(2*A*b + 3*B*a)/(3*a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_531():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**5
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(4*a*x**4) + 3*b*sqrt(a + b*x**2)*(A*b + 4*B*a)/(8*a) - (a + b*x**2)**(sympy.S(3)/2)*(A*b + 4*B*a)/(8*a*x**2) - 3*b*(A*b + 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_532():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**6
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(5*a*x**5) + B*b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - B*b*sqrt(a + b*x**2)/x - B*(a + b*x**2)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_533():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**7
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(6*a*x**6) + b*sqrt(a + b*x**2)*(A*b - 6*B*a)/(16*a*x**2) + (a + b*x**2)**(sympy.S(3)/2)*(A*b - 6*B*a)/(24*a*x**4) + b**2*(A*b - 6*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_534():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**8
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(7*a*x**7) + (a + b*x**2)**(sympy.S(5)/2)*(2*A*b - 7*B*a)/(35*a**2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_535():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**9
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(8*a*x**8) + b*sqrt(a + b*x**2)*(3*A*b - 8*B*a)/(64*a*x**4) + (a + b*x**2)**(sympy.S(3)/2)*(3*A*b - 8*B*a)/(48*a*x**6) + b**2*sqrt(a + b*x**2)*(3*A*b - 8*B*a)/(128*a**2*x**2) - b**3*(3*A*b - 8*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_536():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**10
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(9*a*x**9) + (a + b*x**2)**(sympy.S(5)/2)*(4*A*b - 9*B*a)/(63*a**2*x**7) - 2*b*(a + b*x**2)**(sympy.S(5)/2)*(4*A*b - 9*B*a)/(315*a**3*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_537():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/x**11
    F = -A*(a + b*x**2)**(sympy.S(5)/2)/(10*a*x**10) + b*sqrt(a + b*x**2)*(A*b - 2*B*a)/(32*a*x**6) + (a + b*x**2)**(sympy.S(3)/2)*(A*b - 2*B*a)/(16*a*x**8) + b**2*sqrt(a + b*x**2)*(A*b - 2*B*a)/(128*a**2*x**4) - 3*b**3*sqrt(a + b*x**2)*(A*b - 2*B*a)/(256*a**3*x**2) + 3*b**4*(A*b - 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(256*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_538():
    f = x**5*(A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*(a + b*x**2)**(sympy.S(13)/2)/(13*b**4) + a**2*(a + b*x**2)**(sympy.S(7)/2)*(A*b - B*a)/(7*b**4) - a*(a + b*x**2)**(sympy.S(9)/2)*(2*A*b - 3*B*a)/(9*b**4) + (a + b*x**2)**(sympy.S(11)/2)*(A*b - 3*B*a)/(11*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_539():
    f = x**4*(A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*x**5*(a + b*x**2)**(sympy.S(7)/2)/(12*b) + a**5*(12*A*b - 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(1024*b**(sympy.S(7)/2)) - a**4*x*sqrt(a + b*x**2)*(12*A*b - 5*B*a)/(1024*b**3) + a**3*x**3*sqrt(a + b*x**2)*(12*A*b - 5*B*a)/(1536*b**2) + a**2*x**5*sqrt(a + b*x**2)*(12*A*b - 5*B*a)/(384*b) + a*x**5*(a + b*x**2)**(sympy.S(3)/2)*(12*A*b - 5*B*a)/(192*b) + x**5*(a + b*x**2)**(sympy.S(5)/2)*(12*A*b - 5*B*a)/(120*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_540():
    f = x**3*(A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*(a + b*x**2)**(sympy.S(11)/2)/(11*b**3) - a*(a + b*x**2)**(sympy.S(7)/2)*(A*b - B*a)/(7*b**3) + (a + b*x**2)**(sympy.S(9)/2)*(A*b - 2*B*a)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_541():
    f = x**2*(A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*x**3*(a + b*x**2)**(sympy.S(7)/2)/(10*b) - a**4*(10*A*b - 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(5)/2)) + a**3*x*sqrt(a + b*x**2)*(10*A*b - 3*B*a)/(256*b**2) + a**2*x**3*sqrt(a + b*x**2)*(10*A*b - 3*B*a)/(128*b) + a*x**3*(a + b*x**2)**(sympy.S(3)/2)*(10*A*b - 3*B*a)/(96*b) + x**3*(a + b*x**2)**(sympy.S(5)/2)*(10*A*b - 3*B*a)/(80*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_542():
    f = x*(A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*(a + b*x**2)**(sympy.S(9)/2)/(9*b**2) + (a + b*x**2)**(sympy.S(7)/2)*(A*b - B*a)/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_543():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)
    F = B*x*(a + b*x**2)**(sympy.S(7)/2)/(8*b) + 5*a**3*(8*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(3)/2)) + 5*a**2*x*sqrt(a + b*x**2)*(8*A*b - B*a)/(128*b) + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)*(8*A*b - B*a)/(192*b) + x*(a + b*x**2)**(sympy.S(5)/2)*(8*A*b - B*a)/(48*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_544():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x
    F = -A*a**(sympy.S(5)/2)*atanh(sqrt(a + b*x**2)/sqrt(a)) + A*a**2*sqrt(a + b*x**2) + A*a*(a + b*x**2)**(sympy.S(3)/2)/3 + A*(a + b*x**2)**(sympy.S(5)/2)/5 + B*(a + b*x**2)**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_545():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**2
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(a*x) + 5*a**2*(6*A*b + B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*sqrt(b)) + 5*a*x*sqrt(a + b*x**2)*(6*A*b + B*a)/16 + x*(a + b*x**2)**(sympy.S(3)/2)*(30*A*b + 5*B*a)/24 + x*(a + b*x**2)**(sympy.S(5)/2)*(6*A*b + B*a)/(6*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_546():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**3
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(2*a*x**2) - a**(sympy.S(3)/2)*(5*A*b + 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/2 + a*sqrt(a + b*x**2)*(5*A*b + 2*B*a)/2 + (a + b*x**2)**(sympy.S(3)/2)*(5*A*b + 2*B*a)/6 + (a + b*x**2)**(sympy.S(5)/2)*(5*A*b + 2*B*a)/(10*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_547():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**4
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(3*a*x**3) + 5*a*sqrt(b)*(4*A*b + 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/8 + 5*b*x*sqrt(a + b*x**2)*(4*A*b + 3*B*a)/8 + 5*b*x*(a + b*x**2)**(sympy.S(3)/2)*(4*A*b + 3*B*a)/(12*a) - (a + b*x**2)**(sympy.S(5)/2)*(4*A*b + 3*B*a)/(3*a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_548():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**5
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(4*a*x**4) - 5*sqrt(a)*b*(3*A*b + 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/8 + 5*b*sqrt(a + b*x**2)*(3*A*b + 4*B*a)/8 + 5*b*(a + b*x**2)**(sympy.S(3)/2)*(3*A*b + 4*B*a)/(24*a) - (a + b*x**2)**(sympy.S(5)/2)*(3*A*b + 4*B*a)/(8*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_549():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**6
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(5*a*x**5) + b**(sympy.S(3)/2)*(2*A*b + 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/2 + b**2*x*sqrt(a + b*x**2)*(2*A*b + 5*B*a)/(2*a) - b*(a + b*x**2)**(sympy.S(3)/2)*(2*A*b + 5*B*a)/(3*a*x) - (a + b*x**2)**(sympy.S(5)/2)*(2*A*b + 5*B*a)/(15*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_550():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**7
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(6*a*x**6) + 5*b**2*sqrt(a + b*x**2)*(A*b + 6*B*a)/(16*a) - 5*b*(a + b*x**2)**(sympy.S(3)/2)*(A*b + 6*B*a)/(48*a*x**2) - (a + b*x**2)**(sympy.S(5)/2)*(A*b + 6*B*a)/(24*a*x**4) - 5*b**2*(A*b + 6*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_551():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**8
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(7*a*x**7) + B*b**(sympy.S(5)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2)) - B*b**2*sqrt(a + b*x**2)/x - B*b*(a + b*x**2)**(sympy.S(3)/2)/(3*x**3) - B*(a + b*x**2)**(sympy.S(5)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_552():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**9
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(8*a*x**8) + 5*b**2*sqrt(a + b*x**2)*(A*b - 8*B*a)/(128*a*x**2) + 5*b*(a + b*x**2)**(sympy.S(3)/2)*(A*b - 8*B*a)/(192*a*x**4) + (a + b*x**2)**(sympy.S(5)/2)*(A*b - 8*B*a)/(48*a*x**6) + 5*b**3*(A*b - 8*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(128*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_553():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**10
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(9*a*x**9) + (a + b*x**2)**(sympy.S(7)/2)*(2*A*b - 9*B*a)/(63*a**2*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_554():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(5)/2)/x**11
    F = -A*(a + b*x**2)**(sympy.S(7)/2)/(10*a*x**10) + b**2*sqrt(a + b*x**2)*(3*A*b - 10*B*a)/(128*a*x**4) + b*(a + b*x**2)**(sympy.S(3)/2)*(3*A*b - 10*B*a)/(96*a*x**6) + (a + b*x**2)**(sympy.S(5)/2)*(3*A*b - 10*B*a)/(80*a*x**8) + b**3*sqrt(a + b*x**2)*(3*A*b - 10*B*a)/(256*a**2*x**2) - b**4*(3*A*b - 10*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(256*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_555():
    f = x**5*(A + B*x**2)/sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(7)/2)/(7*b**4) + a**2*sqrt(a + b*x**2)*(A*b - B*a)/b**4 - a*(a + b*x**2)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(3*b**4) + (a + b*x**2)**(sympy.S(5)/2)*(A*b - 3*B*a)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_556():
    f = x**4*(A + B*x**2)/sqrt(a + b*x**2)
    F = B*x**5*sqrt(a + b*x**2)/(6*b) + a**2*(6*A*b - 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(7)/2)) - a*x*sqrt(a + b*x**2)*(6*A*b - 5*B*a)/(16*b**3) + x**3*sqrt(a + b*x**2)*(6*A*b - 5*B*a)/(24*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_557():
    f = x**3*(A + B*x**2)/sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(5)/2)/(5*b**3) - a*sqrt(a + b*x**2)*(A*b - B*a)/b**3 + (a + b*x**2)**(sympy.S(3)/2)*(A*b - 2*B*a)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_558():
    f = x**2*(A + B*x**2)/sqrt(a + b*x**2)
    F = B*x**3*sqrt(a + b*x**2)/(4*b) - a*(4*A*b - 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2)) + x*sqrt(a + b*x**2)*(4*A*b - 3*B*a)/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_559():
    f = x*(A + B*x**2)/sqrt(a + b*x**2)
    F = B*(a + b*x**2)**(sympy.S(3)/2)/(3*b**2) + sqrt(a + b*x**2)*(A*b - B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_560():
    f = (A + B*x**2)/sqrt(a + b*x**2)
    F = B*x*sqrt(a + b*x**2)/(2*b) + (2*A*b - B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_561():
    f = (A + B*x**2)/(x*sqrt(a + b*x**2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/sqrt(a) + B*sqrt(a + b*x**2)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_562():
    f = (A + B*x**2)/(x**2*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(a*x) + B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_563():
    f = (A + B*x**2)/(x**3*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(2*a*x**2) + (A*b - 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_564():
    f = (A + B*x**2)/(x**4*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(3*a*x**3) + sqrt(a + b*x**2)*(2*A*b - 3*B*a)/(3*a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_565():
    f = (A + B*x**2)/(x**5*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(4*a*x**4) + sqrt(a + b*x**2)*(3*A*b - 4*B*a)/(8*a**2*x**2) - b*(3*A*b - 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_566():
    f = (A + B*x**2)/(x**6*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(5*a*x**5) + sqrt(a + b*x**2)*(4*A*b - 5*B*a)/(15*a**2*x**3) - 2*b*sqrt(a + b*x**2)*(4*A*b - 5*B*a)/(15*a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_567():
    f = (A + B*x**2)/(x**7*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(6*a*x**6) + sqrt(a + b*x**2)*(5*A*b - 6*B*a)/(24*a**2*x**4) - b*sqrt(a + b*x**2)*(5*A*b - 6*B*a)/(16*a**3*x**2) + b**2*(5*A*b - 6*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_568():
    f = (A + B*x**2)/(x**8*sqrt(a + b*x**2))
    F = -A*sqrt(a + b*x**2)/(7*a*x**7) + sqrt(a + b*x**2)*(6*A*b - 7*B*a)/(35*a**2*x**5) - 4*b*sqrt(a + b*x**2)*(6*A*b - 7*B*a)/(105*a**3*x**3) + 8*b**2*sqrt(a + b*x**2)*(6*A*b - 7*B*a)/(105*a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_569():
    f = x**6*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*x**7/(6*b*sqrt(a + b*x**2)) + 5*a**2*(6*A*b - 7*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(9)/2)) - 5*a*x*sqrt(a + b*x**2)*(6*A*b - 7*B*a)/(16*b**4) - x**5*(6*A*b - 7*B*a)/(6*b**2*sqrt(a + b*x**2)) + x**3*sqrt(a + b*x**2)*(30*A*b - 35*B*a)/(24*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_570():
    f = x**5*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*(a + b*x**2)**(sympy.S(5)/2)/(5*b**4) - a**2*(A*b - B*a)/(b**4*sqrt(a + b*x**2)) - a*sqrt(a + b*x**2)*(2*A*b - 3*B*a)/b**4 + (a + b*x**2)**(sympy.S(3)/2)*(A*b - 3*B*a)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_571():
    f = x**4*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*x**5/(4*b*sqrt(a + b*x**2)) - 3*a*(4*A*b - 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(7)/2)) - x**3*(4*A*b - 5*B*a)/(4*b**2*sqrt(a + b*x**2)) + x*sqrt(a + b*x**2)*(12*A*b - 15*B*a)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_572():
    f = x**3*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*(a + b*x**2)**(sympy.S(3)/2)/(3*b**3) + a*(A*b - B*a)/(b**3*sqrt(a + b*x**2)) + sqrt(a + b*x**2)*(A*b - 2*B*a)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_573():
    f = x**2*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*x*sqrt(a + b*x**2)/(2*b**2) - x*(A*b - B*a)/(b**2*sqrt(a + b*x**2)) + (2*A*b - 3*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_574():
    f = x*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*sqrt(a + b*x**2)/b**2 - (A*b - B*a)/(b**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_575():
    f = (A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2) + x*(A*b - B*a)/(a*b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_576():
    f = (A + B*x**2)/(x*(a + b*x**2)**(sympy.S(3)/2))
    F = -A*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(3)/2) + (A*b - B*a)/(a*b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_577():
    f = (A + B*x**2)/(x**2*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(a*x*sqrt(a + b*x**2)) - x*(2*A*b - B*a)/(a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_578():
    f = (A + B*x**2)/(x**3*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(2*a*x**2*sqrt(a + b*x**2)) - (3*A*b - 2*B*a)/(2*a**2*sqrt(a + b*x**2)) + (3*A*b - 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_579():
    f = (A + B*x**2)/(x**4*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(3*a*x**3*sqrt(a + b*x**2)) + (4*A*b - 3*B*a)/(3*a**2*x*sqrt(a + b*x**2)) + 2*b*x*(4*A*b - 3*B*a)/(3*a**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_580():
    f = (A + B*x**2)/(x**5*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(4*a*x**4*sqrt(a + b*x**2)) + (5*A*b - 4*B*a)/(8*a**2*x**2*sqrt(a + b*x**2)) + 3*b*(5*A*b - 4*B*a)/(8*a**3*sqrt(a + b*x**2)) - 3*b*(5*A*b - 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_581():
    f = (A + B*x**2)/(x**6*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(5*a*x**5*sqrt(a + b*x**2)) + (6*A*b - 5*B*a)/(15*a**2*x**3*sqrt(a + b*x**2)) - 4*b*(6*A*b - 5*B*a)/(15*a**3*x*sqrt(a + b*x**2)) - 8*b**2*x*(6*A*b - 5*B*a)/(15*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_582():
    f = (A + B*x**2)/(x**7*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(6*a*x**6*sqrt(a + b*x**2)) + (7*A*b - 6*B*a)/(24*a**2*x**4*sqrt(a + b*x**2)) - 5*b*(7*A*b - 6*B*a)/(48*a**3*x**2*sqrt(a + b*x**2)) - 5*b**2*(7*A*b - 6*B*a)/(16*a**4*sqrt(a + b*x**2)) + 5*b**2*(7*A*b - 6*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(16*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_583():
    f = (A + B*x**2)/(x**8*(a + b*x**2)**(sympy.S(3)/2))
    F = -A/(7*a*x**7*sqrt(a + b*x**2)) + (8*A*b - 7*B*a)/(35*a**2*x**5*sqrt(a + b*x**2)) - 2*b*(8*A*b - 7*B*a)/(35*a**3*x**3*sqrt(a + b*x**2)) + 8*b**2*(8*A*b - 7*B*a)/(35*a**4*x*sqrt(a + b*x**2)) + 16*b**3*x*(8*A*b - 7*B*a)/(35*a**5*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_584():
    f = x**7*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*(a + b*x**2)**(sympy.S(5)/2)/(5*b**5) + a**3*(A*b - B*a)/(3*b**5*(a + b*x**2)**(sympy.S(3)/2)) - a**2*(3*A*b - 4*B*a)/(b**5*sqrt(a + b*x**2)) - 3*a*sqrt(a + b*x**2)*(A*b - 2*B*a)/b**5 + (a + b*x**2)**(sympy.S(3)/2)*(A*b - 4*B*a)/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_585():
    f = x**6*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*x**7/(4*b*(a + b*x**2)**(sympy.S(3)/2)) - 5*a*(4*A*b - 7*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(9)/2)) - x**5*(4*A*b - 7*B*a)/(12*b**2*(a + b*x**2)**(sympy.S(3)/2)) - x**3*(20*A*b - 35*B*a)/(12*b**3*sqrt(a + b*x**2)) + x*sqrt(a + b*x**2)*(20*A*b - 35*B*a)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_586():
    f = x**5*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*(a + b*x**2)**(sympy.S(3)/2)/(3*b**4) - a**2*(A*b - B*a)/(3*b**4*(a + b*x**2)**(sympy.S(3)/2)) + a*(2*A*b - 3*B*a)/(b**4*sqrt(a + b*x**2)) + sqrt(a + b*x**2)*(A*b - 3*B*a)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_587():
    f = x**4*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*x*sqrt(a + b*x**2)/(2*b**3) + a*x*(A*b - B*a)/(3*b**3*(a + b*x**2)**(sympy.S(3)/2)) - x*(4*A*b - 7*B*a)/(3*b**3*sqrt(a + b*x**2)) + (2*A*b - 5*B*a)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_588():
    f = x**3*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = B*sqrt(a + b*x**2)/b**3 + a*(A*b - B*a)/(3*b**3*(a + b*x**2)**(sympy.S(3)/2)) - (A*b - 2*B*a)/(b**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_589():
    f = x**2*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = -B*x/(b**2*sqrt(a + b*x**2)) + B*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(5)/2) + x**3*(A*b - B*a)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_590():
    f = x*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = -B/(b**2*sqrt(a + b*x**2)) + (-A*b + B*a)/(3*b**2*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_591():
    f = (A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = 2*A*x/(3*a**2*sqrt(a + b*x**2)) + x*(A + B*x**2)/(3*a*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_592():
    f = (A + B*x**2)/(x*(a + b*x**2)**(sympy.S(5)/2))
    F = A/(a**2*sqrt(a + b*x**2)) - A*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(5)/2) + (A*b - B*a)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_593():
    f = (A + B*x**2)/(x**2*(a + b*x**2)**(sympy.S(5)/2))
    F = -A/(a*x*(a + b*x**2)**(sympy.S(3)/2)) - x*(4*A*b - B*a)/(3*a**2*(a + b*x**2)**(sympy.S(3)/2)) - x*(8*A*b - 2*B*a)/(3*a**3*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_594():
    f = (A + B*x**2)/(x**3*(a + b*x**2)**(sympy.S(5)/2))
    F = -A/(2*a*x**2*(a + b*x**2)**(sympy.S(3)/2)) + (-5*A*b + 2*B*a)/(6*a**2*(a + b*x**2)**(sympy.S(3)/2)) - (5*A*b - 2*B*a)/(2*a**3*sqrt(a + b*x**2)) + (5*A*b - 2*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_595():
    f = (A + B*x**2)/(x**4*(a + b*x**2)**(sympy.S(5)/2))
    F = -A/(3*a*x**3*(a + b*x**2)**(sympy.S(3)/2)) + (2*A*b - B*a)/(a**2*x*(a + b*x**2)**(sympy.S(3)/2)) + 4*b*x*(2*A*b - B*a)/(3*a**3*(a + b*x**2)**(sympy.S(3)/2)) + 8*b*x*(2*A*b - B*a)/(3*a**4*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_596():
    f = (A + B*x**2)/(x**5*(a + b*x**2)**(sympy.S(5)/2))
    F = -A/(4*a*x**4*(a + b*x**2)**(sympy.S(3)/2)) + (7*A*b - 4*B*a)/(8*a**2*x**2*(a + b*x**2)**(sympy.S(3)/2)) + 5*b*(7*A*b - 4*B*a)/(24*a**3*(a + b*x**2)**(sympy.S(3)/2)) + 5*b*(7*A*b - 4*B*a)/(8*a**4*sqrt(a + b*x**2)) - 5*b*(7*A*b - 4*B*a)*atanh(sqrt(a + b*x**2)/sqrt(a))/(8*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_597():
    f = (A + B*x**2)/(x**6*(a + b*x**2)**(sympy.S(5)/2))
    F = -A/(5*a*x**5*(a + b*x**2)**(sympy.S(3)/2)) + (8*A*b - 5*B*a)/(15*a**2*x**3*(a + b*x**2)**(sympy.S(3)/2)) - 2*b*(8*A*b - 5*B*a)/(5*a**3*x*(a + b*x**2)**(sympy.S(3)/2)) - 8*b**2*x*(8*A*b - 5*B*a)/(15*a**4*(a + b*x**2)**(sympy.S(3)/2)) - 16*b**2*x*(8*A*b - 5*B*a)/(15*a**5*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_598():
    f = x**5*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = b**2*(c + d*x**2)**(sympy.S(11)/2)/(11*d**5) - 2*b*(c + d*x**2)**(sympy.S(9)/2)*(-a*d + 2*b*c)/(9*d**5) + c**2*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**5) - 2*c*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*d + 2*b*c)/(5*d**5) + (c + d*x**2)**(sympy.S(7)/2)*(a**2*d**2 - 6*a*b*c*d + 6*b**2*c**2)/(7*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_599():
    f = x**3*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = b**2*(c + d*x**2)**(sympy.S(9)/2)/(9*d**4) - b*(c + d*x**2)**(sympy.S(7)/2)*(-2*a*d + 3*b*c)/(7*d**4) - c*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**4) + (c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*d + 3*b*c)/(5*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_600():
    f = x*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = b**2*(c + d*x**2)**(sympy.S(7)/2)/(7*d**3) - 2*b*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**3) + (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_601():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x
    F = -a**2*sqrt(c)*atanh(sqrt(c + d*x**2)/sqrt(c)) + a**2*sqrt(c + d*x**2) + b**2*(c + d*x**2)**(sympy.S(5)/2)/(5*d**2) - b*(c + d*x**2)**(sympy.S(3)/2)*(-2*a*d + b*c)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_602():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**3
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(2*c*x**2) + a*sqrt(c + d*x**2)*(a*d + 4*b*c)/(2*c) - a*(a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*sqrt(c)) + b**2*(c + d*x**2)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_603():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**5
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(4*c*x**4) - a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 8*b*c)/(8*c**2*x**2) + sqrt(c + d*x**2)*(a*d*(-a*d + 8*b*c) + 8*b**2*c**2)/(8*c**2) - (a*d*(-a*d + 8*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(8*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_604():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**7
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(6*c*x**6) - a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 4*b*c)/(8*c**2*x**4) - sqrt(c + d*x**2)*(a**2*d**2 - 4*a*b*c*d + 8*b**2*c**2)/(16*c**2*x**2) - d*(a**2*d**2 - 4*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(16*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_605():
    f = x**2*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = b**2*x**5*(c + d*x**2)**(sympy.S(3)/2)/(8*d) - b*x**3*(c + d*x**2)**(sympy.S(3)/2)*(-16*a*d + 5*b*c)/(48*d**2) - c**2*(16*a**2*d**2 + b*c*(-16*a*d + 5*b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(128*d**(sympy.S(7)/2)) + c*x*sqrt(c + d*x**2)*(16*a**2*d**2 + b*c*(-16*a*d + 5*b*c))/(128*d**3) + x**3*sqrt(c + d*x**2)*(16*a**2*d**2 + b*c*(-16*a*d + 5*b*c))/(64*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_606():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)
    F = b*x*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)/(6*d) - b*x*(c + d*x**2)**(sympy.S(3)/2)*(-8*a*d + 3*b*c)/(24*d**2) + c*(8*a**2*d**2 - 4*a*b*c*d + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*d**(sympy.S(5)/2)) + x*sqrt(c + d*x**2)*(8*a**2*d**2 - 4*a*b*c*d + b**2*c**2)/(16*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_607():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**2
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(c*x) + b**2*x*(c + d*x**2)**(sympy.S(3)/2)/(4*d) - (-8*a*d*(a*d + b*c) + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(3)/2)) - x*sqrt(c + d*x**2)*(-8*a*d*(a*d + b*c) + b**2*c**2)/(8*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_608():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**4
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(3*c*x**3) - 2*a*b*(c + d*x**2)**(sympy.S(3)/2)/(c*x) + b*(4*a*d + b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*sqrt(d)) + b*x*sqrt(c + d*x**2)*(4*a*d + b*c)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_609():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**6
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(5*c*x**5) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 5*b*c)/(15*c**2*x**3) + b**2*sqrt(d)*atanh(sqrt(d)*x/sqrt(c + d*x**2)) - b**2*sqrt(c + d*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_610():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**8
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(7*c*x**7) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-2*a*d + 7*b*c)/(35*c**2*x**5) - (c + d*x**2)**(sympy.S(3)/2)*(-4*a*d*(-2*a*d + 7*b*c) + 35*b**2*c**2)/(105*c**3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_611():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**10
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(9*c*x**9) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 3*b*c)/(21*c**2*x**7) - (c + d*x**2)**(sympy.S(3)/2)*(-8*a*d*(-a*d + 3*b*c) + 21*b**2*c**2)/(105*c**3*x**5) + 2*d*(c + d*x**2)**(sympy.S(3)/2)*(-8*a*d*(-a*d + 3*b*c) + 21*b**2*c**2)/(315*c**4*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_612():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**12
    F = -a**2*(c + d*x**2)**(sympy.S(3)/2)/(11*c*x**11) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-4*a*d + 11*b*c)/(99*c**2*x**9) - (c + d*x**2)**(sympy.S(3)/2)*(-4*a*d*(-4*a*d + 11*b*c) + 33*b**2*c**2)/(231*c**3*x**7) + 4*d*(c + d*x**2)**(sympy.S(3)/2)*(-4*a*d*(-4*a*d + 11*b*c) + 33*b**2*c**2)/(1155*c**4*x**5) - 8*d**2*(c + d*x**2)**(sympy.S(3)/2)*(-4*a*d*(-4*a*d + 11*b*c) + 33*b**2*c**2)/(3465*c**5*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_613():
    f = x**4*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*x**7*(c + d*x**2)**(sympy.S(5)/2)/(12*d) - b*x**5*(c + d*x**2)**(sympy.S(5)/2)*(-24*a*d + 7*b*c)/(120*d**2) + c**4*(24*a**2*d**2 + b*c*(-24*a*d + 7*b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(1024*d**(sympy.S(9)/2)) - c**3*x*sqrt(c + d*x**2)*(24*a**2*d**2 + b*c*(-24*a*d + 7*b*c))/(1024*d**4) + c**2*x**3*sqrt(c + d*x**2)*(24*a**2*d**2 + b*c*(-24*a*d + 7*b*c))/(1536*d**3) + c*x**5*sqrt(c + d*x**2)*(24*a**2*d**2 + b*c*(-24*a*d + 7*b*c))/(384*d**2) + x**5*(c + d*x**2)**(sympy.S(3)/2)*(24*a**2*d**2 + b*c*(-24*a*d + 7*b*c))/(192*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_614():
    f = x**3*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*(c + d*x**2)**(sympy.S(11)/2)/(11*d**4) - b*(c + d*x**2)**(sympy.S(9)/2)*(-2*a*d + 3*b*c)/(9*d**4) - c*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**4) + (c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c)*(-a*d + 3*b*c)/(7*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_615():
    f = x**2*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*x**5*(c + d*x**2)**(sympy.S(5)/2)/(10*d) - b*x**3*(c + d*x**2)**(sympy.S(5)/2)*(-4*a*d + b*c)/(16*d**2) - c**3*(16*a**2*d**2 + 3*b*c*(-4*a*d + b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(256*d**(sympy.S(7)/2)) + c**2*x*sqrt(c + d*x**2)*(16*a**2*d**2 + 3*b*c*(-4*a*d + b*c))/(256*d**3) + c*x**3*sqrt(c + d*x**2)*(16*a**2*d**2 + 3*b*c*(-4*a*d + b*c))/(128*d**2) + x**3*(c + d*x**2)**(sympy.S(3)/2)*(16*a**2*d**2 + 3*b*c*(-4*a*d + b*c))/(96*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_616():
    f = x*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*(c + d*x**2)**(sympy.S(9)/2)/(9*d**3) - 2*b*(c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c)/(7*d**3) + (c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)**2/(5*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_617():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = b*x*(a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2)/(8*d) - b*x*(c + d*x**2)**(sympy.S(5)/2)*(-10*a*d + 3*b*c)/(48*d**2) + c**2*(48*a**2*d**2 - 16*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(128*d**(sympy.S(5)/2)) + c*x*sqrt(c + d*x**2)*(48*a**2*d**2 - 16*a*b*c*d + 3*b**2*c**2)/(128*d**2) + x*(c + d*x**2)**(sympy.S(3)/2)*(48*a**2*d**2 - 16*a*b*c*d + 3*b**2*c**2)/(192*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_618():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x
    F = -a**2*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**2)/sqrt(c)) + a**2*c*sqrt(c + d*x**2) + a**2*(c + d*x**2)**(sympy.S(3)/2)/3 + b**2*(c + d*x**2)**(sympy.S(7)/2)/(7*d**2) - b*(c + d*x**2)**(sympy.S(5)/2)*(-2*a*d + b*c)/(5*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_619():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**2
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(c*x) + b**2*x*(c + d*x**2)**(sympy.S(5)/2)/(6*d) - c*(-12*a*d*(2*a*d + b*c) + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*d**(sympy.S(3)/2)) - x*sqrt(c + d*x**2)*(-12*a*d*(2*a*d + b*c) + b**2*c**2)/(16*d) - x*(c + d*x**2)**(sympy.S(3)/2)*(-12*a*d*(2*a*d + b*c) + b**2*c**2)/(24*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_620():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**3
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(2*c*x**2) - a*sqrt(c)*(3*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/2 + a*sqrt(c + d*x**2)*(3*a*d + 4*b*c)/2 + a*(c + d*x**2)**(sympy.S(3)/2)*(3*a*d + 4*b*c)/(6*c) + b**2*(c + d*x**2)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_621():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**4
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(3*c*x**3) - 2*a*(c + d*x**2)**(sympy.S(5)/2)*(a*d + 3*b*c)/(3*c**2*x) + (8*a*d*(a*d + 3*b*c) + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*sqrt(d)) + x*sqrt(c + d*x**2)*(8*a*d*(a*d + 3*b*c) + 3*b**2*c**2)/(8*c) + x*(c + d*x**2)**(sympy.S(3)/2)*(8*a*d*(a*d + 3*b*c) + 3*b**2*c**2)/(12*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_622():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**5
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(4*c*x**4) - a*(c + d*x**2)**(sympy.S(5)/2)*(a*d + 8*b*c)/(8*c**2*x**2) + sqrt(c + d*x**2)*(3*a*d*(a*d + 8*b*c) + 8*b**2*c**2)/(8*c) + (c + d*x**2)**(sympy.S(3)/2)*(3*a*d*(a*d + 8*b*c) + 8*b**2*c**2)/(24*c**2) - (3*a*d*(a*d + 8*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(8*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_623():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**6
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(5*c*x**5) - 2*a*b*(c + d*x**2)**(sympy.S(5)/2)/(3*c*x**3) + b*sqrt(d)*(4*a*d + 3*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/2 + b*d*x*sqrt(c + d*x**2)*(4*a*d + 3*b*c)/(2*c) - b*(c + d*x**2)**(sympy.S(3)/2)*(4*a*d + 3*b*c)/(3*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_624():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/x**7
    F = -a**2*(c + d*x**2)**(sympy.S(5)/2)/(6*c*x**6) - a*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + 12*b*c)/(24*c**2*x**4) + d*sqrt(c + d*x**2)*(a*d*(-a*d + 12*b*c) + 24*b**2*c**2)/(16*c**2) - (c + d*x**2)**(sympy.S(3)/2)*(a*d*(-a*d + 12*b*c) + 24*b**2*c**2)/(48*c**2*x**2) - d*(a*d*(-a*d + 12*b*c) + 24*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(16*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_625():
    f = x**3*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*(c + d*x**2)**(sympy.S(13)/2)/(13*d**4) - b*(c + d*x**2)**(sympy.S(11)/2)*(-2*a*d + 3*b*c)/(11*d**4) - c*(c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**4) + (c + d*x**2)**(sympy.S(9)/2)*(-a*d + b*c)*(-a*d + 3*b*c)/(9*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_626():
    f = x**2*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*x**5*(c + d*x**2)**(sympy.S(7)/2)/(12*d) - b*x**3*(c + d*x**2)**(sympy.S(7)/2)*(-24*a*d + 5*b*c)/(120*d**2) - c**4*(40*a**2*d**2 + b*c*(-24*a*d + 5*b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(1024*d**(sympy.S(7)/2)) + c**3*x*sqrt(c + d*x**2)*(40*a**2*d**2 + b*c*(-24*a*d + 5*b*c))/(1024*d**3) + c**2*x**3*sqrt(c + d*x**2)*(40*a**2*d**2 + b*c*(-24*a*d + 5*b*c))/(512*d**2) + c*x**3*(c + d*x**2)**(sympy.S(3)/2)*(40*a**2*d**2 + b*c*(-24*a*d + 5*b*c))/(384*d**2) + x**3*(c + d*x**2)**(sympy.S(5)/2)*(40*a**2*d**2 + b*c*(-24*a*d + 5*b*c))/(320*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_627():
    f = x*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*(c + d*x**2)**(sympy.S(11)/2)/(11*d**3) - 2*b*(c + d*x**2)**(sympy.S(9)/2)*(-a*d + b*c)/(9*d**3) + (c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c)**2/(7*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_628():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)
    F = b*x*(a + b*x**2)*(c + d*x**2)**(sympy.S(7)/2)/(10*d) - 3*b*x*(c + d*x**2)**(sympy.S(7)/2)*(-4*a*d + b*c)/(80*d**2) + c**3*(80*a**2*d**2 - 20*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(256*d**(sympy.S(5)/2)) + c**2*x*sqrt(c + d*x**2)*(80*a**2*d**2 - 20*a*b*c*d + 3*b**2*c**2)/(256*d**2) + c*x*(c + d*x**2)**(sympy.S(3)/2)*(80*a**2*d**2 - 20*a*b*c*d + 3*b**2*c**2)/(384*d**2) + x*(c + d*x**2)**(sympy.S(5)/2)*(80*a**2*d**2 - 20*a*b*c*d + 3*b**2*c**2)/(480*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_629():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x
    F = -a**2*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**2)/sqrt(c)) + a**2*c**2*sqrt(c + d*x**2) + a**2*c*(c + d*x**2)**(sympy.S(3)/2)/3 + a**2*(c + d*x**2)**(sympy.S(5)/2)/5 + b**2*(c + d*x**2)**(sympy.S(9)/2)/(9*d**2) - b*(c + d*x**2)**(sympy.S(7)/2)*(-2*a*d + b*c)/(7*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_630():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**2
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(c*x) + b**2*x*(c + d*x**2)**(sympy.S(7)/2)/(8*d) - 5*c**2*(-16*a*d*(3*a*d + b*c) + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(128*d**(sympy.S(3)/2)) - 5*c*x*sqrt(c + d*x**2)*(-16*a*d*(3*a*d + b*c) + b**2*c**2)/(128*d) - x*(c + d*x**2)**(sympy.S(3)/2)*(-80*a*d*(3*a*d + b*c) + 5*b**2*c**2)/(192*d) - x*(c + d*x**2)**(sympy.S(5)/2)*(-16*a*d*(3*a*d + b*c) + b**2*c**2)/(48*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_631():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**3
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(2*c*x**2) - a*c**(sympy.S(3)/2)*(5*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/2 + a*c*sqrt(c + d*x**2)*(5*a*d + 4*b*c)/2 + a*(c + d*x**2)**(sympy.S(3)/2)*(5*a*d + 4*b*c)/6 + a*(c + d*x**2)**(sympy.S(5)/2)*(5*a*d + 4*b*c)/(10*c) + b**2*(c + d*x**2)**(sympy.S(7)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_632():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**4
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(3*c*x**3) - 2*a*(c + d*x**2)**(sympy.S(7)/2)*(2*a*d + 3*b*c)/(3*c**2*x) + 5*c*(4*a*d*(2*a*d + 3*b*c) + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*sqrt(d)) + x*sqrt(c + d*x**2)*(5*a*d*(2*a*d + 3*b*c)/4 + 5*b**2*c**2/16) + x*(c + d*x**2)**(sympy.S(3)/2)*(20*a*d*(2*a*d + 3*b*c) + 5*b**2*c**2)/(24*c) + x*(c + d*x**2)**(sympy.S(5)/2)*(4*a*d*(2*a*d + 3*b*c) + b**2*c**2)/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_633():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**5
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(4*c*x**4) - a*(c + d*x**2)**(sympy.S(7)/2)*(3*a*d + 8*b*c)/(8*c**2*x**2) - sqrt(c)*(5*a*d*(3*a*d + 8*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/8 + sqrt(c + d*x**2)*(5*a*d*(3*a*d + 8*b*c)/8 + b**2*c**2) + (c + d*x**2)**(sympy.S(3)/2)*(5*a*d*(3*a*d + 8*b*c) + 8*b**2*c**2)/(24*c) + (c + d*x**2)**(sympy.S(5)/2)*(5*a*d*(3*a*d + 8*b*c) + 8*b**2*c**2)/(40*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_634():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**6
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(5*c*x**5) - 2*a*(c + d*x**2)**(sympy.S(7)/2)*(a*d + 5*b*c)/(15*c**2*x**3) + sqrt(d)*(8*a*d*(a*d + 5*b*c) + 15*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/8 + d*x*sqrt(c + d*x**2)*(8*a*d*(a*d + 5*b*c) + 15*b**2*c**2)/(8*c) + d*x*(c + d*x**2)**(sympy.S(3)/2)*(8*a*d*(a*d + 5*b*c) + 15*b**2*c**2)/(12*c**2) - (c + d*x**2)**(sympy.S(5)/2)*(8*a*d*(a*d + 5*b*c) + 15*b**2*c**2)/(15*c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_635():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2)/x**7
    F = -a**2*(c + d*x**2)**(sympy.S(7)/2)/(6*c*x**6) - a*(c + d*x**2)**(sympy.S(7)/2)*(a*d + 12*b*c)/(24*c**2*x**4) + 5*d*sqrt(c + d*x**2)*(a*d*(a*d + 12*b*c) + 8*b**2*c**2)/(16*c) + 5*d*(c + d*x**2)**(sympy.S(3)/2)*(a*d*(a*d + 12*b*c) + 8*b**2*c**2)/(48*c**2) - (c + d*x**2)**(sympy.S(5)/2)*(a*d*(a*d + 12*b*c) + 8*b**2*c**2)/(16*c**2*x**2) - 5*d*(a*d*(a*d + 12*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(16*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_636():
    f = x**4*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = b**2*x**7*sqrt(c + d*x**2)/(8*d) - b*x**5*sqrt(c + d*x**2)*(-16*a*d + 7*b*c)/(48*d**2) + c**2*(48*a**2*d**2 + 5*b*c*(-16*a*d + 7*b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(128*d**(sympy.S(9)/2)) - c*x*sqrt(c + d*x**2)*(48*a**2*d**2 + 5*b*c*(-16*a*d + 7*b*c))/(128*d**4) + x**3*sqrt(c + d*x**2)*(48*a**2*d**2 + 5*b*c*(-16*a*d + 7*b*c))/(192*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_637():
    f = x**3*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = b**2*(c + d*x**2)**(sympy.S(7)/2)/(7*d**4) - b*(c + d*x**2)**(sympy.S(5)/2)*(-2*a*d + 3*b*c)/(5*d**4) - c*sqrt(c + d*x**2)*(-a*d + b*c)**2/d**4 + (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*d + 3*b*c)/(3*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_638():
    f = x**2*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = b**2*x**5*sqrt(c + d*x**2)/(6*d) - b*x**3*sqrt(c + d*x**2)*(-12*a*d + 5*b*c)/(24*d**2) - c*(8*a**2*d**2 + b*c*(-12*a*d + 5*b*c))*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*d**(sympy.S(7)/2)) + x*sqrt(c + d*x**2)*(8*a**2*d**2 + b*c*(-12*a*d + 5*b*c))/(16*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_639():
    f = x*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = b**2*(c + d*x**2)**(sympy.S(5)/2)/(5*d**3) - 2*b*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**3) + sqrt(c + d*x**2)*(-a*d + b*c)**2/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_640():
    f = (a + b*x**2)**2/sqrt(c + d*x**2)
    F = b*x*(a + b*x**2)*sqrt(c + d*x**2)/(4*d) - 3*b*x*sqrt(c + d*x**2)*(-2*a*d + b*c)/(8*d**2) + (8*a**2*d**2 - 8*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_641():
    f = (a + b*x**2)**2/(x*sqrt(c + d*x**2))
    F = -a**2*atanh(sqrt(c + d*x**2)/sqrt(c))/sqrt(c) + b**2*(c + d*x**2)**(sympy.S(3)/2)/(3*d**2) - b*sqrt(c + d*x**2)*(-2*a*d + b*c)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_642():
    f = (a + b*x**2)**2/(x**2*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(c*x) + b**2*x*sqrt(c + d*x**2)/(2*d) - b*(-4*a*d + b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_643():
    f = (a + b*x**2)**2/(x**3*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(2*c*x**2) - a*(-a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*c**(sympy.S(3)/2)) + b**2*sqrt(c + d*x**2)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_644():
    f = (a + b*x**2)**2/(x**4*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(3*c*x**3) - 2*a*sqrt(c + d*x**2)*(-a*d + 3*b*c)/(3*c**2*x) + b**2*atanh(sqrt(d)*x/sqrt(c + d*x**2))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_645():
    f = (a + b*x**2)**2/(x**5*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(4*c*x**4) - a*sqrt(c + d*x**2)*(-3*a*d + 8*b*c)/(8*c**2*x**2) - (3*a**2*d**2 - 8*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_646():
    f = (a + b*x**2)**2/(x**6*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(5*c*x**5) - 2*a*sqrt(c + d*x**2)*(-2*a*d + 5*b*c)/(15*c**2*x**3) - sqrt(c + d*x**2)*(-4*a*d*(-2*a*d + 5*b*c) + 15*b**2*c**2)/(15*c**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_647():
    f = (a + b*x**2)**2/(x**7*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(6*c*x**6) - a*sqrt(c + d*x**2)*(-5*a*d + 12*b*c)/(24*c**2*x**4) - sqrt(c + d*x**2)*(5*a**2*d**2 - 12*a*b*c*d + 8*b**2*c**2)/(16*c**3*x**2) + d*(5*a**2*d**2 - 12*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(16*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_648():
    f = x**4*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*x**5*sqrt(c + d*x**2)/(6*d**2) - c*(24*a**2*d**2 - 60*a*b*c*d + 35*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*d**(sympy.S(9)/2)) + x*sqrt(c + d*x**2)*(24*a**2*d**2 - 60*a*b*c*d + 35*b**2*c**2)/(16*d**4) + x**5*(-a*d + b*c)**2/(c*d**2*sqrt(c + d*x**2)) - x**3*sqrt(c + d*x**2)*(24*a**2*d**2 - 60*a*b*c*d + 35*b**2*c**2)/(24*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_649():
    f = x**3*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*(c + d*x**2)**(sympy.S(5)/2)/(5*d**4) - b*(c + d*x**2)**(sympy.S(3)/2)*(-2*a*d + 3*b*c)/(3*d**4) + c*(-a*d + b*c)**2/(d**4*sqrt(c + d*x**2)) + sqrt(c + d*x**2)*(-a*d + b*c)*(-a*d + 3*b*c)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_650():
    f = x**2*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*x**3*sqrt(c + d*x**2)/(4*d**2) + (8*a**2*d**2 - 24*a*b*c*d + 15*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(7)/2)) + x**3*(-a*d + b*c)**2/(c*d**2*sqrt(c + d*x**2)) - x*sqrt(c + d*x**2)*(8*a**2*d**2 - 24*a*b*c*d + 15*b**2*c**2)/(8*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_651():
    f = x*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = b**2*(c + d*x**2)**(sympy.S(3)/2)/(3*d**3) - 2*b*sqrt(c + d*x**2)*(-a*d + b*c)/d**3 - (-a*d + b*c)**2/(d**3*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_652():
    f = (a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = -b*(-4*a*d + 3*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*d**(sympy.S(5)/2)) + b*x*sqrt(c + d*x**2)*(-2*a*d + 3*b*c)/(2*c*d**2) - x*(a + b*x**2)*(-a*d + b*c)/(c*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_653():
    f = (a + b*x**2)**2/(x*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2*atanh(sqrt(c + d*x**2)/sqrt(c))/c**(sympy.S(3)/2) + b**2*sqrt(c + d*x**2)/d**2 + (-a*d + b*c)**2/(c*d**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_654():
    f = (a + b*x**2)**2/(x**2*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(c*x*sqrt(c + d*x**2)) + b**2*atanh(sqrt(d)*x/sqrt(c + d*x**2))/d**(sympy.S(3)/2) - x*(-2*a*d*(-a*d + b*c) + b**2*c**2)/(c**2*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_655():
    f = (a + b*x**2)**2/(x**3*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(2*c*x**2*sqrt(c + d*x**2)) - a*(-3*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*c**(sympy.S(5)/2)) + (-3*a**2*d/c + 4*a*b - 2*b**2*c/d)/(2*c*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_656():
    f = (a + b*x**2)**2/(x**4*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(3*c*x**3*sqrt(c + d*x**2)) - 2*a*(-2*a*d + 3*b*c)/(3*c**2*x*sqrt(c + d*x**2)) + x*(-4*a*d*(-2*a*d + 3*b*c) + 3*b**2*c**2)/(3*c**3*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_657():
    f = (a + b*x**2)**2/(x**5*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(4*c*x**4*sqrt(c + d*x**2)) - a*(-5*a*d + 8*b*c)/(8*c**2*x**2*sqrt(c + d*x**2)) + (-3*a*d*(-5*a*d + 8*b*c) + 8*b**2*c**2)/(8*c**3*sqrt(c + d*x**2)) - (-3*a*d*(-5*a*d + 8*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(8*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_658():
    f = (a + b*x**2)**2/(x**6*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(5*c*x**5*sqrt(c + d*x**2)) - 2*a*(-3*a*d + 5*b*c)/(15*c**2*x**3*sqrt(c + d*x**2)) - (-8*a*d*(-3*a*d + 5*b*c) + 15*b**2*c**2)/(15*c**3*x*sqrt(c + d*x**2)) - 2*d*x*(-8*a*d*(-3*a*d + 5*b*c) + 15*b**2*c**2)/(15*c**4*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_659():
    f = (a + b*x**2)**2/(x**7*(c + d*x**2)**(sympy.S(3)/2))
    F = -a**2/(6*c*x**6*sqrt(c + d*x**2)) - a*(-7*a*d + 12*b*c)/(24*c**2*x**4*sqrt(c + d*x**2)) - (-5*a*d*(-7*a*d + 12*b*c) + 24*b**2*c**2)/(48*c**3*x**2*sqrt(c + d*x**2)) - d*(-5*a*d*(-7*a*d + 12*b*c) + 24*b**2*c**2)/(16*c**4*sqrt(c + d*x**2)) + d*(-5*a*d*(-7*a*d + 12*b*c) + 24*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(16*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_660():
    f = x**4*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*x**5/(4*d**2*sqrt(c + d*x**2)) + (8*a**2*d**2 - 40*a*b*c*d + 35*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(9)/2)) + x**5*(-a*d + b*c)**2/(3*c*d**2*(c + d*x**2)**(sympy.S(3)/2)) + x**3*(8*a**2*d**2 - 40*a*b*c*d + 35*b**2*c**2)/(12*c*d**3*sqrt(c + d*x**2)) - x*sqrt(c + d*x**2)*(8*a**2*d**2 - 40*a*b*c*d + 35*b**2*c**2)/(8*c*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_661():
    f = x**3*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*(c + d*x**2)**(sympy.S(3)/2)/(3*d**4) - b*sqrt(c + d*x**2)*(-2*a*d + 3*b*c)/d**4 + c*(-a*d + b*c)**2/(3*d**4*(c + d*x**2)**(sympy.S(3)/2)) - (-a*d + b*c)*(-a*d + 3*b*c)/(d**4*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_662():
    f = x**2*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*x*sqrt(c + d*x**2)/(2*d**3) + 2*b*x*(-a*d + b*c)/(d**3*sqrt(c + d*x**2)) - b*(-4*a*d + 5*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*d**(sympy.S(7)/2)) + x**3*(-a*d + b*c)**2/(3*c*d**2*(c + d*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_663():
    f = x*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*sqrt(c + d*x**2)/d**3 + 2*b*(-a*d + b*c)/(d**3*sqrt(c + d*x**2)) - (-a*d + b*c)**2/(3*d**3*(c + d*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_664():
    f = (a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = b**2*atanh(sqrt(d)*x/sqrt(c + d*x**2))/d**(sympy.S(5)/2) - x*(a + b*x**2)*(-a*d + b*c)/(3*c*d*(c + d*x**2)**(sympy.S(3)/2)) - x*(-a*d + b*c)*(2*a*d + 3*b*c)/(3*c**2*d**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_665():
    f = (a + b*x**2)**2/(x*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2*atanh(sqrt(c + d*x**2)/sqrt(c))/c**(sympy.S(5)/2) + (a**2/c**2 - b**2/d**2)/sqrt(c + d*x**2) + (-a*d + b*c)**2/(3*c*d**2*(c + d*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_666():
    f = (a + b*x**2)**2/(x**2*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2/(c*x*(c + d*x**2)**(sympy.S(3)/2)) + 4*a*x*(-2*a*d + b*c)/(3*c**3*sqrt(c + d*x**2)) + x*(2*a*(-2*a*d + b*c) + b**2*c*x**2)/(3*c**2*(c + d*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_667():
    f = (a + b*x**2)**2/(x**3*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2/(2*c*x**2*(c + d*x**2)**(sympy.S(3)/2)) + a*(-5*a*d + 4*b*c)/(2*c**3*sqrt(c + d*x**2)) - a*(-5*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*c**(sympy.S(7)/2)) + (-5*a**2*d/c + 4*a*b - 2*b**2*c/d)/(6*c*(c + d*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_668():
    f = (a + b*x**2)**2/(x**4*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2/(3*c*x**3*(c + d*x**2)**(sympy.S(3)/2)) - 2*a*(-a*d + b*c)/(c**2*x*(c + d*x**2)**(sympy.S(3)/2)) + x*(-8*a*d*(-a*d + b*c) + b**2*c**2)/(3*c**3*(c + d*x**2)**(sympy.S(3)/2)) + x*(-16*a*d*(-a*d + b*c) + 2*b**2*c**2)/(3*c**4*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_669():
    f = (a + b*x**2)**2/(x**5*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2/(4*c*x**4*(c + d*x**2)**(sympy.S(3)/2)) - a*(-7*a*d + 8*b*c)/(8*c**2*x**2*(c + d*x**2)**(sympy.S(3)/2)) + (-5*a*d*(-7*a*d + 8*b*c) + 8*b**2*c**2)/(24*c**3*(c + d*x**2)**(sympy.S(3)/2)) + (-5*a*d*(-7*a*d + 8*b*c) + 8*b**2*c**2)/(8*c**4*sqrt(c + d*x**2)) - (-5*a*d*(-7*a*d + 8*b*c) + 8*b**2*c**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(8*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_670():
    f = (a + b*x**2)**2/(x**6*(c + d*x**2)**(sympy.S(5)/2))
    F = -a**2/(5*c*x**5*(c + d*x**2)**(sympy.S(3)/2)) - 2*a*(-4*a*d + 5*b*c)/(15*c**2*x**3*(c + d*x**2)**(sympy.S(3)/2)) - (-4*a*d*(-4*a*d + 5*b*c) + 5*b**2*c**2)/(5*c**3*x*(c + d*x**2)**(sympy.S(3)/2)) - 4*d*x*(-4*a*d*(-4*a*d + 5*b*c) + 5*b**2*c**2)/(15*c**4*(c + d*x**2)**(sympy.S(3)/2)) - 8*d*x*(-4*a*d*(-4*a*d + 5*b*c) + 5*b**2*c**2)/(15*c**5*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_671():
    f = x**5/(sqrt(d*x**2)*(a + b*x**2))
    F = a**(sympy.S(3)/2)*x*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(5)/2)*sqrt(d*x**2)) - a*x**2/(b**2*sqrt(d*x**2)) + x**4/(3*b*sqrt(d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_672():
    f = x**3/(sqrt(d*x**2)*(a + b*x**2))
    F = -sqrt(a)*x*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(3)/2)*sqrt(d*x**2)) + x**2/(b*sqrt(d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_673():
    f = x/(sqrt(d*x**2)*(a + b*x**2))
    F = x*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b)*sqrt(d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_674():
    f = 1/(x*sqrt(d*x**2)*(a + b*x**2))
    F = -1/(a*sqrt(d*x**2)) - sqrt(b)*x*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_675():
    f = 1/(x**3*sqrt(d*x**2)*(a + b*x**2))
    F = -1/(3*a*x**2*sqrt(d*x**2)) + b/(a**2*sqrt(d*x**2)) + b**(sympy.S(3)/2)*x*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(5)/2)*sqrt(d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_676():
    f = x**4*sqrt(c + d*x**2)/(a + b*x**2)
    F = a**(sympy.S(3)/2)*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**3 + x**3*sqrt(c + d*x**2)/(4*b) + x*sqrt(c + d*x**2)*(-4*a*d + b*c)/(8*b**2*d) - (-8*a**2*d**2 + 4*a*b*c*d + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*b**3*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_677():
    f = x**3*sqrt(c + d*x**2)/(a + b*x**2)
    F = -a*sqrt(c + d*x**2)/b**2 + a*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(5)/2) + (c + d*x**2)**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_678():
    f = x**2*sqrt(c + d*x**2)/(a + b*x**2)
    F = -sqrt(a)*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**2 + x*sqrt(c + d*x**2)/(2*b) + (-2*a*d + b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_679():
    f = x*sqrt(c + d*x**2)/(a + b*x**2)
    F = sqrt(c + d*x**2)/b - sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_680():
    f = sqrt(c + d*x**2)/(a + b*x**2)
    F = sqrt(d)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b + sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_681():
    f = sqrt(c + d*x**2)/(x*(a + b*x**2))
    F = -sqrt(c)*atanh(sqrt(c + d*x**2)/sqrt(c))/a + sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_682():
    f = sqrt(c + d*x**2)/(x**2*(a + b*x**2))
    F = -sqrt(c + d*x**2)/(a*x) - sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_683():
    f = sqrt(c + d*x**2)/(x**3*(a + b*x**2))
    F = -sqrt(c + d*x**2)/(2*a*x**2) - sqrt(b)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/a**2 + (-a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_684():
    f = sqrt(c + d*x**2)/(x**4*(a + b*x**2))
    F = -sqrt(c + d*x**2)/(3*a*x**3) + sqrt(c + d*x**2)*(-a*d + 3*b*c)/(3*a**2*c*x) + b*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_685():
    f = x**4*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**4 + d*x**5*sqrt(c + d*x**2)/(6*b) + x**3*sqrt(c + d*x**2)*(-6*a*d + 7*b*c)/(24*b**2) + x*sqrt(c + d*x**2)*(8*a**2*d**2 - 10*a*b*c*d + b**2*c**2)/(16*b**3*d) - (-2*a*d + b*c)*(-8*a**2*d**2 + 8*a*b*c*d + b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*b**4*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_686():
    f = x**3*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)
    F = -a*(c + d*x**2)**(sympy.S(3)/2)/(3*b**2) - a*sqrt(c + d*x**2)*(-a*d + b*c)/b**3 + a*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(7)/2) + (c + d*x**2)**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_687():
    f = x**2*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)
    F = -sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**3 + d*x**3*sqrt(c + d*x**2)/(4*b) + x*sqrt(c + d*x**2)*(-4*a*d + 5*b*c)/(8*b**2) + (8*a**2*d**2 - 12*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*b**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_688():
    f = x*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)
    F = (c + d*x**2)**(sympy.S(3)/2)/(3*b) + sqrt(c + d*x**2)*(-a*d + b*c)/b**2 - (-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_689():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)
    F = d*x*sqrt(c + d*x**2)/(2*b) + sqrt(d)*(-2*a*d + 3*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**2) + (-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_690():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x*(a + b*x**2))
    F = d*sqrt(c + d*x**2)/b - c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**2)/sqrt(c))/a + (-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_691():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**2*(a + b*x**2))
    F = d**(sympy.S(3)/2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b - c*sqrt(c + d*x**2)/(a*x) - (-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(3)/2)*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_692():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**3*(a + b*x**2))
    F = -c*sqrt(c + d*x**2)/(2*a*x**2) + sqrt(c)*(-3*a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2) - (-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a**2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_693():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**4*(a + b*x**2))
    F = -c*sqrt(c + d*x**2)/(3*a*x**3) + sqrt(c + d*x**2)*(-4*a*d + 3*b*c)/(3*a**2*x) + (-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_694():
    f = x**4*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)
    F = a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**5 + d*x**5*(c + d*x**2)**(sympy.S(3)/2)/(8*b) + d*x**5*sqrt(c + d*x**2)*(-8*a*d + 11*b*c)/(48*b**2) + x**3*sqrt(c + d*x**2)*(48*a**2*d**2 - 104*a*b*c*d + 59*b**2*c**2)/(192*b**3) + x*sqrt(c + d*x**2)*(-64*a**3*d**3 + 144*a**2*b*c*d**2 - 88*a*b**2*c**2*d + 5*b**3*c**3)/(128*b**4*d) - (-128*a**4*d**4 + 320*a**3*b*c*d**3 - 240*a**2*b**2*c**2*d**2 + 40*a*b**3*c**3*d + 5*b**4*c**4)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(128*b**5*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_695():
    f = x**3*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)
    F = -a*(c + d*x**2)**(sympy.S(5)/2)/(5*b**2) - a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*b**3) - a*sqrt(c + d*x**2)*(-a*d + b*c)**2/b**4 + a*(-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(9)/2) + (c + d*x**2)**(sympy.S(7)/2)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_696():
    f = x**2*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)
    F = -sqrt(a)*(-a*d + b*c)**(sympy.S(5)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/b**4 + d*x**3*(c + d*x**2)**(sympy.S(3)/2)/(6*b) + d*x**3*sqrt(c + d*x**2)*(-2*a*d + 3*b*c)/(8*b**2) + x*sqrt(c + d*x**2)*(8*a**2*d**2 - 18*a*b*c*d + 11*b**2*c**2)/(16*b**3) + (-16*a**3*d**3 + 40*a**2*b*c*d**2 - 30*a*b**2*c**2*d + 5*b**3*c**3)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*b**4*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_697():
    f = x*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)
    F = (c + d*x**2)**(sympy.S(5)/2)/(5*b) + (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*b**2) + sqrt(c + d*x**2)*(-a*d + b*c)**2/b**3 - (-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_698():
    f = (c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)
    F = d*x*(c + d*x**2)**(sympy.S(3)/2)/(4*b) + d*x*sqrt(c + d*x**2)*(-4*a*d + 7*b*c)/(8*b**2) + sqrt(d)*(8*a**2*d**2 - 20*a*b*c*d + 15*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*b**3) + (-a*d + b*c)**(sympy.S(5)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_699():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x*(a + b*x**2))
    F = d*(c + d*x**2)**(sympy.S(3)/2)/(3*b) + d*sqrt(c + d*x**2)*(-a*d + 2*b*c)/b**2 - c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**2)/sqrt(c))/a + (-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_700():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**2*(a + b*x**2))
    F = d**(sympy.S(3)/2)*(-2*a*d + 5*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**2) - c*(c + d*x**2)**(sympy.S(3)/2)/(a*x) + d*x*sqrt(c + d*x**2)*(a*d + 2*b*c)/(2*a*b) - (-a*d + b*c)**(sympy.S(5)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(3)/2)*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_701():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**3*(a + b*x**2))
    F = -c*(c + d*x**2)**(sympy.S(3)/2)/(2*a*x**2) + d*sqrt(c + d*x**2)*(2*a*d + b*c)/(2*a*b) + c**(sympy.S(3)/2)*(-5*a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2) - (-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a**2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_702():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**4*(a + b*x**2))
    F = d**(sympy.S(5)/2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b - c*(c + d*x**2)**(sympy.S(3)/2)/(3*a*x**3) + c*sqrt(c + d*x**2)*(-2*a*d + b*c)/(a**2*x) + (-a*d + b*c)**(sympy.S(5)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(5)/2)*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_703():
    f = x**5/((a + b*x**2)*sqrt(c + d*x**2))
    F = -a**2*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(b**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + (c + d*x**2)**(sympy.S(3)/2)/(3*b*d**2) - sqrt(c + d*x**2)*(a*d + b*c)/(b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_704():
    f = x**3/((a + b*x**2)*sqrt(c + d*x**2))
    F = a*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + sqrt(c + d*x**2)/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_705():
    f = x/((a + b*x**2)*sqrt(c + d*x**2))
    F = -atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_706():
    f = 1/(x*(a + b*x**2)*sqrt(c + d*x**2))
    F = sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*sqrt(-a*d + b*c)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_707():
    f = 1/(x**3*(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/(2*a*c*x**2) - b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a**2*sqrt(-a*d + b*c)) + (a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_708():
    f = x**4/((a + b*x**2)*sqrt(c + d*x**2))
    F = a**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(b**2*sqrt(-a*d + b*c)) + x*sqrt(c + d*x**2)/(2*b*d) - (2*a*d + b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_709():
    f = x**2/((a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(a)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(b*sqrt(-a*d + b*c)) + atanh(sqrt(d)*x/sqrt(c + d*x**2))/(b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_710():
    f = 1/((a + b*x**2)*sqrt(c + d*x**2))
    F = atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_711():
    f = 1/(x**2*(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/(a*c*x) - b*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_712():
    f = 1/(x**4*(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/(3*a*c*x**3) + sqrt(c + d*x**2)*(2*a*d + 3*b*c)/(3*a**2*c**2*x) + b**2*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_713():
    f = x**4/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = a**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(b*(-a*d + b*c)**(sympy.S(3)/2)) - c*x/(d*sqrt(c + d*x**2)*(-a*d + b*c)) + atanh(sqrt(d)*x/sqrt(c + d*x**2))/(b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_714():
    f = x**3/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = a*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2)) - c/(d*sqrt(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_715():
    f = x**2/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -sqrt(a)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(-a*d + b*c)**(sympy.S(3)/2) + x/(sqrt(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_716():
    f = x/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(3)/2) + 1/(sqrt(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_717():
    f = 1/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -d*x/(c*sqrt(c + d*x**2)*(-a*d + b*c)) + b*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_718():
    f = 1/(x*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -d/(c*sqrt(c + d*x**2)*(-a*d + b*c)) + b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*(-a*d + b*c)**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_719():
    f = 1/(x**2*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -d/(c*x*sqrt(c + d*x**2)*(-a*d + b*c)) - sqrt(c + d*x**2)*(-2*a*d + b*c)/(a*c**2*x*(-a*d + b*c)) - b**2*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_720():
    f = 1/(x**3*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -1/(2*a*c*x**2*sqrt(c + d*x**2)) - d*(-3*a*d + b*c)/(2*a*c**2*sqrt(c + d*x**2)*(-a*d + b*c)) - b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a**2*(-a*d + b*c)**(sympy.S(3)/2)) + (3*a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_721():
    f = 1/(x**4*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -d/(c*x**3*sqrt(c + d*x**2)*(-a*d + b*c)) - sqrt(c + d*x**2)*(-4*a*d + b*c)/(3*a*c**2*x**3*(-a*d + b*c)) + sqrt(c + d*x**2)*(-4*a*d + 3*b*c)*(2*a*d + b*c)/(3*a**2*c**3*x*(-a*d + b*c)) + b**3*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_722():
    f = x**4/((a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = a**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(-a*d + b*c)**(sympy.S(5)/2) - c*x/(3*d*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + x*(-4*a*d + b*c)/(3*d*sqrt(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_723():
    f = x**3/((a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = a*sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(5)/2) - a/(sqrt(c + d*x**2)*(-a*d + b*c)**2) - c/(3*d*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_724():
    f = x**2/((a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -sqrt(a)*b*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(-a*d + b*c)**(sympy.S(5)/2) + x/((c + d*x**2)**(sympy.S(3)/2)*(-3*a*d + 3*b*c)) + x*(a*d + 2*b*c)/(3*c*sqrt(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_725():
    f = x/((a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(5)/2) + b/(sqrt(c + d*x**2)*(-a*d + b*c)**2) + 1/((c + d*x**2)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_726():
    f = 1/((a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -d*x/(3*c*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*x*(-2*a*d + 5*b*c)/(3*c**2*sqrt(c + d*x**2)*(-a*d + b*c)**2) + b**2*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_727():
    f = 1/(x*(a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -d/(3*c*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(-a*d + 2*b*c)/(c**2*sqrt(c + d*x**2)*(-a*d + b*c)**2) + b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a*(-a*d + b*c)**(sympy.S(5)/2)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_728():
    f = 1/(x**2*(a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -d/(3*c*x*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(-4*a*d + 7*b*c)/(3*c**2*x*sqrt(c + d*x**2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)*(-4*a*d + b*c)*(-2*a*d + 3*b*c)/(3*a*c**3*x*(-a*d + b*c)**2) - b**3*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_729():
    f = 1/(x**3*(a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -1/(2*a*c*x**2*(c + d*x**2)**(sympy.S(3)/2)) - d*(-5*a*d + 3*b*c)/(6*a*c**2*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(5*a**2*d**2 - 8*a*b*c*d + b**2*c**2)/(2*a*c**3*sqrt(c + d*x**2)*(-a*d + b*c)**2) - b**(sympy.S(7)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(a**2*(-a*d + b*c)**(sympy.S(5)/2)) + (5*a*d + 2*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**2*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_730():
    f = 1/(x**4*(a + b*x**2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -d/(3*c*x**3*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(-2*a*d + 3*b*c)/(c**2*x**3*sqrt(c + d*x**2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)*(8*a**2*d**2 - 12*a*b*c*d + b**2*c**2)/(3*a*c**3*x**3*(-a*d + b*c)**2) + sqrt(c + d*x**2)*(-2*a*d + b*c)*(-8*a**2*d**2 + 8*a*b*c*d + 3*b**2*c**2)/(3*a**2*c**4*x*(-a*d + b*c)**2) + b**4*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_731():
    f = x**4*sqrt(c + d*x**2)/(a + b*x**2)**2
    F = -sqrt(a)*(-4*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*b**3*sqrt(-a*d + b*c)) - x**3*sqrt(c + d*x**2)/(2*b*(a + b*x**2)) + x*sqrt(c + d*x**2)/b**2 + (-4*a*d + b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_732():
    f = x**3*sqrt(c + d*x**2)/(a + b*x**2)**2
    F = a*(c + d*x**2)**(sympy.S(3)/2)/(2*b*(a + b*x**2)*(-a*d + b*c)) + sqrt(c + d*x**2)*(-3*a*d + 2*b*c)/(2*b**2*(-a*d + b*c)) - (-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_733():
    f = x**2*sqrt(c + d*x**2)/(a + b*x**2)**2
    F = -x*sqrt(c + d*x**2)/(2*b*(a + b*x**2)) + sqrt(d)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b**2 + (-2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*b**2*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_734():
    f = x*sqrt(c + d*x**2)/(a + b*x**2)**2
    F = -sqrt(c + d*x**2)/(2*b*(a + b*x**2)) - d*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_735():
    f = sqrt(c + d*x**2)/(a + b*x**2)**2
    F = x*sqrt(c + d*x**2)/(2*a*(a + b*x**2)) + c*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_736():
    f = sqrt(c + d*x**2)/(x*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)/(2*a*(a + b*x**2)) - sqrt(c)*atanh(sqrt(c + d*x**2)/sqrt(c))/a**2 + (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_737():
    f = sqrt(c + d*x**2)/(x**2*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)/(2*a*x*(a + b*x**2)) - 3*sqrt(c + d*x**2)/(2*a**2*x) - (-2*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_738():
    f = sqrt(c + d*x**2)/(x**3*(a + b*x**2)**2)
    F = -sqrt(c + d*x**2)/(2*a*x**2*(a + b*x**2)) - b*sqrt(c + d*x**2)/(a**2*(a + b*x**2)) - sqrt(b)*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*sqrt(-a*d + b*c)) + (-a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_739():
    f = sqrt(c + d*x**2)/(x**4*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)/(2*a*x**3*(a + b*x**2)) - 5*sqrt(c + d*x**2)/(6*a**2*x**3) + sqrt(c + d*x**2)*(-2*a*d + 15*b*c)/(6*a**3*c*x) + b*(-4*a*d + 5*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_740():
    f = x**4*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**2
    F = -3*sqrt(a)*(-2*a*d + b*c)*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*b**4) - x**3*(c + d*x**2)**(sympy.S(3)/2)/(2*b*(a + b*x**2)) + 3*d*x**3*sqrt(c + d*x**2)/(4*b**2) + x*sqrt(c + d*x**2)*(-12*a*d + 9*b*c)/(8*b**3) + (24*a**2*d**2 - 24*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*b**4*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_741():
    f = x**3*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**2
    F = a*(c + d*x**2)**(sympy.S(5)/2)/(2*b*(a + b*x**2)*(-a*d + b*c)) + (c + d*x**2)**(sympy.S(3)/2)*(-5*a*d + 2*b*c)/(6*b**2*(-a*d + b*c)) + sqrt(c + d*x**2)*(-5*a*d + 2*b*c)/(2*b**3) - (-5*a*d + 2*b*c)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_742():
    f = x**2*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**2
    F = -x*(c + d*x**2)**(sympy.S(3)/2)/(2*b*(a + b*x**2)) + d*x*sqrt(c + d*x**2)/b**2 + sqrt(d)*(-4*a*d + 3*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**3) + (-4*a*d + b*c)*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_743():
    f = x*(c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**2
    F = -(c + d*x**2)**(sympy.S(3)/2)/(2*b*(a + b*x**2)) + 3*d*sqrt(c + d*x**2)/(2*b**2) - 3*d*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_744():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**2
    F = d**(sympy.S(3)/2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b**2 + x*sqrt(c + d*x**2)*(-a*d + b*c)/(2*a*b*(a + b*x**2)) + sqrt(-a*d + b*c)*(2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_745():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)*(-a*d + b*c)/(2*a*b*(a + b*x**2)) - c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**2)/sqrt(c))/a**2 + sqrt(-a*d + b*c)*(a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_746():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**2*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)*(-a*d + b*c)/(2*a*b*x*(a + b*x**2)) - sqrt(c + d*x**2)*(-a*d + 3*b*c)/(2*a**2*b*x) - 3*c*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_747():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**3*(a + b*x**2)**2)
    F = -c*sqrt(c + d*x**2)/(2*a*x**2*(a + b*x**2)) - sqrt(c + d*x**2)*(-a*d + 2*b*c)/(2*a**2*(a + b*x**2)) + sqrt(c)*(-3*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3) - sqrt(-a*d + b*c)*(-a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_748():
    f = (c + d*x**2)**(sympy.S(3)/2)/(x**4*(a + b*x**2)**2)
    F = sqrt(c + d*x**2)*(-a*d + b*c)/(2*a*b*x**3*(a + b*x**2)) - sqrt(c + d*x**2)*(-3*a*d + 5*b*c)/(6*a**2*b*x**3) + sqrt(c + d*x**2)*(-11*a*d + 15*b*c)/(6*a**3*x) + (-2*a*d + 5*b*c)*sqrt(-a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_749():
    f = x**4*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)**2
    F = -sqrt(a)*(-8*a*d + 3*b*c)*(-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*b**5) - x**3*(c + d*x**2)**(sympy.S(5)/2)/(2*b*(a + b*x**2)) + 2*d*x**3*(c + d*x**2)**(sympy.S(3)/2)/(3*b**2) + d*x**3*sqrt(c + d*x**2)*(-8*a*d + 7*b*c)/(8*b**3) + x*sqrt(c + d*x**2)*(32*a**2*d**2 - 52*a*b*c*d + 19*b**2*c**2)/(16*b**4) + (-64*a**3*d**3 + 120*a**2*b*c*d**2 - 60*a*b**2*c**2*d + 5*b**3*c**3)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*b**5*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_750():
    f = x**3*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)**2
    F = a*(c + d*x**2)**(sympy.S(7)/2)/(2*b*(a + b*x**2)*(-a*d + b*c)) + (c + d*x**2)**(sympy.S(5)/2)*(-7*a*d + 2*b*c)/(10*b**2*(-a*d + b*c)) + (c + d*x**2)**(sympy.S(3)/2)*(-7*a*d + 2*b*c)/(6*b**3) + sqrt(c + d*x**2)*(-7*a*d + 2*b*c)*(-a*d + b*c)/(2*b**4) - (-7*a*d + 2*b*c)*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_751():
    f = x**2*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)**2
    F = -x*(c + d*x**2)**(sympy.S(5)/2)/(2*b*(a + b*x**2)) + 3*d*x*(c + d*x**2)**(sympy.S(3)/2)/(4*b**2) + d*x*sqrt(c + d*x**2)*(-12*a*d + 11*b*c)/(8*b**3) + sqrt(d)*(24*a**2*d**2 - 40*a*b*c*d + 15*b**2*c**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*b**4) + (-6*a*d + b*c)*(-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_752():
    f = x*(c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)**2
    F = -(c + d*x**2)**(sympy.S(5)/2)/(2*b*(a + b*x**2)) + 5*d*(c + d*x**2)**(sympy.S(3)/2)/(6*b**2) + 5*d*sqrt(c + d*x**2)*(-a*d + b*c)/(2*b**3) - 5*d*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_753():
    f = (c + d*x**2)**(sympy.S(5)/2)/(a + b*x**2)**2
    F = d**(sympy.S(3)/2)*(-4*a*d + 5*b*c)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*b**3) + x*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(2*a*b*(a + b*x**2)) - d*x*sqrt(c + d*x**2)*(-2*a*d + b*c)/(2*a*b**2) + (-a*d + b*c)**(sympy.S(3)/2)*(4*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_754():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x*(a + b*x**2)**2)
    F = (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(2*a*b*(a + b*x**2)) - d*sqrt(c + d*x**2)*(-3*a*d + b*c)/(2*a*b**2) - c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**2)/sqrt(c))/a**2 + (-a*d + b*c)**(sympy.S(3)/2)*(3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_755():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**2*(a + b*x**2)**2)
    F = d**(sympy.S(5)/2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/b**2 + (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(2*a*b*x*(a + b*x**2)) - c*sqrt(c + d*x**2)*(-a*d + 3*b*c)/(2*a**2*b*x) - (-a*d + b*c)**(sympy.S(3)/2)*(2*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2)*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_756():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**3*(a + b*x**2)**2)
    F = -c*(c + d*x**2)**(sympy.S(3)/2)/(2*a*x**2*(a + b*x**2)) - sqrt(c + d*x**2)*(-a*d + b*c)*(-a*d + 2*b*c)/(2*a**2*b*(a + b*x**2)) + c**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3) - (-a*d + b*c)**(sympy.S(3)/2)*(a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_757():
    f = (c + d*x**2)**(sympy.S(5)/2)/(x**4*(a + b*x**2)**2)
    F = (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(2*a*b*x**3*(a + b*x**2)) - c*sqrt(c + d*x**2)*(-3*a*d + 5*b*c)/(6*a**2*b*x**3) + sqrt(c + d*x**2)*(3*a**2*d**2 - 20*a*b*c*d + 15*b**2*c**2)/(6*a**3*b*x) + 5*c*(-a*d + b*c)**(sympy.S(3)/2)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_758():
    f = x**4/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = -sqrt(a)*(-2*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*b**2*(-a*d + b*c)**(sympy.S(3)/2)) + a*x*sqrt(c + d*x**2)/(2*b*(a + b*x**2)*(-a*d + b*c)) + atanh(sqrt(d)*x/sqrt(c + d*x**2))/(b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_759():
    f = x**3/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = a*sqrt(c + d*x**2)/(2*b*(a + b*x**2)*(-a*d + b*c)) - (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_760():
    f = x**2/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = -x*sqrt(c + d*x**2)/((a + b*x**2)*(-2*a*d + 2*b*c)) + c*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_761():
    f = x/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/((a + b*x**2)*(-2*a*d + 2*b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_762():
    f = 1/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = b*x*sqrt(c + d*x**2)/(2*a*(a + b*x**2)*(-a*d + b*c)) + (-2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_763():
    f = 1/(x*(a + b*x**2)**2*sqrt(c + d*x**2))
    F = b*sqrt(c + d*x**2)/(2*a*(a + b*x**2)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*(-a*d + b*c)**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_764():
    f = 1/(x**2*(a + b*x**2)**2*sqrt(c + d*x**2))
    F = b*sqrt(c + d*x**2)/(2*a*x*(a + b*x**2)*(-a*d + b*c)) - sqrt(c + d*x**2)*(-2*a*d + 3*b*c)/(2*a**2*c*x*(-a*d + b*c)) - b*(-4*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_765():
    f = 1/(x**3*(a + b*x**2)**2*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/(2*a*c*x**2*(a + b*x**2)) - b*sqrt(c + d*x**2)*(-a*d + 2*b*c)/(2*a**2*c*(a + b*x**2)*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*(-a*d + b*c)**(sympy.S(3)/2)) + (a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_766():
    f = 1/(x**4*(a + b*x**2)**2*sqrt(c + d*x**2))
    F = b*sqrt(c + d*x**2)/(2*a*x**3*(a + b*x**2)*(-a*d + b*c)) - sqrt(c + d*x**2)*(-2*a*d + 5*b*c)/(6*a**2*c*x**3*(-a*d + b*c)) + sqrt(c + d*x**2)*(-4*a**2*d**2 - 8*a*b*c*d + 15*b**2*c**2)/(6*a**3*c**2*x*(-a*d + b*c)) + b**2*(-6*a*d + 5*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_767():
    f = x**4/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = -3*sqrt(a)*c*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*(-a*d + b*c)**(sympy.S(5)/2)) + a*x/(2*b*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + x*(a*d + 2*b*c)/(2*b*sqrt(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_768():
    f = x**3/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = a/(2*b*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + (a*d + 2*b*c)/(2*b*sqrt(c + d*x**2)*(-a*d + b*c)**2) - (a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*sqrt(b)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_769():
    f = x**2/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = -3*d*x/(2*sqrt(c + d*x**2)*(-a*d + b*c)**2) - x/((a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + 2*b*c)) + (2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_770():
    f = x/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = 3*sqrt(b)*d*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*(-a*d + b*c)**(sympy.S(5)/2)) - 3*d/(2*sqrt(c + d*x**2)*(-a*d + b*c)**2) - 1/((a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_771():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = b*x/(2*a*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + d*x*(2*a*d + b*c)/(2*a*c*sqrt(c + d*x**2)*(-a*d + b*c)**2) + b*(-4*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_772():
    f = 1/(x*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = b/(2*a*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + d*(2*a*d + b*c)/(2*a*c*sqrt(c + d*x**2)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*(-5*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*(-a*d + b*c)**(sympy.S(5)/2)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_773():
    f = 1/(x**2*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = b/(2*a*x*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + d*(2*a*d + b*c)/(2*a*c*x*sqrt(c + d*x**2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)*(4*a**2*d**2 - 4*a*b*c*d + 3*b**2*c**2)/(2*a**2*c**2*x*(-a*d + b*c)**2) - 3*b**2*(-2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_774():
    f = 1/(x**3*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = -1/(2*a*c*x**2*(a + b*x**2)*sqrt(c + d*x**2)) - b*(-a*d + 2*b*c)/(2*a**2*c*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) - d*(3*a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)/(2*a**2*c**2*sqrt(c + d*x**2)*(-a*d + b*c)**2) - b**(sympy.S(5)/2)*(-7*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*(-a*d + b*c)**(sympy.S(5)/2)) + (3*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_775():
    f = 1/(x**4*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2))
    F = b/(2*a*x**3*(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + d*(2*a*d + b*c)/(2*a*c*x**3*sqrt(c + d*x**2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)*(8*a**2*d**2 - 4*a*b*c*d + 5*b**2*c**2)/(6*a**2*c**2*x**3*(-a*d + b*c)**2) + sqrt(c + d*x**2)*(16*a**3*d**3 - 8*a**2*b*c*d**2 - 14*a*b**2*c**2*d + 15*b**3*c**3)/(6*a**3*c**3*x*(-a*d + b*c)**2) + b**3*(-8*a*d + 5*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_776():
    f = x**4/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = -sqrt(a)*(2*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*(-a*d + b*c)**(sympy.S(7)/2)) + a*x/(2*b*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + x*(11*a*d + 4*b*c)/(6*sqrt(c + d*x**2)*(-a*d + b*c)**3) + x*(3*a*d + 2*b*c)/(6*b*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_777():
    f = x**3/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = a/(2*b*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt(b)*(3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*(-a*d + b*c)**(sympy.S(7)/2)) + (3*a*d + 2*b*c)/(2*sqrt(c + d*x**2)*(-a*d + b*c)**3) + (3*a*d + 2*b*c)/(6*b*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_778():
    f = x**2/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = -5*d*x/(6*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - x/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)) - d*x*(2*a*d + 13*b*c)/(6*c*sqrt(c + d*x**2)*(-a*d + b*c)**3) + b*(4*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_779():
    f = x/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = 5*b**(sympy.S(3)/2)*d*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*(-a*d + b*c)**(sympy.S(7)/2)) - 5*b*d/(2*sqrt(c + d*x**2)*(-a*d + b*c)**3) - 5*d/(6*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - 1/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_780():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = b*x/(2*a*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + d*x*(2*a*d + 3*b*c)/(6*a*c*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + d*x*(-4*a**2*d**2 + 16*a*b*c*d + 3*b**2*c**2)/(6*a*c**2*sqrt(c + d*x**2)*(-a*d + b*c)**3) + b**2*(-6*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_781():
    f = 1/(x*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = b/(2*a*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + d*(2*a*d + 3*b*c)/(6*a*c*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + d*(-2*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/(2*a*c**2*sqrt(c + d*x**2)*(-a*d + b*c)**3) + b**(sympy.S(5)/2)*(-7*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**2*(-a*d + b*c)**(sympy.S(7)/2)) - atanh(sqrt(c + d*x**2)/sqrt(c))/(a**2*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_782():
    f = 1/(x**2*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = b/(2*a*x*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + d*(2*a*d + 3*b*c)/(6*a*c*x*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + d*(-8*a**2*d**2 + 20*a*b*c*d + 3*b**2*c**2)/(6*a*c**2*x*sqrt(c + d*x**2)*(-a*d + b*c)**3) - sqrt(c + d*x**2)*(-16*a**3*d**3 + 40*a**2*b*c*d**2 - 18*a*b**2*c**2*d + 9*b**3*c**3)/(6*a**2*c**3*x*(-a*d + b*c)**3) - b**3*(-8*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_783():
    f = 1/(x**3*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = -1/(2*a*c*x**2*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)) - b*(-a*d + 2*b*c)/(2*a**2*c*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*(5*a**2*d**2 - 6*a*b*c*d + 6*b**2*c**2)/(6*a**2*c**2*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - d*(-a*d + 2*b*c)*(5*a**2*d**2 - a*b*c*d + b**2*c**2)/(2*a**2*c**3*sqrt(c + d*x**2)*(-a*d + b*c)**3) - b**(sympy.S(7)/2)*(-9*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**2)/sqrt(-a*d + b*c))/(2*a**3*(-a*d + b*c)**(sympy.S(7)/2)) + (5*a*d + 4*b*c)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*a**3*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_784():
    f = 1/(x**4*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(5)/2))
    F = b/(2*a*x**3*(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + d*(2*a*d + 3*b*c)/(6*a*c*x**3*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + d*(-4*a**2*d**2 + 8*a*b*c*d + b**2*c**2)/(2*a*c**2*x**3*sqrt(c + d*x**2)*(-a*d + b*c)**3) - sqrt(c + d*x**2)*(-16*a**3*d**3 + 32*a**2*b*c*d**2 - 6*a*b**2*c**2*d + 5*b**3*c**3)/(6*a**2*c**3*x**3*(-a*d + b*c)**3) + sqrt(c + d*x**2)*(-32*a**4*d**4 + 64*a**3*b*c*d**3 - 12*a**2*b**2*c**2*d**2 - 20*a*b**3*c**3*d + 15*b**4*c**4)/(6*a**3*c**4*x*(-a*d + b*c)**3) + 5*b**4*(-2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_785():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**2)*sqrt(a + b*x**2)
    F = 2*B*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2)/(11*b*e) - 2*a**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(11*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) + 4*a*e*sqrt(e*x)*sqrt(a + b*x**2)*(11*A*b - 5*B*a)/(231*b**2) + (e*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)*(22*A*b - 10*B*a)/(77*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_786():
    f = sqrt(e*x)*(A + B*x**2)*sqrt(a + b*x**2)
    F = 2*B*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)/(9*b*e) - 4*a**(sympy.S(5)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + 2*a**(sympy.S(5)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + 4*a*sqrt(e*x)*sqrt(a + b*x**2)*(3*A*b - B*a)/(15*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(6*A*b - 2*B*a)/(15*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_787():
    f = (A + B*x**2)*sqrt(a + b*x**2)/sqrt(e*x)
    F = 2*B*sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/2)/(7*b*e) + 2*a**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*b**(sympy.S(5)/4)*sqrt(e)*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(14*A*b - 2*B*a)/(21*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_788():
    f = (A + B*x**2)*sqrt(a + b*x**2)/(e*x)**(sympy.S(3)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(a*e*sqrt(e*x)) - 4*a**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b + B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 2*a**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(20*A*b + 4*B*a)/(5*sqrt(b)*e**2*(sqrt(a) + sqrt(b)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(10*A*b + 2*B*a)/(5*a*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_789():
    f = (A + B*x**2)*sqrt(a + b*x**2)/(e*x)**(sympy.S(5)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(3*a*e*(e*x)**(sympy.S(3)/2)) + sqrt(e*x)*sqrt(a + b*x**2)*(2*A*b + 2*B*a)/(3*a*e**3) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(2*A*b + 2*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_790():
    f = (A + B*x**2)*sqrt(a + b*x**2)/(e*x)**(sympy.S(7)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(5*a*e*(e*x)**(sympy.S(5)/2)) + 4*sqrt(b)*sqrt(e*x)*sqrt(a + b*x**2)*(A*b + 5*B*a)/(5*a*e**4*(sqrt(a) + sqrt(b)*x)) - sqrt(a + b*x**2)*(2*A*b + 10*B*a)/(5*a*e**3*sqrt(e*x)) - 4*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + 5*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*a**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 2*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*a**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_791():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**(sympy.S(9)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(7*a*x**(sympy.S(7)/2)) + sqrt(a + b*x**2)*(2*A*b - 14*B*a)/(21*a*x**(sympy.S(3)/2)) - 2*b**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 7*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(21*a**(sympy.S(5)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_792():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**(sympy.S(11)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(9*a*x**(sympy.S(9)/2)) + sqrt(a + b*x**2)*(2*A*b - 6*B*a)/(15*a*x**(sympy.S(5)/2)) - 4*b**(sympy.S(3)/2)*sqrt(x)*sqrt(a + b*x**2)*(A*b - 3*B*a)/(15*a**2*(sqrt(a) + sqrt(b)*x)) + 4*b*sqrt(a + b*x**2)*(A*b - 3*B*a)/(15*a**2*sqrt(x)) + 4*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 3*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(7)/4)*sqrt(a + b*x**2)) - 2*b**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(7)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_793():
    f = (A + B*x**2)*sqrt(a + b*x**2)/x**(sympy.S(13)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(3)/2)/(11*a*x**(sympy.S(11)/2)) + sqrt(a + b*x**2)*(10*A*b - 22*B*a)/(77*a*x**(sympy.S(7)/2)) + 4*b*sqrt(a + b*x**2)*(5*A*b - 11*B*a)/(231*a**2*x**(sympy.S(3)/2)) + 2*b**(sympy.S(7)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 11*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(231*a**(sympy.S(9)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_794():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/2)/(15*b*e) - 4*a**(sympy.S(11)/4)*e**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) + 8*a**2*e*sqrt(e*x)*sqrt(a + b*x**2)*(3*A*b - B*a)/(231*b**2) + 4*a*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)*(3*A*b - B*a)/(77*b*e) + (e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2)*(6*A*b - 2*B*a)/(33*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_795():
    f = sqrt(e*x)*(A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/2)/(13*b*e) - 8*a**(sympy.S(9)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(13*A*b - 3*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + 4*a**(sympy.S(9)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(13*A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + 8*a**2*sqrt(e*x)*sqrt(a + b*x**2)*(13*A*b - 3*B*a)/(195*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + 4*a*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(13*A*b - 3*B*a)/(195*b*e) + (e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)*(26*A*b - 6*B*a)/(117*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_796():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/sqrt(e*x)
    F = 2*B*sqrt(e*x)*(a + b*x**2)**(sympy.S(5)/2)/(11*b*e) + 4*a**(sympy.S(7)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(11*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(77*b**(sympy.S(5)/4)*sqrt(e)*sqrt(a + b*x**2)) + 4*a*sqrt(e*x)*sqrt(a + b*x**2)*(11*A*b - B*a)/(77*b*e) + sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/2)*(22*A*b - 2*B*a)/(77*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_797():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(3)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(5)/2)/(a*e*sqrt(e*x)) - 8*a**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(9*A*b + B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 4*a**(sympy.S(5)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(9*A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + 8*a*sqrt(e*x)*sqrt(a + b*x**2)*(9*A*b + B*a)/(15*sqrt(b)*e**2*(sqrt(a) + sqrt(b)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(36*A*b + 4*B*a)/(15*e**3) + (e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)*(18*A*b + 2*B*a)/(9*a*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_798():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(5)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(5)/2)/(3*a*e*(e*x)**(sympy.S(3)/2)) + 4*a**(sympy.S(3)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b + 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*b**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(28*A*b + 12*B*a)/(21*e**3) + sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/2)*(14*A*b + 6*B*a)/(21*a*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_799():
    f = (A + B*x**2)*(a + b*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(7)/2)
    F = -2*A*(a + b*x**2)**(sympy.S(5)/2)/(5*a*e*(e*x)**(sympy.S(5)/2)) - 24*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*e**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 12*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*e**(sympy.S(7)/2)*sqrt(a + b*x**2)) + 24*sqrt(b)*sqrt(e*x)*sqrt(a + b*x**2)*(A*b + B*a)/(5*e**4*(sqrt(a) + sqrt(b)*x)) + 12*b*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(A*b + B*a)/(5*a*e**5) - (a + b*x**2)**(sympy.S(3)/2)*(2*A*b + 2*B*a)/(a*e**3*sqrt(e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_800():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**2)/sqrt(a + b*x**2)
    F = 2*B*(e*x)**(sympy.S(7)/2)*sqrt(a + b*x**2)/(9*b*e) + 2*a**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(9*A*b - 7*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(a + b*x**2)) - a**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(9*A*b - 7*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(a + b*x**2)) - 2*a*e**2*sqrt(e*x)*sqrt(a + b*x**2)*(9*A*b - 7*B*a)/(15*b**(sympy.S(5)/2)*(sqrt(a) + sqrt(b)*x)) + e*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)*(18*A*b - 14*B*a)/(45*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_801():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**2)/sqrt(a + b*x**2)
    F = 2*B*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)/(7*b*e) - a**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*b**(sympy.S(9)/4)*sqrt(a + b*x**2)) + e*sqrt(e*x)*sqrt(a + b*x**2)*(14*A*b - 10*B*a)/(21*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_802():
    f = sqrt(e*x)*(A + B*x**2)/sqrt(a + b*x**2)
    F = 2*B*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)/(5*b*e) - 2*a**(sympy.S(1)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 3*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + a**(sympy.S(1)/4)*sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(10*A*b - 6*B*a)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_803():
    f = (A + B*x**2)/(sqrt(e*x)*sqrt(a + b*x**2))
    F = 2*B*sqrt(e*x)*sqrt(a + b*x**2)/(3*b*e) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt(e)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_804():
    f = (A + B*x**2)/((e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2))
    F = -2*A*sqrt(a + b*x**2)/(a*e*sqrt(e*x)) + sqrt(e*x)*sqrt(a + b*x**2)*(2*A*b + 2*B*a)/(a*sqrt(b)*e**2*(sqrt(a) + sqrt(b)*x)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(2*A*b + 2*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_805():
    f = (A + B*x**2)/((e*x)**(sympy.S(5)/2)*sqrt(a + b*x**2))
    F = -2*A*sqrt(a + b*x**2)/(3*a*e*(e*x)**(sympy.S(3)/2)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_806():
    f = (A + B*x**2)/((e*x)**(sympy.S(7)/2)*sqrt(a + b*x**2))
    F = -2*A*sqrt(a + b*x**2)/(5*a*e*(e*x)**(sympy.S(5)/2)) - 2*sqrt(b)*sqrt(e*x)*sqrt(a + b*x**2)*(3*A*b - 5*B*a)/(5*a**2*e**4*(sqrt(a) + sqrt(b)*x)) + sqrt(a + b*x**2)*(6*A*b - 10*B*a)/(5*a**2*e**3*sqrt(e*x)) + 2*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - 5*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*a**(sympy.S(7)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2)) - b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*a**(sympy.S(7)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_807():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(sympy.S(9)/2)/(7*b*e*sqrt(a + b*x**2)) - 5*a**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - 9*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(42*b**(sympy.S(13)/4)*sqrt(a + b*x**2)) - e*(e*x)**(sympy.S(5)/2)*(7*A*b - 9*B*a)/(7*b**2*sqrt(a + b*x**2)) + e**3*sqrt(e*x)*sqrt(a + b*x**2)*(35*A*b - 45*B*a)/(21*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_808():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(sympy.S(7)/2)/(5*b*e*sqrt(a + b*x**2)) - 3*a**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 7*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*b**(sympy.S(11)/4)*sqrt(a + b*x**2)) + 3*a**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 7*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*b**(sympy.S(11)/4)*sqrt(a + b*x**2)) - e*(e*x)**(sympy.S(3)/2)*(5*A*b - 7*B*a)/(5*b**2*sqrt(a + b*x**2)) + e**2*sqrt(e*x)*sqrt(a + b*x**2)*(15*A*b - 21*B*a)/(5*b**(sympy.S(5)/2)*(sqrt(a) + sqrt(b)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_809():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(sympy.S(5)/2)/(3*b*e*sqrt(a + b*x**2)) - e*sqrt(e*x)*(3*A*b - 5*B*a)/(3*b**2*sqrt(a + b*x**2)) + e**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(6*a**(sympy.S(1)/4)*b**(sympy.S(9)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_810():
    f = sqrt(e*x)*(A + B*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = (e*x)**(sympy.S(3)/2)*(A*b - B*a)/(a*b*e*sqrt(a + b*x**2)) - sqrt(e*x)*sqrt(a + b*x**2)*(A*b - 3*B*a)/(a*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 3*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_811():
    f = (A + B*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/2))
    F = sqrt(e*x)*(A*b - B*a)/(a*b*e*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(5)/4)*b**(sympy.S(5)/4)*sqrt(e)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_812():
    f = (A + B*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = -2*A/(a*e*sqrt(e*x)*sqrt(a + b*x**2)) - (e*x)**(sympy.S(3)/2)*(3*A*b - B*a)/(a**2*e**3*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(3*A*b - B*a)/(a**2*sqrt(b)*e**2*(sqrt(a) + sqrt(b)*x)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(a**(sympy.S(7)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(3*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_813():
    f = (A + B*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = -2*A/(3*a*e*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**2)) - sqrt(e*x)*(5*A*b - 3*B*a)/(3*a**2*e**3*sqrt(a + b*x**2)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 3*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(6*a**(sympy.S(9)/4)*b**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_814():
    f = (A + B*x**2)/((e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/2))
    F = -2*A/(5*a*e*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**2)) - (7*A*b - 5*B*a)/(5*a**2*e**3*sqrt(e*x)*sqrt(a + b*x**2)) - 3*sqrt(b)*sqrt(e*x)*sqrt(a + b*x**2)*(7*A*b - 5*B*a)/(5*a**3*e**4*(sqrt(a) + sqrt(b)*x)) + sqrt(a + b*x**2)*(21*A*b - 15*B*a)/(5*a**3*e**3*sqrt(e*x)) + 3*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - 5*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*a**(sympy.S(11)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2)) - 3*b**(sympy.S(1)/4)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*a**(sympy.S(11)/4)*e**(sympy.S(7)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_815():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = 2*B*(e*x)**(sympy.S(9)/2)/(3*b*e*(a + b*x**2)**(sympy.S(3)/2)) - e*(e*x)**(sympy.S(5)/2)*(A*b - 3*B*a)/(3*b**2*(a + b*x**2)**(sympy.S(3)/2)) - e**3*sqrt(e*x)*(5*A*b - 15*B*a)/(6*b**3*sqrt(a + b*x**2)) + e**(sympy.S(7)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b - 15*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*a**(sympy.S(1)/4)*b**(sympy.S(13)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_816():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(7)/2)*(A*b - B*a)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/2)) + e*(e*x)**(sympy.S(3)/2)*(A*b - 7*B*a)/(6*a*b**2*sqrt(a + b*x**2)) - e**2*sqrt(e*x)*sqrt(a + b*x**2)*(A*b - 7*B*a)/(2*a*b**(sympy.S(5)/2)*(sqrt(a) + sqrt(b)*x)) + e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 7*B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(11)/4)*sqrt(a + b*x**2)) - e**(sympy.S(5)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b - 7*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4*a**(sympy.S(3)/4)*b**(sympy.S(11)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_817():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(5)/2)*(A*b - B*a)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/2)) - e*sqrt(e*x)*(A*b + 5*B*a)/(6*a*b**2*sqrt(a + b*x**2)) + e**(sympy.S(3)/2)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*a**(sympy.S(5)/4)*b**(sympy.S(9)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_818():
    f = sqrt(e*x)*(A + B*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(3)/2)*(A*b - B*a)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(A*b + B*a)/(2*a**2*b*e*sqrt(a + b*x**2)) - sqrt(e*x)*sqrt(a + b*x**2)*(A*b + B*a)/(2*a**2*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)) + sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2)) - sqrt(e)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_819():
    f = (A + B*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(5)/2))
    F = sqrt(e*x)*(A*b - B*a)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/2)) + sqrt(e*x)*(5*A*b + B*a)/(6*a**2*b*e*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(5*A*b + B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*a**(sympy.S(9)/4)*b**(sympy.S(5)/4)*sqrt(e)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_820():
    f = (A + B*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/2))
    F = -2*A/(a*e*sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(7*A*b - B*a)/(3*a**2*e**3*(a + b*x**2)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(7*A*b - B*a)/(2*a**3*e**3*sqrt(a + b*x**2)) + sqrt(e*x)*sqrt(a + b*x**2)*(7*A*b - B*a)/(2*a**3*sqrt(b)*e**2*(sqrt(a) + sqrt(b)*x)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - B*a)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2)) + sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(7*A*b - B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_821():
    f = (A + B*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/2))
    F = -2*A/(3*a*e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)) - sqrt(e*x)*(3*A*b - B*a)/(3*a**2*e**3*(a + b*x**2)**(sympy.S(3)/2)) - sqrt(e*x)*(15*A*b - 5*B*a)/(6*a**3*e**3*sqrt(a + b*x**2)) - sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*(15*A*b - 5*B*a)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(a**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*a**(sympy.S(13)/4)*b**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_822():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = 2*b**2*(e*x)**(sympy.S(9)/2)*(c + d*x**2)**(sympy.S(3)/2)/(15*d*e**3) - 2*b*(e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2)*(-10*a*d + 3*b*c)/(55*d**2*e) - 2*c**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(11*a**2*d**2 + b*c*(-10*a*d + 3*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*d**(sympy.S(13)/4)*sqrt(c + d*x**2)) + 4*c*e*sqrt(e*x)*sqrt(c + d*x**2)*(11*a**2*d**2 + b*c*(-10*a*d + 3*b*c))/(231*d**3) + (e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(22*a**2*d**2 + 2*b*c*(-10*a*d + 3*b*c))/(77*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_823():
    f = sqrt(e*x)*(a + b*x**2)**2*sqrt(c + d*x**2)
    F = 2*b**2*(e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(3)/2)/(13*d*e**3) - 2*b*(e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)*(-26*a*d + 7*b*c)/(117*d**2*e) - 4*c**(sympy.S(5)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(39*a**2*d**2 + b*c*(-26*a*d + 7*b*c))*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + 2*c**(sympy.S(5)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(39*a**2*d**2 + b*c*(-26*a*d + 7*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + 4*c*sqrt(e*x)*sqrt(c + d*x**2)*(39*a**2*d**2 + b*c*(-26*a*d + 7*b*c))/(195*d**(sympy.S(5)/2)*(sqrt(c) + sqrt(d)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(78*a**2*d**2 + 2*b*c*(-26*a*d + 7*b*c))/(195*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_824():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/sqrt(e*x)
    F = 2*b**2*(e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2)/(11*d*e**3) - 2*b*sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)*(-22*a*d + 5*b*c)/(77*d**2*e) + 2*c**(sympy.S(3)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(77*a**2*d**2 - 22*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*d**(sympy.S(9)/4)*sqrt(e)*sqrt(c + d*x**2)) + sqrt(e*x)*sqrt(c + d*x**2)*(154*a**2*d**2 - 44*a*b*c*d + 10*b**2*c**2)/(231*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_825():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/(e*x)**(sympy.S(3)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(c*e*sqrt(e*x)) + 2*b**2*(e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)/(9*d*e**3) + 4*c**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a*d*(5*a*d + 2*b*c) + b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - 2*c**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a*d*(5*a*d + 2*b*c) + b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(-12*a*d*(5*a*d + 2*b*c) + 4*b**2*c**2)/(15*d**(sympy.S(3)/2)*e**2*(sqrt(c) + sqrt(d)*x)) - (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-6*a*d*(5*a*d + 2*b*c) + 2*b**2*c**2)/(15*c*d*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_826():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/(e*x)**(sympy.S(5)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(3*c*e*(e*x)**(sympy.S(3)/2)) + 2*b**2*sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)/(7*d*e**3) - sqrt(e*x)*sqrt(c + d*x**2)*(-14*a*d*(a*d + 2*b*c) + 2*b**2*c**2)/(21*c*d*e**3) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-14*a*d*(a*d + 2*b*c) + 2*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*c**(sympy.S(1)/4)*d**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_827():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/(e*x)**(sympy.S(7)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(5*c*e*(e*x)**(sympy.S(5)/2)) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(a*d + 10*b*c)/(5*c**2*e**3*sqrt(e*x)) + sqrt(e*x)*sqrt(c + d*x**2)*(4*a*d*(a*d + 10*b*c) + 4*b**2*c**2)/(5*c*sqrt(d)*e**4*(sqrt(c) + sqrt(d)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(2*a*d*(a*d + 10*b*c) + 2*b**2*c**2)/(5*c**2*e**5) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(2*a*d*(a*d + 10*b*c) + 2*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(3)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(4*a*d*(a*d + 10*b*c) + 4*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(3)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_828():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**(sympy.S(9)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(7*c*x**(sympy.S(7)/2)) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 14*b*c)/(21*c**2*x**(sympy.S(3)/2)) + sqrt(x)*sqrt(c + d*x**2)*(2*a*d*(-a*d + 14*b*c) + 14*b**2*c**2)/(21*c**2) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(2*a*d*(-a*d + 14*b*c) + 14*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(5)/4)*d**(sympy.S(1)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_829():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**(sympy.S(11)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(9*c*x**(sympy.S(9)/2)) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + 6*b*c)/(15*c**2*x**(sympy.S(5)/2)) + 4*sqrt(d)*sqrt(x)*sqrt(c + d*x**2)*(a*d*(-a*d + 6*b*c) + 15*b**2*c**2)/(15*c**2*(sqrt(c) + sqrt(d)*x)) - sqrt(c + d*x**2)*(2*a*d*(-a*d + 6*b*c) + 30*b**2*c**2)/(15*c**2*sqrt(x)) - 4*d**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a*d*(-a*d + 6*b*c) + 15*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(c + d*x**2)) + 2*d**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a*d*(-a*d + 6*b*c) + 15*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_830():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**(sympy.S(13)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(11*c*x**(sympy.S(11)/2)) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-5*a*d + 22*b*c)/(77*c**2*x**(sympy.S(7)/2)) - sqrt(c + d*x**2)*(10*a**2*d**2 - 44*a*b*c*d + 154*b**2*c**2)/(231*c**2*x**(sympy.S(3)/2)) + 2*d**(sympy.S(3)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 22*a*b*c*d + 77*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(9)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_831():
    f = (a + b*x**2)**2*sqrt(c + d*x**2)/x**(sympy.S(15)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(3)/2)/(13*c*x**(sympy.S(13)/2)) - 2*a*(c + d*x**2)**(sympy.S(3)/2)*(-7*a*d + 26*b*c)/(117*c**2*x**(sympy.S(9)/2)) - sqrt(c + d*x**2)*(14*a**2*d**2 - 52*a*b*c*d + 78*b**2*c**2)/(195*c**2*x**(sympy.S(5)/2)) + 4*d**(sympy.S(3)/2)*sqrt(x)*sqrt(c + d*x**2)*(7*a**2*d**2 - 26*a*b*c*d + 39*b**2*c**2)/(195*c**3*(sqrt(c) + sqrt(d)*x)) - 4*d*sqrt(c + d*x**2)*(7*a**2*d**2 - 26*a*b*c*d + 39*b**2*c**2)/(195*c**3*sqrt(x)) - 4*d**(sympy.S(5)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(7*a**2*d**2 - 26*a*b*c*d + 39*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(c + d*x**2)) + 2*d**(sympy.S(5)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(7*a**2*d**2 - 26*a*b*c*d + 39*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(x)/c**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_832():
    f = (e*x)**(sympy.S(5)/2)*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(11)/2)*(c + d*x**2)**(sympy.S(5)/2)/(21*d*e**3) - 2*b*(e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(5)/2)*(-42*a*d + 11*b*c)/(357*d**2*e) + 8*c**(sympy.S(13)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(51*a**2*d**2 + b*c*(-42*a*d + 11*b*c))*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3315*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) - 4*c**(sympy.S(13)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(51*a**2*d**2 + b*c*(-42*a*d + 11*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3315*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) - 8*c**3*e**2*sqrt(e*x)*sqrt(c + d*x**2)*(51*a**2*d**2 + b*c*(-42*a*d + 11*b*c))/(3315*d**(sympy.S(7)/2)*(sqrt(c) + sqrt(d)*x)) + 8*c**2*e*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(51*a**2*d**2 + b*c*(-42*a*d + 11*b*c))/(9945*d**3) + 4*c*(e*x)**(sympy.S(7)/2)*sqrt(c + d*x**2)*(51*a**2*d**2 + b*c*(-42*a*d + 11*b*c))/(1989*d**2*e) + (e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(3)/2)*(102*a**2*d**2 + 2*b*c*(-42*a*d + 11*b*c))/(663*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_833():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(9)/2)*(c + d*x**2)**(sympy.S(5)/2)/(19*d*e**3) - 2*b*(e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(5)/2)*(-38*a*d + 9*b*c)/(285*d**2*e) - 4*c**(sympy.S(11)/4)*e**(sympy.S(3)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(57*a**2*d**2 + b*c*(-38*a*d + 9*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4389*d**(sympy.S(13)/4)*sqrt(c + d*x**2)) + 8*c**2*e*sqrt(e*x)*sqrt(c + d*x**2)*(57*a**2*d**2 + b*c*(-38*a*d + 9*b*c))/(4389*d**3) + 4*c*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(57*a**2*d**2 + b*c*(-38*a*d + 9*b*c))/(1463*d**2*e) + (e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2)*(114*a**2*d**2 + 2*b*c*(-38*a*d + 9*b*c))/(627*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_834():
    f = sqrt(e*x)*(a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(5)/2)/(17*d*e**3) - 2*b*(e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(5)/2)*(-34*a*d + 7*b*c)/(221*d**2*e) - 8*c**(sympy.S(9)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(221*a**2*d**2 + 3*b*c*(-34*a*d + 7*b*c))*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3315*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + 4*c**(sympy.S(9)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(221*a**2*d**2 + 3*b*c*(-34*a*d + 7*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3315*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + 8*c**2*sqrt(e*x)*sqrt(c + d*x**2)*(221*a**2*d**2 + 3*b*c*(-34*a*d + 7*b*c))/(3315*d**(sympy.S(5)/2)*(sqrt(c) + sqrt(d)*x)) + 4*c*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(221*a**2*d**2 + 3*b*c*(-34*a*d + 7*b*c))/(3315*d**2*e) + (e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)*(442*a**2*d**2 + 6*b*c*(-34*a*d + 7*b*c))/(1989*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_835():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/sqrt(e*x)
    F = 2*b**2*(e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(5)/2)/(15*d*e**3) - 2*b*sqrt(e*x)*(c + d*x**2)**(sympy.S(5)/2)*(-6*a*d + b*c)/(33*d**2*e) + 4*c**(sympy.S(7)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(33*a**2*d**2 + b*c*(-6*a*d + b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*d**(sympy.S(9)/4)*sqrt(e)*sqrt(c + d*x**2)) + 4*c*sqrt(e*x)*sqrt(c + d*x**2)*(33*a**2*d**2 + b*c*(-6*a*d + b*c))/(231*d**2*e) + sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)*(66*a**2*d**2 + 2*b*c*(-6*a*d + b*c))/(231*d**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_836():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(3)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(5)/2)/(c*e*sqrt(e*x)) + 2*b**2*(e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(5)/2)/(13*d*e**3) + 8*c**(sympy.S(5)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-13*a*d*(9*a*d + 2*b*c) + 3*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - 4*c**(sympy.S(5)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-13*a*d*(9*a*d + 2*b*c) + 3*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - 8*c*sqrt(e*x)*sqrt(c + d*x**2)*(-13*a*d*(9*a*d + 2*b*c) + 3*b**2*c**2)/(195*d**(sympy.S(3)/2)*e**2*(sqrt(c) + sqrt(d)*x)) - (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-52*a*d*(9*a*d + 2*b*c) + 12*b**2*c**2)/(195*d*e**3) - (e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)*(-26*a*d*(9*a*d + 2*b*c) + 6*b**2*c**2)/(117*c*d*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_837():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(5)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(5)/2)/(3*c*e*(e*x)**(sympy.S(3)/2)) + 2*b**2*sqrt(e*x)*(c + d*x**2)**(sympy.S(5)/2)/(11*d*e**3) - 4*c**(sympy.S(3)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-11*a*d*(7*a*d + 6*b*c) + 3*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*d**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(-44*a*d*(7*a*d + 6*b*c) + 12*b**2*c**2)/(231*d*e**3) - sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)*(-22*a*d*(7*a*d + 6*b*c) + 6*b**2*c**2)/(231*c*d*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_838():
    f = (a + b*x**2)**2*(c + d*x**2)**(sympy.S(3)/2)/(e*x)**(sympy.S(7)/2)
    F = -2*a**2*(c + d*x**2)**(sympy.S(5)/2)/(5*c*e*(e*x)**(sympy.S(5)/2)) - 2*a*(c + d*x**2)**(sympy.S(5)/2)*(a*d + 2*b*c)/(c**2*e**3*sqrt(e*x)) - 8*c**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(9*a*d*(a*d + 2*b*c) + b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) + 4*c**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(9*a*d*(a*d + 2*b*c) + b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) + sqrt(e*x)*sqrt(c + d*x**2)*(72*a*d*(a*d + 2*b*c) + 8*b**2*c**2)/(15*sqrt(d)*e**4*(sqrt(c) + sqrt(d)*x)) + (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(36*a*d*(a*d + 2*b*c) + 4*b**2*c**2)/(15*c*e**5) + (e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)*(18*a*d*(a*d + 2*b*c) + 2*b**2*c**2)/(9*c**2*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_839():
    f = (e*x)**(sympy.S(5)/2)*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = 2*b**2*(e*x)**(sympy.S(11)/2)*sqrt(c + d*x**2)/(13*d*e**3) - 2*b*(e*x)**(sympy.S(7)/2)*sqrt(c + d*x**2)*(-26*a*d + 11*b*c)/(117*d**2*e) + 2*c**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(117*a**2*d**2 + 7*b*c*(-26*a*d + 11*b*c))*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) - c**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(117*a**2*d**2 + 7*b*c*(-26*a*d + 11*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(195*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) - 2*c*e**2*sqrt(e*x)*sqrt(c + d*x**2)*(117*a**2*d**2 + 7*b*c*(-26*a*d + 11*b*c))/(195*d**(sympy.S(7)/2)*(sqrt(c) + sqrt(d)*x)) + e*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(234*a**2*d**2 + 14*b*c*(-26*a*d + 11*b*c))/(585*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_840():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = 2*b**2*(e*x)**(sympy.S(9)/2)*sqrt(c + d*x**2)/(11*d*e**3) - 2*b*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(-22*a*d + 9*b*c)/(77*d**2*e) - c**(sympy.S(3)/4)*e**(sympy.S(3)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(77*a**2*d**2 + 5*b*c*(-22*a*d + 9*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*d**(sympy.S(13)/4)*sqrt(c + d*x**2)) + e*sqrt(e*x)*sqrt(c + d*x**2)*(154*a**2*d**2 + 10*b*c*(-22*a*d + 9*b*c))/(231*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_841():
    f = sqrt(e*x)*(a + b*x**2)**2/sqrt(c + d*x**2)
    F = 2*b**2*(e*x)**(sympy.S(7)/2)*sqrt(c + d*x**2)/(9*d*e**3) - 2*b*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-18*a*d + 7*b*c)/(45*d**2*e) - 2*c**(sympy.S(1)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(15*a**2*d**2 + b*c*(-18*a*d + 7*b*c))*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + c**(sympy.S(1)/4)*sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(15*a**2*d**2 + b*c*(-18*a*d + 7*b*c))*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + sqrt(e*x)*sqrt(c + d*x**2)*(30*a**2*d**2 + 2*b*c*(-18*a*d + 7*b*c))/(15*d**(sympy.S(5)/2)*(sqrt(c) + sqrt(d)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_842():
    f = (a + b*x**2)**2/(sqrt(e*x)*sqrt(c + d*x**2))
    F = 2*b**2*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(7*d*e**3) - 2*b*sqrt(e*x)*sqrt(c + d*x**2)*(-14*a*d + 5*b*c)/(21*d**2*e) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(21*a**2*d**2 - 14*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*c**(sympy.S(1)/4)*d**(sympy.S(9)/4)*sqrt(e)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_843():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(c*e*sqrt(e*x)) + 2*b**2*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(5*d*e**3) - sqrt(e*x)*sqrt(c + d*x**2)*(-10*a*d*(a*d + 2*b*c) + 6*b**2*c**2)/(5*c*d**(sympy.S(3)/2)*e**2*(sqrt(c) + sqrt(d)*x)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-10*a*d*(a*d + 2*b*c) + 6*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(3)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-5*a*d*(a*d + 2*b*c) + 3*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(3)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_844():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(3*c*e*(e*x)**(sympy.S(3)/2)) + 2*b**2*sqrt(e*x)*sqrt(c + d*x**2)/(3*d*e**3) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a**2*d**2 - 6*a*b*c*d + b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(3*c**(sympy.S(5)/4)*d**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_845():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(7)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(5*c*e*(e*x)**(sympy.S(5)/2)) - 2*a*sqrt(c + d*x**2)*(-3*a*d + 10*b*c)/(5*c**2*e**3*sqrt(e*x)) + sqrt(e*x)*sqrt(c + d*x**2)*(-6*a**2*d**2 + 20*a*b*c*d + 10*b**2*c**2)/(5*c**2*sqrt(d)*e**4*(sqrt(c) + sqrt(d)*x)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-6*a**2*d**2 + 20*a*b*c*d + 10*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(7)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a**2*d**2 + 10*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(7)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_846():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(9)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(7*c*e*(e*x)**(sympy.S(7)/2)) - 2*a*sqrt(c + d*x**2)*(-5*a*d + 14*b*c)/(21*c**2*e**3*(e*x)**(sympy.S(3)/2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 14*a*b*c*d + 21*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(21*c**(sympy.S(9)/4)*d**(sympy.S(1)/4)*e**(sympy.S(9)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_847():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(11)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(9*c*e*(e*x)**(sympy.S(9)/2)) - 2*a*sqrt(c + d*x**2)*(-7*a*d + 18*b*c)/(45*c**2*e**3*(e*x)**(sympy.S(5)/2)) + 2*sqrt(d)*sqrt(e*x)*sqrt(c + d*x**2)*(7*a**2*d**2 - 18*a*b*c*d + 15*b**2*c**2)/(15*c**3*e**6*(sqrt(c) + sqrt(d)*x)) - sqrt(c + d*x**2)*(14*a**2*d**2 - 36*a*b*c*d + 30*b**2*c**2)/(15*c**3*e**5*sqrt(e*x)) - 2*d**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(7*a**2*d**2 - 18*a*b*c*d + 15*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*c**(sympy.S(11)/4)*e**(sympy.S(11)/2)*sqrt(c + d*x**2)) + d**(sympy.S(1)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(7*a**2*d**2 - 18*a*b*c*d + 15*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*c**(sympy.S(11)/4)*e**(sympy.S(11)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_848():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(13)/2)*sqrt(c + d*x**2))
    F = -2*a**2*sqrt(c + d*x**2)/(11*c*e*(e*x)**(sympy.S(11)/2)) - 2*a*sqrt(c + d*x**2)*(-9*a*d + 22*b*c)/(77*c**2*e**3*(e*x)**(sympy.S(7)/2)) - sqrt(c + d*x**2)*(-10*a*d*(-9*a*d + 22*b*c) + 154*b**2*c**2)/(231*c**3*e**5*(e*x)**(sympy.S(3)/2)) - d**(sympy.S(3)/4)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-5*a*d*(-9*a*d + 22*b*c) + 77*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(231*c**(sympy.S(13)/4)*e**(sympy.S(13)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_849():
    f = (e*x)**(sympy.S(7)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(9)/2)*sqrt(c + d*x**2)/(11*d**2*e) - 5*c**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(77*a**2*d**2 - 198*a*b*c*d + 117*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(462*d**(sympy.S(17)/4)*sqrt(c + d*x**2)) + e**3*sqrt(e*x)*sqrt(c + d*x**2)*(385*a**2*d**2 - 990*a*b*c*d + 585*b**2*c**2)/(231*d**4) + (e*x)**(sympy.S(9)/2)*(-a*d + b*c)**2/(c*d**2*e*sqrt(c + d*x**2)) - e*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(77*a**2*d**2 - 198*a*b*c*d + 117*b**2*c**2)/(77*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_850():
    f = (e*x)**(sympy.S(5)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(7)/2)*sqrt(c + d*x**2)/(9*d**2*e) - c**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(45*a**2*d**2 - 126*a*b*c*d + 77*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(15*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) + c**(sympy.S(1)/4)*e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(45*a**2*d**2 - 126*a*b*c*d + 77*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(30*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) + e**2*sqrt(e*x)*sqrt(c + d*x**2)*(45*a**2*d**2 - 126*a*b*c*d + 77*b**2*c**2)/(15*d**(sympy.S(7)/2)*(sqrt(c) + sqrt(d)*x)) + (e*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(c*d**2*e*sqrt(c + d*x**2)) - e*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(45*a**2*d**2 - 126*a*b*c*d + 77*b**2*c**2)/(45*c*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_851():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(7*d**2*e) + (e*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/(c*d**2*e*sqrt(c + d*x**2)) - e*sqrt(e*x)*sqrt(c + d*x**2)*(21*a**2*d**2 - 70*a*b*c*d + 45*b**2*c**2)/(21*c*d**3) + e**(sympy.S(3)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(21*a**2*d**2 - 70*a*b*c*d + 45*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(42*c**(sympy.S(1)/4)*d**(sympy.S(13)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_852():
    f = sqrt(e*x)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(3)/2)
    F = 2*b**2*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(5*d**2*e) + (e*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/(c*d**2*e*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(5*a**2*d**2 - 30*a*b*c*d + 21*b**2*c**2)/(5*c*d**(sympy.S(5)/2)*(sqrt(c) + sqrt(d)*x)) + sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 30*a*b*c*d + 21*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(3)/4)*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) - sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 30*a*b*c*d + 21*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*c**(sympy.S(3)/4)*d**(sympy.S(11)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_853():
    f = (a + b*x**2)**2/(sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2))
    F = 2*b**2*sqrt(e*x)*sqrt(c + d*x**2)/(3*d**2*e) + sqrt(e*x)*(-a*d + b*c)**2/(c*d**2*e*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a**2*d**2 - 6*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(6*c**(sympy.S(5)/4)*d**(sympy.S(9)/4)*sqrt(e)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_854():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -2*a**2/(c*e*sqrt(e*x)*sqrt(c + d*x**2)) - (e*x)**(sympy.S(3)/2)*(3*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(c**2*d*e**3*sqrt(c + d*x**2)) + sqrt(e*x)*sqrt(c + d*x**2)*(3*a**2*d**2 - 2*a*b*c*d + 3*b**2*c**2)/(c**2*d**(sympy.S(3)/2)*e**2*(sqrt(c) + sqrt(d)*x)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(3*a**2*d**2 - 2*a*b*c*d + 3*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(c**(sympy.S(7)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(3*a**2*d**2 - 2*a*b*c*d + 3*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*c**(sympy.S(7)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_855():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -2*a**2/(3*c*e*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**2)) - sqrt(e*x)*(5*a**2*d**2 - 6*a*b*c*d + 3*b**2*c**2)/(3*c**2*d*e**3*sqrt(c + d*x**2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a*d*(-5*a*d + 6*b*c) + 3*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(6*c**(sympy.S(9)/4)*d**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_856():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -2*a**2/(5*c*e*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**2)) - 2*a*(-7*a*d + 10*b*c)/(5*c**2*e**3*sqrt(e*x)*sqrt(c + d*x**2)) + (e*x)**(sympy.S(3)/2)*(-3*a*d*(-7*a*d + 10*b*c) + 5*b**2*c**2)/(5*c**3*e**5*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(-3*a*d*(-7*a*d + 10*b*c) + 5*b**2*c**2)/(5*c**3*sqrt(d)*e**4*(sqrt(c) + sqrt(d)*x)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a*d*(-7*a*d + 10*b*c) + 5*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(5*c**(sympy.S(11)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-3*a*d*(-7*a*d + 10*b*c) + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*c**(sympy.S(11)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_857():
    f = (e*x)**(sympy.S(7)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = 2*b**2*(e*x)**(sympy.S(9)/2)/(7*d**2*e*sqrt(c + d*x**2)) + (e*x)**(sympy.S(9)/2)*(-a*d + b*c)**2/(3*c*d**2*e*(c + d*x**2)**(sympy.S(3)/2)) + e*(e*x)**(sympy.S(5)/2)*(7*a**2*d**2 - 42*a*b*c*d + 39*b**2*c**2)/(14*c*d**3*sqrt(c + d*x**2)) - e**3*sqrt(e*x)*sqrt(c + d*x**2)*(35*a**2*d**2 - 210*a*b*c*d + 195*b**2*c**2)/(42*c*d**4) + e**(sympy.S(7)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(35*a**2*d**2 - 210*a*b*c*d + 195*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(84*c**(sympy.S(1)/4)*d**(sympy.S(17)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_858():
    f = (e*x)**(sympy.S(5)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = 2*b**2*(e*x)**(sympy.S(7)/2)/(5*d**2*e*sqrt(c + d*x**2)) + (e*x)**(sympy.S(7)/2)*(-a*d + b*c)**2/(3*c*d**2*e*(c + d*x**2)**(sympy.S(3)/2)) + e*(e*x)**(sympy.S(3)/2)*(5*a**2*d**2 - 70*a*b*c*d + 77*b**2*c**2)/(30*c*d**3*sqrt(c + d*x**2)) - e**2*sqrt(e*x)*sqrt(c + d*x**2)*(5*a**2*d**2 - 70*a*b*c*d + 77*b**2*c**2)/(10*c*d**(sympy.S(7)/2)*(sqrt(c) + sqrt(d)*x)) + e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 70*a*b*c*d + 77*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*c**(sympy.S(3)/4)*d**(sympy.S(15)/4)*sqrt(c + d*x**2)) - e**(sympy.S(5)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 - 70*a*b*c*d + 77*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(20*c**(sympy.S(3)/4)*d**(sympy.S(15)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_859():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = 2*b**2*(e*x)**(sympy.S(5)/2)/(3*d**2*e*sqrt(c + d*x**2)) + (e*x)**(sympy.S(5)/2)*(-a*d + b*c)**2/(3*c*d**2*e*(c + d*x**2)**(sympy.S(3)/2)) + e*sqrt(e*x)*(-a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)/(6*c*d**3*sqrt(c + d*x**2)) - e**(sympy.S(3)/2)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*c**(sympy.S(5)/4)*d**(sympy.S(13)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_860():
    f = sqrt(e*x)*(a + b*x**2)**2/(c + d*x**2)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*c*d**2*e*(c + d*x**2)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(-a*d + b*c)*(a*d + 3*b*c)/(2*c**2*d**2*e*sqrt(c + d*x**2)) + sqrt(e*x)*sqrt(c + d*x**2)*(-a**2*d**2 - 2*a*b*c*d + 7*b**2*c**2)/(2*c**2*d**(sympy.S(5)/2)*(sqrt(c) + sqrt(d)*x)) - sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-a**2*d**2 - 2*a*b*c*d + 7*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*c**(sympy.S(7)/4)*d**(sympy.S(11)/4)*sqrt(c + d*x**2)) + sqrt(e)*sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(-a**2*d**2 - 2*a*b*c*d + 7*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4*c**(sympy.S(7)/4)*d**(sympy.S(11)/4)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_861():
    f = (a + b*x**2)**2/(sqrt(e*x)*(c + d*x**2)**(sympy.S(5)/2))
    F = sqrt(e*x)*(-a*d + b*c)**2/(3*c*d**2*e*(c + d*x**2)**(sympy.S(3)/2)) - sqrt(e*x)*(-a*d + b*c)*(5*a*d + 7*b*c)/(6*c**2*d**2*e*sqrt(c + d*x**2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a**2*d**2 + 2*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*c**(sympy.S(9)/4)*d**(sympy.S(9)/4)*sqrt(e)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_862():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -2*a**2/(c*e*sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(7*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(3*c**2*d*e**3*(c + d*x**2)**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(a*d*(-7*a*d + 2*b*c) + b**2*c**2)/(2*c**3*d*e**3*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(a*d*(-7*a*d + 2*b*c) + b**2*c**2)/(2*c**3*d**(sympy.S(3)/2)*e**2*(sqrt(c) + sqrt(d)*x)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a*d*(-7*a*d + 2*b*c) + b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(2*c**(sympy.S(11)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(a*d*(-7*a*d + 2*b*c) + b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(4*c**(sympy.S(11)/4)*d**(sympy.S(7)/4)*e**(sympy.S(3)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_863():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -2*a**2/(3*c*e*(e*x)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)) - sqrt(e*x)*(3*a**2*d**2 - 2*a*b*c*d + b**2*c**2)/(3*c**2*d*e**3*(c + d*x**2)**(sympy.S(3)/2)) + sqrt(e*x)*(5*a*d*(-3*a*d + 2*b*c) + b**2*c**2)/(6*c**3*d*e**3*sqrt(c + d*x**2)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(5*a*d*(-3*a*d + 2*b*c) + b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(12*c**(sympy.S(13)/4)*d**(sympy.S(5)/4)*e**(sympy.S(5)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_864():
    f = (a + b*x**2)**2/((e*x)**(sympy.S(7)/2)*(c + d*x**2)**(sympy.S(5)/2))
    F = -2*a**2/(5*c*e*(e*x)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2)) - 2*a*(-11*a*d + 10*b*c)/(5*c**2*e**3*sqrt(e*x)*(c + d*x**2)**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(77*a**2*d**2 - 70*a*b*c*d + 5*b**2*c**2)/(15*c**3*e**5*(c + d*x**2)**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(77*a**2*d**2 - 70*a*b*c*d + 5*b**2*c**2)/(10*c**4*e**5*sqrt(c + d*x**2)) - sqrt(e*x)*sqrt(c + d*x**2)*(77*a**2*d**2 - 70*a*b*c*d + 5*b**2*c**2)/(10*c**4*sqrt(d)*e**4*(sqrt(c) + sqrt(d)*x)) + sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(77*a**2*d**2 - 70*a*b*c*d + 5*b**2*c**2)*elliptic_e(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(10*c**(sympy.S(15)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2)) - sqrt((c + d*x**2)/(sqrt(c) + sqrt(d)*x)**2)*(sqrt(c) + sqrt(d)*x)*(77*a**2*d**2 - 70*a*b*c*d + 5*b**2*c**2)*elliptic_f(2*atan(d**(sympy.S(1)/4)*sqrt(e*x)/(c**(sympy.S(1)/4)*sqrt(e))), sympy.S.Half)/(20*c**(sympy.S(15)/4)*d**(sympy.S(3)/4)*e**(sympy.S(7)/2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_865():
    f = (e*x)**(sympy.S(7)/2)*sqrt(c - d*x**2)/(a - b*x**2)
    F = ((Integer(2) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(21) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(7) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(14) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(21) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(21) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_866():
    f = (e*x)**(sympy.S(5)/2)*sqrt(c - d*x**2)/(a - b*x**2)
    F = ((Integer(-2) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_867():
    f = (e*x)**(sympy.S(3)/2)*sqrt(c - d*x**2)/(a - b*x**2)
    F = ((Integer(-2) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_868():
    f = sqrt(e*x)*sqrt(c - d*x**2)/(a - b*x**2)
    F = ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_869():
    f = sqrt(c - d*x**2)/(sqrt(e*x)*(a - b*x**2))
    F = ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_870():
    f = sqrt(c - d*x**2)/((e*x)**(sympy.S(3)/2)*(a - b*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_871():
    f = sqrt(c - d*x**2)/((e*x)**(sympy.S(5)/2)*(a - b*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_872():
    f = sqrt(c - d*x**2)/((e*x)**(sympy.S(7)/2)*(a - b*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_873():
    f = (e*x)**(sympy.S(5)/2)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)
    F = ((Integer(-2) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(45) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * ((Symbol('e') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(9) * Symbol('b') * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(21) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(15) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(21) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(15) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_874():
    f = (e*x)**(sympy.S(3)/2)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)
    F = ((Integer(-2) * ((Integer(9) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(21) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(7) * Symbol('b') * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(35) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(21) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(21) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_875():
    f = sqrt(e*x)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)
    F = ((Integer(2) * Symbol('d') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('b') * Symbol('e')))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_876():
    f = (c - d*x**2)**(sympy.S(3)/2)/(sqrt(e*x)*(a - b*x**2))
    F = ((Integer(2) * Symbol('d') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b') * Symbol('e')))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_877():
    f = (c - d*x**2)**(sympy.S(3)/2)/((e*x)**(sympy.S(3)/2)*(a - b*x**2))
    F = ((Integer(-2) * Symbol('c') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_878():
    f = (c - d*x**2)**(sympy.S(3)/2)/((e*x)**(sympy.S(5)/2)*(a - b*x**2))
    F = ((Integer(-2) * Symbol('c') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * Symbol('a') * Symbol('b') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_879():
    f = (c - d*x**2)**(sympy.S(3)/2)/((e*x)**(sympy.S(7)/2)*(a - b*x**2))
    F = ((Integer(-2) * Symbol('c') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_880():
    f = (e*x)**(sympy.S(7)/2)/((a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(2) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_881():
    f = (e*x)**(sympy.S(5)/2)/((a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(-2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_882():
    f = (e*x)**(sympy.S(3)/2)/((a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(-2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_883():
    f = sqrt(e*x)/((a - b*x**2)*sqrt(c - d*x**2))
    F = (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_884():
    f = 1/(sqrt(e*x)*(a - b*x**2)*sqrt(c - d*x**2))
    F = (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_885():
    f = 1/((e*x)**(sympy.S(3)/2)*(a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('a') * Symbol('c') * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_886():
    f = 1/((e*x)**(sympy.S(5)/2)*(a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * Symbol('c') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_887():
    f = 1/((e*x)**(sympy.S(7)/2)*(a - b*x**2)*sqrt(c - d*x**2))
    F = ((Integer(-2) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * Symbol('c') * Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_888():
    f = (e*x)**(sympy.S(9)/2)/((a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('c') * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_889():
    f = (e*x)**(sympy.S(7)/2)/((a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('c') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_890():
    f = (e*x)**(sympy.S(5)/2)/((a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_891():
    f = (e*x)**(sympy.S(3)/2)/((a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_892():
    f = sqrt(e*x)/((a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('d') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('a')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_893():
    f = 1/(sqrt(e*x)*(a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('e') * x))) * ((Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_894():
    f = 1/((e*x)**(sympy.S(3)/2)*(a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * (Symbol('d') * ((Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_895():
    f = 1/((e*x)**(sympy.S(5)/2)*(a - b*x**2)*(c - d*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * (Symbol('d') * ((Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * Symbol('a') * (Symbol('c'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_896():
    f = (e*x)**(sympy.S(7)/2)*sqrt(c - d*x**2)/(a - b*x**2)**2
    F = ((Integer(7) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(8) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(21) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_897():
    f = (e*x)**(sympy.S(5)/2)*sqrt(c - d*x**2)/(a - b*x**2)**2
    F = ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_898():
    f = (e*x)**(sympy.S(3)/2)*sqrt(c - d*x**2)/(a - b*x**2)**2
    F = ((Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_899():
    f = sqrt(e*x)*sqrt(c - d*x**2)/(a - b*x**2)**2
    F = ((((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_900():
    f = sqrt(c - d*x**2)/(sqrt(e*x)*(a - b*x**2)**2)
    F = ((sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_901():
    f = sqrt(c - d*x**2)/((e*x)**(sympy.S(3)/2)*(a - b*x**2)**2)
    F = ((Integer(-5) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + (sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_902():
    f = sqrt(c - d*x**2)/((e*x)**(sympy.S(5)/2)*(a - b*x**2)**2)
    F = ((Integer(-7) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('a') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(7) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_903():
    f = (e*x)**(sympy.S(7)/2)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)**2
    F = ((((Integer(57) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(77) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(42) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(11) * Symbol('d') * Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(14) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * ((Symbol('e') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(48) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(259) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(231) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(42) * (Symbol('b'))**(Integer(4)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_904():
    f = (e*x)**(sympy.S(5)/2)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)**2
    F = ((Integer(-9) * Symbol('d') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(10) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(15) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(10) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(15) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(10) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_905():
    f = (e*x)**(sympy.S(3)/2)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)**2
    F = ((Integer(-7) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * sympy.sqrt((Symbol('e') * x)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(17) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(21) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_906():
    f = sqrt(e*x)*(c - d*x**2)**(sympy.S(3)/2)/(a - b*x**2)**2
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_907():
    f = (c - d*x**2)**(sympy.S(3)/2)/(sqrt(e*x)*(a - b*x**2)**2)
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_908():
    f = (c - d*x**2)**(sympy.S(3)/2)/((e*x)**(sympy.S(3)/2)*(a - b*x**2)**2)
    F = ((Integer(-1) * (((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_909():
    f = (c - d*x**2)**(sympy.S(3)/2)/((e*x)**(sympy.S(5)/2)*(a - b*x**2)**2)
    F = ((Integer(-1) * (((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_910():
    f = (e*x)**(sympy.S(9)/2)/((a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('a') * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_911():
    f = (e*x)**(sympy.S(7)/2)/((a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('a') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_912():
    f = (e*x)**(sympy.S(5)/2)/((a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_913():
    f = (e*x)**(sympy.S(3)/2)/((a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_914():
    f = sqrt(e*x)/((a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('b') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_915():
    f = 1/(sqrt(e*x)*(a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Symbol('b') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_916():
    f = 1/((e*x)**(sympy.S(3)/2)*(a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Integer(-1) * (((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_917():
    f = 1/((e*x)**(sympy.S(5)/2)*(a - b*x**2)**2*sqrt(c - d*x**2))
    F = ((Integer(-1) * (((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_918():
    f = (e*x)**(sympy.S(9)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_919():
    f = (e*x)**(sympy.S(7)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_920():
    f = (e*x)**(sympy.S(5)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Integer(3) * Symbol('d') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_921():
    f = (e*x)**(sympy.S(3)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Integer(3) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(5) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(5) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_922():
    f = sqrt(e*x)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_923():
    f = 1/(sqrt(e*x)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_924():
    f = 1/((e*x)**(sympy.S(3)/2)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Symbol('b') * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_925():
    f = 1/((e*x)**(sympy.S(5)/2)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(3)/2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Symbol('b') * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(7) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_926():
    f = (e*x)**(sympy.S(9)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(6) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(4) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(4) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(4) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_927():
    f = (e*x)**(sympy.S(7)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Integer(5) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_928():
    f = (e*x)**(sympy.S(5)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Integer(5) * Symbol('d') * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(6) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * ((Integer(4) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(7) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(7) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_929():
    f = (e*x)**(sympy.S(3)/2)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Integer(5) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * ((Integer(14) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * Symbol('e') * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(14) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(9) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(9) * Symbol('a') * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_930():
    f = sqrt(e*x)/((a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(6) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(5) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(5) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(5) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_931():
    f = 1/(sqrt(e*x)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('e') * x))) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(17) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Symbol('e') * x))) * ((Integer(6) * Symbol('a') * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(17) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * Symbol('a') * (Symbol('c'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_932():
    f = 1/((e*x)**(sympy.S(3)/2)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d')))) * ((Integer(6) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Symbol('b') * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(19) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))))) * ((Integer(6) * Symbol('a') * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * sympy.sqrt((Symbol('e') * x)) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(19) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * sympy.sqrt((Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(19) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(19) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_933():
    f = 1/((e*x)**(sympy.S(5)/2)*(a - b*x**2)**2*(c - d*x**2)**(sympy.S(5)/2))
    F = ((Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d')))) * ((Integer(6) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Symbol('b') * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2))))) * ((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(7) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(7) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(35) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * ((Symbol('e') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(35) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(11) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(17) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(17) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('e') * x))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_934():
    f = x**5*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = x**2*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(6*b*d) - (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(3*a*d + 5*b*c)/(24*b**2*d**2) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(a**2*d**2 + 2*a*b*c*d + 5*b**2*c**2)/(16*b**2*d**3) - (-a*d + b*c)*(a**2*d**2 + 2*a*b*c*d + 5*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(16*b**(sympy.S(5)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_935():
    f = x**3*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*b*d) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(a*d + 3*b*c)/(8*b*d**2) + (-a*d + b*c)*(a*d + 3*b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(8*b**(sympy.S(3)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_936():
    f = x*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*d) - (-a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(2*sqrt(b)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_937():
    f = sqrt(a + b*x**2)/(x*sqrt(c + d*x**2))
    F = -sqrt(a)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/sqrt(c) + sqrt(b)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_938():
    f = sqrt(a + b*x**2)/(x**3*sqrt(c + d*x**2))
    F = -sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*c*x**2) - (-a*d + b*c)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(2*sqrt(a)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_939():
    f = sqrt(a + b*x**2)/(x**5*sqrt(c + d*x**2))
    F = -(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*a*c*x**4) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(3*a*d + b*c)/(8*a*c**2*x**2) + (-a*d + b*c)*(3*a*d + b*c)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(8*a**(sympy.S(3)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_940():
    f = x**4*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(5*d) + c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-a*d + 4*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + 4*b*c)/(15*b*d**2) - sqrt(c)*sqrt(a + b*x**2)*(-2*a**2*d**2 - 3*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b**2*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-2*a**2*d**2 - 3*a*b*c*d + 8*b**2*c**2)/(15*b**2*d**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_941():
    f = x**2*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = -c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*d) + sqrt(c)*sqrt(a + b*x**2)*(-a*d + 2*b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*b*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(3*b*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_942():
    f = sqrt(a + b*x**2)/(x**2*sqrt(c + d*x**2))
    F = d*x*sqrt(a + b*x**2)/(c*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*sqrt(c + d*x**2)/(c*x) - sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + b*sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_943():
    f = sqrt(a + b*x**2)/(x**4*sqrt(c + d*x**2))
    F = -sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*c*x**3) - b*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a*sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + d*x*sqrt(a + b*x**2)*(-2*a*d + b*c)/(3*a*c**2*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + b*c)/(3*a*c**2*x) - sqrt(d)*sqrt(a + b*x**2)*(-2*a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a*c**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_944():
    f = x**5*(a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = x**2*(a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(8*b*d) - (a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(3*a*d + 7*b*c)/(48*b**2*d**2) + (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)/(192*b**2*d**3) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)/(128*b**2*d**4) + (-a*d + b*c)**2*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(128*b**(sympy.S(5)/2)*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_945():
    f = x**3*(a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = (a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(6*b*d) - (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(a*d + 5*b*c)/(24*b*d**2) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)*(a*d + 5*b*c)/(16*b*d**3) - (-a*d + b*c)**2*(a*d + 5*b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(16*b**(sympy.S(3)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_946():
    f = x*(a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*d) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-3*a*d + 3*b*c)/(8*d**2) + 3*(-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(8*sqrt(b)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_947():
    f = (a + b*x**2)**(sympy.S(3)/2)/(x*sqrt(c + d*x**2))
    F = -a**(sympy.S(3)/2)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/sqrt(c) - sqrt(b)*(-3*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(2*d**(sympy.S(3)/2)) + b*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_948():
    f = (a + b*x**2)**(sympy.S(3)/2)/(x**3*sqrt(c + d*x**2))
    F = -sqrt(a)*(-a*d + 3*b*c)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(2*c**(sympy.S(3)/2)) - a*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*c*x**2) + b**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_949():
    f = (a + b*x**2)**(sympy.S(3)/2)/(x**5*sqrt(c + d*x**2))
    F = -(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*c*x**4) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-3*a*d + 3*b*c)/(8*c**2*x**2) - 3*(-a*d + b*c)**2*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(8*sqrt(a)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_950():
    f = x**4*(a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = b*x**5*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(7*d) - x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-8*a*d + 6*b*c)/(35*d**2) - c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(a**2*d**2 - 11*a*b*c*d + 8*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(35*b*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(a**2*d**2 - 11*a*b*c*d + 8*b**2*c**2)/(35*b*d**3) + 2*sqrt(c)*sqrt(a + b*x**2)*(-a*d + 2*b*c)*(-a**2*d**2 - 4*a*b*c*d + 4*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(35*b**2*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-2*a*d + 4*b*c)*(-a**2*d**2 - 4*a*b*c*d + 4*b**2*c**2)/(35*b**2*d**3*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_951():
    f = x**2*(a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = b*x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(5*d) + 2*c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-3*a*d + 2*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-3*a**2*d/b + 13*a*c - 8*b*c**2/d)/(15*d*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-6*a*d + 4*b*c)/(15*d**2) - sqrt(c)*sqrt(a + b*x**2)*(3*a**2*d**2 - 13*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_952():
    f = (a + b*x**2)**(sympy.S(3)/2)/(x**2*sqrt(c + d*x**2))
    F = -a*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(c*x) + 2*b*sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(a*d + b*c)/(c*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*(a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(c)*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_953():
    f = (a + b*x**2)**(sympy.S(3)/2)/(x**4*sqrt(c + d*x**2))
    F = -a*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*c*x**3) + 2*d*x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(3*c**2*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + 4*b*c)/(3*c**2*x) - 2*sqrt(d)*sqrt(a + b*x**2)*(-a*d + 2*b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*c**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + b*sqrt(a + b*x**2)*(-a*d + 3*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a*sqrt(c)*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_954():
    f = x**5*(a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = x**2*(a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2)/(10*b*d) - (a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2)*(3*a*d + 9*b*c)/(80*b**2*d**2) + (a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(480*b**2*d**3) - (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-a*d + b*c)*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(384*b**2*d**4) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)**2*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(256*b**2*d**5) - (-a*d + b*c)**3*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(256*b**(sympy.S(5)/2)*d**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_955():
    f = x**3*(a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = (a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2)/(8*b*d) - (a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)*(a*d + 7*b*c)/(48*b*d**2) + (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-5*a*d + 5*b*c)*(a*d + 7*b*c)/(192*b*d**3) - 5*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)**2*(a*d + 7*b*c)/(128*b*d**4) + 5*(-a*d + b*c)**3*(a*d + 7*b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(128*b**(sympy.S(3)/2)*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_956():
    f = x*(a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = (a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(6*d) - (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-5*a*d + 5*b*c)/(24*d**2) + 5*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)**2/(16*d**3) - 5*(-a*d + b*c)**3*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(16*sqrt(b)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_957():
    f = (a + b*x**2)**(sympy.S(5)/2)/(x*sqrt(c + d*x**2))
    F = -a**(sympy.S(5)/2)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/sqrt(c) + sqrt(b)*(15*a**2*d**2 - 10*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(8*d**(sympy.S(5)/2)) + b*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*d) - b*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-7*a*d + 3*b*c)/(8*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_958():
    f = (a + b*x**2)**(sympy.S(5)/2)/(x**3*sqrt(c + d*x**2))
    F = -a**(sympy.S(3)/2)*(-a*d + 5*b*c)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(2*c**(sympy.S(3)/2)) - a*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(2*c*x**2) - b**(sympy.S(3)/2)*(-5*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(2*d**(sympy.S(3)/2)) + b*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(a*d + b*c)/(2*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_959():
    f = (a + b*x**2)**(sympy.S(5)/2)/(x**5*sqrt(c + d*x**2))
    F = -sqrt(a)*(3*a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(8*c**(sympy.S(5)/2)) - a*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(4*c*x**4) - a*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-3*a*d + 7*b*c)/(8*c**2*x**2) + b**(sympy.S(5)/2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_960():
    f = x**4*(a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = b*x**5*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(9*d) - 4*b*x**5*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-3*a*d + 2*b*c)/(63*d**2) + x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(75*a**2*d**2 - 115*a*b*c*d + 48*b**2*c**2)/(315*d**3) + c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-5*a**3*d**3 + 105*a**2*b*c*d**2 - 156*a*b**2*c**2*d + 64*b**3*c**3)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(315*b*d**(sympy.S(9)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-5*a**3*d**3 + 105*a**2*b*c*d**2 - 156*a*b**2*c**2*d + 64*b**3*c**3)/(315*b*d**4) - sqrt(c)*sqrt(a + b*x**2)*(-10*a**4*d**4 - 25*a**3*b*c*d**3 + 243*a**2*b**2*c**2*d**2 - 328*a*b**3*c**3*d + 128*b**4*c**4)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(315*b**2*d**(sympy.S(9)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-10*a**4*d**4 - 25*a**3*b*c*d**3 + 243*a**2*b**2*c**2*d**2 - 328*a*b**3*c**3*d + 128*b**4*c**4)/(315*b**2*d**4*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_961():
    f = x**2*(a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = b*x**3*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(7*d) - 2*b*x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-5*a*d + 3*b*c)/(35*d**2) - c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(45*a**2*d**2 - 61*a*b*c*d + 24*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(105*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(45*a**2*d**2 - 61*a*b*c*d + 24*b**2*c**2)/(105*d**3) + sqrt(c)*sqrt(a + b*x**2)*(-15*a**3*d**3 + 103*a**2*b*c*d**2 - 128*a*b**2*c**2*d + 48*b**3*c**3)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(105*b*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-15*a**3*d**3 + 103*a**2*b*c*d**2 - 128*a*b**2*c**2*d + 48*b**3*c**3)/(105*b*d**3*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_962():
    f = (a + b*x**2)**(sympy.S(5)/2)/(x**2*sqrt(c + d*x**2))
    F = -a*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(c*x) - b*sqrt(c)*sqrt(a + b*x**2)*(-9*a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(3*a*d + b*c)/(3*c*d) + x*sqrt(a + b*x**2)*(3*a**2*d/c + 7*a*b - 2*b**2*c/d)/(3*sqrt(c + d*x**2)) + sqrt(a + b*x**2)*(-3*a**2*d**2 - 7*a*b*c*d + 2*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*sqrt(c)*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_963():
    f = (a + b*x**2)**(sympy.S(5)/2)/(x**4*sqrt(c + d*x**2))
    F = -a*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(3*c*x**3) - 2*a*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + 3*b*c)/(3*c**2*x) + b*sqrt(a + b*x**2)*(-a*d + 9*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*sqrt(c)*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-2*a**2*d**2 + 7*a*b*c*d + 3*b**2*c**2)/(3*c**2*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*(-2*a**2*d**2 + 7*a*b*c*d + 3*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*c**(sympy.S(3)/2)*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_964():
    f = x**4*sqrt(3*x**2 - 1)/sqrt(2 - 3*x**2)
    F = -x**3*sqrt(2 - 3*x**2)*sqrt(3*x**2 - 1)/15 - 7*x*sqrt(2 - 3*x**2)*sqrt(3*x**2 - 1)/135 - 8*sqrt(3)*elliptic_e(acos(sqrt(6)*x/2), 2)/135 - 2*sqrt(3)*elliptic_f(acos(sqrt(6)*x/2), 2)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_965():
    f = x**3*sqrt(3*x**2 - 1)/sqrt(2 - 3*x**2)
    F = -sqrt(2 - 3*x**2)*(3*x**2 - 1)**(sympy.S(3)/2)/36 - 7*sqrt(2 - 3*x**2)*sqrt(3*x**2 - 1)/72 + 7*asin(6*x**2 - 3)/144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_966():
    f = x**2*sqrt(3*x**2 - 1)/sqrt(2 - 3*x**2)
    F = -x*sqrt(2 - 3*x**2)*sqrt(3*x**2 - 1)/9 - sqrt(3)*elliptic_e(acos(sqrt(6)*x/2), 2)/9 - sqrt(3)*elliptic_f(acos(sqrt(6)*x/2), 2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_967():
    f = x*sqrt(3*x**2 - 1)/sqrt(2 - 3*x**2)
    F = -sqrt(2 - 3*x**2)*sqrt(3*x**2 - 1)/6 + asin(6*x**2 - 3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_968():
    f = x**2*sqrt(b*x**2 + 2)/sqrt(d*x**2 + 3)
    F = x*sqrt(b*x**2 + 2)*sqrt(d*x**2 + 3)/(3*d) - sqrt(2)*sqrt(b*x**2 + 2)*elliptic_f(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(d**(sympy.S(3)/2)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3)) - x*(6*b - 2*d)*sqrt(b*x**2 + 2)/(3*b*d*sqrt(d*x**2 + 3)) + 2*sqrt(2)*(3*b - d)*sqrt(b*x**2 + 2)*elliptic_e(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(3*b*d**(sympy.S(3)/2)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_969():
    f = x**5/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = x**2*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(4*b*d) - sqrt(a + b*x**2)*sqrt(c + d*x**2)*(3*a*d + 3*b*c)/(8*b**2*d**2) - (4*a*b*c*d - 3*(a*d + b*c)**2)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(8*b**(sympy.S(5)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_970():
    f = x**3/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*b*d) - (a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(2*b**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_971():
    f = x/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_972():
    f = 1/(x*sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_973():
    f = 1/(x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*a*c*x**2) + (a*d + b*c)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_974():
    f = 1/(x**5*sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(a + b*x**2)*sqrt(c + d*x**2)/(4*a*c*x**4) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(3*a*d + 3*b*c)/(8*a**2*c**2*x**2) - (3*a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)*atanh(sqrt(c)*sqrt(a + b*x**2)/(sqrt(a)*sqrt(c + d*x**2)))/(8*a**(sympy.S(5)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_975():
    f = x**6/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = x**3*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(5*b*d) + 4*c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b**2*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(4*a*d + 4*b*c)/(15*b**2*d**2) - sqrt(c)*sqrt(a + b*x**2)*(8*a**2*d**2 + 7*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b**3*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(8*a**2*d**2 + 7*a*b*c*d + 8*b**2*c**2)/(15*b**3*d**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_976():
    f = x**4/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*b*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*b*d) + 2*sqrt(c)*sqrt(a + b*x**2)*(a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*b**2*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(2*a*d + 2*b*c)/(3*b**2*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_977():
    f = x**2/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(c)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)/(b*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_978():
    f = 1/(x**2*sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = d*x*sqrt(a + b*x**2)/(a*c*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*sqrt(c + d*x**2)/(a*c*x) - sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_979():
    f = 1/(x**4*sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = -sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*a*c*x**3) - b*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - 2*d*x*sqrt(a + b*x**2)*(a*d + b*c)/(3*a**2*c**2*sqrt(c + d*x**2)) + sqrt(a + b*x**2)*sqrt(c + d*x**2)*(2*a*d + 2*b*c)/(3*a**2*c**2*x) + 2*sqrt(d)*sqrt(a + b*x**2)*(a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*c**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_980():
    f = x**5/((a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(b**2*sqrt(a + b*x**2)*(-a*d + b*c)) + sqrt(a + b*x**2)*sqrt(c + d*x**2)/(2*b**2*d) - (3*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(2*b**(sympy.S(5)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_981():
    f = x**3/((a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2))
    F = a*sqrt(c + d*x**2)/(b*sqrt(a + b*x**2)*(-a*d + b*c)) + atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(b**(sympy.S(3)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_982():
    f = x/((a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2))
    F = -sqrt(c + d*x**2)/(sqrt(a + b*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_983():
    f = x**5/((a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(3*b**2*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + 2*a*sqrt(c + d*x**2)*(-2*a*d + 3*b*c)/(3*b**2*sqrt(a + b*x**2)*(-a*d + b*c)**2) + atanh(sqrt(d)*sqrt(a + b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(b**(sympy.S(5)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_984():
    f = x**3/((a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2))
    F = a*sqrt(c + d*x**2)/(3*b*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt(c + d*x**2)*(-a*d + 3*b*c)/(3*b*sqrt(a + b*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_985():
    f = x/((a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2))
    F = 2*d*sqrt(c + d*x**2)/(3*sqrt(a + b*x**2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)/((a + b*x**2)**(sympy.S(3)/2)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_986():
    f = x**5/((a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(5*b**2*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + b*c)) + 2*a*sqrt(c + d*x**2)*(-3*a*d + 5*b*c)/(15*b**2*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)*(3*a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)/(15*b**2*sqrt(a + b*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_987():
    f = x**3/((a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2))
    F = a*sqrt(c + d*x**2)/(5*b*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + b*c)) + 2*d*sqrt(c + d*x**2)*(-a*d + 5*b*c)/(15*b*sqrt(a + b*x**2)*(-a*d + b*c)**3) - sqrt(c + d*x**2)*(-a*d + 5*b*c)/(15*b*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_988():
    f = x/((a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2))
    F = -8*d**2*sqrt(c + d*x**2)/(15*sqrt(a + b*x**2)*(-a*d + b*c)**3) + 4*d*sqrt(c + d*x**2)/(15*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - sqrt(c + d*x**2)/((a + b*x**2)**(sympy.S(5)/2)*(-5*a*d + 5*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_989():
    f = x**5/((a + b*x**2)**(sympy.S(9)/2)*sqrt(c + d*x**2))
    F = -a**2*sqrt(c + d*x**2)/(7*b**2*(a + b*x**2)**(sympy.S(7)/2)*(-a*d + b*c)) + 2*a*sqrt(c + d*x**2)*(-4*a*d + 7*b*c)/(35*b**2*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + b*c)**2) + 2*d*sqrt(c + d*x**2)*(3*a**2*d**2 - 14*a*b*c*d + 35*b**2*c**2)/(105*b**2*sqrt(a + b*x**2)*(-a*d + b*c)**4) - sqrt(c + d*x**2)*(3*a**2*d**2 - 14*a*b*c*d + 35*b**2*c**2)/(105*b**2*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_990():
    f = x/(sqrt(a - b*x**2)*sqrt(c + d*x**2))
    F = -atan(sqrt(d)*sqrt(a - b*x**2)/(sqrt(b)*sqrt(c + d*x**2)))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_991():
    f = x/(sqrt(a - b*x**2)*sqrt(c - d*x**2))
    F = -atanh(sqrt(d)*sqrt(a - b*x**2)/(sqrt(b)*sqrt(c - d*x**2)))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_992():
    f = x**2/(sqrt(b*x**2 + 2)*sqrt(d*x**2 + 3))
    F = x*sqrt(b*x**2 + 2)/(b*sqrt(d*x**2 + 3)) - sqrt(2)*sqrt(b*x**2 + 2)*elliptic_e(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(b*sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_993():
    f = x**2/(sqrt(4 - x**2)*sqrt(c + d*x**2))
    F = -c*sqrt(1 + d*x**2/c)*elliptic_f(asin(x/2), -4*d/c)/(d*sqrt(c + d*x**2)) + sqrt(c + d*x**2)*elliptic_e(asin(x/2), -4*d/c)/(d*sqrt(1 + d*x**2/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_994():
    f = x**2/(sqrt(c + d*x**2)*sqrt(x**2 + 4))
    F = x*sqrt(c + d*x**2)/(d*sqrt(x**2 + 4)) - sqrt(c + d*x**2)*elliptic_e(atan(x/2), 1 - 4*d/c)/(d*sqrt((c + d*x**2)/(c*(x**2 + 4)))*sqrt(x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_995():
    f = x**2/(sqrt(1 - x**2)*sqrt(3*x**2 + 2))
    F = sqrt(2)*elliptic_e(asin(x), sympy.S(-3)/2)/3 - sqrt(2)*elliptic_f(asin(x), sympy.S(-3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_996():
    f = x**2/(sqrt(1 - x**2)*sqrt(2 - 3*x**2))
    F = -sqrt(2)*elliptic_e(asin(x), sympy.S(3)/2)/3 + sqrt(2)*elliptic_f(asin(x), sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_997():
    f = x**2/(sqrt(4 - x**2)*sqrt(3*x**2 + 2))
    F = sqrt(2)*elliptic_e(asin(x/2), -6)/3 - sqrt(2)*elliptic_f(asin(x/2), -6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_998():
    f = x**2/(sqrt(2 - 3*x**2)*sqrt(4 - x**2))
    F = -sqrt(2)*elliptic_e(asin(x/2), 6)/3 + sqrt(2)*elliptic_f(asin(x/2), 6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_999():
    f = x**2/(sqrt(1 - 4*x**2)*sqrt(3*x**2 + 2))
    F = sqrt(2)*elliptic_e(asin(2*x), sympy.S(-3)/8)/6 - sqrt(2)*elliptic_f(asin(2*x), sympy.S(-3)/8)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1000():
    f = x**2/(sqrt(1 - 4*x**2)*sqrt(2 - 3*x**2))
    F = -sqrt(2)*elliptic_e(asin(2*x), sympy.S(3)/8)/6 + sqrt(2)*elliptic_f(asin(2*x), sympy.S(3)/8)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1001():
    f = x**2/(sqrt(2 - 3*x**2)*sqrt(x**2 + 1))
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-2)/3)/3 - sqrt(3)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1002():
    f = x**2/(sqrt(2 - 3*x**2)*sqrt(x**2 + 4))
    F = 2*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-1)/6)/3 - 2*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-1)/6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1003():
    f = x**2/(sqrt(2 - 3*x**2)*sqrt(4*x**2 + 1))
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-8)/3)/12 - sqrt(3)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-8)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1004():
    f = x**2/(sqrt(x**2 + 1)*sqrt(3*x**2 + 2))
    F = x*sqrt(3*x**2 + 2)/(3*sqrt(x**2 + 1)) - sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(x), sympy.S(-1)/2)/(3*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1005():
    f = x**2/(sqrt(x**2 + 4)*sqrt(3*x**2 + 2))
    F = x*sqrt(3*x**2 + 2)/(3*sqrt(x**2 + 4)) - sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(x/2), -5)/(3*sqrt((3*x**2 + 2)/(x**2 + 4))*sqrt(x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1006():
    f = x**2/(sqrt(3*x**2 + 2)*sqrt(4*x**2 + 1))
    F = x*sqrt(3*x**2 + 2)/(3*sqrt(4*x**2 + 1)) - sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(2*x), sympy.S(5)/8)/(6*sqrt((3*x**2 + 2)/(4*x**2 + 1))*sqrt(4*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1007():
    f = x**2/(sqrt(1 - x**2)*sqrt(2*x**2 - 1))
    F = -elliptic_e(acos(x), 2)/2 - elliptic_f(acos(x), 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1008():
    f = x**5/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 3*(1 - x**2)**(sympy.S(5)/3)/10 + 3*(1 - x**2)**(sympy.S(2)/3)/2 - 9*2**(sympy.S(1)/3)*log(x**2 + 3)/8 + 27*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 + 9*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1009():
    f = x**3/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -3*(1 - x**2)**(sympy.S(2)/3)/4 + 3*2**(sympy.S(1)/3)*log(x**2 + 3)/8 - 9*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 - 3*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1010():
    f = x/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -2**(sympy.S(1)/3)*log(x**2 + 3)/8 + 3*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1011():
    f = 1/(x*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -log(x)/6 + log(1 - (1 - x**2)**(sympy.S(1)/3))/4 + 2**(sympy.S(1)/3)*log(x**2 + 3)/24 - 2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/6 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1012():
    f = 1/(x**3*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -2**(sympy.S(1)/3)*log(x**2 + 3)/72 + 2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/24 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/36 - (1 - x**2)**(sympy.S(2)/3)/(6*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1013():
    f = 1/(x**5*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -log(x)/27 + log(1 - (1 - x**2)**(sympy.S(1)/3))/18 + 2**(sympy.S(1)/3)*log(x**2 + 3)/216 - 2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/72 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/27 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/108 - (1 - x**2)**(sympy.S(2)/3)/(18*x**2) - (1 - x**2)**(sympy.S(2)/3)/(12*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1014():
    f = x**4/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -3*x*(1 - x**2)**(sympy.S(2)/3)/7 + 54*x/(-7*(1 - x**2)**(sympy.S(1)/3) - 7*sqrt(3) + 7) + 3*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/4 + 3*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/4 - 3*2**(sympy.S(1)/3)*atanh(x)/4 + 9*2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4 + 27*3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(7*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - 18*sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(7*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1015():
    f = x**2/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -3*x/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1) - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/4 + 2**(sympy.S(1)/3)*atanh(x)/4 - 3*2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4 - 3*3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1016():
    f = 1/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/12 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/12 - 2**(sympy.S(1)/3)*atanh(x)/12 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1017():
    f = 1/(x**2*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = x/(-3*(1 - x**2)**(sympy.S(1)/3) - 3*sqrt(3) + 3) - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/36 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/36 + 2**(sympy.S(1)/3)*atanh(x)/36 - 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/12 + 3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(6*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - (1 - x**2)**(sympy.S(2)/3)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1018():
    f = 1/(x**4*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2*x/(-27*(1 - x**2)**(sympy.S(1)/3) - 27*sqrt(3) + 27) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/108 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/108 - 2**(sympy.S(1)/3)*atanh(x)/108 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/36 + 3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(27*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(81*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - 2*(1 - x**2)**(sympy.S(2)/3)/(27*x) - (1 - x**2)**(sympy.S(2)/3)/(9*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1019():
    f = x**7/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -3*x**4*(1 - x**2)**(sympy.S(2)/3)/(10*x**2 + 30) + 9*(1 - x**2)**(sympy.S(2)/3)*(14*x**2 + 69)/(40*x**2 + 120) - 99*2**(sympy.S(1)/3)*log(x**2 + 3)/32 + 297*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/32 + 99*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1020():
    f = x**5/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -3*(1 - x**2)**(sympy.S(2)/3)/4 - 9*(1 - x**2)**(sympy.S(2)/3)/(8*x**2 + 24) + 21*2**(sympy.S(1)/3)*log(x**2 + 3)/32 - 63*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/32 - 21*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1021():
    f = x**3/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = 3*(1 - x**2)**(sympy.S(2)/3)/(8*x**2 + 24) - 3*2**(sympy.S(1)/3)*log(x**2 + 3)/32 + 9*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/32 + 3*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1022():
    f = x/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -(1 - x**2)**(sympy.S(2)/3)/(8*x**2 + 24) - 2**(sympy.S(1)/3)*log(x**2 + 3)/96 + 2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/32 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/48
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1023():
    f = 1/(x*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = (1 - x**2)**(sympy.S(2)/3)/(24*x**2 + 72) - log(x)/18 + log(1 - (1 - x**2)**(sympy.S(1)/3))/12 + 5*2**(sympy.S(1)/3)*log(x**2 + 3)/288 - 5*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/96 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/18 - 5*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1024():
    f = 1/(x**3*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -5*(1 - x**2)**(sympy.S(2)/3)/(72*x**2 + 216) + log(x)/54 - log(1 - (1 - x**2)**(sympy.S(1)/3))/36 - 2**(sympy.S(1)/3)*log(x**2 + 3)/96 + 2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/32 - sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/54 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/48 - (1 - x**2)**(sympy.S(2)/3)/(6*x**2*(x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1025():
    f = 1/(x**5*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = (1 - x**2)**(sympy.S(2)/3)/(216*x**2 + 648) - log(x)/54 + log(1 - (1 - x**2)**(sympy.S(1)/3))/36 + 13*2**(sympy.S(1)/3)*log(x**2 + 3)/2592 - 13*2**(sympy.S(1)/3)*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/864 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/54 - 13*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/1296 - (1 - x**2)**(sympy.S(2)/3)/(36*x**2*(x**2 + 3)) - (1 - x**2)**(sympy.S(2)/3)/(12*x**4*(x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1026():
    f = x**4/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = 3*x*(1 - x**2)**(sympy.S(2)/3)/(8*x**2 + 24) - 27*x/(-8*(1 - x**2)**(sympy.S(1)/3) - 8*sqrt(3) + 8) - 5*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/16 - 5*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/16 + 5*2**(sympy.S(1)/3)*atanh(x)/16 - 15*2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/16 - 27*3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(16*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 9*sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(8*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1027():
    f = x**2/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -x*(1 - x**2)**(sympy.S(2)/3)/(8*x**2 + 24) + x/(-8*(1 - x**2)**(sympy.S(1)/3) - 8*sqrt(3) + 8) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/48 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/48 - 2**(sympy.S(1)/3)*atanh(x)/48 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/16 + 3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(16*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(24*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1028():
    f = 1/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = x*(1 - x**2)**(sympy.S(2)/3)/(24*x**2 + 72) - x/(-24*(1 - x**2)**(sympy.S(1)/3) - 24*sqrt(3) + 24) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/48 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/48 - 2**(sympy.S(1)/3)*atanh(x)/48 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/16 - 3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(48*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(72*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1029():
    f = 1/(x**2*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = x/(-8*(1 - x**2)**(sympy.S(1)/3) - 8*sqrt(3) + 8) - 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/432 - 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/432 + 7*2**(sympy.S(1)/3)*atanh(x)/432 - 7*2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/144 + 3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(16*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(24*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) - (1 - x**2)**(sympy.S(2)/3)/(8*x) + (1 - x**2)**(sympy.S(2)/3)/(24*x*(x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1030():
    f = 1/(x**4*(1 - x**2)**(sympy.S(1)/3)*(x**2 + 3)**2)
    F = -11*x/(-648*(1 - x**2)**(sympy.S(1)/3) - 648*sqrt(3) + 648) + 11*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/1296 + 11*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/1296 - 11*2**(sympy.S(1)/3)*atanh(x)/1296 + 11*2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/432 - 11*3**(sympy.S(1)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(1 - (1 - x**2)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_e(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(1296*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt(((1 - x**2)**(sympy.S(2)/3) + (1 - x**2)**(sympy.S(1)/3) + 1)/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(11 - 11*(1 - x**2)**(sympy.S(1)/3))*elliptic_f(asin((-(1 - x**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(1944*x*sqrt(-(1 - (1 - x**2)**(sympy.S(1)/3))/(-(1 - x**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)) + 11*(1 - x**2)**(sympy.S(2)/3)/(648*x) - 11*(1 - x**2)**(sympy.S(2)/3)/(216*x**3) + (1 - x**2)**(sympy.S(2)/3)/(24*x**3*(x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1031():
    f = x**7/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2*(2 - 3*x**2)**(sympy.S(11)/4)/891 - 16*(2 - 3*x**2)**(sympy.S(7)/4)/567 + 56*(2 - 3*x**2)**(sympy.S(3)/4)/243 + 32*2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/81 + 32*2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1032():
    f = x**5/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = -2*(2 - 3*x**2)**(sympy.S(7)/4)/189 + 4*(2 - 3*x**2)**(sympy.S(3)/4)/27 + 8*2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/27 + 8*2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1033():
    f = x**3/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2*(2 - 3*x**2)**(sympy.S(3)/4)/27 + 2*2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/9 + 2*2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1034():
    f = x/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/6 + 2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1035():
    f = 1/(x*(2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/8 + 2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/8 - 2**(sympy.S(3)/4)*atanh(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/8 + 2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1036():
    f = 1/(x**3*(2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 9*2**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/64 + 3*2**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*(-sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/32 - 9*2**(sympy.S(3)/4)*atanh(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/64 + 3*2**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*(sqrt(2 - 3*x**2) + sqrt(2))/(2*(2 - 3*x**2)**(sympy.S(1)/4)))/32 - (2 - 3*x**2)**(sympy.S(3)/4)/(16*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1037():
    f = x**4/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2*x*(2 - 3*x**2)**(sympy.S(3)/4)/45 + 4*2**(sympy.S(1)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/27 + 4*2**(sympy.S(1)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/27 - 16*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/45
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1038():
    f = x**2/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(1)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/9 + 2**(sympy.S(1)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/9 - 2*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1039():
    f = 1/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(1)/4)*sqrt(3)*atan(2**(sympy.S(3)/4)*sqrt(3)*(-sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/12 + 2**(sympy.S(1)/4)*sqrt(3)*atanh(2**(sympy.S(3)/4)*sqrt(3)*(sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1040():
    f = 1/(x**2*(2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(1)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/16 + 2**(sympy.S(1)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/16 - 2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/8 - (2 - 3*x**2)**(sympy.S(3)/4)/(8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1041():
    f = 1/(x**4*(2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 3*2**(sympy.S(1)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/64 + 3*2**(sympy.S(1)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/64 - 3*2**(sympy.S(1)/4)*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2)/2, 2)/16 - 3*(2 - 3*x**2)**(sympy.S(3)/4)/(16*x) - (2 - 3*x**2)**(sympy.S(3)/4)/(24*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1042():
    f = x**7/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(11)/4)/891 + 8*(3*x**2 - 1)**(sympy.S(7)/4)/567 + 14*(3*x**2 - 1)**(sympy.S(3)/4)/243 + 8*atan((3*x**2 - 1)**(sympy.S(1)/4))/81 - 8*atanh((3*x**2 - 1)**(sympy.S(1)/4))/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1043():
    f = x**5/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(7)/4)/189 + 2*(3*x**2 - 1)**(sympy.S(3)/4)/27 + 4*atan((3*x**2 - 1)**(sympy.S(1)/4))/27 - 4*atanh((3*x**2 - 1)**(sympy.S(1)/4))/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1044():
    f = x**3/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(3)/4)/27 + 2*atan((3*x**2 - 1)**(sympy.S(1)/4))/9 - 2*atanh((3*x**2 - 1)**(sympy.S(1)/4))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1045():
    f = x/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = atan((3*x**2 - 1)**(sympy.S(1)/4))/3 - atanh((3*x**2 - 1)**(sympy.S(1)/4))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1046():
    f = 1/(x*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(2)*log(-sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/8 + sqrt(2)*log(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/8 + atan((3*x**2 - 1)**(sympy.S(1)/4))/2 - sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) - 1)/4 - sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + 1)/4 - atanh((3*x**2 - 1)**(sympy.S(1)/4))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1047():
    f = 1/(x**3*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = -9*sqrt(2)*log(-sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/32 + 9*sqrt(2)*log(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/32 + 3*atan((3*x**2 - 1)**(sympy.S(1)/4))/4 - 9*sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) - 1)/16 - 9*sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + 1)/16 - 3*atanh((3*x**2 - 1)**(sympy.S(1)/4))/4 - (3*x**2 - 1)**(sympy.S(3)/4)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1048():
    f = x**4/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 2*x*(3*x**2 - 1)**(sympy.S(3)/4)/45 + 8*x*(3*x**2 - 1)**(sympy.S(1)/4)/(15*sqrt(3*x**2 - 1) + 15) - sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/27 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/27 - 8*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_e(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(45*x) + 4*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(45*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1049():
    f = x**2/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 2*x*(3*x**2 - 1)**(sympy.S(1)/4)/(3*sqrt(3*x**2 - 1) + 3) - sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18 - 2*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_e(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(9*x) + sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1050():
    f = 1/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1051():
    f = 1/(x**2*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 3*x*(3*x**2 - 1)**(sympy.S(1)/4)/(2*sqrt(3*x**2 - 1) + 2) - sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/8 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/8 - sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_e(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(2*x) + sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(4*x) - (3*x**2 - 1)**(sympy.S(3)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1052():
    f = 1/(x**4*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = 9*x*(3*x**2 - 1)**(sympy.S(1)/4)/(2*sqrt(3*x**2 - 1) + 2) - 3*sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/16 - 3*sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/16 - 3*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_e(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(2*x) + 3*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(4*x) - 3*(3*x**2 - 1)**(sympy.S(3)/4)/(2*x) - (3*x**2 - 1)**(sympy.S(3)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1053():
    f = x**2/((3*x**2 + 2)**(sympy.S(3)/4)*(3*x**2 + 4))
    F = -2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(2*2**(sympy.S(1)/4)*sqrt(3*x**2 + 2) + 2*2**(sympy.S(3)/4))/(6*x*(3*x**2 + 2)**(sympy.S(1)/4)))/18 + 2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(-2*2**(sympy.S(1)/4)*sqrt(3*x**2 + 2) + 2*2**(sympy.S(3)/4))/(6*x*(3*x**2 + 2)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1054():
    f = x**2/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*sqrt(3)*atan(2**(sympy.S(3)/4)*sqrt(3)*(-sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/18 - 2**(sympy.S(3)/4)*sqrt(3)*atanh(2**(sympy.S(3)/4)*sqrt(3)*(sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1055():
    f = x**2/((b*x**2 + 2)**(sympy.S(3)/4)*(b*x**2 + 4))
    F = -2**(sympy.S(3)/4)*atan((2*2**(sympy.S(1)/4)*sqrt(b*x**2 + 2) + 2*2**(sympy.S(3)/4))/(2*sqrt(b)*x*(b*x**2 + 2)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2)) + 2**(sympy.S(3)/4)*atanh((-2*2**(sympy.S(1)/4)*sqrt(b*x**2 + 2) + 2*2**(sympy.S(3)/4))/(2*sqrt(b)*x*(b*x**2 + 2)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1056():
    f = x**2/((-b*x**2 + 2)**(sympy.S(3)/4)*(-b*x**2 + 4))
    F = 2**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*(-sqrt(2)*sqrt(-b*x**2 + 2) + 2)/(2*sqrt(b)*x*(-b*x**2 + 2)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2)) - 2**(sympy.S(3)/4)*atanh(2**(sympy.S(3)/4)*(sqrt(2)*sqrt(-b*x**2 + 2) + 2)/(2*sqrt(b)*x*(-b*x**2 + 2)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1057():
    f = x**2/((a + 3*x**2)**(sympy.S(3)/4)*(2*a + 3*x**2))
    F = -sqrt(3)*atan(sqrt(3)*a**(sympy.S(3)/4)*(1 + sqrt(a + 3*x**2)/sqrt(a))/(3*x*(a + 3*x**2)**(sympy.S(1)/4)))/(9*a**(sympy.S(1)/4)) + sqrt(3)*atanh(sqrt(3)*a**(sympy.S(3)/4)*(1 - sqrt(a + 3*x**2)/sqrt(a))/(3*x*(a + 3*x**2)**(sympy.S(1)/4)))/(9*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1058():
    f = x**2/((a - 3*x**2)**(sympy.S(3)/4)*(2*a - 3*x**2))
    F = sqrt(3)*atan(sqrt(3)*a**(sympy.S(3)/4)*(1 - sqrt(a - 3*x**2)/sqrt(a))/(3*x*(a - 3*x**2)**(sympy.S(1)/4)))/(9*a**(sympy.S(1)/4)) - sqrt(3)*atanh(sqrt(3)*a**(sympy.S(3)/4)*(1 + sqrt(a - 3*x**2)/sqrt(a))/(3*x*(a - 3*x**2)**(sympy.S(1)/4)))/(9*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1059():
    f = x**2/((a + b*x**2)**(sympy.S(3)/4)*(2*a + b*x**2))
    F = -atan(a**(sympy.S(3)/4)*(1 + sqrt(a + b*x**2)/sqrt(a))/(sqrt(b)*x*(a + b*x**2)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/2)) + atanh(a**(sympy.S(3)/4)*(1 - sqrt(a + b*x**2)/sqrt(a))/(sqrt(b)*x*(a + b*x**2)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1060():
    f = x**2/((a - b*x**2)**(sympy.S(3)/4)*(2*a - b*x**2))
    F = atan(a**(sympy.S(3)/4)*(1 - sqrt(a - b*x**2)/sqrt(a))/(sqrt(b)*x*(a - b*x**2)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/2)) - atanh(a**(sympy.S(3)/4)*(1 + sqrt(a - b*x**2)/sqrt(a))/(sqrt(b)*x*(a - b*x**2)**(sympy.S(1)/4)))/(a**(sympy.S(1)/4)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1061():
    f = x**7/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2*(2 - 3*x**2)**(sympy.S(9)/4)/729 - 16*(2 - 3*x**2)**(sympy.S(5)/4)/405 + 56*(2 - 3*x**2)**(sympy.S(1)/4)/81 + 8*2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/81 - 8*2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/81 - 16*2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/81 - 16*2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1062():
    f = x**5/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = -2*(2 - 3*x**2)**(sympy.S(5)/4)/135 + 4*(2 - 3*x**2)**(sympy.S(1)/4)/9 + 2*2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/27 - 2*2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/27 - 4*2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/27 - 4*2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1063():
    f = x**3/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2*(2 - 3*x**2)**(sympy.S(1)/4)/9 + 2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/18 - 2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/18 - 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/9 - 2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1064():
    f = x/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/24 - 2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/24 - 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/12 - 2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1065():
    f = 1/(x*(2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/32 - 2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/32 - 2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/8 - 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/16 - 2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/16 - 2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1066():
    f = 1/(x**3*(2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 3*2**(sympy.S(3)/4)*log(-2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/128 - 3*2**(sympy.S(3)/4)*log(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4) + sqrt(2 - 3*x**2) + sqrt(2))/128 - 15*2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/64 - 3*2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*(2 - 3*x**2)**(sympy.S(1)/4) - 1)/64 - 3*2**(sympy.S(3)/4)*atan((4 - 6*x**2)**(sympy.S(1)/4) + 1)/64 - 15*2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*(2 - 3*x**2)**(sympy.S(1)/4)/2)/64 - (2 - 3*x**2)**(sympy.S(1)/4)/(16*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1067():
    f = x**6/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2*x**3*(2 - 3*x**2)**(sympy.S(1)/4)/63 + 80*x*(2 - 3*x**2)**(sympy.S(1)/4)/567 + 8*2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/81 - 8*2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/81 - 160*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/1701
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1068():
    f = x**4/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2*x*(2 - 3*x**2)**(sympy.S(1)/4)/27 + 2*2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/27 - 2*2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/27 - 4*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1069():
    f = x**2/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*sqrt(3)*atan(2**(sympy.S(3)/4)*sqrt(3)*(-sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/18 - 2**(sympy.S(3)/4)*sqrt(3)*atanh(2**(sympy.S(3)/4)*sqrt(3)*(sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1070():
    f = 1/((2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/24 - 2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/24 + 2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1071():
    f = 1/(x**2*(2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/32 - 2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/32 + 2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/8 - (2 - 3*x**2)**(sympy.S(1)/4)/(8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1072():
    f = 1/(x**4*(2 - 3*x**2)**(sympy.S(3)/4)*(4 - 3*x**2))
    F = 3*2**(sympy.S(3)/4)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/128 - 3*2**(sympy.S(3)/4)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/4)*sqrt(2 - 3*x**2) + 2**(sympy.S(3)/4))/(3*x*(2 - 3*x**2)**(sympy.S(1)/4)))/128 + 11*2**(sympy.S(3)/4)*sqrt(3)*elliptic_f(asin(sqrt(6)*x/2)/2, 2)/64 - (2 - 3*x**2)**(sympy.S(1)/4)/(4*x) - (2 - 3*x**2)**(sympy.S(1)/4)/(24*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1073():
    f = x**2/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1074():
    f = x**2/((-3*x**2 - 2)*(-3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*(-3*x**2 - 1)**(sympy.S(1)/4)))/18 - sqrt(6)*atanh(sqrt(6)*x/(2*(-3*x**2 - 1)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1075():
    f = x**2/((b*x**2 - 2)*(b*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*(b*x**2 - 1)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*(b*x**2 - 1)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1076():
    f = x**2/((-b*x**2 - 2)*(-b*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*(-b*x**2 - 1)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*(-b*x**2 - 1)**(sympy.S(1)/4)))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1077():
    f = x**2/((-2*a + 3*x**2)*(-a + 3*x**2)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a + 3*x**2)**(sympy.S(1)/4)))/(18*a**(sympy.S(1)/4)) - sqrt(6)*atanh(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a + 3*x**2)**(sympy.S(1)/4)))/(18*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1078():
    f = x**2/((-2*a - 3*x**2)*(-a - 3*x**2)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a - 3*x**2)**(sympy.S(1)/4)))/(18*a**(sympy.S(1)/4)) - sqrt(6)*atanh(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a - 3*x**2)**(sympy.S(1)/4)))/(18*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1079():
    f = x**2/((-2*a + b*x**2)*(-a + b*x**2)**(sympy.S(3)/4))
    F = sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a + b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a + b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1080():
    f = x**2/((-2*a - b*x**2)*(-a - b*x**2)**(sympy.S(3)/4))
    F = sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a - b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a - b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1081():
    f = x**7/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(9)/4)/729 + 8*(3*x**2 - 1)**(sympy.S(5)/4)/405 + 14*(3*x**2 - 1)**(sympy.S(1)/4)/81 - 8*atan((3*x**2 - 1)**(sympy.S(1)/4))/81 - 8*atanh((3*x**2 - 1)**(sympy.S(1)/4))/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1082():
    f = x**5/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(5)/4)/135 + 2*(3*x**2 - 1)**(sympy.S(1)/4)/9 - 4*atan((3*x**2 - 1)**(sympy.S(1)/4))/27 - 4*atanh((3*x**2 - 1)**(sympy.S(1)/4))/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1083():
    f = x**3/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 2*(3*x**2 - 1)**(sympy.S(1)/4)/9 - 2*atan((3*x**2 - 1)**(sympy.S(1)/4))/9 - 2*atanh((3*x**2 - 1)**(sympy.S(1)/4))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1084():
    f = x/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = -atan((3*x**2 - 1)**(sympy.S(1)/4))/3 - atanh((3*x**2 - 1)**(sympy.S(1)/4))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1085():
    f = 1/(x*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(2)*log(-sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/8 - sqrt(2)*log(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/8 - atan((3*x**2 - 1)**(sympy.S(1)/4))/2 - sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) - 1)/4 - sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + 1)/4 - atanh((3*x**2 - 1)**(sympy.S(1)/4))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1086():
    f = 1/(x**3*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 15*sqrt(2)*log(-sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/32 - 15*sqrt(2)*log(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + sqrt(3*x**2 - 1) + 1)/32 - 3*atan((3*x**2 - 1)**(sympy.S(1)/4))/4 - 15*sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) - 1)/16 - 15*sqrt(2)*atan(sqrt(2)*(3*x**2 - 1)**(sympy.S(1)/4) + 1)/16 - 3*atanh((3*x**2 - 1)**(sympy.S(1)/4))/4 - (3*x**2 - 1)**(sympy.S(1)/4)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1087():
    f = x**6/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 2*x**3*(3*x**2 - 1)**(sympy.S(1)/4)/63 + 40*x*(3*x**2 - 1)**(sympy.S(1)/4)/567 + 2*sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/81 - 2*sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/81 + 40*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(1701*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1088():
    f = x**4/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 2*x*(3*x**2 - 1)**(sympy.S(1)/4)/27 + sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/27 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/27 + 2*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(81*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1089():
    f = x**2/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1090():
    f = 1/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12 - sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1091():
    f = 1/(x**2*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/8 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/8 - sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(2*x) - (3*x**2 - 1)**(sympy.S(1)/4)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1092():
    f = 1/(x**4*(3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(3)/4))
    F = 3*sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/16 - 3*sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/16 - 11*sqrt(3)*sqrt(x**2/(sqrt(3*x**2 - 1) + 1)**2)*(sqrt(3*x**2 - 1) + 1)*elliptic_f(2*atan((3*x**2 - 1)**(sympy.S(1)/4)), sympy.S.Half)/(8*x) - 2*(3*x**2 - 1)**(sympy.S(1)/4)/x - (3*x**2 - 1)**(sympy.S(1)/4)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1093():
    f = (e*x)**(sympy.S(5)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(3)/4)
    F = 3*a*e**(sympy.S(5)/2)*(-7*a*d + 8*b*c)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(32*b**(sympy.S(11)/4)) - 3*a*e**(sympy.S(5)/2)*(-7*a*d + 8*b*c)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(32*b**(sympy.S(11)/4)) + d*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)/(4*b*e) + e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)*(-7*a*d + 8*b*c)/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1094():
    f = sqrt(e*x)*(c + d*x**2)/(a + b*x**2)**(sympy.S(3)/4)
    F = d*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)/(2*b*e) - sqrt(e)*(-3*a*d + 4*b*c)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(7)/4)) + sqrt(e)*(-3*a*d + 4*b*c)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1095():
    f = (c + d*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -d*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(3)/4)*e**(sympy.S(3)/2)) + d*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(3)/4)*e**(sympy.S(3)/2)) - 2*c*(a + b*x**2)**(sympy.S(1)/4)/(a*e*sqrt(e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1096():
    f = (c + d*x**2)/((e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(5*a*e*(e*x)**(sympy.S(5)/2)) + (a + b*x**2)**(sympy.S(1)/4)*(-10*a*d + 8*b*c)/(5*a**2*e**3*sqrt(e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1097():
    f = (c + d*x**2)/((e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(9*a*e*(e*x)**(sympy.S(9)/2)) + (a + b*x**2)**(sympy.S(1)/4)*(-18*a*d + 16*b*c)/(9*a**2*e**3*(e*x)**(sympy.S(5)/2)) - (a + b*x**2)**(sympy.S(5)/4)*(-72*a*d + 64*b*c)/(45*a**3*e**3*(e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1098():
    f = (c + d*x**2)/((e*x)**(sympy.S(15)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(13*a*e*(e*x)**(sympy.S(13)/2)) + (a + b*x**2)**(sympy.S(1)/4)*(-26*a*d + 24*b*c)/(13*a**2*e**3*(e*x)**(sympy.S(9)/2)) - (a + b*x**2)**(sympy.S(5)/4)*(-208*a*d + 192*b*c)/(65*a**3*e**3*(e*x)**(sympy.S(9)/2)) + (a + b*x**2)**(sympy.S(9)/4)*(-832*a*d + 768*b*c)/(585*a**4*e**3*(e*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1099():
    f = (e*x)**(sympy.S(7)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(3)/4)
    F = -a**(sympy.S(3)/2)*e**2*(e*x)**(sympy.S(3)/2)*(-9*a*d + 10*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(12*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) - a*e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)*(-9*a*d + 10*b*c)/(12*b**3) + d*(e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4)/(5*b*e) + e*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)*(-9*a*d + 10*b*c)/(30*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1100():
    f = (e*x)**(sympy.S(3)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(3)/4)
    F = sqrt(a)*(e*x)**(sympy.S(3)/2)*(-5*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(6*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) + d*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)/(3*b*e) + e*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)*(-5*a*d + 6*b*c)/(6*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1101():
    f = (c + d*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/4))
    F = d*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)/(b*e) - (e*x)**(sympy.S(3)/2)*(-a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*sqrt(b)*e**2*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1102():
    f = (c + d*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(3*a*e*(e*x)**(sympy.S(3)/2)) + 2*sqrt(b)*(e*x)**(sympy.S(3)/2)*(-3*a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(3)/2)*e**4*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1103():
    f = (c + d*x**2)/((e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(7*a*e*(e*x)**(sympy.S(7)/2)) + (a + b*x**2)**(sympy.S(1)/4)*(-14*a*d + 12*b*c)/(21*a**2*e**3*(e*x)**(sympy.S(3)/2)) - 4*b**(sympy.S(3)/2)*(e*x)**(sympy.S(3)/2)*(-7*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(5)/2)*e**6*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1104():
    f = (c + d*x**2)/((e*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(3)/4))
    F = -2*c*(a + b*x**2)**(sympy.S(1)/4)/(11*a*e*(e*x)**(sympy.S(11)/2)) + (a + b*x**2)**(sympy.S(1)/4)*(-22*a*d + 20*b*c)/(77*a**2*e**3*(e*x)**(sympy.S(7)/2)) - 4*b*(a + b*x**2)**(sympy.S(1)/4)*(-11*a*d + 10*b*c)/(77*a**3*e**5*(e*x)**(sympy.S(3)/2)) + 8*b**(sympy.S(5)/2)*(e*x)**(sympy.S(3)/2)*(-11*a*d + 10*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(77*a**(sympy.S(7)/2)*e**8*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1105():
    f = (e*x)**(sympy.S(3)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(5)/4)
    F = d*(e*x)**(sympy.S(5)/2)/(2*b*e*(a + b*x**2)**(sympy.S(1)/4)) - e*sqrt(e*x)*(-5*a*d + 4*b*c)/(2*b**2*(a + b*x**2)**(sympy.S(1)/4)) + e**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(9)/4)) + e**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1106():
    f = (c + d*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(5)/4))
    F = d*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(5)/4)*sqrt(e)) + d*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(b**(sympy.S(5)/4)*sqrt(e)) + sqrt(e*x)*(-2*a*d + 2*b*c)/(a*b*e*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1107():
    f = (c + d*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(3*a*e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) - sqrt(e*x)*(-6*a*d + 8*b*c)/(3*a**2*e**3*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1108():
    f = (c + d*x**2)/((e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(7*a*e*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) - (-14*a*d + 16*b*c)/(7*a**2*e**3*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (a + b*x**2)**(sympy.S(3)/4)*(-56*a*d + 64*b*c)/(21*a**3*e**3*(e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1109():
    f = (c + d*x**2)/((e*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(11*a*e*(e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(1)/4)) - (-22*a*d + 24*b*c)/(11*a**2*e**3*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (a + b*x**2)**(sympy.S(3)/4)*(-176*a*d + 192*b*c)/(33*a**3*e**3*(e*x)**(sympy.S(7)/2)) - (a + b*x**2)**(sympy.S(7)/4)*(-704*a*d + 768*b*c)/(231*a**4*e**3*(e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1110():
    f = (e*x)**(sympy.S(9)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(5)/4)
    F = -7*a**(sympy.S(3)/2)*e**4*sqrt(e*x)*(-11*a*d + 10*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(20*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 7*a*e**3*(e*x)**(sympy.S(3)/2)*(-11*a*d + 10*b*c)/(60*b**3*(a + b*x**2)**(sympy.S(1)/4)) + d*(e*x)**(sympy.S(11)/2)/(5*b*e*(a + b*x**2)**(sympy.S(1)/4)) + e*(e*x)**(sympy.S(7)/2)*(-11*a*d + 10*b*c)/(30*b**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1111():
    f = (e*x)**(sympy.S(5)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(5)/4)
    F = sqrt(a)*e**2*sqrt(e*x)*(-7*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(2*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + d*(e*x)**(sympy.S(7)/2)/(3*b*e*(a + b*x**2)**(sympy.S(1)/4)) + e*(e*x)**(sympy.S(3)/2)*(-7*a*d + 6*b*c)/(6*b**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1112():
    f = sqrt(e*x)*(c + d*x**2)/(a + b*x**2)**(sympy.S(5)/4)
    F = d*(e*x)**(sympy.S(3)/2)/(b*e*(a + b*x**2)**(sympy.S(1)/4)) - sqrt(e*x)*(-3*a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(sqrt(a)*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1113():
    f = (c + d*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(a*e*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)) + sqrt(e*x)*(-2*a*d + 4*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(a**(sympy.S(3)/2)*sqrt(b)*e**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1114():
    f = (c + d*x**2)/((e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(5*a*e*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (-10*a*d + 12*b*c)/(5*a**2*e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)) - 4*sqrt(b)*sqrt(e*x)*(-5*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(5)/2)*e**4*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1115():
    f = (c + d*x**2)/((e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(5)/4))
    F = -2*c/(9*a*e*(e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (-18*a*d + 20*b*c)/(45*a**2*e**3*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 4*b*(-9*a*d + 10*b*c)/(15*a**3*e**5*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)) + 8*b**(sympy.S(3)/2)*sqrt(e*x)*(-9*a*d + 10*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(15*a**(sympy.S(7)/2)*e**6*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1116():
    f = (e*x)**(sympy.S(5)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(7)/4)
    F = -e**(sympy.S(5)/2)*(-7*a*d + 4*b*c)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(11)/4)) + e**(sympy.S(5)/2)*(-7*a*d + 4*b*c)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(11)/4)) + (e*x)**(sympy.S(7)/2)*(-2*a*d + 2*b*c)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/4)) - e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)*(-7*a*d + 4*b*c)/(6*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1117():
    f = sqrt(e*x)*(c + d*x**2)/(a + b*x**2)**(sympy.S(7)/4)
    F = -d*sqrt(e)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(7)/4) + d*sqrt(e)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(7)/4) + (e*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1118():
    f = (c + d*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(7)/4))
    F = -2*c/(a*e*sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/4)) - (e*x)**(sympy.S(3)/2)*(-2*a*d + 8*b*c)/(3*a**2*e**3*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1119():
    f = (c + d*x**2)/((e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(7)/4))
    F = -2*c/(5*a*e*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) - (-10*a*d + 16*b*c)/(15*a**2*e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(3)/4)) + (a + b*x**2)**(sympy.S(1)/4)*(-40*a*d + 64*b*c)/(15*a**3*e**3*sqrt(e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1120():
    f = (c + d*x**2)/((e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(7)/4))
    F = -2*c/(9*a*e*(e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(3)/4)) - (-6*a*d + 8*b*c)/(9*a**2*e**3*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) + (a + b*x**2)**(sympy.S(1)/4)*(-48*a*d + 64*b*c)/(9*a**3*e**3*(e*x)**(sympy.S(5)/2)) - (a + b*x**2)**(sympy.S(5)/4)*(-192*a*d + 256*b*c)/(45*a**4*e**3*(e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1121():
    f = (e*x)**(sympy.S(7)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(7)/4)
    F = 5*sqrt(a)*e**2*(e*x)**(sympy.S(3)/2)*(-3*a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(6*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(3)/4)) + e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)*(-15*a*d + 10*b*c)/(6*b**3) + (e*x)**(sympy.S(9)/2)*(-2*a*d + 2*b*c)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/4)) - e*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)*(-3*a*d + 2*b*c)/(3*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1122():
    f = (e*x)**(sympy.S(3)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(7)/4)
    F = (e*x)**(sympy.S(5)/2)*(-2*a*d + 2*b*c)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/4)) - e*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)*(-5*a*d + 2*b*c)/(3*a*b**2) - (e*x)**(sympy.S(3)/2)*(-5*a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*sqrt(a)*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1123():
    f = (c + d*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(7)/4))
    F = sqrt(e*x)*(-2*a*d + 2*b*c)/(3*a*b*e*(a + b*x**2)**(sympy.S(3)/4)) - (e*x)**(sympy.S(3)/2)*(2*a*d + 4*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(3)/2)*sqrt(b)*e**2*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1124():
    f = (c + d*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(7)/4))
    F = -2*c/(3*a*e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) - sqrt(e*x)*(-2*a*d + 4*b*c)/(3*a**2*e**3*(a + b*x**2)**(sympy.S(3)/4)) + 4*sqrt(b)*(e*x)**(sympy.S(3)/2)*(-a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(3*a**(sympy.S(5)/2)*e**4*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1125():
    f = (c + d*x**2)/((e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(7)/4))
    F = -2*c/(7*a*e*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(3)/4)) - (-14*a*d + 20*b*c)/(21*a**2*e**3*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/4)) + (a + b*x**2)**(sympy.S(1)/4)*(-28*a*d + 40*b*c)/(21*a**3*e**3*(e*x)**(sympy.S(3)/2)) - 8*b**(sympy.S(3)/2)*(e*x)**(sympy.S(3)/2)*(-7*a*d + 10*b*c)*(a/(b*x**2) + 1)**(sympy.S(3)/4)*elliptic_f(acot(sqrt(b)*x/sqrt(a))/2, 2)/(21*a**(sympy.S(7)/2)*e**6*(a + b*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1126():
    f = (e*x)**(sympy.S(7)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = -e**3*sqrt(e*x)*(-9*a*d + 4*b*c)/(2*b**3*(a + b*x**2)**(sympy.S(1)/4)) + e**(sympy.S(7)/2)*(-9*a*d + 4*b*c)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(13)/4)) + e**(sympy.S(7)/2)*(-9*a*d + 4*b*c)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/(4*b**(sympy.S(13)/4)) + (e*x)**(sympy.S(9)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) - e*(e*x)**(sympy.S(5)/2)*(-9*a*d + 4*b*c)/(10*a*b**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1127():
    f = (e*x)**(sympy.S(3)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = -2*d*e*sqrt(e*x)/(b**2*(a + b*x**2)**(sympy.S(1)/4)) + d*e**(sympy.S(3)/2)*atan(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(9)/4) + d*e**(sympy.S(3)/2)*atanh(b**(sympy.S(1)/4)*sqrt(e*x)/(sqrt(e)*(a + b*x**2)**(sympy.S(1)/4)))/b**(sympy.S(9)/4) + (e*x)**(sympy.S(5)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1128():
    f = (c + d*x**2)/(sqrt(e*x)*(a + b*x**2)**(sympy.S(9)/4))
    F = sqrt(e*x)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) + sqrt(e*x)*(2*a*d + 8*b*c)/(5*a**2*b*e*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1129():
    f = (c + d*x**2)/((e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(3*a*e*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/4)) - sqrt(e*x)*(-6*a*d + 16*b*c)/(15*a**2*e**3*(a + b*x**2)**(sympy.S(5)/4)) - sqrt(e*x)*(-24*a*d + 64*b*c)/(15*a**3*e**3*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1130():
    f = (c + d*x**2)/((e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(7*a*e*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-14*a*d + 24*b*c)/(35*a**2*e**3*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-112*a*d + 192*b*c)/(35*a**3*e**3*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (a + b*x**2)**(sympy.S(3)/4)*(-448*a*d + 768*b*c)/(105*a**4*e**3*(e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1131():
    f = (c + d*x**2)/((e*x)**(sympy.S(13)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(11*a*e*(e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-22*a*d + 32*b*c)/(55*a**2*e**3*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-264*a*d + 384*b*c)/(55*a**3*e**3*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) + (a + b*x**2)**(sympy.S(3)/4)*(-704*a*d + 1024*b*c)/(55*a**4*e**3*(e*x)**(sympy.S(7)/2)) - (a + b*x**2)**(sympy.S(7)/4)*(-2816*a*d + 4096*b*c)/(385*a**5*e**3*(e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1132():
    f = (e*x)**(sympy.S(13)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = -77*a**(sympy.S(3)/2)*e**6*sqrt(e*x)*(-3*a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(20*b**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 77*a*e**5*(e*x)**(sympy.S(3)/2)*(-3*a*d + 2*b*c)/(60*b**4*(a + b*x**2)**(sympy.S(1)/4)) + e**3*(e*x)**(sympy.S(7)/2)*(-33*a*d + 22*b*c)/(30*b**3*(a + b*x**2)**(sympy.S(1)/4)) + (e*x)**(sympy.S(15)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) - e*(e*x)**(sympy.S(11)/2)*(-3*a*d + 2*b*c)/(5*a*b**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1133():
    f = (e*x)**(sympy.S(9)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = 7*sqrt(a)*e**4*sqrt(e*x)*(-11*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(10*b**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(1)/4)) + e**3*(e*x)**(sympy.S(3)/2)*(-77*a*d + 42*b*c)/(30*b**3*(a + b*x**2)**(sympy.S(1)/4)) + (e*x)**(sympy.S(11)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) - e*(e*x)**(sympy.S(7)/2)*(-11*a*d + 6*b*c)/(15*a*b**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1134():
    f = (e*x)**(sympy.S(5)/2)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = (e*x)**(sympy.S(7)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) - e*(e*x)**(sympy.S(3)/2)*(-7*a*d + 2*b*c)/(5*a*b**2*(a + b*x**2)**(sympy.S(1)/4)) - e**2*sqrt(e*x)*(-21*a*d + 6*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*sqrt(a)*b**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1135():
    f = sqrt(e*x)*(c + d*x**2)/(a + b*x**2)**(sympy.S(9)/4)
    F = (e*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)/(5*a*b*e*(a + b*x**2)**(sympy.S(5)/4)) - sqrt(e*x)*(6*a*d + 4*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1136():
    f = (c + d*x**2)/((e*x)**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(a*e*sqrt(e*x)*(a + b*x**2)**(sympy.S(5)/4)) - (e*x)**(sympy.S(3)/2)*(-2*a*d + 12*b*c)/(5*a**2*e**3*(a + b*x**2)**(sympy.S(5)/4)) + sqrt(e*x)*(-4*a*d + 24*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(5)/2)*sqrt(b)*e**2*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1137():
    f = (c + d*x**2)/((e*x)**(sympy.S(7)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(5*a*e*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-2*a*d + 4*b*c)/(5*a**2*e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(5)/4)) + (-12*a*d + 24*b*c)/(5*a**3*e**3*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)) - 24*sqrt(b)*sqrt(e*x)*(-a*d + 2*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(5*a**(sympy.S(7)/2)*e**4*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1138():
    f = (c + d*x**2)/((e*x)**(sympy.S(11)/2)*(a + b*x**2)**(sympy.S(9)/4))
    F = -2*c/(9*a*e*(e*x)**(sympy.S(9)/2)*(a + b*x**2)**(sympy.S(5)/4)) - (-18*a*d + 28*b*c)/(45*a**2*e**3*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(5)/4)) + (-36*a*d + 56*b*c)/(45*a**3*e**3*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**(sympy.S(1)/4)) - 8*b*(-9*a*d + 14*b*c)/(15*a**4*e**5*sqrt(e*x)*(a + b*x**2)**(sympy.S(1)/4)) + 16*b**(sympy.S(3)/2)*sqrt(e*x)*(-9*a*d + 14*b*c)*(a/(b*x**2) + 1)**(sympy.S(1)/4)*elliptic_e(acot(sqrt(b)*x/sqrt(a))/2, 2)/(15*a**(sympy.S(9)/2)*e**6*(a + b*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1139():
    f = (e*x)**m*(a + b*x**2)**p*(c + d*x**2)**q
    F = (e*x)**(m + 1)*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(m/2 + sympy.S.Half, -p, -q, m/2 + sympy.S(3)/2, -b*x**2/a, -d*x**2/c)/(e*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1140():
    f = x**4*(a + b*x**2)**p*(c + d*x**2)**q
    F = x**5*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(5)/2, -p, -q, sympy.S(7)/2, -b*x**2/a, -d*x**2/c)/(5*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1141():
    f = x**2*(a + b*x**2)**p*(c + d*x**2)**q
    F = x**3*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(3)/2, -p, -q, sympy.S(5)/2, -b*x**2/a, -d*x**2/c)/(3*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1142():
    f = (a + b*x**2)**p*(c + d*x**2)**q
    F = x*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -b*x**2/a, -d*x**2/c)/((1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1143():
    f = (a + b*x**2)**p*(c + d*x**2)**q/x**2
    F = -(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(-1)/2, -p, -q, sympy.S.Half, -b*x**2/a, -d*x**2/c)/(x*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1144():
    f = (a + b*x**2)**p*(c + d*x**2)**q/x**4
    F = -(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(-3)/2, -p, -q, sympy.S(-1)/2, -b*x**2/a, -d*x**2/c)/(3*x**3*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1145():
    f = x**5*(a + b*x**2)**p*(c + d*x**2)**q
    F = x**2*(a + b*x**2)**(p + 1)*(c + d*x**2)**(q + 1)/(2*b*d*(p + q + 3)) - (a + b*x**2)**(p + 1)*(c + d*x**2)**(q + 1)*(a*d*(q + 2) + b*c*(p + 2))/(2*b**2*d**2*(p + q + 2)*(p + q + 3)) + (a + b*x**2)**(p + 1)*(c + d*x**2)**q*(a**2*d**2*(q**2 + 3*q + 2) + 2*a*b*c*d*(p + 1)*(q + 1) + b**2*c**2*(p**2 + 3*p + 2))*hyper((-q, p + 1), (p + 2,), -d*(a + b*x**2)/(-a*d + b*c))/(2*b**3*d**2*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1)*(p + q + 2)*(p + q + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1146():
    f = x**3*(a + b*x**2)**p*(c + d*x**2)**q
    F = (a + b*x**2)**(p + 1)*(c + d*x**2)**(q + 1)/(2*b*d*(p + q + 2)) - (a + b*x**2)**(p + 1)*(c + d*x**2)**q*(a*d*(q + 1) + b*c*(p + 1))*hyper((-q, p + 1), (p + 2,), -d*(a + b*x**2)/(-a*d + b*c))/(2*b**2*d*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1)*(p + q + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1147():
    f = x*(a + b*x**2)**p*(c + d*x**2)**q
    F = (a + b*x**2)**(p + 1)*(c + d*x**2)**q*hyper((-q, p + 1), (p + 2,), -d*(a + b*x**2)/(-a*d + b*c))/(2*b*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1148():
    f = (a + b*x**2)**p*(c + d*x**2)**q/x
    F = -(a + b*x**2)**(p + 1)*(c + d*x**2)**q*appellf1(p + 1, 1, -q, p + 2, (a + b*x**2)/a, -d*(a + b*x**2)/(-a*d + b*c))/(2*a*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1149():
    f = (a + b*x**2)**p*(c + d*x**2)**q/x**3
    F = b*(a + b*x**2)**(p + 1)*(c + d*x**2)**q*appellf1(p + 1, 2, -q, p + 2, (a + b*x**2)/a, -d*(a + b*x**2)/(-a*d + b*c))/(2*a**2*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1150():
    f = (a + b*x**2)**p*(c + d*x**2)**q/x**5
    F = -b**2*(a + b*x**2)**(p + 1)*(c + d*x**2)**q*appellf1(p + 1, 3, -q, p + 2, (a + b*x**2)/a, -d*(a + b*x**2)/(-a*d + b*c))/(2*a**3*(b*(c + d*x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1151():
    f = (e*x)**(sympy.S(5)/2)*(a + b*x**2)**p*(c + d*x**2)**q
    F = 2*(e*x)**(sympy.S(7)/2)*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(7)/4, -p, -q, sympy.S(11)/4, -b*x**2/a, -d*x**2/c)/(7*e*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1152():
    f = (e*x)**(sympy.S(3)/2)*(a + b*x**2)**p*(c + d*x**2)**q
    F = 2*(e*x)**(sympy.S(5)/2)*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(5)/4, -p, -q, sympy.S(9)/4, -b*x**2/a, -d*x**2/c)/(5*e*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1153():
    f = sqrt(e*x)*(a + b*x**2)**p*(c + d*x**2)**q
    F = 2*(e*x)**(sympy.S(3)/2)*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(3)/4, -p, -q, sympy.S(7)/4, -b*x**2/a, -d*x**2/c)/(3*e*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1154():
    f = (a + b*x**2)**p*(c + d*x**2)**q/sqrt(e*x)
    F = 2*sqrt(e*x)*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(1)/4, -p, -q, sympy.S(5)/4, -b*x**2/a, -d*x**2/c)/(e*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1155():
    f = (a + b*x**2)**p*(c + d*x**2)**q/(e*x)**(sympy.S(3)/2)
    F = -2*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(-1)/4, -p, -q, sympy.S(3)/4, -b*x**2/a, -d*x**2/c)/(e*sqrt(e*x)*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_4_e_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1156():
    f = (a + b*x**2)**p*(c + d*x**2)**q/(e*x)**(sympy.S(5)/2)
    F = -2*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S(-3)/4, -p, -q, sympy.S(1)/4, -b*x**2/a, -d*x**2/c)/(3*e*(e*x)**(sympy.S(3)/2)*(1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F

