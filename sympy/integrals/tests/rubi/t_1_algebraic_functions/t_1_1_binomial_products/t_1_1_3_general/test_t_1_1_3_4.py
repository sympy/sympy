"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.4 (e x)^m (a+b x^n)^p (c+d x^n)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, m, n, p, q = symbols('A B a b c d e m n p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_1():
    f = x**2*(A + B*x**3)*(a + b*x**3)
    F = A*a*x**3/3 + B*b*x**9/9 + x**6*(A*b + B*a)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_2():
    f = x*(A + B*x**3)*(a + b*x**3)
    F = A*a*x**2/2 + B*b*x**8/8 + x**5*(A*b + B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_3():
    f = (A + B*x**3)*(a + b*x**3)
    F = A*a*x + B*b*x**7/7 + x**4*(A*b + B*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_4():
    f = (A + B*x**3)*(a + b*x**3)/x
    F = A*a*log(x) + B*b*x**6/6 + x**3*(A*b + B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_5():
    f = (A + B*x**3)*(a + b*x**3)/x**2
    F = -A*a/x + B*b*x**5/5 + x**2*(A*b + B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_6():
    f = (A + B*x**3)*(a + b*x**3)/x**3
    F = -A*a/(2*x**2) + B*b*x**4/4 + x*(A*b + B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_7():
    f = (A + B*x**3)*(a + b*x**3)/x**4
    F = -A*a/(3*x**3) + B*b*x**3/3 + (A*b + B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_8():
    f = (A + B*x**3)*(a + b*x**3)/x**5
    F = -A*a/(4*x**4) + B*b*x**2/2 - (A*b + B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_9():
    f = (A + B*x**3)*(a + b*x**3)/x**6
    F = -A*a/(5*x**5) + B*b*x - (A*b + B*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_10():
    f = (A + B*x**3)*(a + b*x**3)/x**7
    F = -A*a/(6*x**6) + B*b*log(x) - (A*b + B*a)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_11():
    f = x**2*(A + B*x**3)*(a + b*x**3)**2
    F = B*(a + b*x**3)**4/(12*b**2) + (a + b*x**3)**3*(A*b - B*a)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_12():
    f = x*(A + B*x**3)*(a + b*x**3)**2
    F = A*a**2*x**2/2 + B*b**2*x**11/11 + a*x**5*(2*A*b + B*a)/5 + b*x**8*(A*b + 2*B*a)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_13():
    f = (A + B*x**3)*(a + b*x**3)**2
    F = A*a**2*x + B*b**2*x**10/10 + a*x**4*(2*A*b + B*a)/4 + b*x**7*(A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_14():
    f = (A + B*x**3)*(a + b*x**3)**2/x
    F = A*a**2*log(x) + 2*A*a*b*x**3/3 + A*b**2*x**6/6 + B*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_15():
    f = (A + B*x**3)*(a + b*x**3)**2/x**2
    F = -A*a**2/x + B*b**2*x**8/8 + a*x**2*(2*A*b + B*a)/2 + b*x**5*(A*b + 2*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_16():
    f = (A + B*x**3)*(a + b*x**3)**2/x**3
    F = -A*a**2/(2*x**2) + B*b**2*x**7/7 + a*x*(2*A*b + B*a) + b*x**4*(A*b + 2*B*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_17():
    f = (A + B*x**3)*(a + b*x**3)**2/x**4
    F = -A*a**2/(3*x**3) + B*b**2*x**6/6 + a*(2*A*b + B*a)*log(x) + b*x**3*(A*b + 2*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_18():
    f = (A + B*x**3)*(a + b*x**3)**2/x**5
    F = -A*a**2/(4*x**4) + B*b**2*x**5/5 - a*(2*A*b + B*a)/x + b*x**2*(A*b + 2*B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_19():
    f = (A + B*x**3)*(a + b*x**3)**2/x**6
    F = -A*a**2/(5*x**5) + B*b**2*x**4/4 - a*(2*A*b + B*a)/(2*x**2) + b*x*(A*b + 2*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_20():
    f = (A + B*x**3)*(a + b*x**3)**2/x**7
    F = -A*a**2/(6*x**6) + B*b**2*x**3/3 - a*(2*A*b + B*a)/(3*x**3) + b*(A*b + 2*B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_21():
    f = (A + B*x**3)*(a + b*x**3)**2/x**8
    F = -A*a**2/(7*x**7) + B*b**2*x**2/2 - a*(2*A*b + B*a)/(4*x**4) - b*(A*b + 2*B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_22():
    f = (A + B*x**3)*(a + b*x**3)**2/x**9
    F = -A*a**2/(8*x**8) + B*b**2*x - a*(2*A*b + B*a)/(5*x**5) - b*(A*b + 2*B*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_23():
    f = x**9*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**10/10 + B*b**5*x**28/28 + a**4*x**13*(5*A*b + B*a)/13 + 5*a**3*b*x**16*(2*A*b + B*a)/16 + 10*a**2*b**2*x**19*(A*b + B*a)/19 + 5*a*b**3*x**22*(A*b + 2*B*a)/22 + b**4*x**25*(A*b + 5*B*a)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_24():
    f = x**8*(A + B*x**3)*(a + b*x**3)**5
    F = B*(a + b*x**3)**9/(27*b**4) + a**2*(a + b*x**3)**6*(A*b - B*a)/(18*b**4) - a*(a + b*x**3)**7*(2*A*b - 3*B*a)/(21*b**4) + (a + b*x**3)**8*(A*b - 3*B*a)/(24*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_25():
    f = x**7*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**8/8 + B*b**5*x**26/26 + a**4*x**11*(5*A*b + B*a)/11 + 5*a**3*b*x**14*(2*A*b + B*a)/14 + 10*a**2*b**2*x**17*(A*b + B*a)/17 + a*b**3*x**20*(A*b + 2*B*a)/4 + b**4*x**23*(A*b + 5*B*a)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_26():
    f = x**6*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**7/7 + B*b**5*x**25/25 + a**4*x**10*(5*A*b + B*a)/10 + 5*a**3*b*x**13*(2*A*b + B*a)/13 + 5*a**2*b**2*x**16*(A*b + B*a)/8 + 5*a*b**3*x**19*(A*b + 2*B*a)/19 + b**4*x**22*(A*b + 5*B*a)/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_27():
    f = x**5*(A + B*x**3)*(a + b*x**3)**5
    F = B*(a + b*x**3)**8/(24*b**3) - a*(a + b*x**3)**6*(A*b - B*a)/(18*b**3) + (a + b*x**3)**7*(A*b - 2*B*a)/(21*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_28():
    f = x**4*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**5/5 + B*b**5*x**23/23 + a**4*x**8*(5*A*b + B*a)/8 + 5*a**3*b*x**11*(2*A*b + B*a)/11 + 5*a**2*b**2*x**14*(A*b + B*a)/7 + 5*a*b**3*x**17*(A*b + 2*B*a)/17 + b**4*x**20*(A*b + 5*B*a)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_29():
    f = x**3*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**4/4 + B*b**5*x**22/22 + a**4*x**7*(5*A*b + B*a)/7 + a**3*b*x**10*(2*A*b + B*a)/2 + 10*a**2*b**2*x**13*(A*b + B*a)/13 + 5*a*b**3*x**16*(A*b + 2*B*a)/16 + b**4*x**19*(A*b + 5*B*a)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_30():
    f = x**2*(A + B*x**3)*(a + b*x**3)**5
    F = B*(a + b*x**3)**7/(21*b**2) + (a + b*x**3)**6*(A*b - B*a)/(18*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_31():
    f = x*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**2/2 + B*b**5*x**20/20 + a**4*x**5*(5*A*b + B*a)/5 + 5*a**3*b*x**8*(2*A*b + B*a)/8 + 10*a**2*b**2*x**11*(A*b + B*a)/11 + 5*a*b**3*x**14*(A*b + 2*B*a)/14 + b**4*x**17*(A*b + 5*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_32():
    f = (A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x + B*b**5*x**19/19 + a**4*x**4*(5*A*b + B*a)/4 + 5*a**3*b*x**7*(2*A*b + B*a)/7 + a**2*b**2*x**10*(A*b + B*a) + 5*a*b**3*x**13*(A*b + 2*B*a)/13 + b**4*x**16*(A*b + 5*B*a)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_33():
    f = (A + B*x**3)*(a + b*x**3)**5/x
    F = A*a**5*log(x) + 5*A*a**4*b*x**3/3 + 5*A*a**3*b**2*x**6/3 + 10*A*a**2*b**3*x**9/9 + 5*A*a*b**4*x**12/12 + A*b**5*x**15/15 + B*(a + b*x**3)**6/(18*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_34():
    f = (A + B*x**3)*(a + b*x**3)**5/x**2
    F = -A*a**5/x + B*b**5*x**17/17 + a**4*x**2*(5*A*b + B*a)/2 + a**3*b*x**5*(2*A*b + B*a) + 5*a**2*b**2*x**8*(A*b + B*a)/4 + 5*a*b**3*x**11*(A*b + 2*B*a)/11 + b**4*x**14*(A*b + 5*B*a)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_35():
    f = (A + B*x**3)*(a + b*x**3)**5/x**3
    F = -A*a**5/(2*x**2) + B*b**5*x**16/16 + a**4*x*(5*A*b + B*a) + 5*a**3*b*x**4*(2*A*b + B*a)/4 + 10*a**2*b**2*x**7*(A*b + B*a)/7 + a*b**3*x**10*(A*b + 2*B*a)/2 + b**4*x**13*(A*b + 5*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_36():
    f = (A + B*x**3)*(a + b*x**3)**5/x**4
    F = -A*a**5/(3*x**3) + B*b**5*x**15/15 + a**4*(5*A*b + B*a)*log(x) + 5*a**3*b*x**3*(2*A*b + B*a)/3 + 5*a**2*b**2*x**6*(A*b + B*a)/3 + 5*a*b**3*x**9*(A*b + 2*B*a)/9 + b**4*x**12*(A*b + 5*B*a)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_37():
    f = (A + B*x**3)*(a + b*x**3)**5/x**5
    F = -A*a**5/(4*x**4) + B*b**5*x**14/14 - a**4*(5*A*b + B*a)/x + 5*a**3*b*x**2*(2*A*b + B*a)/2 + 2*a**2*b**2*x**5*(A*b + B*a) + 5*a*b**3*x**8*(A*b + 2*B*a)/8 + b**4*x**11*(A*b + 5*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_38():
    f = (A + B*x**3)*(a + b*x**3)**5/x**6
    F = -A*a**5/(5*x**5) + B*b**5*x**13/13 - a**4*(5*A*b + B*a)/(2*x**2) + 5*a**3*b*x*(2*A*b + B*a) + 5*a**2*b**2*x**4*(A*b + B*a)/2 + 5*a*b**3*x**7*(A*b + 2*B*a)/7 + b**4*x**10*(A*b + 5*B*a)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_39():
    f = (A + B*x**3)*(a + b*x**3)**5/x**7
    F = -A*a**5/(6*x**6) + B*b**5*x**12/12 - a**4*(5*A*b + B*a)/(3*x**3) + 5*a**3*b*(2*A*b + B*a)*log(x) + 10*a**2*b**2*x**3*(A*b + B*a)/3 + 5*a*b**3*x**6*(A*b + 2*B*a)/6 + b**4*x**9*(A*b + 5*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_40():
    f = (A + B*x**3)*(a + b*x**3)**5/x**8
    F = -A*a**5/(7*x**7) + B*b**5*x**11/11 - a**4*(5*A*b + B*a)/(4*x**4) - 5*a**3*b*(2*A*b + B*a)/x + 5*a**2*b**2*x**2*(A*b + B*a) + a*b**3*x**5*(A*b + 2*B*a) + b**4*x**8*(A*b + 5*B*a)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_41():
    f = (A + B*x**3)*(a + b*x**3)**5/x**9
    F = -A*a**5/(8*x**8) + B*b**5*x**10/10 - a**4*(5*A*b + B*a)/(5*x**5) - 5*a**3*b*(2*A*b + B*a)/(2*x**2) + 10*a**2*b**2*x*(A*b + B*a) + 5*a*b**3*x**4*(A*b + 2*B*a)/4 + b**4*x**7*(A*b + 5*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_42():
    f = (A + B*x**3)*(a + b*x**3)**5/x**10
    F = -A*a**5/(9*x**9) + B*b**5*x**9/9 - a**4*(5*A*b + B*a)/(6*x**6) - 5*a**3*b*(2*A*b + B*a)/(3*x**3) + 10*a**2*b**2*(A*b + B*a)*log(x) + 5*a*b**3*x**3*(A*b + 2*B*a)/3 + b**4*x**6*(A*b + 5*B*a)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_43():
    f = (A + B*x**3)*(a + b*x**3)**5/x**11
    F = -A*a**5/(10*x**10) + B*b**5*x**8/8 - a**4*(5*A*b + B*a)/(7*x**7) - 5*a**3*b*(2*A*b + B*a)/(4*x**4) - 10*a**2*b**2*(A*b + B*a)/x + 5*a*b**3*x**2*(A*b + 2*B*a)/2 + b**4*x**5*(A*b + 5*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_44():
    f = (A + B*x**3)*(a + b*x**3)**5/x**12
    F = -A*a**5/(11*x**11) + B*b**5*x**7/7 - a**4*(5*A*b + B*a)/(8*x**8) - a**3*b*(2*A*b + B*a)/x**5 - 5*a**2*b**2*(A*b + B*a)/x**2 + 5*a*b**3*x*(A*b + 2*B*a) + b**4*x**4*(A*b + 5*B*a)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_45():
    f = (A + B*x**3)*(a + b*x**3)**5/x**13
    F = -A*a**5/(12*x**12) + B*b**5*x**6/6 - a**4*(5*A*b + B*a)/(9*x**9) - 5*a**3*b*(2*A*b + B*a)/(6*x**6) - 10*a**2*b**2*(A*b + B*a)/(3*x**3) + 5*a*b**3*(A*b + 2*B*a)*log(x) + b**4*x**3*(A*b + 5*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_46():
    f = (A + B*x**3)*(a + b*x**3)**5/x**14
    F = -A*a**5/(13*x**13) + B*b**5*x**5/5 - a**4*(5*A*b + B*a)/(10*x**10) - 5*a**3*b*(2*A*b + B*a)/(7*x**7) - 5*a**2*b**2*(A*b + B*a)/(2*x**4) - 5*a*b**3*(A*b + 2*B*a)/x + b**4*x**2*(A*b + 5*B*a)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_47():
    f = (A + B*x**3)*(a + b*x**3)**5/x**15
    F = -A*a**5/(14*x**14) + B*b**5*x**4/4 - a**4*(5*A*b + B*a)/(11*x**11) - 5*a**3*b*(2*A*b + B*a)/(8*x**8) - 2*a**2*b**2*(A*b + B*a)/x**5 - 5*a*b**3*(A*b + 2*B*a)/(2*x**2) + b**4*x*(A*b + 5*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_48():
    f = (A + B*x**3)*(a + b*x**3)**5/x**16
    F = -A*a**5/(15*x**15) + B*b**5*x**3/3 - a**4*(5*A*b + B*a)/(12*x**12) - 5*a**3*b*(2*A*b + B*a)/(9*x**9) - 5*a**2*b**2*(A*b + B*a)/(3*x**6) - 5*a*b**3*(A*b + 2*B*a)/(3*x**3) + b**4*(A*b + 5*B*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_49():
    f = (A + B*x**3)*(a + b*x**3)**5/x**17
    F = -A*a**5/(16*x**16) + B*b**5*x**2/2 - a**4*(5*A*b + B*a)/(13*x**13) - a**3*b*(2*A*b + B*a)/(2*x**10) - 10*a**2*b**2*(A*b + B*a)/(7*x**7) - 5*a*b**3*(A*b + 2*B*a)/(4*x**4) - b**4*(A*b + 5*B*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_50():
    f = (A + B*x**3)*(a + b*x**3)**5/x**18
    F = -A*a**5/(17*x**17) + B*b**5*x - a**4*(5*A*b + B*a)/(14*x**14) - 5*a**3*b*(2*A*b + B*a)/(11*x**11) - 5*a**2*b**2*(A*b + B*a)/(4*x**8) - a*b**3*(A*b + 2*B*a)/x**5 - b**4*(A*b + 5*B*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_51():
    f = (A + B*x**3)*(a + b*x**3)**5/x**19
    F = -A*(a + b*x**3)**6/(18*a*x**18) - B*a**5/(15*x**15) - 5*B*a**4*b/(12*x**12) - 10*B*a**3*b**2/(9*x**9) - 5*B*a**2*b**3/(3*x**6) - 5*B*a*b**4/(3*x**3) + B*b**5*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_52():
    f = (A + B*x**3)*(a + b*x**3)**5/x**20
    F = -A*a**5/(19*x**19) - B*b**5/x - a**4*(5*A*b + B*a)/(16*x**16) - 5*a**3*b*(2*A*b + B*a)/(13*x**13) - a**2*b**2*(A*b + B*a)/x**10 - 5*a*b**3*(A*b + 2*B*a)/(7*x**7) - b**4*(A*b + 5*B*a)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_53():
    f = (A + B*x**3)*(a + b*x**3)**5/x**21
    F = -A*a**5/(20*x**20) - B*b**5/(2*x**2) - a**4*(5*A*b + B*a)/(17*x**17) - 5*a**3*b*(2*A*b + B*a)/(14*x**14) - 10*a**2*b**2*(A*b + B*a)/(11*x**11) - 5*a*b**3*(A*b + 2*B*a)/(8*x**8) - b**4*(A*b + 5*B*a)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_54():
    f = (A + B*x**3)*(a + b*x**3)**5/x**22
    F = -A*(a + b*x**3)**6/(21*a*x**21) + (a + b*x**3)**6*(A*b - 7*B*a)/(126*a**2*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_55():
    f = (A + B*x**3)*(a + b*x**3)**5/x**23
    F = -A*a**5/(22*x**22) - B*b**5/(4*x**4) - a**4*(5*A*b + B*a)/(19*x**19) - 5*a**3*b*(2*A*b + B*a)/(16*x**16) - 10*a**2*b**2*(A*b + B*a)/(13*x**13) - a*b**3*(A*b + 2*B*a)/(2*x**10) - b**4*(A*b + 5*B*a)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_56():
    f = x**6*(A + B*x**3)/(a + b*x**3)
    F = B*x**7/(7*b) + a**(sympy.S(4)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(10)/3)) - a**(sympy.S(4)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(10)/3)) - sqrt(3)*a**(sympy.S(4)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(10)/3)) - a*x*(A*b - B*a)/b**3 + x**4*(A*b - B*a)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_57():
    f = x**5*(A + B*x**3)/(a + b*x**3)
    F = B*x**6/(6*b) - a*(A*b - B*a)*log(a + b*x**3)/(3*b**3) + x**3*(A*b - B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_58():
    f = x**4*(A + B*x**3)/(a + b*x**3)
    F = B*x**5/(5*b) + a**(sympy.S(2)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(8)/3)) - a**(sympy.S(2)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(8)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(8)/3)) + x**2*(A*b - B*a)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_59():
    f = x**3*(A + B*x**3)/(a + b*x**3)
    F = B*x**4/(4*b) - a**(sympy.S(1)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(7)/3)) + a**(sympy.S(1)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(7)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(7)/3)) + x*(A*b - B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_60():
    f = x**2*(A + B*x**3)/(a + b*x**3)
    F = B*x**3/(3*b) + (A*b - B*a)*log(a + b*x**3)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_61():
    f = x*(A + B*x**3)/(a + b*x**3)
    F = B*x**2/(2*b) - (A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3)) + (A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_62():
    f = (A + B*x**3)/(a + b*x**3)
    F = B*x/b + (A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - (A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_63():
    f = (A + B*x**3)/(x*(a + b*x**3))
    F = A*log(x)/a - (A*b - B*a)*log(a + b*x**3)/(3*a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_64():
    f = (A + B*x**3)/(x**2*(a + b*x**3))
    F = -A/(a*x) + (A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)) - (A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_65():
    f = (A + B*x**3)/(x**3*(a + b*x**3))
    F = -A/(2*a*x**2) - (A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)) + (A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_66():
    f = (A + B*x**3)/(x**4*(a + b*x**3))
    F = -A/(3*a*x**3) - (A*b - B*a)*log(x)/a**2 + (A*b - B*a)*log(a + b*x**3)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_67():
    f = (A + B*x**3)/(x**5*(a + b*x**3))
    F = -A/(4*a*x**4) + (A*b - B*a)/(a**2*x) - b**(sympy.S(1)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(7)/3)) + b**(sympy.S(1)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(7)/3)) - sqrt(3)*b**(sympy.S(1)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_68():
    f = (A + B*x**3)/(x**6*(a + b*x**3))
    F = -A/(5*a*x**5) + (A*b - B*a)/(2*a**2*x**2) + b**(sympy.S(2)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(8)/3)) - b**(sympy.S(2)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(8)/3)) - sqrt(3)*b**(sympy.S(2)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_69():
    f = (A + B*x**3)/(x**7*(a + b*x**3))
    F = -A/(6*a*x**6) + (A*b - B*a)/(3*a**2*x**3) + b*(A*b - B*a)*log(x)/a**3 - b*(A*b - B*a)*log(a + b*x**3)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_70():
    f = (A + B*x**3)/(x**8*(a + b*x**3))
    F = -A/(7*a*x**7) + (A*b - B*a)/(4*a**2*x**4) - b*(A*b - B*a)/(a**3*x) + b**(sympy.S(4)/3)*(A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(10)/3)) - b**(sympy.S(4)/3)*(A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(10)/3)) + sqrt(3)*b**(sympy.S(4)/3)*(A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_71():
    f = x**9*(A + B*x**3)/(a + b*x**3)**2
    F = a**(sympy.S(4)/3)*(7*A*b - 10*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(13)/3)) - a**(sympy.S(4)/3)*(7*A*b - 10*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(13)/3)) - sqrt(3)*a**(sympy.S(4)/3)*(7*A*b - 10*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(13)/3)) - a*x*(7*A*b - 10*B*a)/(3*b**4) + x**4*(7*A*b - 10*B*a)/(12*b**3) + x**10*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x**7*(7*A*b - 10*B*a)/(21*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_72():
    f = x**8*(A + B*x**3)/(a + b*x**3)**2
    F = B*x**6/(6*b**2) - a**2*(A*b - B*a)/(3*b**4*(a + b*x**3)) - a*(2*A*b - 3*B*a)*log(a + b*x**3)/(3*b**4) + x**3*(A*b - 2*B*a)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_73():
    f = x**7*(A + B*x**3)/(a + b*x**3)**2
    F = a**(sympy.S(2)/3)*(5*A*b - 8*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(11)/3)) - a**(sympy.S(2)/3)*(5*A*b - 8*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(11)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(5*A*b - 8*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(11)/3)) + x**2*(5*A*b - 8*B*a)/(6*b**3) + x**8*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x**5*(5*A*b - 8*B*a)/(15*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_74():
    f = x**6*(A + B*x**3)/(a + b*x**3)**2
    F = -a**(sympy.S(1)/3)*(4*A*b - 7*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(10)/3)) + a**(sympy.S(1)/3)*(4*A*b - 7*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(10)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(4*A*b - 7*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(10)/3)) + x*(4*A*b - 7*B*a)/(3*b**3) + x**7*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x**4*(4*A*b - 7*B*a)/(12*a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_75():
    f = x**5*(A + B*x**3)/(a + b*x**3)**2
    F = B*x**3/(3*b**2) + a*(A*b - B*a)/(3*b**3*(a + b*x**3)) + (A*b - 2*B*a)*log(a + b*x**3)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_76():
    f = x**4*(A + B*x**3)/(a + b*x**3)**2
    F = x**5*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x**2*(2*A*b - 5*B*a)/(6*a*b**2) - (2*A*b - 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(1)/3)*b**(sympy.S(8)/3)) + (2*A*b - 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(1)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(2*A*b - 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(1)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_77():
    f = x**3*(A + B*x**3)/(a + b*x**3)**2
    F = x**4*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x*(A*b - 4*B*a)/(3*a*b**2) + (A*b - 4*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - (A*b - 4*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(A*b - 4*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_78():
    f = x**2*(A + B*x**3)/(a + b*x**3)**2
    F = B*log(a + b*x**3)/(3*b**2) + (-A*b + B*a)/(3*b**2*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_79():
    f = x*(A + B*x**3)/(a + b*x**3)**2
    F = x**2*(A*b - B*a)/(3*a*b*(a + b*x**3)) - (A*b + 2*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)) + (A*b + 2*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(A*b + 2*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_80():
    f = (A + B*x**3)/(a + b*x**3)**2
    F = x*(A*b - B*a)/(3*a*b*(a + b*x**3)) + (2*A*b + B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - (2*A*b + B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(2*A*b + B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_81():
    f = (A + B*x**3)/(x*(a + b*x**3)**2)
    F = A*log(x)/a**2 - A*log(a + b*x**3)/(3*a**2) + (A*b - B*a)/(3*a*b*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_82():
    f = (A + B*x**3)/(x**2*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x*(a + b*x**3)) + (-4*A*b + B*a)/(3*a**2*b*x) + (4*A*b - B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)) - (4*A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(4*A*b - B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_83():
    f = (A + B*x**3)/(x**3*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x**2*(a + b*x**3)) + (-5*A*b + 2*B*a)/(6*a**2*b*x**2) - (5*A*b - 2*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)) + (5*A*b - 2*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(5*A*b - 2*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_84():
    f = (A + B*x**3)/(x**4*(a + b*x**3)**2)
    F = -A/(3*a**2*x**3) - (A*b - B*a)/(3*a**2*(a + b*x**3)) - (2*A*b - B*a)*log(x)/a**3 + (2*A*b - B*a)*log(a + b*x**3)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_85():
    f = (A + B*x**3)/(x**5*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x**4*(a + b*x**3)) + (-7*A*b + 4*B*a)/(12*a**2*b*x**4) + (7*A*b - 4*B*a)/(3*a**3*x) - b**(sympy.S(1)/3)*(7*A*b - 4*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(10)/3)) + b**(sympy.S(1)/3)*(7*A*b - 4*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(10)/3)) - sqrt(3)*b**(sympy.S(1)/3)*(7*A*b - 4*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_86():
    f = (A + B*x**3)/(x**6*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x**5*(a + b*x**3)) + (-8*A*b + 5*B*a)/(15*a**2*b*x**5) + (8*A*b - 5*B*a)/(6*a**3*x**2) + b**(sympy.S(2)/3)*(8*A*b - 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(11)/3)) - b**(sympy.S(2)/3)*(8*A*b - 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(11)/3)) - sqrt(3)*b**(sympy.S(2)/3)*(8*A*b - 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_87():
    f = (A + B*x**3)/(x**7*(a + b*x**3)**2)
    F = -A/(6*a**2*x**6) + b*(A*b - B*a)/(3*a**3*(a + b*x**3)) + (2*A*b - B*a)/(3*a**3*x**3) + b*(3*A*b - 2*B*a)*log(x)/a**4 - b*(3*A*b - 2*B*a)*log(a + b*x**3)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_88():
    f = x**11*(A + B*x**3)/(a + b*x**3)**3
    F = B*x**6/(6*b**3) + a**3*(A*b - B*a)/(6*b**5*(a + b*x**3)**2) - a**2*(3*A*b - 4*B*a)/(3*b**5*(a + b*x**3)) - a*(A*b - 2*B*a)*log(a + b*x**3)/b**5 + x**3*(A*b - 3*B*a)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_89():
    f = x**8*(A + B*x**3)/(a + b*x**3)**3
    F = B*x**3/(3*b**3) - a**2*(A*b - B*a)/(6*b**4*(a + b*x**3)**2) + a*(2*A*b - 3*B*a)/(3*b**4*(a + b*x**3)) + (A*b - 3*B*a)*log(a + b*x**3)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_90():
    f = x**5*(A + B*x**3)/(a + b*x**3)**3
    F = B*log(a + b*x**3)/(3*b**3) + a*(A*b - B*a)/(6*b**3*(a + b*x**3)**2) - (A*b - 2*B*a)/(3*b**3*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_91():
    f = x**2*(A + B*x**3)/(a + b*x**3)**3
    F = -(A + B*x**3)**2/((a + b*x**3)**2*(6*A*b - 6*B*a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_92():
    f = (A + B*x**3)/(x*(a + b*x**3)**3)
    F = A/(3*a**2*(a + b*x**3)) + A*log(x)/a**3 - A*log(a + b*x**3)/(3*a**3) + (A*b - B*a)/(6*a*b*(a + b*x**3)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_93():
    f = (A + B*x**3)/(x**4*(a + b*x**3)**3)
    F = -A/(3*a**3*x**3) - (A*b - B*a)/(6*a**2*(a + b*x**3)**2) - (2*A*b - B*a)/(3*a**3*(a + b*x**3)) - (3*A*b - B*a)*log(x)/a**4 + (3*A*b - B*a)*log(a + b*x**3)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_94():
    f = (A + B*x**3)/(x**7*(a + b*x**3)**3)
    F = -A/(6*a**3*x**6) + b*(A*b - B*a)/(6*a**3*(a + b*x**3)**2) + b*(3*A*b - 2*B*a)/(3*a**4*(a + b*x**3)) + (3*A*b - B*a)/(3*a**4*x**3) + 3*b*(2*A*b - B*a)*log(x)/a**5 - b*(2*A*b - B*a)*log(a + b*x**3)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_95():
    f = x**10*(A + B*x**3)/(a + b*x**3)**3
    F = 4*a**(sympy.S(2)/3)*(5*A*b - 11*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*b**(sympy.S(14)/3)) - 2*a**(sympy.S(2)/3)*(5*A*b - 11*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*b**(sympy.S(14)/3)) + 4*sqrt(3)*a**(sympy.S(2)/3)*(5*A*b - 11*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*b**(sympy.S(14)/3)) + x**2*(10*A*b - 22*B*a)/(9*b**4) + x**11*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**8*(5*A*b - 11*B*a)/(18*a*b**2*(a + b*x**3)) - x**5*(20*A*b - 44*B*a)/(45*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_96():
    f = x**9*(A + B*x**3)/(a + b*x**3)**3
    F = -7*a**(sympy.S(1)/3)*(2*A*b - 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*b**(sympy.S(13)/3)) + 7*a**(sympy.S(1)/3)*(2*A*b - 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*b**(sympy.S(13)/3)) + 7*sqrt(3)*a**(sympy.S(1)/3)*(2*A*b - 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*b**(sympy.S(13)/3)) + x*(14*A*b - 35*B*a)/(9*b**4) + x**10*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**7*(2*A*b - 5*B*a)/(9*a*b**2*(a + b*x**3)) - x**4*(14*A*b - 35*B*a)/(36*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_97():
    f = x**7*(A + B*x**3)/(a + b*x**3)**3
    F = x**8*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**5*(A*b - 4*B*a)/(9*a*b**2*(a + b*x**3)) + x**2*(-5*A*b + 20*B*a)/(18*a*b**3) - (5*A*b - 20*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(1)/3)*b**(sympy.S(11)/3)) + (5*A*b - 20*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(1)/3)*b**(sympy.S(11)/3)) - sqrt(3)*(5*A*b - 20*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(1)/3)*b**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_98():
    f = x**6*(A + B*x**3)/(a + b*x**3)**3
    F = x**7*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**4*(A*b - 7*B*a)/(18*a*b**2*(a + b*x**3)) + x*(-2*A*b + 14*B*a)/(9*a*b**3) - (A*b - 7*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) + (2*A*b - 14*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(2*A*b - 14*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(2)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_99():
    f = x**4*(A + B*x**3)/(a + b*x**3)**3
    F = x**5*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) - x**2*(A*b + 5*B*a)/(18*a*b**2*(a + b*x**3)) - (A*b + 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(4)/3)*b**(sympy.S(8)/3)) + (A*b + 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(4)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(A*b + 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(4)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_100():
    f = x**3*(A + B*x**3)/(a + b*x**3)**3
    F = x**4*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) - x*(A*b + 2*B*a)/(9*a*b**2*(a + b*x**3)) + (A*b + 2*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) - (A*b + 2*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(A*b + 2*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_101():
    f = x*(A + B*x**3)/(a + b*x**3)**3
    F = x**2*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**2*(2*A*b + B*a)/(9*a**2*b*(a + b*x**3)) - (2*A*b + B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(7)/3)*b**(sympy.S(5)/3)) + (2*A*b + B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(7)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(2*A*b + B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(7)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_102():
    f = (A + B*x**3)/(a + b*x**3)**3
    F = x*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x*(5*A*b + B*a)/(18*a**2*b*(a + b*x**3)) + (5*A*b + B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(4)/3)) - (5*A*b + B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(5*A*b + B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_103():
    f = (A + B*x**3)/(x**2*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x*(a + b*x**3)**2) + (7*A*b - B*a)/(18*a**2*b*x*(a + b*x**3)) + (-14*A*b + 2*B*a)/(9*a**3*b*x) - (7*A*b - B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(10)/3)*b**(sympy.S(2)/3)) + (14*A*b - 2*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(10)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(14*A*b - 2*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(10)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_104():
    f = (A + B*x**3)/(x**3*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x**2*(a + b*x**3)**2) + (4*A*b - B*a)/(9*a**2*b*x**2*(a + b*x**3)) + (-20*A*b + 5*B*a)/(18*a**3*b*x**2) - (20*A*b - 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)*b**(sympy.S(1)/3)) + (20*A*b - 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(11)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(20*A*b - 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_105():
    f = (A + B*x**3)/(x**5*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x**4*(a + b*x**3)**2) + (5*A*b - 2*B*a)/(9*a**2*b*x**4*(a + b*x**3)) + (-35*A*b + 14*B*a)/(36*a**3*b*x**4) + (35*A*b - 14*B*a)/(9*a**4*x) - 7*b**(sympy.S(1)/3)*(5*A*b - 2*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(13)/3)) + 7*b**(sympy.S(1)/3)*(5*A*b - 2*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(13)/3)) - 7*sqrt(3)*b**(sympy.S(1)/3)*(5*A*b - 2*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_106():
    f = (A + B*x**3)/(x**6*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x**5*(a + b*x**3)**2) + (11*A*b - 5*B*a)/(18*a**2*b*x**5*(a + b*x**3)) + (-44*A*b + 20*B*a)/(45*a**3*b*x**5) + (22*A*b - 10*B*a)/(9*a**4*x**2) + 4*b**(sympy.S(2)/3)*(11*A*b - 5*B*a)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(14)/3)) - 2*b**(sympy.S(2)/3)*(11*A*b - 5*B*a)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(14)/3)) - 4*sqrt(3)*b**(sympy.S(2)/3)*(11*A*b - 5*B*a)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_107():
    f = x**8/((a + b*x**3)*(c + d*x**3))
    F = a**2*log(a + b*x**3)/(3*b**2*(-a*d + b*c)) - c**2*log(c + d*x**3)/(3*d**2*(-a*d + b*c)) + x**3/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_108():
    f = x**7/((a + b*x**3)*(c + d*x**3))
    F = -a**(sympy.S(5)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(5)/3)*(-a*d + b*c)) + a**(sympy.S(5)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(5)/3)*(-a*d + b*c)) - sqrt(3)*a**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(5)/3)*(-a*d + b*c)) + c**(sympy.S(5)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*d**(sympy.S(5)/3)*(-a*d + b*c)) - c**(sympy.S(5)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(5)/3)*(-a*d + b*c)) + sqrt(3)*c**(sympy.S(5)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*d**(sympy.S(5)/3)*(-a*d + b*c)) + x**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_109():
    f = x**6/((a + b*x**3)*(c + d*x**3))
    F = a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(4)/3)*(-a*d + b*c)) - a**(sympy.S(4)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(4)/3)*(-a*d + b*c)) - sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(4)/3)*(-a*d + b*c)) - c**(sympy.S(4)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*d**(sympy.S(4)/3)*(-a*d + b*c)) + c**(sympy.S(4)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(4)/3)*(-a*d + b*c)) + sqrt(3)*c**(sympy.S(4)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*d**(sympy.S(4)/3)*(-a*d + b*c)) + x/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_110():
    f = x**5/((a + b*x**3)*(c + d*x**3))
    F = -a*log(a + b*x**3)/(3*b*(-a*d + b*c)) + c*log(c + d*x**3)/(3*d*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_111():
    f = x**4/((a + b*x**3)*(c + d*x**3))
    F = a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(2)/3)*(-a*d + b*c)) - a**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(2)/3)*(-a*d + b*c)) + sqrt(3)*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(2)/3)*(-a*d + b*c)) - c**(sympy.S(2)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*(-a*d + b*c)) + c**(sympy.S(2)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*(-a*d + b*c)) - sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_112():
    f = x**3/((a + b*x**3)*(c + d*x**3))
    F = -a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(1)/3)*(-a*d + b*c)) + a**(sympy.S(1)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(1)/3)*(-a*d + b*c)) + sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(1)/3)*(-a*d + b*c)) + c**(sympy.S(1)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)*(-a*d + b*c)) - c**(sympy.S(1)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(1)/3)*(-a*d + b*c)) - sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*d**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_113():
    f = x**2/((a + b*x**3)*(c + d*x**3))
    F = log(a + b*x**3)/(-3*a*d + 3*b*c) - log(c + d*x**3)/(-3*a*d + 3*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_114():
    f = x/((a + b*x**3)*(c + d*x**3))
    F = d**(sympy.S(1)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)*(-a*d + b*c)) - d**(sympy.S(1)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(1)/3)*(-a*d + b*c)) + sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(1)/3)*(-a*d + b*c)) - b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*(-a*d + b*c)) + b**(sympy.S(1)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*(-a*d + b*c)) - sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_115():
    f = 1/((a + b*x**3)*(c + d*x**3))
    F = -d**(sympy.S(2)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*(-a*d + b*c)) + d**(sympy.S(2)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*(-a*d + b*c)) + sqrt(3)*d**(sympy.S(2)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*(-a*d + b*c)) + b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*(-a*d + b*c)) - b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*(-a*d + b*c)) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_116():
    f = 1/(x*(a + b*x**3)*(c + d*x**3))
    F = d*log(c + d*x**3)/(3*c*(-a*d + b*c)) - b*log(a + b*x**3)/(3*a*(-a*d + b*c)) + log(x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_117():
    f = 1/(x**2*(a + b*x**3)*(c + d*x**3))
    F = -d**(sympy.S(4)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(4)/3)*(-a*d + b*c)) + d**(sympy.S(4)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(4)/3)*(-a*d + b*c)) - sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(4)/3)*(-a*d + b*c)) - 1/(a*c*x) + b**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*(-a*d + b*c)) - b**(sympy.S(4)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*(-a*d + b*c)) + sqrt(3)*b**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_118():
    f = 1/(x**3*(a + b*x**3)*(c + d*x**3))
    F = d**(sympy.S(5)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(5)/3)*(-a*d + b*c)) - d**(sympy.S(5)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(5)/3)*(-a*d + b*c)) - sqrt(3)*d**(sympy.S(5)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(5)/3)*(-a*d + b*c)) - 1/(2*a*c*x**2) - b**(sympy.S(5)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*(-a*d + b*c)) + b**(sympy.S(5)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*(-a*d + b*c)) + sqrt(3)*b**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_119():
    f = 1/(x**4*(a + b*x**3)*(c + d*x**3))
    F = -d**2*log(c + d*x**3)/(3*c**2*(-a*d + b*c)) - 1/(3*a*c*x**3) + b**2*log(a + b*x**3)/(3*a**2*(-a*d + b*c)) - (a*d + b*c)*log(x)/(a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_120():
    f = 1/(x**5*(a + b*x**3)*(c + d*x**3))
    F = d**(sympy.S(7)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(7)/3)*(-a*d + b*c)) - d**(sympy.S(7)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(7)/3)*(-a*d + b*c)) + sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(7)/3)*(-a*d + b*c)) - 1/(4*a*c*x**4) + (a*d + b*c)/(a**2*c**2*x) - b**(sympy.S(7)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(7)/3)*(-a*d + b*c)) + b**(sympy.S(7)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(7)/3)*(-a*d + b*c)) - sqrt(3)*b**(sympy.S(7)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_121():
    f = 1/(x**6*(a + b*x**3)*(c + d*x**3))
    F = -d**(sympy.S(8)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(8)/3)*(-a*d + b*c)) + d**(sympy.S(8)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(8)/3)*(-a*d + b*c)) + sqrt(3)*d**(sympy.S(8)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(8)/3)*(-a*d + b*c)) - 1/(5*a*c*x**5) + (a*d + b*c)/(2*a**2*c**2*x**2) + b**(sympy.S(8)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(8)/3)*(-a*d + b*c)) - b**(sympy.S(8)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(8)/3)*(-a*d + b*c)) - sqrt(3)*b**(sympy.S(8)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(8)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_122():
    f = 1/(x**7*(a + b*x**3)*(c + d*x**3))
    F = d**3*log(c + d*x**3)/(3*c**3*(-a*d + b*c)) - 1/(6*a*c*x**6) + (a*d + b*c)/(3*a**2*c**2*x**3) - b**3*log(a + b*x**3)/(3*a**3*(-a*d + b*c)) + (a**2*d**2 + a*b*c*d + b**2*c**2)*log(x)/(a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_123():
    f = 1/(x**8*(a + b*x**3)*(c + d*x**3))
    F = -d**(sympy.S(10)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(10)/3)*(-a*d + b*c)) + d**(sympy.S(10)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(10)/3)*(-a*d + b*c)) - sqrt(3)*d**(sympy.S(10)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(10)/3)*(-a*d + b*c)) - 1/(7*a*c*x**7) + (a*d + b*c)/(4*a**2*c**2*x**4) - (a**2*d**2 + a*b*c*d + b**2*c**2)/(a**3*c**3*x) + b**(sympy.S(10)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(10)/3)*(-a*d + b*c)) - b**(sympy.S(10)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(10)/3)*(-a*d + b*c)) + sqrt(3)*b**(sympy.S(10)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(10)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_124():
    f = x**m*(A + B*x**3)*(a + b*x**3)**5
    F = A*a**5*x**(m + 1)/(m + 1) + B*b**5*x**(m + 19)/(m + 19) + a**4*x**(m + 4)*(5*A*b + B*a)/(m + 4) + 5*a**3*b*x**(m + 7)*(2*A*b + B*a)/(m + 7) + 10*a**2*b**2*x**(m + 10)*(A*b + B*a)/(m + 10) + 5*a*b**3*x**(m + 13)*(A*b + 2*B*a)/(m + 13) + b**4*x**(m + 16)*(A*b + 5*B*a)/(m + 16)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_125():
    f = x**m*(A + B*x**3)*(a + b*x**3)**2
    F = A*a**2*x**(m + 1)/(m + 1) + B*b**2*x**(m + 10)/(m + 10) + a*x**(m + 4)*(2*A*b + B*a)/(m + 4) + b*x**(m + 7)*(A*b + 2*B*a)/(m + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_126():
    f = x**m*(A + B*x**3)*(a + b*x**3)
    F = A*a*x**(m + 1)/(m + 1) + B*b*x**(m + 7)/(m + 7) + x**(m + 4)*(A*b + B*a)/(m + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_127():
    f = x**m*(A + B*x**3)/(a + b*x**3)
    F = B*x**(m + 1)/(b*(m + 1)) + x**(m + 1)*(A*b - B*a)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(a*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_128():
    f = x**m*(A + B*x**3)/(a + b*x**3)**2
    F = x**(m + 1)*(A*b - B*a)/(3*a*b*(a + b*x**3)) + x**(m + 1)*(A*b*(2 - m) + B*a*(m + 1))*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(3*a**2*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_129():
    f = x**m*(A + B*x**3)/(a + b*x**3)**3
    F = x**(m + 1)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**(m + 1)*(A*b*(5 - m) + B*a*(m + 1))*hyper((2, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(6*a**3*b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_130():
    f = (e*x)**m/((a + b*x**3)*(c + d*x**3))
    F = -d*(e*x)**(m + 1)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -d*x**3/c)/(c*e*(m + 1)*(-a*d + b*c)) + b*(e*x)**(m + 1)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(a*e*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_131():
    f = x**(sympy.S(7)/2)*(A + B*x**3)*(a + b*x**3)
    F = 2*A*a*x**(sympy.S(9)/2)/9 + 2*B*b*x**(sympy.S(21)/2)/21 + x**(sympy.S(15)/2)*(2*A*b + 2*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_132():
    f = x**(sympy.S(5)/2)*(A + B*x**3)*(a + b*x**3)
    F = 2*A*a*x**(sympy.S(7)/2)/7 + 2*B*b*x**(sympy.S(19)/2)/19 + x**(sympy.S(13)/2)*(2*A*b + 2*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_133():
    f = x**(sympy.S(3)/2)*(A + B*x**3)*(a + b*x**3)
    F = 2*A*a*x**(sympy.S(5)/2)/5 + 2*B*b*x**(sympy.S(17)/2)/17 + x**(sympy.S(11)/2)*(2*A*b + 2*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_134():
    f = sqrt(x)*(A + B*x**3)*(a + b*x**3)
    F = 2*A*a*x**(sympy.S(3)/2)/3 + 2*B*b*x**(sympy.S(15)/2)/15 + x**(sympy.S(9)/2)*(2*A*b + 2*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_135():
    f = (A + B*x**3)*(a + b*x**3)/sqrt(x)
    F = 2*A*a*sqrt(x) + 2*B*b*x**(sympy.S(13)/2)/13 + x**(sympy.S(7)/2)*(2*A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_136():
    f = (A + B*x**3)*(a + b*x**3)/x**(sympy.S(3)/2)
    F = -2*A*a/sqrt(x) + 2*B*b*x**(sympy.S(11)/2)/11 + x**(sympy.S(5)/2)*(2*A*b + 2*B*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_137():
    f = (A + B*x**3)*(a + b*x**3)/x**(sympy.S(5)/2)
    F = -2*A*a/(3*x**(sympy.S(3)/2)) + 2*B*b*x**(sympy.S(9)/2)/9 + x**(sympy.S(3)/2)*(2*A*b + 2*B*a)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_138():
    f = (A + B*x**3)*(a + b*x**3)/x**(sympy.S(7)/2)
    F = -2*A*a/(5*x**(sympy.S(5)/2)) + 2*B*b*x**(sympy.S(7)/2)/7 + sqrt(x)*(2*A*b + 2*B*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_139():
    f = x**(sympy.S(7)/2)*(A + B*x**3)*(a + b*x**3)**2
    F = 2*A*a**2*x**(sympy.S(9)/2)/9 + 2*B*b**2*x**(sympy.S(27)/2)/27 + 2*a*x**(sympy.S(15)/2)*(2*A*b + B*a)/15 + 2*b*x**(sympy.S(21)/2)*(A*b + 2*B*a)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_140():
    f = x**(sympy.S(5)/2)*(A + B*x**3)*(a + b*x**3)**2
    F = 2*A*a**2*x**(sympy.S(7)/2)/7 + 2*B*b**2*x**(sympy.S(25)/2)/25 + 2*a*x**(sympy.S(13)/2)*(2*A*b + B*a)/13 + 2*b*x**(sympy.S(19)/2)*(A*b + 2*B*a)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_141():
    f = x**(sympy.S(3)/2)*(A + B*x**3)*(a + b*x**3)**2
    F = 2*A*a**2*x**(sympy.S(5)/2)/5 + 2*B*b**2*x**(sympy.S(23)/2)/23 + 2*a*x**(sympy.S(11)/2)*(2*A*b + B*a)/11 + 2*b*x**(sympy.S(17)/2)*(A*b + 2*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_142():
    f = sqrt(x)*(A + B*x**3)*(a + b*x**3)**2
    F = 2*A*a**2*x**(sympy.S(3)/2)/3 + 2*B*b**2*x**(sympy.S(21)/2)/21 + 2*a*x**(sympy.S(9)/2)*(2*A*b + B*a)/9 + 2*b*x**(sympy.S(15)/2)*(A*b + 2*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_143():
    f = (A + B*x**3)*(a + b*x**3)**2/sqrt(x)
    F = 2*A*a**2*sqrt(x) + 2*B*b**2*x**(sympy.S(19)/2)/19 + 2*a*x**(sympy.S(7)/2)*(2*A*b + B*a)/7 + 2*b*x**(sympy.S(13)/2)*(A*b + 2*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_144():
    f = (A + B*x**3)*(a + b*x**3)**2/x**(sympy.S(3)/2)
    F = -2*A*a**2/sqrt(x) + 2*B*b**2*x**(sympy.S(17)/2)/17 + 2*a*x**(sympy.S(5)/2)*(2*A*b + B*a)/5 + 2*b*x**(sympy.S(11)/2)*(A*b + 2*B*a)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_145():
    f = (A + B*x**3)*(a + b*x**3)**2/x**(sympy.S(5)/2)
    F = -2*A*a**2/(3*x**(sympy.S(3)/2)) + 2*B*b**2*x**(sympy.S(15)/2)/15 + 2*a*x**(sympy.S(3)/2)*(2*A*b + B*a)/3 + 2*b*x**(sympy.S(9)/2)*(A*b + 2*B*a)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_146():
    f = (A + B*x**3)*(a + b*x**3)**2/x**(sympy.S(7)/2)
    F = -2*A*a**2/(5*x**(sympy.S(5)/2)) + 2*B*b**2*x**(sympy.S(13)/2)/13 + 2*a*sqrt(x)*(2*A*b + B*a) + 2*b*x**(sympy.S(7)/2)*(A*b + 2*B*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_147():
    f = x**(sympy.S(7)/2)*(A + B*x**3)*(a + b*x**3)**3
    F = 2*A*a**3*x**(sympy.S(9)/2)/9 + 2*B*b**3*x**(sympy.S(33)/2)/33 + 2*a**2*x**(sympy.S(15)/2)*(3*A*b + B*a)/15 + 2*a*b*x**(sympy.S(21)/2)*(A*b + B*a)/7 + 2*b**2*x**(sympy.S(27)/2)*(A*b + 3*B*a)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_148():
    f = x**(sympy.S(5)/2)*(A + B*x**3)*(a + b*x**3)**3
    F = 2*A*a**3*x**(sympy.S(7)/2)/7 + 2*B*b**3*x**(sympy.S(31)/2)/31 + 2*a**2*x**(sympy.S(13)/2)*(3*A*b + B*a)/13 + 6*a*b*x**(sympy.S(19)/2)*(A*b + B*a)/19 + 2*b**2*x**(sympy.S(25)/2)*(A*b + 3*B*a)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_149():
    f = x**(sympy.S(3)/2)*(A + B*x**3)*(a + b*x**3)**3
    F = 2*A*a**3*x**(sympy.S(5)/2)/5 + 2*B*b**3*x**(sympy.S(29)/2)/29 + 2*a**2*x**(sympy.S(11)/2)*(3*A*b + B*a)/11 + 6*a*b*x**(sympy.S(17)/2)*(A*b + B*a)/17 + 2*b**2*x**(sympy.S(23)/2)*(A*b + 3*B*a)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_150():
    f = sqrt(x)*(A + B*x**3)*(a + b*x**3)**3
    F = 2*A*a**3*x**(sympy.S(3)/2)/3 + 2*B*b**3*x**(sympy.S(27)/2)/27 + 2*a**2*x**(sympy.S(9)/2)*(3*A*b + B*a)/9 + 2*a*b*x**(sympy.S(15)/2)*(A*b + B*a)/5 + 2*b**2*x**(sympy.S(21)/2)*(A*b + 3*B*a)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_151():
    f = (A + B*x**3)*(a + b*x**3)**3/sqrt(x)
    F = 2*A*a**3*sqrt(x) + 2*B*b**3*x**(sympy.S(25)/2)/25 + 2*a**2*x**(sympy.S(7)/2)*(3*A*b + B*a)/7 + 6*a*b*x**(sympy.S(13)/2)*(A*b + B*a)/13 + 2*b**2*x**(sympy.S(19)/2)*(A*b + 3*B*a)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_152():
    f = (A + B*x**3)*(a + b*x**3)**3/x**(sympy.S(3)/2)
    F = -2*A*a**3/sqrt(x) + 2*B*b**3*x**(sympy.S(23)/2)/23 + 2*a**2*x**(sympy.S(5)/2)*(3*A*b + B*a)/5 + 6*a*b*x**(sympy.S(11)/2)*(A*b + B*a)/11 + 2*b**2*x**(sympy.S(17)/2)*(A*b + 3*B*a)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_153():
    f = (A + B*x**3)*(a + b*x**3)**3/x**(sympy.S(5)/2)
    F = -2*A*a**3/(3*x**(sympy.S(3)/2)) + 2*B*b**3*x**(sympy.S(21)/2)/21 + 2*a**2*x**(sympy.S(3)/2)*(3*A*b + B*a)/3 + 2*a*b*x**(sympy.S(9)/2)*(A*b + B*a)/3 + 2*b**2*x**(sympy.S(15)/2)*(A*b + 3*B*a)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_154():
    f = (A + B*x**3)*(a + b*x**3)**3/x**(sympy.S(7)/2)
    F = -2*A*a**3/(5*x**(sympy.S(5)/2)) + 2*B*b**3*x**(sympy.S(19)/2)/19 + 2*a**2*sqrt(x)*(3*A*b + B*a) + 6*a*b*x**(sympy.S(7)/2)*(A*b + B*a)/7 + 2*b**2*x**(sympy.S(13)/2)*(A*b + 3*B*a)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_155():
    f = x**(sympy.S(7)/2)*(A + B*x**3)/(a + b*x**3)
    F = 2*B*x**(sympy.S(9)/2)/(9*b) - 2*sqrt(a)*(A*b - B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*b**(sympy.S(5)/2)) + x**(sympy.S(3)/2)*(2*A*b - 2*B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_156():
    f = x**(sympy.S(5)/2)*(A + B*x**3)/(a + b*x**3)
    F = 2*B*x**(sympy.S(7)/2)/(7*b) + sqrt(3)*a**(sympy.S(1)/6)*(A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*b**(sympy.S(13)/6)) - sqrt(3)*a**(sympy.S(1)/6)*(A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*b**(sympy.S(13)/6)) - 2*a**(sympy.S(1)/6)*(A*b - B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*b**(sympy.S(13)/6)) + a**(sympy.S(1)/6)*(A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*b**(sympy.S(13)/6)) - a**(sympy.S(1)/6)*(A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*b**(sympy.S(13)/6)) + sqrt(x)*(2*A*b - 2*B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_157():
    f = x**(sympy.S(3)/2)*(A + B*x**3)/(a + b*x**3)
    F = 2*B*x**(sympy.S(5)/2)/(5*b) + sqrt(3)*(A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(1)/6)*b**(sympy.S(11)/6)) - sqrt(3)*(A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(1)/6)*b**(sympy.S(11)/6)) - (A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(1)/6)*b**(sympy.S(11)/6)) + (A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(1)/6)*b**(sympy.S(11)/6)) + (2*A*b - 2*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(1)/6)*b**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_158():
    f = sqrt(x)*(A + B*x**3)/(a + b*x**3)
    F = 2*B*x**(sympy.S(3)/2)/(3*b) + (2*A*b - 2*B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_159():
    f = (A + B*x**3)/(sqrt(x)*(a + b*x**3))
    F = 2*B*sqrt(x)/b - sqrt(3)*(A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(5)/6)*b**(sympy.S(7)/6)) + sqrt(3)*(A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(5)/6)*b**(sympy.S(7)/6)) - (A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*b**(sympy.S(7)/6)) + (A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*b**(sympy.S(7)/6)) + (2*A*b - 2*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*b**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_160():
    f = (A + B*x**3)/(x**(sympy.S(3)/2)*(a + b*x**3))
    F = -2*A/(a*sqrt(x)) - sqrt(3)*(A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(7)/6)*b**(sympy.S(5)/6)) + sqrt(3)*(A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(7)/6)*b**(sympy.S(5)/6)) + (A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(7)/6)*b**(sympy.S(5)/6)) - (A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(7)/6)*b**(sympy.S(5)/6)) - (2*A*b - 2*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(7)/6)*b**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_161():
    f = (A + B*x**3)/(x**(sympy.S(5)/2)*(a + b*x**3))
    F = -2*A/(3*a*x**(sympy.S(3)/2)) - (2*A*b - 2*B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_162():
    f = (A + B*x**3)/(x**(sympy.S(7)/2)*(a + b*x**3))
    F = -2*A/(5*a*x**(sympy.S(5)/2)) + sqrt(3)*(A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(11)/6)*b**(sympy.S(1)/6)) - sqrt(3)*(A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(6*a**(sympy.S(11)/6)*b**(sympy.S(1)/6)) + (A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(11)/6)*b**(sympy.S(1)/6)) - (A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(11)/6)*b**(sympy.S(1)/6)) - (2*A*b - 2*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(11)/6)*b**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_163():
    f = x**(sympy.S(7)/2)*(A + B*x**3)/(a + b*x**3)**2
    F = x**(sympy.S(9)/2)*(A*b - B*a)/(3*a*b*(a + b*x**3)) - x**(sympy.S(3)/2)*(A*b - 3*B*a)/(3*a*b**2) + (A*b - 3*B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_164():
    f = x**(sympy.S(5)/2)*(A + B*x**3)/(a + b*x**3)**2
    F = x**(sympy.S(7)/2)*(A*b - B*a)/(3*a*b*(a + b*x**3)) - sqrt(x)*(A*b - 7*B*a)/(3*a*b**2) - sqrt(3)*(A*b - 7*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(5)/6)*b**(sympy.S(13)/6)) + sqrt(3)*(A*b - 7*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(5)/6)*b**(sympy.S(13)/6)) + (A*b - 7*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(9*a**(sympy.S(5)/6)*b**(sympy.S(13)/6)) - (A*b - 7*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(5)/6)*b**(sympy.S(13)/6)) + (A*b - 7*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(5)/6)*b**(sympy.S(13)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_165():
    f = x**(sympy.S(3)/2)*(A + B*x**3)/(a + b*x**3)**2
    F = x**(sympy.S(5)/2)*(A*b - B*a)/(3*a*b*(a + b*x**3)) + sqrt(3)*(A*b + 5*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(7)/6)*b**(sympy.S(11)/6)) - sqrt(3)*(A*b + 5*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(7)/6)*b**(sympy.S(11)/6)) + (A*b + 5*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(9*a**(sympy.S(7)/6)*b**(sympy.S(11)/6)) - (A*b + 5*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(7)/6)*b**(sympy.S(11)/6)) + (A*b + 5*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(7)/6)*b**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_166():
    f = sqrt(x)*(A + B*x**3)/(a + b*x**3)**2
    F = x**(sympy.S(3)/2)*(A*b - B*a)/(3*a*b*(a + b*x**3)) + (A*b + B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_167():
    f = (A + B*x**3)/(sqrt(x)*(a + b*x**3)**2)
    F = sqrt(x)*(A*b - B*a)/(3*a*b*(a + b*x**3)) - sqrt(3)*(5*A*b + B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(11)/6)*b**(sympy.S(7)/6)) + sqrt(3)*(5*A*b + B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(11)/6)*b**(sympy.S(7)/6)) + (5*A*b + B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(9*a**(sympy.S(11)/6)*b**(sympy.S(7)/6)) - (5*A*b + B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(11)/6)*b**(sympy.S(7)/6)) + (5*A*b + B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(11)/6)*b**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_168():
    f = (A + B*x**3)/(x**(sympy.S(3)/2)*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*sqrt(x)*(a + b*x**3)) - (7*A*b - B*a)/(3*a**2*b*sqrt(x)) - sqrt(3)*(7*A*b - B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(13)/6)*b**(sympy.S(5)/6)) + sqrt(3)*(7*A*b - B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(13)/6)*b**(sympy.S(5)/6)) - (7*A*b - B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(9*a**(sympy.S(13)/6)*b**(sympy.S(5)/6)) + (7*A*b - B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(13)/6)*b**(sympy.S(5)/6)) - (7*A*b - B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(13)/6)*b**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_169():
    f = (A + B*x**3)/(x**(sympy.S(5)/2)*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x**(sympy.S(3)/2)*(a + b*x**3)) + (-3*A*b + B*a)/(3*a**2*b*x**(sympy.S(3)/2)) - (3*A*b - B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(3*a**(sympy.S(5)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_170():
    f = (A + B*x**3)/(x**(sympy.S(7)/2)*(a + b*x**3)**2)
    F = (A*b - B*a)/(3*a*b*x**(sympy.S(5)/2)*(a + b*x**3)) - (11*A*b - 5*B*a)/(15*a**2*b*x**(sympy.S(5)/2)) + sqrt(3)*(11*A*b - 5*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(17)/6)*b**(sympy.S(1)/6)) - sqrt(3)*(11*A*b - 5*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(36*a**(sympy.S(17)/6)*b**(sympy.S(1)/6)) - (11*A*b - 5*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(9*a**(sympy.S(17)/6)*b**(sympy.S(1)/6)) + (11*A*b - 5*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(17)/6)*b**(sympy.S(1)/6)) - (11*A*b - 5*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(18*a**(sympy.S(17)/6)*b**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_171():
    f = x**(sympy.S(7)/2)*(A + B*x**3)/(a + b*x**3)**3
    F = x**(sympy.S(9)/2)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) - x**(sympy.S(3)/2)*(A*b + 3*B*a)/(12*a*b**2*(a + b*x**3)) + (A*b + 3*B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(12*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_172():
    f = x**(sympy.S(5)/2)*(A + B*x**3)/(a + b*x**3)**3
    F = x**(sympy.S(7)/2)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) - sqrt(x)*(5*A*b + 7*B*a)/(36*a*b**2*(a + b*x**3)) - sqrt(3)*(5*A*b + 7*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(11)/6)*b**(sympy.S(13)/6)) + sqrt(3)*(5*A*b + 7*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(11)/6)*b**(sympy.S(13)/6)) + (5*A*b + 7*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(108*a**(sympy.S(11)/6)*b**(sympy.S(13)/6)) - (5*A*b + 7*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(11)/6)*b**(sympy.S(13)/6)) + (5*A*b + 7*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(11)/6)*b**(sympy.S(13)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_173():
    f = x**(sympy.S(3)/2)*(A + B*x**3)/(a + b*x**3)**3
    F = x**(sympy.S(5)/2)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**(sympy.S(5)/2)*(7*A*b + 5*B*a)/(36*a**2*b*(a + b*x**3)) + sqrt(3)*(7*A*b + 5*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(13)/6)*b**(sympy.S(11)/6)) - sqrt(3)*(7*A*b + 5*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(13)/6)*b**(sympy.S(11)/6)) + (7*A*b + 5*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(108*a**(sympy.S(13)/6)*b**(sympy.S(11)/6)) - (7*A*b + 5*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(13)/6)*b**(sympy.S(11)/6)) + (7*A*b + 5*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(13)/6)*b**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_174():
    f = sqrt(x)*(A + B*x**3)/(a + b*x**3)**3
    F = x**(sympy.S(3)/2)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + x**(sympy.S(3)/2)*(3*A*b + B*a)/(12*a**2*b*(a + b*x**3)) + (3*A*b + B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(12*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_175():
    f = (A + B*x**3)/(sqrt(x)*(a + b*x**3)**3)
    F = sqrt(x)*(A*b - B*a)/(6*a*b*(a + b*x**3)**2) + sqrt(x)*(11*A*b + B*a)/(36*a**2*b*(a + b*x**3)) - sqrt(3)*(55*A*b + 5*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(17)/6)*b**(sympy.S(7)/6)) + sqrt(3)*(55*A*b + 5*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(17)/6)*b**(sympy.S(7)/6)) + (55*A*b + 5*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(108*a**(sympy.S(17)/6)*b**(sympy.S(7)/6)) - (55*A*b + 5*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(17)/6)*b**(sympy.S(7)/6)) + (55*A*b + 5*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(17)/6)*b**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_176():
    f = (A + B*x**3)/(x**(sympy.S(3)/2)*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*sqrt(x)*(a + b*x**3)**2) + (13*A*b - B*a)/(36*a**2*b*sqrt(x)*(a + b*x**3)) - (91*A*b - 7*B*a)/(36*a**3*b*sqrt(x)) - sqrt(3)*(91*A*b - 7*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(19)/6)*b**(sympy.S(5)/6)) + sqrt(3)*(91*A*b - 7*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(19)/6)*b**(sympy.S(5)/6)) - (91*A*b - 7*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(108*a**(sympy.S(19)/6)*b**(sympy.S(5)/6)) + (91*A*b - 7*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(19)/6)*b**(sympy.S(5)/6)) - (91*A*b - 7*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(19)/6)*b**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_177():
    f = (A + B*x**3)/(x**(sympy.S(5)/2)*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x**(sympy.S(3)/2)*(a + b*x**3)**2) + (5*A*b - B*a)/(12*a**2*b*x**(sympy.S(3)/2)*(a + b*x**3)) + (-5*A*b + B*a)/(4*a**3*b*x**(sympy.S(3)/2)) - (5*A*b - B*a)*atan(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a))/(4*a**(sympy.S(7)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_178():
    f = (A + B*x**3)/(x**(sympy.S(7)/2)*(a + b*x**3)**3)
    F = (A*b - B*a)/(6*a*b*x**(sympy.S(5)/2)*(a + b*x**3)**2) + (17*A*b - 5*B*a)/(36*a**2*b*x**(sympy.S(5)/2)*(a + b*x**3)) - (187*A*b - 55*B*a)/(180*a**3*b*x**(sympy.S(5)/2)) + sqrt(3)*(187*A*b - 55*B*a)*log(-sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(23)/6)*b**(sympy.S(1)/6)) - sqrt(3)*(187*A*b - 55*B*a)*log(sqrt(3)*a**(sympy.S(1)/6)*b**(sympy.S(1)/6)*sqrt(x) + a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(432*a**(sympy.S(23)/6)*b**(sympy.S(1)/6)) - (187*A*b - 55*B*a)*atan(b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(108*a**(sympy.S(23)/6)*b**(sympy.S(1)/6)) + (187*A*b - 55*B*a)*atan(sqrt(3) - 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(23)/6)*b**(sympy.S(1)/6)) - (187*A*b - 55*B*a)*atan(sqrt(3) + 2*b**(sympy.S(1)/6)*sqrt(x)/a**(sympy.S(1)/6))/(216*a**(sympy.S(23)/6)*b**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_179():
    f = x**8*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(9)/2)/(27*b**4) + 2*a**2*(a + b*x**3)**(sympy.S(3)/2)*(A*b - B*a)/(9*b**4) - 2*a*(a + b*x**3)**(sympy.S(5)/2)*(2*A*b - 3*B*a)/(15*b**4) + (a + b*x**3)**(sympy.S(7)/2)*(2*A*b - 6*B*a)/(21*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_180():
    f = x**5*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(7)/2)/(21*b**3) - 2*a*(a + b*x**3)**(sympy.S(3)/2)*(A*b - B*a)/(9*b**3) + (a + b*x**3)**(sympy.S(5)/2)*(2*A*b - 4*B*a)/(15*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_181():
    f = x**2*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(5)/2)/(15*b**2) + (a + b*x**3)**(sympy.S(3)/2)*(2*A*b - 2*B*a)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_182():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x
    F = -2*A*sqrt(a)*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*A*sqrt(a + b*x**3)/3 + 2*B*(a + b*x**3)**(sympy.S(3)/2)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_183():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**4
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(3*a*x**3) + sqrt(a + b*x**3)*(A*b + 2*B*a)/(3*a) - (A*b + 2*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_184():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**7
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(6*a*x**6) + sqrt(a + b*x**3)*(A*b - 4*B*a)/(12*a*x**3) + b*(A*b - 4*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(12*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_185():
    f = x**3*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*x**4*(a + b*x**3)**(sympy.S(3)/2)/(17*b) - 4*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(935*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*x*sqrt(a + b*x**3)*(17*A*b - 8*B*a)/(935*b**2) + x**4*sqrt(a + b*x**3)*(34*A*b - 16*B*a)/(187*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_186():
    f = (A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*x*(a + b*x**3)**(sympy.S(3)/2)/(11*b) + 2*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(55*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + x*sqrt(a + b*x**3)*(22*A*b - 4*B*a)/(55*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_187():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**3
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(2*a*x**2) + 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b + 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(10*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + x*sqrt(a + b*x**3)*(5*A*b + 4*B*a)/(10*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_188():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**6
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(5*a*x**5) - 3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(20*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(A*b - 10*B*a)/(20*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_189():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**9
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(8*a*x**8) + sqrt(a + b*x**3)*(7*A*b - 16*B*a)/(80*a*x**5) + 3**(sympy.S(3)/4)*b**(sympy.S(5)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 16*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(320*a**2*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3*b*sqrt(a + b*x**3)*(7*A*b - 16*B*a)/(320*a**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_190():
    f = x**4*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*x**5*(a + b*x**3)**(sympy.S(3)/2)/(19*b) + 12*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(19*A*b - 10*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(19*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 24*a**2*sqrt(a + b*x**3)*(19*A*b - 10*B*a)/(1729*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 6*a*x**2*sqrt(a + b*x**3)*(19*A*b - 10*B*a)/(1729*b**2) + x**5*sqrt(a + b*x**3)*(38*A*b - 20*B*a)/(247*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_191():
    f = x*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*x**2*(a + b*x**3)**(sympy.S(3)/2)/(13*b) - 3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*sqrt(a + b*x**3)*(13*A*b - 4*B*a)/(91*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + x**2*sqrt(a + b*x**3)*(26*A*b - 8*B*a)/(91*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_192():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**2
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(a*x) - 3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b + 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(21*A*b + 6*B*a)/(7*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + x**2*sqrt(a + b*x**3)*(7*A*b + 2*B*a)/(7*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_193():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**5
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(4*a*x**4) + 3*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(A*b + 8*B*a)/(8*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - sqrt(a + b*x**3)*(A*b + 8*B*a)/(8*a*x) - 3*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 8*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_194():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**8
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(7*a*x**7) + sqrt(a + b*x**3)*(5*A*b - 14*B*a)/(56*a*x**4) - 3*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(5*A*b - 14*B*a)/(112*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3*b*sqrt(a + b*x**3)*(5*A*b - 14*B*a)/(112*a**2*x) + 3*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 14*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_195():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**11
    F = -A*(a + b*x**3)**(sympy.S(3)/2)/(10*a*x**10) + sqrt(a + b*x**3)*(11*A*b - 20*B*a)/(140*a*x**7) + 3*b*sqrt(a + b*x**3)*(11*A*b - 20*B*a)/(1120*a**2*x**4) + 3*b**(sympy.S(7)/3)*sqrt(a + b*x**3)*(11*A*b - 20*B*a)/(448*a**3*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3*b**2*sqrt(a + b*x**3)*(11*A*b - 20*B*a)/(448*a**3*x) - 3*3**(sympy.S(1)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 20*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(896*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 20*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(448*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_196():
    f = x**8*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(11)/2)/(33*b**4) + 2*a**2*(a + b*x**3)**(sympy.S(5)/2)*(A*b - B*a)/(15*b**4) - 2*a*(a + b*x**3)**(sympy.S(7)/2)*(2*A*b - 3*B*a)/(21*b**4) + (a + b*x**3)**(sympy.S(9)/2)*(2*A*b - 6*B*a)/(27*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_197():
    f = x**5*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(9)/2)/(27*b**3) - 2*a*(a + b*x**3)**(sympy.S(5)/2)*(A*b - B*a)/(15*b**3) + (a + b*x**3)**(sympy.S(7)/2)*(2*A*b - 4*B*a)/(21*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_198():
    f = x**2*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(7)/2)/(21*b**2) + (a + b*x**3)**(sympy.S(5)/2)*(2*A*b - 2*B*a)/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_199():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x
    F = -2*A*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*A*a*sqrt(a + b*x**3)/3 + 2*A*(a + b*x**3)**(sympy.S(3)/2)/9 + 2*B*(a + b*x**3)**(sympy.S(5)/2)/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_200():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**4
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(3*a*x**3) - sqrt(a)*(3*A*b + 2*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + sqrt(a + b*x**3)*(3*A*b + 2*B*a)/3 + (a + b*x**3)**(sympy.S(3)/2)*(3*A*b + 2*B*a)/(9*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_201():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**7
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(6*a*x**6) + b*sqrt(a + b*x**3)*(A*b + 4*B*a)/(4*a) - (a + b*x**3)**(sympy.S(3)/2)*(A*b + 4*B*a)/(12*a*x**3) - b*(A*b + 4*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_202():
    f = x**3*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**4*(a + b*x**3)**(sympy.S(5)/2)/(23*b) - 36*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(23*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21505*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*x*sqrt(a + b*x**3)*(23*A*b - 8*B*a)/(21505*b**2) + 18*a*x**4*sqrt(a + b*x**3)*(23*A*b - 8*B*a)/(4301*b) + x**4*(a + b*x**3)**(sympy.S(3)/2)*(46*A*b - 16*B*a)/(391*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_203():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x*(a + b*x**3)**(sympy.S(5)/2)/(17*b) + 18*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(935*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 18*a*x*sqrt(a + b*x**3)*(17*A*b - 2*B*a)/(935*b) + x*(a + b*x**3)**(sympy.S(3)/2)*(34*A*b - 4*B*a)/(187*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_204():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**3
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(2*a*x**2) + 9*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b + 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(110*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + x*sqrt(a + b*x**3)*(9*A*b/10 + 18*B*a/55) + x*(a + b*x**3)**(sympy.S(3)/2)*(11*A*b + 4*B*a)/(22*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_205():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**6
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(5*a*x**5) + 9*3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(20*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*b*x*sqrt(a + b*x**3)*(A*b + 2*B*a)/(20*a) - (a + b*x**3)**(sympy.S(3)/2)*(A*b + 2*B*a)/(4*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_206():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**9
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(8*a*x**8) - 9*3**(sympy.S(3)/4)*b**(sympy.S(5)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 16*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(320*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*b*sqrt(a + b*x**3)*(A*b - 16*B*a)/(320*a*x**2) + (a + b*x**3)**(sympy.S(3)/2)*(A*b - 16*B*a)/(80*a*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_207():
    f = x**4*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**5*(a + b*x**3)**(sympy.S(5)/2)/(25*b) + 108*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8645*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 72*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8645*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 216*a**3*sqrt(a + b*x**3)*(5*A*b - 2*B*a)/(8645*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 54*a**2*x**2*sqrt(a + b*x**3)*(5*A*b - 2*B*a)/(8645*b**2) + 18*a*x**5*sqrt(a + b*x**3)*(5*A*b - 2*B*a)/(1235*b) + x**5*(a + b*x**3)**(sympy.S(3)/2)*(10*A*b - 4*B*a)/(95*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_208():
    f = x*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**2*(a + b*x**3)**(sympy.S(5)/2)/(19*b) - 27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(19*A*b - 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 18*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(19*A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*sqrt(a + b*x**3)*(19*A*b - 4*B*a)/(1729*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 18*a*x**2*sqrt(a + b*x**3)*(19*A*b - 4*B*a)/(1729*b) + x**2*(a + b*x**3)**(sympy.S(3)/2)*(38*A*b - 8*B*a)/(247*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_209():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**2
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(a*x) - 27*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b + 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(182*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*a*sqrt(a + b*x**3)*(13*A*b + 2*B*a)/(91*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + x**2*sqrt(a + b*x**3)*(9*A*b/7 + 18*B*a/91) + x**2*(a + b*x**3)**(sympy.S(3)/2)*(13*A*b + 2*B*a)/(13*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_210():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**5
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(4*a*x**4) - 27*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b + 8*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b + 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(56*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(7*A*b + 8*B*a)/(56*a**(sympy.S(1)/3)*(1 + sqrt(3)) + 56*b**(sympy.S(1)/3)*x) + 9*b*x**2*sqrt(a + b*x**3)*(7*A*b + 8*B*a)/(56*a) - (a + b*x**3)**(sympy.S(3)/2)*(7*A*b + 8*B*a)/(8*a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_211():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**8
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(7*a*x**7) + 27*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(A*b + 14*B*a)/(112*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 9*b*sqrt(a + b*x**3)*(A*b + 14*B*a)/(112*a*x) - (a + b*x**3)**(sympy.S(3)/2)*(A*b + 14*B*a)/(56*a*x**4) - 27*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 14*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_212():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/x**11
    F = -A*(a + b*x**3)**(sympy.S(5)/2)/(10*a*x**10) + 9*b*sqrt(a + b*x**3)*(A*b - 4*B*a)/(224*a*x**4) + (a + b*x**3)**(sympy.S(3)/2)*(A*b - 4*B*a)/(28*a*x**7) - 27*b**(sympy.S(7)/3)*sqrt(a + b*x**3)*(A*b - 4*B*a)/(448*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 27*b**2*sqrt(a + b*x**3)*(A*b - 4*B*a)/(448*a**2*x) + 27*3**(sympy.S(1)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(896*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 9*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(448*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_213():
    f = x**8*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(7)/2)/(21*b**4) + 2*a**2*sqrt(a + b*x**3)*(A*b - B*a)/(3*b**4) - 2*a*(a + b*x**3)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(9*b**4) + (a + b*x**3)**(sympy.S(5)/2)*(2*A*b - 6*B*a)/(15*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_214():
    f = x**5*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(5)/2)/(15*b**3) - 2*a*sqrt(a + b*x**3)*(A*b - B*a)/(3*b**3) + (a + b*x**3)**(sympy.S(3)/2)*(2*A*b - 4*B*a)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_215():
    f = x**2*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*(a + b*x**3)**(sympy.S(3)/2)/(9*b**2) + sqrt(a + b*x**3)*(2*A*b - 2*B*a)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_216():
    f = (A + B*x**3)/(x*sqrt(a + b*x**3))
    F = -2*A*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)) + 2*B*sqrt(a + b*x**3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_217():
    f = (A + B*x**3)/(x**4*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(3*a*x**3) + (A*b - 2*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_218():
    f = (A + B*x**3)/(x**7*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(6*a*x**6) + sqrt(a + b*x**3)*(3*A*b - 4*B*a)/(12*a**2*x**3) - b*(3*A*b - 4*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(12*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_219():
    f = x**3*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*x**4*sqrt(a + b*x**3)/(11*b) - 4*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(165*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + x*sqrt(a + b*x**3)*(22*A*b - 16*B*a)/(55*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_220():
    f = (A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*x*sqrt(a + b*x**3)/(5*b) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(15*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_221():
    f = (A + B*x**3)/(x**3*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(2*a*x**2) - 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*a*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_222():
    f = (A + B*x**3)/(x**6*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(5*a*x**5) + 3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(60*a**2*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(7*A*b - 10*B*a)/(20*a**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_223():
    f = x**4*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*x**5*sqrt(a + b*x**3)/(13*b) + 4*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 10*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(273*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*a*sqrt(a + b*x**3)*(13*A*b - 10*B*a)/(91*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + x**2*sqrt(a + b*x**3)*(26*A*b - 20*B*a)/(91*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_224():
    f = x*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*x**2*sqrt(a + b*x**3)/(7*b) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(14*A*b - 8*B*a)/(7*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_225():
    f = (A + B*x**3)/(x**2*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(a*x) + sqrt(a + b*x**3)*(A*b + 2*B*a)/(a*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_226():
    f = (A + B*x**3)/(x**5*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(4*a*x**4) - b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(5*A*b - 8*B*a)/(8*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + sqrt(a + b*x**3)*(5*A*b - 8*B*a)/(8*a**2*x) + 3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 8*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(24*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_227():
    f = (A + B*x**3)/(x**8*sqrt(a + b*x**3))
    F = -A*sqrt(a + b*x**3)/(7*a*x**7) + sqrt(a + b*x**3)*(11*A*b - 14*B*a)/(56*a**2*x**4) + 5*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(11*A*b - 14*B*a)/(112*a**3*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 5*b*sqrt(a + b*x**3)*(11*A*b - 14*B*a)/(112*a**3*x) - 5*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 14*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 5*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(336*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_228():
    f = x**8*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(5)/2)/(15*b**4) - 2*a**2*(A*b - B*a)/(3*b**4*sqrt(a + b*x**3)) - 2*a*sqrt(a + b*x**3)*(2*A*b - 3*B*a)/(3*b**4) + (a + b*x**3)**(sympy.S(3)/2)*(2*A*b - 6*B*a)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_229():
    f = x**5*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(3)/2)/(9*b**3) + 2*a*(A*b - B*a)/(3*b**3*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(2*A*b - 4*B*a)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_230():
    f = x**2*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*sqrt(a + b*x**3)/(3*b**2) + (-2*A*b + 2*B*a)/(3*b**2*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_231():
    f = (A + B*x**3)/(x*(a + b*x**3)**(sympy.S(3)/2))
    F = -2*A*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(3)/2)) + (2*A*b - 2*B*a)/(3*a*b*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_232():
    f = (A + B*x**3)/(x**4*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(3*a*x**3*sqrt(a + b*x**3)) + (-3*A*b + 2*B*a)/(3*a**2*sqrt(a + b*x**3)) + (3*A*b - 2*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_233():
    f = (A + B*x**3)/(x**7*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(6*a*x**6*sqrt(a + b*x**3)) + (5*A*b - 4*B*a)/(12*a**2*x**3*sqrt(a + b*x**3)) + b*(5*A*b - 4*B*a)/(4*a**3*sqrt(a + b*x**3)) - b*(5*A*b - 4*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_234():
    f = x**6*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**7/(11*b*sqrt(a + b*x**3)) - 32*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(495*b**(sympy.S(10)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - x**4*(22*A*b - 28*B*a)/(33*b**2*sqrt(a + b*x**3)) + x*sqrt(a + b*x**3)*(176*A*b - 224*B*a)/(165*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_235():
    f = x**3*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**4/(5*b*sqrt(a + b*x**3)) - x*(10*A*b - 16*B*a)/(15*b**2*sqrt(a + b*x**3)) + 4*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(45*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_236():
    f = (A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = x*(2*A*b - 2*B*a)/(3*a*b*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_237():
    f = (A + B*x**3)/(x**3*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(2*a*x**2*sqrt(a + b*x**3)) - x*(7*A*b - 4*B*a)/(6*a**2*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(18*a**2*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_238():
    f = (A + B*x**3)/(x**6*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(5*a*x**5*sqrt(a + b*x**3)) - (13*A*b - 10*B*a)/(15*a**2*x**2*sqrt(a + b*x**3)) + 7*3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(180*a**3*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(91*A*b - 70*B*a)/(60*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_239():
    f = x**4*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*x**5/(7*b*sqrt(a + b*x**3)) - 4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 10*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 8*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(63*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - x**2*(14*A*b - 20*B*a)/(21*b**2*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(56*A*b - 80*B*a)/(21*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_240():
    f = x*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = x**2*(2*A*b - 2*B*a)/(3*a*b*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(2*A*b - 8*B*a)/(3*a*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_241():
    f = (A + B*x**3)/(x**2*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(a*x*sqrt(a + b*x**3)) - x**2*(5*A*b - 2*B*a)/(3*a**2*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(5*A*b - 2*B*a)/(3*a**2*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_242():
    f = (A + B*x**3)/(x**5*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(4*a*x**4*sqrt(a + b*x**3)) - (11*A*b - 8*B*a)/(12*a**2*x*sqrt(a + b*x**3)) - 5*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(11*A*b - 8*B*a)/(24*a**3*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + sqrt(a + b*x**3)*(55*A*b - 40*B*a)/(24*a**3*x) + 5*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 8*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(48*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 5*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(72*a**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_243():
    f = (A + B*x**3)/(x**8*(a + b*x**3)**(sympy.S(3)/2))
    F = -A/(7*a*x**7*sqrt(a + b*x**3)) - (17*A*b - 14*B*a)/(21*a**2*x**4*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(187*A*b - 154*B*a)/(168*a**3*x**4) + 55*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(17*A*b - 14*B*a)/(336*a**4*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 55*b*sqrt(a + b*x**3)*(17*A*b - 14*B*a)/(336*a**4*x) - 55*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 14*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(672*a**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 55*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1008*a**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_244():
    f = x**8*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*B*(a + b*x**3)**(sympy.S(3)/2)/(9*b**4) - 2*a**2*(A*b - B*a)/(9*b**4*(a + b*x**3)**(sympy.S(3)/2)) + 2*a*(2*A*b - 3*B*a)/(3*b**4*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(2*A*b - 6*B*a)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_245():
    f = x**5*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*B*sqrt(a + b*x**3)/(3*b**3) + 2*a*(A*b - B*a)/(9*b**3*(a + b*x**3)**(sympy.S(3)/2)) - (2*A*b - 4*B*a)/(3*b**3*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_246():
    f = x**2*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = -2*B/(3*b**2*sqrt(a + b*x**3)) + (-2*A*b + 2*B*a)/(9*b**2*(a + b*x**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_247():
    f = (A + B*x**3)/(x*(a + b*x**3)**(sympy.S(5)/2))
    F = 2*A/(3*a**2*sqrt(a + b*x**3)) - 2*A*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(5)/2)) + (2*A*b - 2*B*a)/(9*a*b*(a + b*x**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_248():
    f = (A + B*x**3)/(x**4*(a + b*x**3)**(sympy.S(5)/2))
    F = -A/(3*a*x**3*(a + b*x**3)**(sympy.S(3)/2)) + (-5*A*b + 2*B*a)/(9*a**2*(a + b*x**3)**(sympy.S(3)/2)) - (5*A*b - 2*B*a)/(3*a**3*sqrt(a + b*x**3)) + (5*A*b - 2*B*a)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_249():
    f = x**6*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*B*x**7/(5*b*(a + b*x**3)**(sympy.S(3)/2)) - x**4*(10*A*b - 28*B*a)/(45*b**2*(a + b*x**3)**(sympy.S(3)/2)) - x*(80*A*b - 224*B*a)/(135*b**3*sqrt(a + b*x**3)) + 32*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b - 14*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(405*b**(sympy.S(10)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_250():
    f = x**3*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = x**4*(2*A*b - 2*B*a)/(9*a*b*(a + b*x**3)**(sympy.S(3)/2)) - x*(2*A*b + 16*B*a)/(27*a*b**2*sqrt(a + b*x**3)) + 4*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b + 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_251():
    f = (A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = x*(2*A*b - 2*B*a)/(9*a*b*(a + b*x**3)**(sympy.S(3)/2)) + x*(14*A*b + 4*B*a)/(27*a**2*b*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b + 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**2*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_252():
    f = (A + B*x**3)/(x**3*(a + b*x**3)**(sympy.S(5)/2))
    F = -A/(2*a*x**2*(a + b*x**3)**(sympy.S(3)/2)) - x*(13*A*b - 4*B*a)/(18*a**2*(a + b*x**3)**(sympy.S(3)/2)) - x*(91*A*b - 28*B*a)/(54*a**3*sqrt(a + b*x**3)) - 7*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(13*A*b - 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(162*a**3*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_253():
    f = (A + B*x**3)/(x**6*(a + b*x**3)**(sympy.S(5)/2))
    F = -A/(5*a*x**5*(a + b*x**3)**(sympy.S(3)/2)) - (19*A*b - 10*B*a)/(45*a**2*x**2*(a + b*x**3)**(sympy.S(3)/2)) - (247*A*b - 130*B*a)/(135*a**3*x**2*sqrt(a + b*x**3)) + 91*3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(19*A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1620*a**4*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(1729*A*b - 910*B*a)/(540*a**4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_254():
    f = x**7*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*B*x**8/(7*b*(a + b*x**3)**(sympy.S(3)/2)) - 40*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 16*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(189*b**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 80*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*A*b - 16*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(567*b**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - x**5*(14*A*b - 32*B*a)/(63*b**2*(a + b*x**3)**(sympy.S(3)/2)) - x**2*(140*A*b - 320*B*a)/(189*b**3*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(560*A*b - 1280*B*a)/(189*b**(sympy.S(11)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_255():
    f = x**4*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = x**5*(2*A*b - 2*B*a)/(9*a*b*(a + b*x**3)**(sympy.S(3)/2)) + x**2*(2*A*b - 20*B*a)/(27*a*b**2*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(8*A*b - 80*B*a)/(27*a*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 4*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 10*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(A*b - 10*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_256():
    f = x*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = x**2*(2*A*b - 2*B*a)/(9*a*b*(a + b*x**3)**(sympy.S(3)/2)) + x**2*(10*A*b + 8*B*a)/(27*a**2*b*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(10*A*b + 8*B*a)/(27*a**2*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b + 4*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(5*A*b + 4*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_257():
    f = (A + B*x**3)/(x**2*(a + b*x**3)**(sympy.S(5)/2))
    F = -A/(a*x*(a + b*x**3)**(sympy.S(3)/2)) - x**2*(11*A*b - 2*B*a)/(9*a**2*(a + b*x**3)**(sympy.S(3)/2)) - x**2*(55*A*b - 10*B*a)/(27*a**3*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(55*A*b - 10*B*a)/(27*a**3*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 5*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 2*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(54*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 5*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(11*A*b - 2*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_258():
    f = (A + B*x**3)/(x**5*(a + b*x**3)**(sympy.S(5)/2))
    F = -A/(4*a*x**4*(a + b*x**3)**(sympy.S(3)/2)) - (17*A*b - 8*B*a)/(36*a**2*x*(a + b*x**3)**(sympy.S(3)/2)) - (187*A*b - 88*B*a)/(108*a**3*x*sqrt(a + b*x**3)) - 55*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(17*A*b - 8*B*a)/(216*a**4*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + sqrt(a + b*x**3)*(935*A*b - 440*B*a)/(216*a**4*x) + 55*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 8*B*a)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(432*a**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 55*sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(17*A*b - 8*B*a)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(648*a**(sympy.S(11)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_259():
    f = x**8*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = -32*sqrt(3)*c**(sympy.S(5)/2)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(3*d**3) + 32*c**2*sqrt(c + d*x**3)/(3*d**3) - 10*c*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3) + 2*(c + d*x**3)**(sympy.S(5)/2)/(15*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_260():
    f = x**5*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = 8*sqrt(3)*c**(sympy.S(3)/2)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(3*d**2) - 8*c*sqrt(c + d*x**3)/(3*d**2) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_261():
    f = x**2*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = -2*sqrt(3)*sqrt(c)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(3*d) + 2*sqrt(c + d*x**3)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_262():
    f = sqrt(c + d*x**3)/(x*(4*c + d*x**3))
    F = sqrt(3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(6*sqrt(c)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(6*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_263():
    f = sqrt(c + d*x**3)/(x**4*(4*c + d*x**3))
    F = -sqrt(c + d*x**3)/(12*c*x**3) - sqrt(3)*d*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(24*c**(sympy.S(3)/2)) - d*atanh(sqrt(c + d*x**3)/sqrt(c))/(24*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_264():
    f = x**4*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = 2*2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(3*d**(sympy.S(5)/3)) - 2*2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(3*d**(sympy.S(5)/3)) + 2*2**(sympy.S(1)/3)*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*d**(sympy.S(5)/3)) - 2*2**(sympy.S(1)/3)*c**(sympy.S(7)/6)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(5)/3) + 25*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 50*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 50*c*sqrt(c + d*x**3)/(7*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 2*x**2*sqrt(c + d*x**3)/(7*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_265():
    f = x*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = -2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(6*d**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(6*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/sqrt(c))/(6*d**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*c**(sympy.S(1)/6)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(2*d**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(c + d*x**3)/(d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_266():
    f = sqrt(c + d*x**3)/(x**2*(4*c + d*x**3))
    F = d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(4*c*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(4*c*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(12*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2**(sympy.S(1)/3)*sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(24*c**(sympy.S(5)/6)) - 2**(sympy.S(1)/3)*sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(24*c**(sympy.S(5)/6)) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/sqrt(c))/(24*c**(sympy.S(5)/6)) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(8*c**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_267():
    f = x**3*sqrt(c + d*x**3)/(4*c + d*x**3)
    F = x**4*sqrt(c + d*x**3)*appellf1(sympy.S(4)/3, sympy.S(-1)/2, 1, sympy.S(7)/3, -d*x**3/c, -d*x**3/(4*c))/(16*c*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_268():
    f = sqrt(c + d*x**3)/(4*c + d*x**3)
    F = x*sqrt(c + d*x**3)*appellf1(sympy.S(1)/3, sympy.S(-1)/2, 1, sympy.S(4)/3, -d*x**3/c, -d*x**3/(4*c))/(4*c*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_269():
    f = sqrt(c + d*x**3)/(x**3*(4*c + d*x**3))
    F = -sqrt(c + d*x**3)*appellf1(sympy.S(-2)/3, sympy.S(-1)/2, 1, sympy.S(1)/3, -d*x**3/c, -d*x**3/(4*c))/(8*c*x**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_270():
    f = x**8/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = 32*sqrt(3)*c**(sympy.S(3)/2)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**3) - 10*c*sqrt(c + d*x**3)/(3*d**3) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_271():
    f = x**5/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = -8*sqrt(3)*sqrt(c)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**2) + 2*sqrt(c + d*x**3)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_272():
    f = x**2/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = 2*sqrt(3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(9*sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_273():
    f = 1/(x*sqrt(c + d*x**3)*(4*c + d*x**3))
    F = -sqrt(3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(18*c**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(6*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_274():
    f = 1/(x**4*sqrt(c + d*x**3)*(4*c + d*x**3))
    F = -sqrt(c + d*x**3)/(12*c**2*x**3) + sqrt(3)*d*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(72*c**(sympy.S(5)/2)) + d*atanh(sqrt(c + d*x**3)/sqrt(c))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_275():
    f = x**4/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = -2*2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**(sympy.S(5)/3)) + 2*2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9*d**(sympy.S(5)/3)) - 2*2**(sympy.S(1)/3)*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/sqrt(c))/(9*d**(sympy.S(5)/3)) + 2*2**(sympy.S(1)/3)*c**(sympy.S(1)/6)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(3*d**(sympy.S(5)/3)) - 3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(c + d*x**3)/(d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_276():
    f = x/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/sqrt(c))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(6*c**(sympy.S(5)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_277():
    f = 1/(x**2*sqrt(c + d*x**3)*(4*c + d*x**3))
    F = d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(4*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(4*c**2*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(12*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2**(sympy.S(1)/3)*sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*sqrt(c + d*x**3)/(3*sqrt(c)))/(72*c**(sympy.S(11)/6)) + 2**(sympy.S(1)/3)*sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(72*c**(sympy.S(11)/6)) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/sqrt(c))/(72*c**(sympy.S(11)/6)) + 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(24*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_278():
    f = x**3/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -d*x**3/c, -d*x**3/(4*c))/(16*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_279():
    f = 1/(sqrt(c + d*x**3)*(4*c + d*x**3))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -d*x**3/c, -d*x**3/(4*c))/(4*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_280():
    f = 1/(x**3*sqrt(c + d*x**3)*(4*c + d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, -d*x**3/c, -d*x**3/(4*c))/(8*c*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_281():
    f = x/(sqrt(1 - x**3)*(4 - x**3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(1 - x**3)/3)/18 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/18 - 2**(sympy.S(1)/3)*atanh((2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/6 + 2**(sympy.S(1)/3)*atanh(sqrt(1 - x**3))/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_282():
    f = x**11*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = 1024*c**(sympy.S(7)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**4 - 1024*c**3*sqrt(c + d*x**3)/(3*d**4) - 38*c**2*(c + d*x**3)**(sympy.S(3)/2)/(3*d**4) - 4*c*(c + d*x**3)**(sympy.S(5)/2)/(5*d**4) - 2*(c + d*x**3)**(sympy.S(7)/2)/(21*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_283():
    f = x**8*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = 128*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**3 - 128*c**2*sqrt(c + d*x**3)/(3*d**3) - 14*c*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3) - 2*(c + d*x**3)**(sympy.S(5)/2)/(15*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_284():
    f = x**5*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = 16*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**2 - 16*c*sqrt(c + d*x**3)/(3*d**2) - 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_285():
    f = x**2*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = 2*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d - 2*sqrt(c + d*x**3)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_286():
    f = sqrt(c + d*x**3)/(x*(8*c - d*x**3))
    F = atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(4*sqrt(c)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(12*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_287():
    f = sqrt(c + d*x**3)/(x**4*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(24*c*x**3) + d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(32*c**(sympy.S(3)/2)) - 5*d*atanh(sqrt(c + d*x**3)/sqrt(c))/(96*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_288():
    f = sqrt(c + d*x**3)/(x**7*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(48*c*x**6) - d*sqrt(c + d*x**3)/(64*c**2*x**3) + d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(256*c**(sympy.S(5)/2)) + d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(256*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_289():
    f = x**7*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = -32*sqrt(3)*c**(sympy.S(13)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(8)/3) - 32*c**(sympy.S(13)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(8)/3) + 32*c**(sympy.S(13)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(8)/3) + 6124*3**(sympy.S(1)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 12248*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(273*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 12248*c**2*sqrt(c + d*x**3)/(91*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 214*c*x**2*sqrt(c + d*x**3)/(91*d**2) - 2*x**5*sqrt(c + d*x**3)/(13*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_290():
    f = x**4*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = -4*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(5)/3) - 4*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(5)/3) + 4*c**(sympy.S(7)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(5)/3) + 59*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 118*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 118*c*sqrt(c + d*x**3)/(7*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 2*x**2*sqrt(c + d*x**3)/(7*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_291():
    f = x*sqrt(c + d*x**3)/(8*c - d*x**3)
    F = -sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(2*d**(sympy.S(2)/3)) - c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2*d**(sympy.S(2)/3)) + c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(2*d**(sympy.S(2)/3)) + 3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2*sqrt(c + d*x**3)/(d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_292():
    f = sqrt(c + d*x**3)/(x**2*(8*c - d*x**3))
    F = d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(8*c*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(8*c*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(24*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(16*c**(sympy.S(5)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(16*c**(sympy.S(5)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(16*c**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_293():
    f = sqrt(c + d*x**3)/(x**5*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(32*c*x**4) + d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(16*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - d*sqrt(c + d*x**3)/(16*c**2*x) - 3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(32*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(48*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(128*c**(sympy.S(11)/6)) - d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(128*c**(sympy.S(11)/6)) + d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(128*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_294():
    f = sqrt(c + d*x**3)/(x**8*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(56*c*x**7) - 19*d*sqrt(c + d*x**3)/(1792*c**2*x**4) - d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(112*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + d**2*sqrt(c + d*x**3)/(112*c**3*x) + 3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(336*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1024*c**(sympy.S(17)/6)) - d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1024*c**(sympy.S(17)/6)) + d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1024*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_295():
    f = x**11*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = 9216*c**(sympy.S(9)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**4 - 3072*c**4*sqrt(c + d*x**3)/d**4 - 1024*c**3*(c + d*x**3)**(sympy.S(3)/2)/(9*d**4) - 38*c**2*(c + d*x**3)**(sympy.S(5)/2)/(5*d**4) - 4*c*(c + d*x**3)**(sympy.S(7)/2)/(7*d**4) - 2*(c + d*x**3)**(sympy.S(9)/2)/(27*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_296():
    f = x**8*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = 1152*c**(sympy.S(7)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**3 - 384*c**3*sqrt(c + d*x**3)/d**3 - 128*c**2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3) - 14*c*(c + d*x**3)**(sympy.S(5)/2)/(15*d**3) - 2*(c + d*x**3)**(sympy.S(7)/2)/(21*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_297():
    f = x**5*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = 144*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**2 - 48*c**2*sqrt(c + d*x**3)/d**2 - 16*c*(c + d*x**3)**(sympy.S(3)/2)/(9*d**2) - 2*(c + d*x**3)**(sympy.S(5)/2)/(15*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_298():
    f = x**2*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = 18*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d - 6*c*sqrt(c + d*x**3)/d - 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_299():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x*(8*c - d*x**3))
    F = 9*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/4 - sqrt(c)*atanh(sqrt(c + d*x**3)/sqrt(c))/12 - 2*sqrt(c + d*x**3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_300():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**4*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(24*x**3) + 9*d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(32*sqrt(c)) - 13*d*atanh(sqrt(c + d*x**3)/sqrt(c))/(96*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_301():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**7*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(48*x**6) - 11*d*sqrt(c + d*x**3)/(192*c*x**3) + 9*d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(256*c**(sympy.S(3)/2)) - 37*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(768*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_302():
    f = x**7*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = -288*sqrt(3)*c**(sympy.S(19)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(8)/3) - 288*c**(sympy.S(19)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(8)/3) + 288*c**(sympy.S(19)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(8)/3) + 1047324*3**(sympy.S(1)/4)*c**(sympy.S(10)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 698216*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(10)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2094648*c**3*sqrt(c + d*x**3)/(1729*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 36534*c**2*x**2*sqrt(c + d*x**3)/(1729*d**2) - 348*c*x**5*sqrt(c + d*x**3)/(247*d) - 2*x**8*sqrt(c + d*x**3)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_303():
    f = x**4*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = -36*sqrt(3)*c**(sympy.S(13)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(5)/3) - 36*c**(sympy.S(13)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(5)/3) + 36*c**(sympy.S(13)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(5)/3) + 6891*3**(sympy.S(1)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 4594*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 13782*c**2*sqrt(c + d*x**3)/(91*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 240*c*x**2*sqrt(c + d*x**3)/(91*d) - 2*x**5*sqrt(c + d*x**3)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_304():
    f = x*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)
    F = -9*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(2*d**(sympy.S(2)/3)) - 9*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2*d**(sympy.S(2)/3)) + 9*c**(sympy.S(7)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(2*d**(sympy.S(2)/3)) + 66*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 44*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 132*c*sqrt(c + d*x**3)/(7*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 2*x**2*sqrt(c + d*x**3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_305():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**2*(8*c - d*x**3))
    F = -9*sqrt(3)*c**(sympy.S(1)/6)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/16 - 9*c**(sympy.S(1)/6)*d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/16 + 9*c**(sympy.S(1)/6)*d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/16 + 15*3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 5*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 15*d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(8*c**(sympy.S(1)/3)*(1 + sqrt(3)) + 8*d**(sympy.S(1)/3)*x) - sqrt(c + d*x**3)/(8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_306():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**5*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(32*x**4) + 3*d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(16*c*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3*d*sqrt(c + d*x**3)/(16*c*x) - 3*3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(32*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 9*sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(128*c**(sympy.S(5)/6)) - 9*d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(128*c**(sympy.S(5)/6)) + 9*d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(128*c**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_307():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**8*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(56*x**7) - 75*d*sqrt(c + d*x**3)/(1792*c*x**4) + 3*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(56*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3*d**2*sqrt(c + d*x**3)/(56*c**2*x) - 3*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(56*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 9*sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1024*c**(sympy.S(11)/6)) - 9*d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1024*c**(sympy.S(11)/6)) + 9*d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1024*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_308():
    f = x**11/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = 1024*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**4) - 38*c**2*sqrt(c + d*x**3)/d**4 - 4*c*(c + d*x**3)**(sympy.S(3)/2)/(3*d**4) - 2*(c + d*x**3)**(sympy.S(5)/2)/(15*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_309():
    f = x**8/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = 128*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**3) - 14*c*sqrt(c + d*x**3)/(3*d**3) - 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_310():
    f = x**5/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = 16*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**2) - 2*sqrt(c + d*x**3)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_311():
    f = x**2/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = 2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_312():
    f = 1/(x*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(36*c**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(12*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_313():
    f = 1/(x**4*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(24*c**2*x**3) + d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(288*c**(sympy.S(5)/2)) + d*atanh(sqrt(c + d*x**3)/sqrt(c))/(32*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_314():
    f = 1/(x**7*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(48*c**2*x**6) + 5*d*sqrt(c + d*x**3)/(192*c**3*x**3) + d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2304*c**(sympy.S(7)/2)) - 7*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(256*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_315():
    f = x**7/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -32*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9*d**(sympy.S(8)/3)) - 32*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**(sympy.S(8)/3)) + 32*c**(sympy.S(7)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9*d**(sympy.S(8)/3)) + 52*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 104*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 104*c*sqrt(c + d*x**3)/(7*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 2*x**2*sqrt(c + d*x**3)/(7*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_316():
    f = x**4/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -4*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9*d**(sympy.S(5)/3)) - 4*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**(sympy.S(5)/3)) + 4*c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9*d**(sympy.S(5)/3)) + 3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2*sqrt(c + d*x**3)/(d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_317():
    f = x/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) - atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) + atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(18*c**(sympy.S(5)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_318():
    f = 1/(x**2*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(8*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(8*c**2*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(24*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(144*c**(sympy.S(11)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(144*c**(sympy.S(11)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(144*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_319():
    f = 1/(x**5*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(32*c**2*x**4) - d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(16*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + d*sqrt(c + d*x**3)/(16*c**3*x) + 3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(32*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(48*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1152*c**(sympy.S(17)/6)) - d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1152*c**(sympy.S(17)/6)) + d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1152*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_320():
    f = 1/(x**8*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(c + d*x**3)/(56*c**2*x**7) + 37*d*sqrt(c + d*x**3)/(1792*c**3*x**4) + 3*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(56*c**4*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3*d**2*sqrt(c + d*x**3)/(56*c**4*x) - 3*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(56*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9216*c**(sympy.S(23)/6)) - d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9216*c**(sympy.S(23)/6)) + d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9216*c**(sympy.S(23)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_321():
    f = x**3/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -d*x**3/c, d*x**3/(8*c))/(32*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_322():
    f = 1/(sqrt(c + d*x**3)*(8*c - d*x**3))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -d*x**3/c, d*x**3/(8*c))/(8*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_323():
    f = 1/(x**3*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, -d*x**3/c, d*x**3/(8*c))/(16*c*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_324():
    f = 1/(x**6*sqrt(c + d*x**3)*(8*c - d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-5)/3, sympy.S.Half, 1, sympy.S(-2)/3, -d*x**3/c, d*x**3/(8*c))/(40*c*x**5*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_325():
    f = x**11/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 1024*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**4) + 2*c**2/(27*d**4*sqrt(c + d*x**3)) - 4*c*sqrt(c + d*x**3)/d**4 - 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_326():
    f = x**8/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 128*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**3) - 2*c/(27*d**3*sqrt(c + d*x**3)) - 2*sqrt(c + d*x**3)/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_327():
    f = x**5/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2/(27*d**2*sqrt(c + d*x**3)) + 16*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*sqrt(c)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_328():
    f = x**2/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -2/(27*c*d*sqrt(c + d*x**3)) + 2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*c**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_329():
    f = 1/(x*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2/(27*c**2*sqrt(c + d*x**3)) + atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(324*c**(sympy.S(5)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(12*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_330():
    f = 1/(x**4*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -1/(24*c**2*x**3*sqrt(c + d*x**3)) - 25*d/(216*c**3*sqrt(c + d*x**3)) + d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2592*c**(sympy.S(7)/2)) + 11*d*atanh(sqrt(c + d*x**3)/sqrt(c))/(96*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_331():
    f = 1/(x**7*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -1/(48*c**2*x**6*sqrt(c + d*x**3)) + 3*d/(64*c**3*x**3*sqrt(c + d*x**3)) + 245*d**2/(1728*c**4*sqrt(c + d*x**3)) + d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(20736*c**(sympy.S(9)/2)) - 109*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(768*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_332():
    f = x**7/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -32*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(81*d**(sympy.S(8)/3)) - 32*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**(sympy.S(8)/3)) + 32*c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(81*d**(sympy.S(8)/3)) + 28*3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 56*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*x**2/(27*d**2*sqrt(c + d*x**3)) - 56*sqrt(c + d*x**3)/(27*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_333():
    f = x**4/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -2*x**2/(27*c*d*sqrt(c + d*x**3)) + 2*sqrt(c + d*x**3)/(27*c*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 4*sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3)) - 4*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3)) + 4*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_334():
    f = x/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2*x**2/(27*c**2*sqrt(c + d*x**3)) - 2*sqrt(c + d*x**3)/(27*c**2*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*c**(sympy.S(5)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*c**(sympy.S(5)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(162*c**(sympy.S(11)/6)*d**(sympy.S(2)/3)) - atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(162*c**(sympy.S(11)/6)*d**(sympy.S(2)/3)) + atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(162*c**(sympy.S(11)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_335():
    f = 1/(x**2*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2/(27*c**2*x*sqrt(c + d*x**3)) + 43*d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(216*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 43*sqrt(c + d*x**3)/(216*c**3*x) - 43*3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(432*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 43*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(648*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1296*c**(sympy.S(17)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1296*c**(sympy.S(17)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1296*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_336():
    f = 1/(x**5*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2/(27*c**2*x**4*sqrt(c + d*x**3)) - 91*sqrt(c + d*x**3)/(864*c**3*x**4) - 113*d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(432*c**4*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 113*d*sqrt(c + d*x**3)/(432*c**4*x) + 113*3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(864*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 113*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1296*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(10368*c**(sympy.S(23)/6)) - d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(10368*c**(sympy.S(23)/6)) + d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(10368*c**(sympy.S(23)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_337():
    f = 1/(x**8*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = 2/(27*c**2*x**7*sqrt(c + d*x**3)) - 139*sqrt(c + d*x**3)/(1512*c**3*x**7) + 6095*d*sqrt(c + d*x**3)/(48384*c**4*x**4) + 953*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(3024*c**5*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 953*d**2*sqrt(c + d*x**3)/(3024*c**5*x) - 953*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6048*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 953*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9072*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(82944*c**(sympy.S(29)/6)) - d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(82944*c**(sympy.S(29)/6)) + d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(82944*c**(sympy.S(29)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_338():
    f = x**3/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, 1, sympy.S(3)/2, sympy.S(7)/3, d*x**3/(8*c), -d*x**3/c)/(32*c**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_339():
    f = 1/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, 1, sympy.S(3)/2, sympy.S(4)/3, d*x**3/(8*c), -d*x**3/c)/(8*c**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_340():
    f = 1/(x**3*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, 1, sympy.S(3)/2, sympy.S(1)/3, d*x**3/(8*c), -d*x**3/c)/(16*c**2*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_341():
    f = 1/(x**6*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-5)/3, 1, sympy.S(3)/2, sympy.S(-2)/3, d*x**3/(8*c), -d*x**3/c)/(40*c**2*x**5*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_342():
    f = x*sqrt(a + b*x**3)/(a*(10 + 6*sqrt(3)) + b*x**3)
    F = sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(a + b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_343():
    f = x*sqrt(a - b*x**3)/(a*(10 + 6*sqrt(3)) - b*x**3)
    F = sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(a - b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 2*sqrt(a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_344():
    f = x*sqrt(-a + b*x**3)/(a*(-6*sqrt(3) - 10) + b*x**3)
    F = sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(4*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(-a + b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(4*b**(sympy.S(2)/3)) + 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) - 2*sqrt(-a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_345():
    f = x*sqrt(-a - b*x**3)/(a*(-6*sqrt(3) - 10) - b*x**3)
    F = sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(4*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(-a - b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(4*b**(sympy.S(2)/3)) + 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*sqrt(-a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_346():
    f = x*sqrt(a + b*x**3)/(a*(10 - 6*sqrt(3)) + b*x**3)
    F = -sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(2*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(a + b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_347():
    f = x*sqrt(a - b*x**3)/(a*(10 - 6*sqrt(3)) - b*x**3)
    F = -sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(2*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(a - b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 2*sqrt(a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_348():
    f = x*sqrt(-a + b*x**3)/(a*(10 - 6*sqrt(3)) - b*x**3)
    F = sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(-a + b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 2*sqrt(-a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_349():
    f = x*sqrt(-a - b*x**3)/(a*(10 - 6*sqrt(3)) + b*x**3)
    F = sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(-a - b*x**3)/(6*sqrt(a)))/(6*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/6)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(4*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(2*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(4*b**(sympy.S(2)/3)) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) + 2*sqrt(-a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_350():
    f = x/(sqrt(a + b*x**3)*(a*(10 + 6*sqrt(3)) + b*x**3))
    F = -sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(a + b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_351():
    f = x/(sqrt(a - b*x**3)*(a*(10 + 6*sqrt(3)) - b*x**3))
    F = -sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(a - b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_352():
    f = x/(sqrt(-a + b*x**3)*(a*(-6*sqrt(3) - 10) + b*x**3))
    F = sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(-a + b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_353():
    f = x/(sqrt(-a - b*x**3)*(a*(-6*sqrt(3) - 10) - b*x**3))
    F = sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(-a - b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_354():
    f = x/(sqrt(a + b*x**3)*(a*(10 - 6*sqrt(3)) + b*x**3))
    F = -sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(a + b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(a + b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_355():
    f = x/(sqrt(a - b*x**3)*(a*(10 - 6*sqrt(3)) - b*x**3))
    F = -sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(a - b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(a - b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_356():
    f = x/(sqrt(-a + b*x**3)*(a*(10 - 6*sqrt(3)) - b*x**3))
    F = -sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(-a + b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(2*sqrt(-a + b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_357():
    f = x/(sqrt(-a - b*x**3)*(a*(10 - 6*sqrt(3)) + b*x**3))
    F = -sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(-a - b*x**3)/(6*sqrt(a)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(12*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - 2*b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(18*a**(sympy.S(5)/6)*b**(sympy.S(2)/3)) - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*a**(sympy.S(1)/6)*(1 + sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(2*sqrt(-a - b*x**3)))/(36*a**(sympy.S(5)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_358():
    f = x**8*sqrt(c + d*x**3)/(a + b*x**3)
    F = 2*a**2*sqrt(c + d*x**3)/(3*b**3) - 2*a**2*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(7)/2)) + 2*(c + d*x**3)**(sympy.S(5)/2)/(15*b*d**2) - (c + d*x**3)**(sympy.S(3)/2)*(2*a*d + 2*b*c)/(9*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_359():
    f = x**5*sqrt(c + d*x**3)/(a + b*x**3)
    F = -2*a*sqrt(c + d*x**3)/(3*b**2) + 2*a*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2)) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_360():
    f = x**2*sqrt(c + d*x**3)/(a + b*x**3)
    F = 2*sqrt(c + d*x**3)/(3*b) - 2*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_361():
    f = sqrt(c + d*x**3)/(x*(a + b*x**3))
    F = -2*sqrt(c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a) + 2*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_362():
    f = sqrt(c + d*x**3)/(x**4*(a + b*x**3))
    F = -sqrt(c + d*x**3)/(3*a*x**3) - 2*sqrt(b)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2) + (-a*d + 2*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_363():
    f = x**3*sqrt(c + d*x**3)/(a + b*x**3)
    F = x**4*sqrt(c + d*x**3)*appellf1(sympy.S(4)/3, sympy.S(-1)/2, 1, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_364():
    f = x*sqrt(c + d*x**3)/(a + b*x**3)
    F = x**2*sqrt(c + d*x**3)*appellf1(sympy.S(2)/3, sympy.S(-1)/2, 1, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_365():
    f = sqrt(c + d*x**3)/(a + b*x**3)
    F = x*sqrt(c + d*x**3)*appellf1(sympy.S(1)/3, sympy.S(-1)/2, 1, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_366():
    f = sqrt(c + d*x**3)/(x**2*(a + b*x**3))
    F = -sqrt(c + d*x**3)*appellf1(sympy.S(-1)/3, sympy.S(-1)/2, 1, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a*x*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_367():
    f = sqrt(c + d*x**3)/(x**3*(a + b*x**3))
    F = -sqrt(c + d*x**3)*appellf1(sympy.S(-2)/3, sympy.S(-1)/2, 1, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a*x**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_368():
    f = x**8*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = 2*a**2*(c + d*x**3)**(sympy.S(3)/2)/(9*b**3) + 2*a**2*sqrt(c + d*x**3)*(-a*d + b*c)/(3*b**4) - 2*a**2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(9)/2)) + 2*(c + d*x**3)**(sympy.S(7)/2)/(21*b*d**2) - (c + d*x**3)**(sympy.S(5)/2)*(2*a*d + 2*b*c)/(15*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_369():
    f = x**5*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = -2*a*(c + d*x**3)**(sympy.S(3)/2)/(9*b**2) - 2*a*sqrt(c + d*x**3)*(-a*d + b*c)/(3*b**3) + 2*a*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(7)/2)) + 2*(c + d*x**3)**(sympy.S(5)/2)/(15*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_370():
    f = x**2*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = 2*(c + d*x**3)**(sympy.S(3)/2)/(9*b) + sqrt(c + d*x**3)*(-2*a*d + 2*b*c)/(3*b**2) - 2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_371():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x*(a + b*x**3))
    F = 2*d*sqrt(c + d*x**3)/(3*b) - 2*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a) + 2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_372():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**4*(a + b*x**3))
    F = -c*sqrt(c + d*x**3)/(3*a*x**3) + sqrt(c)*(-3*a*d + 2*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2) - 2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_373():
    f = x**3*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = c*x**4*sqrt(c + d*x**3)*appellf1(sympy.S(4)/3, sympy.S(-3)/2, 1, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_374():
    f = x*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = c*x**2*sqrt(c + d*x**3)*appellf1(sympy.S(2)/3, sympy.S(-3)/2, 1, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_375():
    f = (c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)
    F = c*x*sqrt(c + d*x**3)*appellf1(sympy.S(1)/3, sympy.S(-3)/2, 1, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_376():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**2*(a + b*x**3))
    F = -c*sqrt(c + d*x**3)*appellf1(sympy.S(-1)/3, sympy.S(-3)/2, 1, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a*x*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_377():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**3*(a + b*x**3))
    F = -c*sqrt(c + d*x**3)*appellf1(sympy.S(-2)/3, sympy.S(-3)/2, 1, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a*x**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_378():
    f = x**8/((a + b*x**3)*sqrt(c + d*x**3))
    F = -2*a**2*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*b*d**2) - sqrt(c + d*x**3)*(2*a*d + 2*b*c)/(3*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_379():
    f = x**5/((a + b*x**3)*sqrt(c + d*x**3))
    F = 2*a*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + 2*sqrt(c + d*x**3)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_380():
    f = x**2/((a + b*x**3)*sqrt(c + d*x**3))
    F = -2*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_381():
    f = 1/(x*(a + b*x**3)*sqrt(c + d*x**3))
    F = 2*sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a*sqrt(-a*d + b*c)) - 2*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_382():
    f = 1/(x**4*(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(c + d*x**3)/(3*a*c*x**3) - 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*sqrt(-a*d + b*c)) + (a*d + 2*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_383():
    f = x**3/((a + b*x**3)*sqrt(c + d*x**3))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_384():
    f = x/((a + b*x**3)*sqrt(c + d*x**3))
    F = x**2*sqrt(1 + d*x**3/c)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_385():
    f = 1/((a + b*x**3)*sqrt(c + d*x**3))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_386():
    f = 1/(x**2*(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a*x*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_387():
    f = 1/(x**3*(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_388():
    f = x**8/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -2*a**2*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + 2*c**2/(3*d**2*sqrt(c + d*x**3)*(-a*d + b*c)) + 2*sqrt(c + d*x**3)/(3*b*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_389():
    f = x**5/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = 2*a*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2)) - 2*c/(3*d*sqrt(c + d*x**3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_390():
    f = x**2/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -2*sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*(-a*d + b*c)**(sympy.S(3)/2)) + 2/(sqrt(c + d*x**3)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_391():
    f = 1/(x*(a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -2*d/(3*c*sqrt(c + d*x**3)*(-a*d + b*c)) + 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a*(-a*d + b*c)**(sympy.S(3)/2)) - 2*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_392():
    f = 1/(x**4*(a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -1/(3*a*c*x**3*sqrt(c + d*x**3)) - d*(-3*a*d + b*c)/(3*a*c**2*sqrt(c + d*x**3)*(-a*d + b*c)) - 2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*(-a*d + b*c)**(sympy.S(3)/2)) + (3*a*d + 2*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_393():
    f = x**3/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, 1, sympy.S(3)/2, sympy.S(7)/3, -b*x**3/a, -d*x**3/c)/(4*a*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_394():
    f = x/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = x**2*sqrt(1 + d*x**3/c)*appellf1(sympy.S(2)/3, 1, sympy.S(3)/2, sympy.S(5)/3, -b*x**3/a, -d*x**3/c)/(2*a*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_395():
    f = 1/((a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, 1, sympy.S(3)/2, sympy.S(4)/3, -b*x**3/a, -d*x**3/c)/(a*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_396():
    f = 1/(x**2*(a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-1)/3, 1, sympy.S(3)/2, sympy.S(2)/3, -b*x**3/a, -d*x**3/c)/(a*c*x*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_397():
    f = 1/(x**3*(a + b*x**3)*(c + d*x**3)**(sympy.S(3)/2))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, 1, sympy.S(3)/2, sympy.S(1)/3, -b*x**3/a, -d*x**3/c)/(2*a*c*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_398():
    f = x**11*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = -3968*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**4) + 2*c*sqrt(c + d*x**3)*(1146*c + 47*d*x**3)/(15*d**4) + x**9*sqrt(c + d*x**3)/(3*d*(8*c - d*x**3)) + 7*x**6*sqrt(c + d*x**3)/(15*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_399():
    f = x**8*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = -352*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**3) + 64*c*(c + d*x**3)**(sympy.S(3)/2)/(27*d**3*(8*c - d*x**3)) + 352*c*sqrt(c + d*x**3)/(27*d**3) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_400():
    f = x**5*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = -26*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**2) + 8*(c + d*x**3)**(sympy.S(3)/2)/(27*d**2*(8*c - d*x**3)) + 26*sqrt(c + d*x**3)/(27*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_401():
    f = x**2*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = sqrt(c + d*x**3)/(3*d*(8*c - d*x**3)) - atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_402():
    f = sqrt(c + d*x**3)/(x*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(24*c*(8*c - d*x**3)) + 5*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(288*c**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(96*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_403():
    f = sqrt(c + d*x**3)/(x**4*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(24*c*x**3*(8*c - d*x**3)) + d*sqrt(c + d*x**3)/(96*c**2*(8*c - d*x**3)) + 7*d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1152*c**(sympy.S(5)/2)) - d*atanh(sqrt(c + d*x**3)/sqrt(c))/(128*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_404():
    f = sqrt(c + d*x**3)/(x**7*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(48*c*x**6*(8*c - d*x**3)) - 7*d*sqrt(c + d*x**3)/(384*c**2*x**3*(8*c - d*x**3)) + 5*d**2*sqrt(c + d*x**3)/(1536*c**3*(8*c - d*x**3)) + 23*d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(18432*c**(sympy.S(7)/2)) - d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(2048*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_405():
    f = x**7*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = 76*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9*d**(sympy.S(8)/3)) + 76*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**(sympy.S(8)/3)) - 76*c**(sympy.S(7)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9*d**(sympy.S(8)/3)) - 373*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 746*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(63*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 746*c*sqrt(c + d*x**3)/(21*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + x**5*sqrt(c + d*x**3)/(3*d*(8*c - d*x**3)) + 13*x**2*sqrt(c + d*x**3)/(21*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_406():
    f = x**4*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = 5*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9*d**(sympy.S(5)/3)) + 5*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9*d**(sympy.S(5)/3)) - 5*c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9*d**(sympy.S(5)/3)) - 7*3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 7*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + x**2*sqrt(c + d*x**3)/(3*d*(8*c - d*x**3)) + 7*sqrt(c + d*x**3)/(3*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_407():
    f = x*sqrt(c + d*x**3)/(8*c - d*x**3)**2
    F = x**2*sqrt(c + d*x**3)/(24*c*(8*c - d*x**3)) + sqrt(c + d*x**3)/(24*c*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(48*c**(sympy.S(2)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(72*c**(sympy.S(2)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(144*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) + atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(144*c**(sympy.S(5)/6)*d**(sympy.S(2)/3)) - atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(144*c**(sympy.S(5)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_408():
    f = sqrt(c + d*x**3)/(x**2*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(24*c*x*(8*c - d*x**3)) + d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(48*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(48*c**2*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(96*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(144*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(144*c**(sympy.S(11)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(144*c**(sympy.S(11)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(144*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_409():
    f = sqrt(c + d*x**3)/(x**5*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(24*c*x**4*(8*c - d*x**3)) - 7*sqrt(c + d*x**3)/(768*c**2*x**4) + d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(96*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - d*sqrt(c + d*x**3)/(96*c**3*x) - 3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(192*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(288*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 17*sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(9216*c**(sympy.S(17)/6)) - 17*d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(9216*c**(sympy.S(17)/6)) + 17*d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(9216*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_410():
    f = sqrt(c + d*x**3)/(x**8*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(24*c*x**7*(8*c - d*x**3)) - 5*sqrt(c + d*x**3)/(672*c**2*x**7) - 53*d*sqrt(c + d*x**3)/(21504*c**3*x**4) + d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(5376*c**4*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - d**2*sqrt(c + d*x**3)/(5376*c**4*x) - 3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(10752*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16128*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 13*sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(36864*c**(sympy.S(23)/6)) - 13*d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(36864*c**(sympy.S(23)/6)) + 13*d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(36864*c**(sympy.S(23)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_411():
    f = x**11*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = -4992*c**(sympy.S(7)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**4 + 1664*c**3*sqrt(c + d*x**3)/d**4 + 2*c*(c + d*x**3)**(sympy.S(3)/2)*(694*c + 51*d*x**3)/(21*d**4) + x**9*(c + d*x**3)**(sympy.S(3)/2)/(3*d*(8*c - d*x**3)) + 3*x**6*(c + d*x**3)**(sympy.S(3)/2)/(7*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_412():
    f = x**8*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = -480*c**(sympy.S(5)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**3 + 160*c**2*sqrt(c + d*x**3)/d**3 + 64*c*(c + d*x**3)**(sympy.S(5)/2)/(27*d**3*(8*c - d*x**3)) + 160*c*(c + d*x**3)**(sympy.S(3)/2)/(27*d**3) + 2*(c + d*x**3)**(sympy.S(5)/2)/(15*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_413():
    f = x**5*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = -42*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**2 + 14*c*sqrt(c + d*x**3)/d**2 + 8*(c + d*x**3)**(sympy.S(5)/2)/(27*d**2*(8*c - d*x**3)) + 14*(c + d*x**3)**(sympy.S(3)/2)/(27*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_414():
    f = x**2*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = -3*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d + (c + d*x**3)**(sympy.S(3)/2)/(3*d*(8*c - d*x**3)) + sqrt(c + d*x**3)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_415():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x*(8*c - d*x**3)**2)
    F = 3*sqrt(c + d*x**3)/(64*c - 8*d*x**3) - 3*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(32*sqrt(c)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(96*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_416():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**4*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(24*x**3*(8*c - d*x**3)) + 5*d*sqrt(c + d*x**3)/(96*c*(8*c - d*x**3)) + 3*d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(128*c**(sympy.S(3)/2)) - 7*d*atanh(sqrt(c + d*x**3)/sqrt(c))/(384*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_417():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**7*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(48*x**6*(8*c - d*x**3)) - 23*d*sqrt(c + d*x**3)/(384*c*x**3*(8*c - d*x**3)) + 7*d**2*sqrt(c + d*x**3)/(512*c**2*(8*c - d*x**3)) + 15*d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2048*c**(sympy.S(5)/2)) - 17*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(2048*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_418():
    f = x**7*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = 108*sqrt(3)*c**(sympy.S(13)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(8)/3) + 108*c**(sympy.S(13)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(8)/3) - 108*c**(sympy.S(13)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(8)/3) - 2953*3**(sympy.S(1)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(13*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 5906*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(39*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 5906*c**2*sqrt(c + d*x**3)/(13*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 103*c*x**2*sqrt(c + d*x**3)/(13*d**2) + x**5*(c + d*x**3)**(sympy.S(3)/2)/(3*d*(8*c - d*x**3)) + 19*x**5*sqrt(c + d*x**3)/(39*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_419():
    f = x**4*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = 9*sqrt(3)*c**(sympy.S(7)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/d**(sympy.S(5)/3) + 9*c**(sympy.S(7)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/d**(sympy.S(5)/3) - 9*c**(sympy.S(7)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/d**(sympy.S(5)/3) - 265*3**(sympy.S(1)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(14*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 265*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 265*c*sqrt(c + d*x**3)/(7*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + x**2*(c + d*x**3)**(sympy.S(3)/2)/(3*d*(8*c - d*x**3)) + 13*x**2*sqrt(c + d*x**3)/(21*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_420():
    f = x*(c + d*x**3)**(sympy.S(3)/2)/(8*c - d*x**3)**2
    F = 9*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(16*d**(sympy.S(2)/3)) + 9*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(16*d**(sympy.S(2)/3)) - 9*c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(16*d**(sympy.S(2)/3)) - 19*3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 19*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(24*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 3*x**2*sqrt(c + d*x**3)/(64*c - 8*d*x**3) + 19*sqrt(c + d*x**3)/(8*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_421():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**2*(8*c - d*x**3)**2)
    F = 3*sqrt(c + d*x**3)/(8*x*(8*c - d*x**3)) + d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(16*c*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - sqrt(c + d*x**3)/(16*c*x) - 3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(32*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(48*c**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_422():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**5*(8*c - d*x**3)**2)
    F = 3*sqrt(c + d*x**3)/(8*x**4*(8*c - d*x**3)) - 13*sqrt(c + d*x**3)/(256*c*x**4) + d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(32*c**2*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - d*sqrt(c + d*x**3)/(32*c**2*x) - 3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(64*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(96*c**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 9*sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1024*c**(sympy.S(11)/6)) - 9*d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1024*c**(sympy.S(11)/6)) + 9*d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1024*c**(sympy.S(11)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_423():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**8*(8*c - d*x**3)**2)
    F = 3*sqrt(c + d*x**3)/(8*x**7*(8*c - d*x**3)) - 11*sqrt(c + d*x**3)/(224*c*x**7) - 83*d*sqrt(c + d*x**3)/(7168*c**2*x**4) + 19*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(1792*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 19*d**2*sqrt(c + d*x**3)/(1792*c**3*x) - 19*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3584*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 19*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5376*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 9*sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(4096*c**(sympy.S(17)/6)) - 9*d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(4096*c**(sympy.S(17)/6)) + 9*d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(4096*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_424():
    f = x**11/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -2944*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**4) + 8*x**6*sqrt(c + d*x**3)/(27*d**2*(8*c - d*x**3)) + 2*sqrt(c + d*x**3)*(170*c + 7*d*x**3)/(27*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_425():
    f = x**8/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -224*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**3) + 64*c*sqrt(c + d*x**3)/(27*d**3*(8*c - d*x**3)) + 2*sqrt(c + d*x**3)/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_426():
    f = x**5/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = 8*sqrt(c + d*x**3)/(27*d**2*(8*c - d*x**3)) - 10*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*sqrt(c)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_427():
    f = x**2/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(27*c*d*(8*c - d*x**3)) + atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*c**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_428():
    f = 1/(x*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(216*c**2*(8*c - d*x**3)) + 13*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(2592*c**(sympy.S(5)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(96*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_429():
    f = 1/(x**4*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(24*c**2*x**3*(8*c - d*x**3)) + 5*d*sqrt(c + d*x**3)/(864*c**3*(8*c - d*x**3)) + 11*d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(10368*c**(sympy.S(7)/2)) + d*atanh(sqrt(c + d*x**3)/sqrt(c))/(384*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_430():
    f = 1/(x**7*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -sqrt(c + d*x**3)/(48*c**2*x**6*(8*c - d*x**3)) + 3*d*sqrt(c + d*x**3)/(128*c**3*x**3*(8*c - d*x**3)) - 35*d**2*sqrt(c + d*x**3)/(13824*c**4*(8*c - d*x**3)) + 31*d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(165888*c**(sympy.S(9)/2)) - 19*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(6144*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_431():
    f = x**7/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = 44*sqrt(3)*c**(sympy.S(1)/6)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(81*d**(sympy.S(8)/3)) + 44*c**(sympy.S(1)/6)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*d**(sympy.S(8)/3)) - 44*c**(sympy.S(1)/6)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(81*d**(sympy.S(8)/3)) - 31*3**(sympy.S(1)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 62*sqrt(2)*3**(sympy.S(3)/4)*c**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 8*x**2*sqrt(c + d*x**3)/(27*d**2*(8*c - d*x**3)) + 62*sqrt(c + d*x**3)/(27*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_432():
    f = x**4/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = x**2*sqrt(c + d*x**3)/(27*c*d*(8*c - d*x**3)) + sqrt(c + d*x**3)/(27*c*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(54*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3)) + atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3)) - atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(81*c**(sympy.S(5)/6)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_433():
    f = x/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = x**2*sqrt(c + d*x**3)/(216*c**2*(8*c - d*x**3)) + sqrt(c + d*x**3)/(216*c**2*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(432*c**(sympy.S(5)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(648*c**(sympy.S(5)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 7*sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(1296*c**(sympy.S(11)/6)*d**(sympy.S(2)/3)) - 7*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(1296*c**(sympy.S(11)/6)*d**(sympy.S(2)/3)) + 7*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(1296*c**(sympy.S(11)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_434():
    f = 1/(x**2*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(216*c**2*x*(8*c - d*x**3)) + 7*d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(432*c**3*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 7*sqrt(c + d*x**3)/(432*c**3*x) - 7*3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(864*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 7*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1296*c**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(648*c**(sympy.S(17)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(648*c**(sympy.S(17)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(648*c**(sympy.S(17)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_435():
    f = 1/(x**5*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(216*c**2*x**4*(8*c - d*x**3)) - 31*sqrt(c + d*x**3)/(6912*c**3*x**4) - 5*d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(864*c**4*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 5*d*sqrt(c + d*x**3)/(864*c**4*x) + 5*3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1728*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 5*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2592*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 25*sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(82944*c**(sympy.S(23)/6)) - 25*d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(82944*c**(sympy.S(23)/6)) + 25*d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(82944*c**(sympy.S(23)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_436():
    f = 1/(x**8*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = sqrt(c + d*x**3)/(216*c**2*x**7*(8*c - d*x**3)) - 17*sqrt(c + d*x**3)/(6048*c**3*x**7) + 391*d*sqrt(c + d*x**3)/(193536*c**4*x**4) + 289*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(48384*c**5*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 289*d**2*sqrt(c + d*x**3)/(48384*c**5*x) - 289*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(96768*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 289*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(145152*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 17*sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(331776*c**(sympy.S(29)/6)) - 17*d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(331776*c**(sympy.S(29)/6)) + 17*d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(331776*c**(sympy.S(29)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_437():
    f = x**6/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = x**7*sqrt(1 + d*x**3/c)*appellf1(sympy.S(7)/3, sympy.S.Half, 2, sympy.S(10)/3, -d*x**3/c, d*x**3/(8*c))/(448*c**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_438():
    f = x**3/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, 2, sympy.S(7)/3, -d*x**3/c, d*x**3/(8*c))/(256*c**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_439():
    f = 1/(sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 2, sympy.S(4)/3, -d*x**3/c, d*x**3/(8*c))/(64*c**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_440():
    f = 1/(x**3*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 2, sympy.S(1)/3, -d*x**3/c, d*x**3/(8*c))/(128*c**2*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_441():
    f = 1/(x**6*sqrt(c + d*x**3)*(8*c - d*x**3)**2)
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-5)/3, sympy.S.Half, 2, sympy.S(-2)/3, -d*x**3/c, d*x**3/(8*c))/(320*c**2*x**5*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_442():
    f = x**11/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = -640*sqrt(c)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*d**4) + 8*x**6/(27*d**2*sqrt(c + d*x**3)*(8*c - d*x**3)) + (76*c + 78*d*x**3)/(81*d**4*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_443():
    f = x**8/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 64*c/(27*d**3*sqrt(c + d*x**3)*(8*c - d*x**3)) - 22/(81*d**3*sqrt(c + d*x**3)) - 32*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*sqrt(c)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_444():
    f = x**5/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 8/(27*d**2*sqrt(c + d*x**3)*(8*c - d*x**3)) - 2/(81*c*d**2*sqrt(c + d*x**3)) + 2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*c**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_445():
    f = x**2/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 1/(27*c*d*sqrt(c + d*x**3)*(8*c - d*x**3)) - 1/(81*c**2*d*sqrt(c + d*x**3)) + atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*c**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_446():
    f = 1/(x*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 1/(216*c**2*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5/(648*c**3*sqrt(c + d*x**3)) + 7*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(7776*c**(sympy.S(7)/2)) - atanh(sqrt(c + d*x**3)/sqrt(c))/(96*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_447():
    f = 1/(x**4*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = -1/(24*c**2*x**3*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5*d/(864*c**3*sqrt(c + d*x**3)*(8*c - d*x**3)) - 35*d/(2592*c**4*sqrt(c + d*x**3)) + 5*d*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(31104*c**(sympy.S(9)/2)) + 5*d*atanh(sqrt(c + d*x**3)/sqrt(c))/(384*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_448():
    f = 1/(x**7*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = -1/(48*c**2*x**6*sqrt(c + d*x**3)*(8*c - d*x**3)) + 17*d/(384*c**3*x**3*sqrt(c + d*x**3)*(8*c - d*x**3)) - 71*d**2/(13824*c**4*sqrt(c + d*x**3)*(8*c - d*x**3)) + 665*d**2/(41472*c**5*sqrt(c + d*x**3)) + 13*d**2*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(497664*c**(sympy.S(11)/2)) - 33*d**2*atanh(sqrt(c + d*x**3)/sqrt(c))/(2048*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_449():
    f = x**7/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 8*x**2/(27*d**2*sqrt(c + d*x**3)*(8*c - d*x**3)) - 2*x**2/(81*c*d**2*sqrt(c + d*x**3)) + 2*sqrt(c + d*x**3)/(81*c*d**(sympy.S(8)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*c**(sympy.S(2)/3)*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(243*c**(sympy.S(2)/3)*d**(sympy.S(8)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 4*sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(243*c**(sympy.S(5)/6)*d**(sympy.S(8)/3)) + 4*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*c**(sympy.S(5)/6)*d**(sympy.S(8)/3)) - 4*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(243*c**(sympy.S(5)/6)*d**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_450():
    f = x**4/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = x**2/(27*c*d*sqrt(c + d*x**3)*(8*c - d*x**3)) - x**2/(81*c**2*d*sqrt(c + d*x**3)) + sqrt(c + d*x**3)/(81*c**2*d**(sympy.S(5)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(162*c**(sympy.S(5)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(243*c**(sympy.S(5)/3)*d**(sympy.S(5)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(243*c**(sympy.S(11)/6)*d**(sympy.S(5)/3)) - atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(243*c**(sympy.S(11)/6)*d**(sympy.S(5)/3)) + atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(243*c**(sympy.S(11)/6)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_451():
    f = x/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = x**2/(216*c**2*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5*x**2/(648*c**3*sqrt(c + d*x**3)) - 5*sqrt(c + d*x**3)/(648*c**3*d**(sympy.S(2)/3)*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 5*3**(sympy.S(1)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1296*c**(sympy.S(8)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(5*c**(sympy.S(1)/3) + 5*d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1944*c**(sympy.S(8)/3)*d**(sympy.S(2)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 5*sqrt(3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(3888*c**(sympy.S(17)/6)*d**(sympy.S(2)/3)) - 5*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(3888*c**(sympy.S(17)/6)*d**(sympy.S(2)/3)) + 5*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(3888*c**(sympy.S(17)/6)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_452():
    f = 1/(x**2*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 1/(216*c**2*x*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5/(648*c**3*x*sqrt(c + d*x**3)) + 31*d**(sympy.S(1)/3)*sqrt(c + d*x**3)/(1296*c**4*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 31*sqrt(c + d*x**3)/(1296*c**4*x) - 31*3**(sympy.S(1)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2592*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 31*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(1)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3888*c**(sympy.S(11)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - sqrt(3)*d**(sympy.S(1)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(3888*c**(sympy.S(23)/6)) - d**(sympy.S(1)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(3888*c**(sympy.S(23)/6)) + d**(sympy.S(1)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(3888*c**(sympy.S(23)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_453():
    f = 1/(x**5*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 1/(216*c**2*x**4*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5/(648*c**3*x**4*sqrt(c + d*x**3)) - 253*sqrt(c + d*x**3)/(20736*c**4*x**4) - 77*d**(sympy.S(4)/3)*sqrt(c + d*x**3)/(2592*c**5*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) + 77*d*sqrt(c + d*x**3)/(2592*c**5*x) + 77*3**(sympy.S(1)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5184*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 77*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(4)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7776*c**(sympy.S(14)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 11*sqrt(3)*d**(sympy.S(4)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(248832*c**(sympy.S(29)/6)) - 11*d**(sympy.S(4)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(248832*c**(sympy.S(29)/6)) + 11*d**(sympy.S(4)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(248832*c**(sympy.S(29)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_454():
    f = 1/(x**8*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 1/(216*c**2*x**7*sqrt(c + d*x**3)*(8*c - d*x**3)) + 5/(648*c**3*x**7*sqrt(c + d*x**3)) - 191*sqrt(c + d*x**3)/(18144*c**4*x**7) + 8257*d*sqrt(c + d*x**3)/(580608*c**5*x**4) + 5179*d**(sympy.S(7)/3)*sqrt(c + d*x**3)/(145152*c**6*(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)) - 5179*d**2*sqrt(c + d*x**3)/(145152*c**6*x) - 5179*3**(sympy.S(1)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_e(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(290304*c**(sympy.S(17)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) + 5179*sqrt(2)*3**(sympy.S(3)/4)*d**(sympy.S(7)/3)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(435456*c**(sympy.S(17)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3)) - 7*sqrt(3)*d**(sympy.S(7)/3)*atan(sqrt(3)*c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/sqrt(c + d*x**3))/(995328*c**(sympy.S(35)/6)) - 7*d**(sympy.S(7)/3)*atanh(sqrt(c + d*x**3)/(3*sqrt(c)))/(995328*c**(sympy.S(35)/6)) + 7*d**(sympy.S(7)/3)*atanh((c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)**2/(3*c**(sympy.S(1)/6)*sqrt(c + d*x**3)))/(995328*c**(sympy.S(35)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_455():
    f = x**6/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = 2*x*(4*c + d*x**3)/(81*c*d**2*sqrt(c + d*x**3)*(8*c - d*x**3)) - 2*3**(sympy.S(3)/4)*sqrt((c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)*elliptic_f(asin((c**(sympy.S(1)/3)*(1 - sqrt(3)) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(243*c*d**(sympy.S(7)/3)*sqrt(c**(sympy.S(1)/3)*(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(c**(sympy.S(1)/3)*(1 + sqrt(3)) + d**(sympy.S(1)/3)*x)**2)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_456():
    f = x**3/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S(3)/2, 2, sympy.S(7)/3, -d*x**3/c, d*x**3/(8*c))/(256*c**3*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_457():
    f = 1/((c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S(3)/2, 2, sympy.S(4)/3, -d*x**3/c, d*x**3/(8*c))/(64*c**3*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_458():
    f = 1/(x**3*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S(3)/2, 2, sympy.S(1)/3, -d*x**3/c, d*x**3/(8*c))/(128*c**3*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_459():
    f = 1/(x**6*(c + d*x**3)**(sympy.S(3)/2)*(8*c - d*x**3)**2)
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-5)/3, sympy.S(3)/2, 2, sympy.S(-2)/3, -d*x**3/c, d*x**3/(8*c))/(320*c**3*x**5*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_460():
    f = x**8*sqrt(c + d*x**3)/(a + b*x**3)**2
    F = -a**2*(c + d*x**3)**(sympy.S(3)/2)/(3*b**2*(a + b*x**3)*(-a*d + b*c)) - a*sqrt(c + d*x**3)*(-5*a*d + 4*b*c)/(3*b**3*(-a*d + b*c)) + a*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(7)/2)*sqrt(-a*d + b*c)) + 2*(c + d*x**3)**(sympy.S(3)/2)/(9*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_461():
    f = x**5*sqrt(c + d*x**3)/(a + b*x**3)**2
    F = a*(c + d*x**3)**(sympy.S(3)/2)/(3*b*(a + b*x**3)*(-a*d + b*c)) + sqrt(c + d*x**3)*(-3*a*d + 2*b*c)/(3*b**2*(-a*d + b*c)) - (-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_462():
    f = x**2*sqrt(c + d*x**3)/(a + b*x**3)**2
    F = -sqrt(c + d*x**3)/(3*b*(a + b*x**3)) - d*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_463():
    f = sqrt(c + d*x**3)/(x*(a + b*x**3)**2)
    F = sqrt(c + d*x**3)/(3*a*(a + b*x**3)) - 2*sqrt(c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2) + (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_464():
    f = sqrt(c + d*x**3)/(x**4*(a + b*x**3)**2)
    F = -sqrt(c + d*x**3)/(3*a*x**3*(a + b*x**3)) - 2*b*sqrt(c + d*x**3)/(3*a**2*(a + b*x**3)) - sqrt(b)*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**3*sqrt(-a*d + b*c)) + (-a*d + 4*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**3*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_465():
    f = x**3*sqrt(c + d*x**3)/(a + b*x**3)**2
    F = x**4*sqrt(c + d*x**3)*appellf1(sympy.S(4)/3, sympy.S(-1)/2, 2, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_466():
    f = x*sqrt(c + d*x**3)/(a + b*x**3)**2
    F = x**2*sqrt(c + d*x**3)*appellf1(sympy.S(2)/3, sympy.S(-1)/2, 2, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_467():
    f = sqrt(c + d*x**3)/(a + b*x**3)**2
    F = x*sqrt(c + d*x**3)*appellf1(sympy.S(1)/3, sympy.S(-1)/2, 2, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_468():
    f = sqrt(c + d*x**3)/(x**2*(a + b*x**3)**2)
    F = -sqrt(c + d*x**3)*appellf1(sympy.S(-1)/3, sympy.S(-1)/2, 2, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a**2*x*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_469():
    f = sqrt(c + d*x**3)/(x**3*(a + b*x**3)**2)
    F = -sqrt(c + d*x**3)*appellf1(sympy.S(-2)/3, sympy.S(-1)/2, 2, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*x**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_470():
    f = x**8*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = -a**2*(c + d*x**3)**(sympy.S(5)/2)/(3*b**2*(a + b*x**3)*(-a*d + b*c)) - a*(c + d*x**3)**(sympy.S(3)/2)*(-7*a*d + 4*b*c)/(9*b**3*(-a*d + b*c)) - a*sqrt(c + d*x**3)*(-7*a*d + 4*b*c)/(3*b**4) + a*(-7*a*d + 4*b*c)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(9)/2)) + 2*(c + d*x**3)**(sympy.S(5)/2)/(15*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_471():
    f = x**5*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = a*(c + d*x**3)**(sympy.S(5)/2)/(3*b*(a + b*x**3)*(-a*d + b*c)) + (c + d*x**3)**(sympy.S(3)/2)*(-5*a*d + 2*b*c)/(9*b**2*(-a*d + b*c)) + sqrt(c + d*x**3)*(-5*a*d + 2*b*c)/(3*b**3) - (-5*a*d + 2*b*c)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_472():
    f = x**2*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = -(c + d*x**3)**(sympy.S(3)/2)/(3*b*(a + b*x**3)) + d*sqrt(c + d*x**3)/b**2 - d*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_473():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x*(a + b*x**3)**2)
    F = sqrt(c + d*x**3)*(-a*d + b*c)/(3*a*b*(a + b*x**3)) - 2*c**(sympy.S(3)/2)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2) + sqrt(-a*d + b*c)*(a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_474():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**4*(a + b*x**3)**2)
    F = -c*sqrt(c + d*x**3)/(3*a*x**3*(a + b*x**3)) - sqrt(c + d*x**3)*(-a*d + 2*b*c)/(3*a**2*(a + b*x**3)) + sqrt(c)*(-3*a*d + 4*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**3) - sqrt(-a*d + b*c)*(-a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_475():
    f = x**3*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = c*x**4*sqrt(c + d*x**3)*appellf1(sympy.S(4)/3, sympy.S(-3)/2, 2, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_476():
    f = x*(c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = c*x**2*sqrt(c + d*x**3)*appellf1(sympy.S(2)/3, sympy.S(-3)/2, 2, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_477():
    f = (c + d*x**3)**(sympy.S(3)/2)/(a + b*x**3)**2
    F = c*x*sqrt(c + d*x**3)*appellf1(sympy.S(1)/3, sympy.S(-3)/2, 2, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_478():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**2*(a + b*x**3)**2)
    F = -c*sqrt(c + d*x**3)*appellf1(sympy.S(-1)/3, sympy.S(-3)/2, 2, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a**2*x*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_479():
    f = (c + d*x**3)**(sympy.S(3)/2)/(x**3*(a + b*x**3)**2)
    F = -c*sqrt(c + d*x**3)*appellf1(sympy.S(-2)/3, sympy.S(-3)/2, 2, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*x**2*sqrt(1 + d*x**3/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_480():
    f = x**8/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = -a**2*sqrt(c + d*x**3)/(3*b**2*(a + b*x**3)*(-a*d + b*c)) + a*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + 2*sqrt(c + d*x**3)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_481():
    f = x**5/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = a*sqrt(c + d*x**3)/(3*b*(a + b*x**3)*(-a*d + b*c)) - (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_482():
    f = x**2/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = -sqrt(c + d*x**3)/((a + b*x**3)*(-3*a*d + 3*b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_483():
    f = 1/(x*(a + b*x**3)**2*sqrt(c + d*x**3))
    F = b*sqrt(c + d*x**3)/(3*a*(a + b*x**3)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*(-a*d + b*c)**(sympy.S(3)/2)) - 2*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_484():
    f = 1/(x**4*(a + b*x**3)**2*sqrt(c + d*x**3))
    F = -sqrt(c + d*x**3)/(3*a*c*x**3*(a + b*x**3)) - b*sqrt(c + d*x**3)*(-a*d + 2*b*c)/(3*a**2*c*(a + b*x**3)*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**3*(-a*d + b*c)**(sympy.S(3)/2)) + (a*d + 4*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_485():
    f = x**3/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, 2, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_486():
    f = x/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = x**2*sqrt(1 + d*x**3/c)*appellf1(sympy.S(2)/3, sympy.S.Half, 2, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_487():
    f = 1/((a + b*x**3)**2*sqrt(c + d*x**3))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 2, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_488():
    f = 1/(x**2*(a + b*x**3)**2*sqrt(c + d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-1)/3, sympy.S.Half, 2, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a**2*x*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_489():
    f = 1/(x**3*(a + b*x**3)**2*sqrt(c + d*x**3))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 2, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_490():
    f = x**8/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = -a**2/(3*b**2*(a + b*x**3)*sqrt(c + d*x**3)*(-a*d + b*c)) + a*(-a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2)) + (-a**2*d**2 - 2*b**2*c**2)/(3*b**2*d*sqrt(c + d*x**3)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_491():
    f = x**5/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = a/(3*b*(a + b*x**3)*sqrt(c + d*x**3)*(-a*d + b*c)) + (a*d + 2*b*c)/(3*b*sqrt(c + d*x**3)*(-a*d + b*c)**2) - (a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*sqrt(b)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_492():
    f = x**2/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = sqrt(b)*d*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(-a*d + b*c)**(sympy.S(5)/2) - d/(sqrt(c + d*x**3)*(-a*d + b*c)**2) - 1/((a + b*x**3)*sqrt(c + d*x**3)*(-3*a*d + 3*b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_493():
    f = 1/(x*(a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = b/(3*a*(a + b*x**3)*sqrt(c + d*x**3)*(-a*d + b*c)) + d*(2*a*d + b*c)/(3*a*c*sqrt(c + d*x**3)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*(-5*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**2*(-a*d + b*c)**(sympy.S(5)/2)) - 2*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_494():
    f = 1/(x**4*(a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = -1/(3*a*c*x**3*(a + b*x**3)*sqrt(c + d*x**3)) - b*(-a*d + 2*b*c)/(3*a**2*c*(a + b*x**3)*sqrt(c + d*x**3)*(-a*d + b*c)) - d*(3*a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)/(3*a**2*c**2*sqrt(c + d*x**3)*(-a*d + b*c)**2) - b**(sympy.S(5)/2)*(-7*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**3)/sqrt(-a*d + b*c))/(3*a**3*(-a*d + b*c)**(sympy.S(5)/2)) + (3*a*d + 4*b*c)*atanh(sqrt(c + d*x**3)/sqrt(c))/(3*a**3*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_495():
    f = x**3/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = x**4*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S(3)/2, 2, sympy.S(7)/3, -d*x**3/c, -b*x**3/a)/(4*a**2*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_496():
    f = x/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = x**2*sqrt(1 + d*x**3/c)*appellf1(sympy.S(2)/3, sympy.S(3)/2, 2, sympy.S(5)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_497():
    f = 1/((a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = x*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S(3)/2, 2, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(a**2*c*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_498():
    f = 1/(x**2*(a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-1)/3, sympy.S(3)/2, 2, sympy.S(2)/3, -d*x**3/c, -b*x**3/a)/(a**2*c*x*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_499():
    f = 1/(x**3*(a + b*x**3)**2*(c + d*x**3)**(sympy.S(3)/2))
    F = -sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S(3)/2, 2, sympy.S(1)/3, -d*x**3/c, -b*x**3/a)/(2*a**2*c*x**2*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_500():
    f = (e*x)**m*(A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)
    F = 2*B*(e*x)**(m + 1)*(a + b*x**3)**(sympy.S(7)/2)/(b*e*(2*m + 23)) - a**2*(e*x)**(m + 1)*sqrt(a + b*x**3)*(-A*b*(2*m + 23) + 2*B*a*(m + 1))*hyper((sympy.S(-5)/2, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(b*e*sqrt(1 + b*x**3/a)*(m + 1)*(2*m + 23))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_501():
    f = (e*x)**m*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*(e*x)**(m + 1)*(a + b*x**3)**(sympy.S(5)/2)/(b*e*(2*m + 17)) - a*(e*x)**(m + 1)*sqrt(a + b*x**3)*(-A*b*(2*m + 17) + 2*B*a*(m + 1))*hyper((sympy.S(-3)/2, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(b*e*sqrt(1 + b*x**3/a)*(m + 1)*(2*m + 17))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_502():
    f = (e*x)**m*(A + B*x**3)*sqrt(a + b*x**3)
    F = 2*B*(e*x)**(m + 1)*(a + b*x**3)**(sympy.S(3)/2)/(b*e*(2*m + 11)) - (e*x)**(m + 1)*sqrt(a + b*x**3)*(-A*b*(2*m + 11) + 2*B*a*(m + 1))*hyper((sympy.S(-1)/2, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(b*e*sqrt(1 + b*x**3/a)*(m + 1)*(2*m + 11))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_503():
    f = (e*x)**m*(A + B*x**3)/sqrt(a + b*x**3)
    F = 2*B*(e*x)**(m + 1)*sqrt(a + b*x**3)/(b*e*(2*m + 5)) - (e*x)**(m + 1)*sqrt(1 + b*x**3/a)*(-A*b*(2*m + 5) + 2*B*a*(m + 1))*hyper((sympy.S.Half, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(b*e*sqrt(a + b*x**3)*(m + 1)*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_504():
    f = (e*x)**m*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = (e*x)**(m + 1)*sqrt(1 + b*x**3/a)*(A*(-2*b*m + b) + 2*B*a*(m + 1))*hyper((sympy.S.Half, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(3*a*b*e*sqrt(a + b*x**3)*(m + 1)) + (e*x)**(m + 1)*(2*A*b - 2*B*a)/(3*a*b*e*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_505():
    f = (e*x)**m*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = (e*x)**(m + 1)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2)) + (e*x)**(m + 1)*sqrt(1 + b*x**3/a)*(A*b*(7 - 2*m) + 2*B*a*(m + 1))*hyper((sympy.S(3)/2, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(9*a**2*b*e*sqrt(a + b*x**3)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_506():
    f = x**5/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = sqrt(a + b*x**3)*sqrt(c + d*x**3)/(3*b*d) - (a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**3)/(sqrt(b)*sqrt(c + d*x**3)))/(3*b**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_507():
    f = x**2/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = 2*atanh(sqrt(d)*sqrt(a + b*x**3)/(sqrt(b)*sqrt(c + d*x**3)))/(3*sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_508():
    f = 1/(x*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = -2*atanh(sqrt(c)*sqrt(a + b*x**3)/(sqrt(a)*sqrt(c + d*x**3)))/(3*sqrt(a)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_509():
    f = 1/(x**4*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(a + b*x**3)*sqrt(c + d*x**3)/(3*a*c*x**3) + (a*d + b*c)*atanh(sqrt(c)*sqrt(a + b*x**3)/(sqrt(a)*sqrt(c + d*x**3)))/(3*a**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_510():
    f = x**4/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = x**5*sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(5)/3, sympy.S.Half, sympy.S.Half, sympy.S(8)/3, -b*x**3/a, -d*x**3/c)/(5*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_511():
    f = x**3/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = x**4*sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(4)/3, sympy.S.Half, sympy.S.Half, sympy.S(7)/3, -b*x**3/a, -d*x**3/c)/(4*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_512():
    f = x/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = x**2*sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(2)/3, sympy.S.Half, sympy.S.Half, sympy.S(5)/3, -b*x**3/a, -d*x**3/c)/(2*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_513():
    f = 1/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = x*sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(1)/3, sympy.S.Half, sympy.S.Half, sympy.S(4)/3, -b*x**3/a, -d*x**3/c)/(sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_514():
    f = 1/(x**2*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(-1)/3, sympy.S.Half, sympy.S.Half, sympy.S(2)/3, -b*x**3/a, -d*x**3/c)/(x*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_515():
    f = 1/(x**3*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    F = -sqrt(1 + b*x**3/a)*sqrt(1 + d*x**3/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, sympy.S.Half, sympy.S(1)/3, -b*x**3/a, -d*x**3/c)/(2*x**2*sqrt(a + b*x**3)*sqrt(c + d*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_516():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)*sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(3)/2)/(9*b*e) - a**2*e**(sympy.S(7)/2)*(2*A*b - B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(24*b**(sympy.S(5)/2)) + a*e**2*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(2*A*b - B*a)/(24*b**2) + (e*x)**(sympy.S(9)/2)*sqrt(a + b*x**3)*(2*A*b - B*a)/(12*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_517():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)*sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(3)/2)/(8*b*e) - 3**(sympy.S(3)/4)*a**(sympy.S(5)/3)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(16*A*b - 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(640*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 3*a*e**2*sqrt(e*x)*sqrt(a + b*x**3)*(16*A*b - 7*B*a)/(320*b**2) + (e*x)**(sympy.S(7)/2)*sqrt(a + b*x**3)*(16*A*b - 7*B*a)/(80*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_518():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)*sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)/(7*b*e) - 3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*A*b - 5*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(112*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*A*b - 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(224*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + a*e*sqrt(e*x)*(3 + 3*sqrt(3))*sqrt(a + b*x**3)*(14*A*b - 5*B*a)/(112*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + (e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(14*A*b - 5*B*a)/(56*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_519():
    f = sqrt(e*x)*(A + B*x**3)*sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)/(6*b*e) + a*sqrt(e)*(4*A*b - B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(12*b**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(4*A*b - B*a)/(12*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_520():
    f = (A + B*x**3)*sqrt(a + b*x**3)/sqrt(e*x)
    F = B*sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)/(5*b*e) + 3**(sympy.S(3)/4)*a**(sympy.S(2)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(10*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(40*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + sqrt(e*x)*sqrt(a + b*x**3)*(10*A*b - B*a)/(20*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_521():
    f = (A + B*x**3)*sqrt(a + b*x**3)/(e*x)**(sympy.S(3)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(a*e*sqrt(e*x)) - 3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b + B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(16*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + sqrt(e*x)*(3 + 3*sqrt(3))*sqrt(a + b*x**3)*(8*A*b + B*a)/(8*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + (e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(8*A*b + B*a)/(4*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_522():
    f = (A + B*x**3)*sqrt(a + b*x**3)/(e*x)**(sympy.S(5)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(3*a*e*(e*x)**(sympy.S(3)/2)) + (2*A*b + B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*sqrt(b)*e**(sympy.S(5)/2)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(2*A*b + B*a)/(3*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_523():
    f = (A + B*x**3)*sqrt(a + b*x**3)/(e*x)**(sympy.S(7)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(5*a*e*(e*x)**(sympy.S(5)/2)) + sqrt(e*x)*sqrt(a + b*x**3)*(4*A*b + 5*B*a)/(10*a*e**4) + 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b + 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(20*a**(sympy.S(1)/3)*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_524():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**(sympy.S(9)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(7*a*x**(sympy.S(7)/2)) + b**(sympy.S(1)/3)*sqrt(x)*(3 + 3*sqrt(3))*sqrt(a + b*x**3)*(2*A*b + 7*B*a)/(7*a*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) - sqrt(a + b*x**3)*(4*A*b + 14*B*a)/(7*a*sqrt(x)) - 3*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt(x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + 7*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(7*a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt(x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(14*a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_525():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**(sympy.S(11)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(9*a*x**(sympy.S(9)/2)) + 2*B*sqrt(b)*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a + b*x**3))/3 - 2*B*sqrt(a + b*x**3)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_526():
    f = (A + B*x**3)*sqrt(a + b*x**3)/x**(sympy.S(13)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(3)/2)/(11*a*x**(sympy.S(11)/2)) + sqrt(a + b*x**3)*(4*A*b - 22*B*a)/(55*a*x**(sympy.S(5)/2)) - 3**(sympy.S(3)/4)*b*sqrt(x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b - 11*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(55*a**(sympy.S(4)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_527():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(5)/2)/(12*b*e) - a**3*e**(sympy.S(7)/2)*(8*A*b - 3*B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(192*b**(sympy.S(5)/2)) + a**2*e**2*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(8*A*b - 3*B*a)/(192*b**2) + a*(e*x)**(sympy.S(9)/2)*sqrt(a + b*x**3)*(8*A*b - 3*B*a)/(96*b*e) + (e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(3)/2)*(8*A*b - 3*B*a)/(72*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_528():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(5)/2)/(11*b*e) - 9*3**(sympy.S(3)/4)*a**(sympy.S(8)/3)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(22*A*b - 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(14080*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 27*a**2*e**2*sqrt(e*x)*sqrt(a + b*x**3)*(22*A*b - 7*B*a)/(7040*b**2) + 9*a*(e*x)**(sympy.S(7)/2)*sqrt(a + b*x**3)*(22*A*b - 7*B*a)/(1760*b*e) + (e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(3)/2)*(22*A*b - 7*B*a)/(176*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_529():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(5)/2)/(10*b*e) - 27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(448*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 9*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(896*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + a**2*e*sqrt(e*x)*(27 + 27*sqrt(3))*sqrt(a + b*x**3)*(4*A*b - B*a)/(448*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + 9*a*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(4*A*b - B*a)/(224*b*e) + (e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)*(4*A*b - B*a)/(28*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_530():
    f = sqrt(e*x)*(A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(5)/2)/(9*b*e) + a**2*sqrt(e)*(6*A*b - B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(24*b**(sympy.S(3)/2)) + a*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(6*A*b - B*a)/(24*b*e) + (e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)*(6*A*b - B*a)/(36*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_531():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/sqrt(e*x)
    F = B*sqrt(e*x)*(a + b*x**3)**(sympy.S(5)/2)/(8*b*e) + 9*3**(sympy.S(3)/4)*a**(sympy.S(5)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(16*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(640*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 9*a*sqrt(e*x)*sqrt(a + b*x**3)*(16*A*b - B*a)/(320*b*e) + sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)*(16*A*b - B*a)/(80*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_532():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/(e*x)**(sympy.S(3)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(5)/2)/(a*e*sqrt(e*x)) - 27*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*A*b + B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(112*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 9*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(224*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + a*sqrt(e*x)*(27 + 27*sqrt(3))*sqrt(a + b*x**3)*(14*A*b + B*a)/(112*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + (e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(126*A*b + 9*B*a)/(56*e**4) + (e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)*(14*A*b + B*a)/(7*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_533():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/(e*x)**(sympy.S(5)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(5)/2)/(3*a*e*(e*x)**(sympy.S(3)/2)) + a*(4*A*b + B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(4*sqrt(b)*e**(sympy.S(5)/2)) + (e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(4*A*b + B*a)/(4*e**4) + (e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)*(4*A*b + B*a)/(6*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_534():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(3)/2)/(e*x)**(sympy.S(7)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(5)/2)/(5*a*e*(e*x)**(sympy.S(5)/2)) + 9*3**(sympy.S(3)/4)*a**(sympy.S(2)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(40*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + sqrt(e*x)*sqrt(a + b*x**3)*(18*A*b + 9*B*a)/(20*e**4) + sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)*(2*A*b + B*a)/(5*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_535():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)
    F = B*(e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(7)/2)/(15*b*e) - a**4*e**(sympy.S(7)/2)*(10*A*b - 3*B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(384*b**(sympy.S(5)/2)) + a**3*e**2*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(10*A*b - 3*B*a)/(384*b**2) + a**2*(e*x)**(sympy.S(9)/2)*sqrt(a + b*x**3)*(10*A*b - 3*B*a)/(192*b*e) + a*(e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(3)/2)*(10*A*b - 3*B*a)/(144*b*e) + (e*x)**(sympy.S(9)/2)*(a + b*x**3)**(sympy.S(5)/2)*(10*A*b - 3*B*a)/(120*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_536():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)
    F = B*(e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(7)/2)/(14*b*e) - 27*3**(sympy.S(3)/4)*a**(sympy.S(11)/3)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(11264*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 81*a**3*e**2*sqrt(e*x)*sqrt(a + b*x**3)*(4*A*b - B*a)/(5632*b**2) + 27*a**2*(e*x)**(sympy.S(7)/2)*sqrt(a + b*x**3)*(4*A*b - B*a)/(1408*b*e) + 15*a*(e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(3)/2)*(4*A*b - B*a)/(704*b*e) + (e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(5)/2)*(4*A*b - B*a)/(44*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_537():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)
    F = B*(e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(7)/2)/(13*b*e) - 81*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(26*A*b - 5*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(11648*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 27*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(26*A*b - 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(23296*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + a**3*e*sqrt(e*x)*(81 + 81*sqrt(3))*sqrt(a + b*x**3)*(26*A*b - 5*B*a)/(11648*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + 27*a**2*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(26*A*b - 5*B*a)/(5824*b*e) + 3*a*(e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)*(26*A*b - 5*B*a)/(728*b*e) + (e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(5)/2)*(26*A*b - 5*B*a)/(260*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_538():
    f = sqrt(e*x)*(A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)
    F = B*(e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(7)/2)/(12*b*e) + 5*a**3*sqrt(e)*(8*A*b - B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(192*b**(sympy.S(3)/2)) + 5*a**2*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(8*A*b - B*a)/(192*b*e) + 5*a*(e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)*(8*A*b - B*a)/(288*b*e) + (e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(5)/2)*(8*A*b - B*a)/(72*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_539():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)/sqrt(e*x)
    F = B*sqrt(e*x)*(a + b*x**3)**(sympy.S(7)/2)/(11*b*e) + 27*3**(sympy.S(3)/4)*a**(sympy.S(8)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(22*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(2816*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 27*a**2*sqrt(e*x)*sqrt(a + b*x**3)*(22*A*b - B*a)/(1408*b*e) + 3*a*sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)*(22*A*b - B*a)/(352*b*e) + sqrt(e*x)*(a + b*x**3)**(sympy.S(5)/2)*(22*A*b - B*a)/(176*b*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_540():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)/(e*x)**(sympy.S(3)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(7)/2)/(a*e*sqrt(e*x)) - 81*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(20*A*b + B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(448*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 27*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(20*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(896*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + a**2*sqrt(e*x)*(81 + 81*sqrt(3))*sqrt(a + b*x**3)*(20*A*b + B*a)/(448*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + 27*a*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)*(20*A*b + B*a)/(224*e**4) + (e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)*(60*A*b + 3*B*a)/(28*e**4) + (e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(5)/2)*(20*A*b + B*a)/(10*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_541():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)/(e*x)**(sympy.S(5)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(7)/2)/(3*a*e*(e*x)**(sympy.S(3)/2)) + 5*a**2*(6*A*b + B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(24*sqrt(b)*e**(sympy.S(5)/2)) + 5*a*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(6*A*b + B*a)/(24*e**4) + (e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)*(30*A*b + 5*B*a)/(36*e**4) + (e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(5)/2)*(6*A*b + B*a)/(9*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_542():
    f = (A + B*x**3)*(a + b*x**3)**(sympy.S(5)/2)/(e*x)**(sympy.S(7)/2)
    F = -2*A*(a + b*x**3)**(sympy.S(7)/2)/(5*a*e*(e*x)**(sympy.S(5)/2)) + 27*3**(sympy.S(3)/4)*a**(sympy.S(5)/3)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(16*A*b + 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(640*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 27*a*sqrt(e*x)*sqrt(a + b*x**3)*(16*A*b + 5*B*a)/(320*e**4) + sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)*(48*A*b + 15*B*a)/(80*e**4) + sqrt(e*x)*(a + b*x**3)**(sympy.S(5)/2)*(16*A*b + 5*B*a)/(40*a*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_543():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)/sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(9)/2)*sqrt(a + b*x**3)/(6*b*e) - a*e**(sympy.S(7)/2)*(4*A*b - 3*B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(12*b**(sympy.S(5)/2)) + e**2*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)*(4*A*b - 3*B*a)/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_544():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)/sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(7)/2)*sqrt(a + b*x**3)/(5*b*e) - 3**(sympy.S(3)/4)*a**(sympy.S(2)/3)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(10*A*b - 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(120*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + e**2*sqrt(e*x)*sqrt(a + b*x**3)*(10*A*b - 7*B*a)/(20*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_545():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)/sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)/(4*b*e) - 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b - 5*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b - 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(48*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + e*sqrt(e*x)*(1 + sqrt(3))*sqrt(a + b*x**3)*(8*A*b - 5*B*a)/(8*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_546():
    f = sqrt(e*x)*(A + B*x**3)/sqrt(a + b*x**3)
    F = B*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)/(3*b*e) + sqrt(e)*(2*A*b - B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_547():
    f = (A + B*x**3)/(sqrt(e*x)*sqrt(a + b*x**3))
    F = B*sqrt(e*x)*sqrt(a + b*x**3)/(2*b*e) + 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(12*a**(sympy.S(1)/3)*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_548():
    f = (A + B*x**3)/((e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3))
    F = -2*A*sqrt(a + b*x**3)/(a*e*sqrt(e*x)) + sqrt(e*x)*(1 + sqrt(3))*sqrt(a + b*x**3)*(2*A*b + B*a)/(a*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) - 3**(sympy.S(1)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(6*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_549():
    f = (A + B*x**3)/((e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3))
    F = -2*A*sqrt(a + b*x**3)/(3*a*e*(e*x)**(sympy.S(3)/2)) + 2*B*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*sqrt(b)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_550():
    f = (A + B*x**3)/((e*x)**(sympy.S(7)/2)*sqrt(a + b*x**3))
    F = -2*A*sqrt(a + b*x**3)/(5*a*e*(e*x)**(sympy.S(5)/2)) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b - 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(15*a**(sympy.S(4)/3)*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_551():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(9)/2)/(3*b*e*sqrt(a + b*x**3)) - e**2*(e*x)**(sympy.S(3)/2)*(2*A*b - 3*B*a)/(3*b**2*sqrt(a + b*x**3)) + e**(sympy.S(7)/2)*(2*A*b - 3*B*a)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_552():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = B*(e*x)**(sympy.S(7)/2)/(2*b*e*sqrt(a + b*x**3)) - e**2*sqrt(e*x)*(4*A*b - 7*B*a)/(6*b**2*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(36*a**(sympy.S(1)/3)*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_553():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = (e*x)**(sympy.S(5)/2)*(2*A*b - 2*B*a)/(3*a*b*e*sqrt(a + b*x**3)) - e*sqrt(e*x)*(1 + sqrt(3))*sqrt(a + b*x**3)*(2*A*b - 5*B*a)/(3*a*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + 3**(sympy.S(1)/4)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b - 5*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b - 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(18*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_554():
    f = sqrt(e*x)*(A + B*x**3)/(a + b*x**3)**(sympy.S(3)/2)
    F = 2*B*sqrt(e)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*b**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(2*A*b - 2*B*a)/(3*a*b*e*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_555():
    f = (A + B*x**3)/(sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2))
    F = sqrt(e*x)*(2*A*b - 2*B*a)/(3*a*b*e*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(9*a**(sympy.S(4)/3)*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_556():
    f = (A + B*x**3)/((e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2))
    F = -2*A/(a*e*sqrt(e*x)*sqrt(a + b*x**3)) - (e*x)**(sympy.S(5)/2)*(8*A*b - 2*B*a)/(3*a**2*e**4*sqrt(a + b*x**3)) + sqrt(e*x)*(2 + 2*sqrt(3))*sqrt(a + b*x**3)*(4*A*b - B*a)/(3*a**2*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(1)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b - 2*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_557():
    f = (A + B*x**3)/((e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2))
    F = -2*A/(3*a*e*(e*x)**(sympy.S(3)/2)*sqrt(a + b*x**3)) - (e*x)**(sympy.S(3)/2)*(4*A*b - 2*B*a)/(3*a**2*e**4*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_558():
    f = (A + B*x**3)/((e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(3)/2))
    F = -2*A/(5*a*e*(e*x)**(sympy.S(5)/2)*sqrt(a + b*x**3)) - sqrt(e*x)*(16*A*b - 10*B*a)/(15*a**2*e**4*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(16*A*b - 10*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(45*a**(sympy.S(7)/3)*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_559():
    f = (e*x)**(sympy.S(7)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = -2*B*e**2*(e*x)**(sympy.S(3)/2)/(3*b**2*sqrt(a + b*x**3)) + 2*B*e**(sympy.S(7)/2)*atanh(sqrt(b)*(e*x)**(sympy.S(3)/2)/(e**(sympy.S(3)/2)*sqrt(a + b*x**3)))/(3*b**(sympy.S(5)/2)) + (e*x)**(sympy.S(9)/2)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_560():
    f = (e*x)**(sympy.S(5)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(7)/2)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2)) - e**2*sqrt(e*x)*(4*A*b + 14*B*a)/(27*a*b**2*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*e**2*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*A*b + 7*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(81*a**(sympy.S(4)/3)*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_561():
    f = (e*x)**(sympy.S(3)/2)*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(5)/2)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2)) + (e*x)**(sympy.S(5)/2)*(8*A*b + 10*B*a)/(27*a**2*b*e*sqrt(a + b*x**3)) - e*sqrt(e*x)*(2 + 2*sqrt(3))*sqrt(a + b*x**3)*(4*A*b + 5*B*a)/(27*a**2*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) + 3**(sympy.S(3)/4)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*A*b + 5*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(81*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) + 3**(sympy.S(1)/4)*e*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*A*b + 10*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(27*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_562():
    f = sqrt(e*x)*(A + B*x**3)/(a + b*x**3)**(sympy.S(5)/2)
    F = (e*x)**(sympy.S(3)/2)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2)) + (e*x)**(sympy.S(3)/2)*(4*A*b + 2*B*a)/(9*a**2*b*e*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_563():
    f = (A + B*x**3)/(sqrt(e*x)*(a + b*x**3)**(sympy.S(5)/2))
    F = sqrt(e*x)*(2*A*b - 2*B*a)/(9*a*b*e*(a + b*x**3)**(sympy.S(3)/2)) + sqrt(e*x)*(16*A*b + 2*B*a)/(27*a**2*b*e*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(16*A*b + 2*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(81*a**(sympy.S(7)/3)*b*e*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_564():
    f = (A + B*x**3)/((e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(5)/2))
    F = -2*A/(a*e*sqrt(e*x)*(a + b*x**3)**(sympy.S(3)/2)) - (e*x)**(sympy.S(5)/2)*(20*A*b - 2*B*a)/(9*a**2*e**4*(a + b*x**3)**(sympy.S(3)/2)) - (e*x)**(sympy.S(5)/2)*(80*A*b - 8*B*a)/(27*a**3*e**4*sqrt(a + b*x**3)) + sqrt(e*x)*(8 + 8*sqrt(3))*sqrt(a + b*x**3)*(10*A*b - B*a)/(27*a**3*b**(sympy.S(2)/3)*e**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(4 - 4*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(10*A*b - B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(81*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3)) - 3**(sympy.S(1)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(80*A*b - 8*B*a)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)*e**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_565():
    f = (A + B*x**3)/((e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(5)/2))
    F = -2*A/(3*a*e*(e*x)**(sympy.S(3)/2)*(a + b*x**3)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(8*A*b - 2*B*a)/(9*a**2*e**4*(a + b*x**3)**(sympy.S(3)/2)) - (e*x)**(sympy.S(3)/2)*(16*A*b - 4*B*a)/(9*a**3*e**4*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_566():
    f = (A + B*x**3)/((e*x)**(sympy.S(7)/2)*(a + b*x**3)**(sympy.S(5)/2))
    F = -2*A/(5*a*e*(e*x)**(sympy.S(5)/2)*(a + b*x**3)**(sympy.S(3)/2)) - sqrt(e*x)*(28*A*b - 10*B*a)/(45*a**2*e**4*(a + b*x**3)**(sympy.S(3)/2)) - sqrt(e*x)*(224*A*b - 80*B*a)/(135*a**3*e**4*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt(e*x)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(224*A*b - 80*B*a)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(405*a**(sympy.S(10)/3)*e**4*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_567():
    f = x**14/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = (1 - x**3)**(sympy.S(11)/3)/11 - (1 - x**3)**(sympy.S(8)/3)/4 + 2*(1 - x**3)**(sympy.S(5)/3)/5 - 2**(sympy.S(2)/3)*log(x**3 + 1)/12 + 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_568():
    f = x**11/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -(1 - x**3)**(sympy.S(8)/3)/8 + (1 - x**3)**(sympy.S(5)/3)/5 - (1 - x**3)**(sympy.S(2)/3)/2 + 2**(sympy.S(2)/3)*log(x**3 + 1)/12 - 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_569():
    f = x**8/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = (1 - x**3)**(sympy.S(5)/3)/5 - 2**(sympy.S(2)/3)*log(x**3 + 1)/12 + 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_570():
    f = x**5/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -(1 - x**3)**(sympy.S(2)/3)/2 + 2**(sympy.S(2)/3)*log(x**3 + 1)/12 - 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_571():
    f = x**2/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -2**(sympy.S(2)/3)*log(x**3 + 1)/12 + 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_572():
    f = 1/(x*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -log(x)/2 + log(1 - (1 - x**3)**(sympy.S(1)/3))/2 + 2**(sympy.S(2)/3)*log(x**3 + 1)/12 - 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_573():
    f = 1/(x**4*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = log(x)/3 - log(1 - (1 - x**3)**(sympy.S(1)/3))/3 - 2**(sympy.S(2)/3)*log(x**3 + 1)/12 + 2**(sympy.S(2)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - 2*sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/9 - (1 - x**3)**(sympy.S(2)/3)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_574():
    f = x**6/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -x*(1 - x**3)**(sympy.S(2)/3)/3 - 2*log(x/(1 - x**3)**(sympy.S(1)/3) + 1)/9 + 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + log(x**2/(1 - x**3)**(sympy.S(2)/3) - x/(1 - x**3)**(sympy.S(1)/3) + 1)/9 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2*sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/9 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_575():
    f = x**3/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = log(x/(1 - x**3)**(sympy.S(1)/3) + 1)/3 - 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - log(x**2/(1 - x**3)**(sympy.S(2)/3) - x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_576():
    f = 1/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_577():
    f = 1/(x**3*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - (1 - x**3)**(sympy.S(2)/3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_578():
    f = 1/(x**6*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - (1 - x**3)**(sympy.S(5)/3)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_579():
    f = 1/(x**9*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - (1 - x**3)**(sympy.S(2)/3)/(2*x**2) - (1 - x**3)**(sympy.S(5)/3)/(5*x**5) - (1 - x**3)**(sympy.S(8)/3)/(8*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_580():
    f = x**7/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -x**2*(1 - x**3)**(sympy.S(2)/3)/4 - x**2*hyper((sympy.S(1)/3, sympy.S(2)/3), (sympy.S(5)/3,), x**3)/4 - 2**(sympy.S(2)/3)*log(-(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/12 - 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log((1 - x)**2/(1 - x**3)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2*2**(sympy.S(1)/3))/24 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_581():
    f = x**4/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = x**2*hyper((sympy.S(1)/3, sympy.S(2)/3), (sympy.S(5)/3,), x**3)/2 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3) + (x - 1)/(1 - x**3)**(sympy.S(1)/3))/12 + 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(2)/3)*log((1 - x)**2/(1 - x**3)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2*2**(sympy.S(1)/3))/24 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_582():
    f = x/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3) + (x - 1)/(1 - x**3)**(sympy.S(1)/3))/12 - 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log((1 - x)**2/(1 - x**3)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2*2**(sympy.S(1)/3))/24 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_583():
    f = 1/(x**2*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -x**2*hyper((sympy.S(1)/3, sympy.S(2)/3), (sympy.S(5)/3,), x**3)/2 + 2**(sympy.S(2)/3)*log(-(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/12 + 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(2)/3)*log((1 - x)**2/(1 - x**3)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2*2**(sympy.S(1)/3))/24 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/12 - (1 - x**3)**(sympy.S(2)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_584():
    f = 1/(x**5*(1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = x**2*hyper((sympy.S(1)/3, sympy.S(2)/3), (sympy.S(5)/3,), x**3)/4 - 2**(sympy.S(2)/3)*log(-(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/12 - 2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(2)/3)*log((1 - x)**2/(1 - x**3)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 2*2**(sympy.S(1)/3))/24 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/12 + (1 - x**3)**(sympy.S(2)/3)/(2*x) - (1 - x**3)**(sympy.S(2)/3)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_585():
    f = x**11/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -(1 - x**3)**(sympy.S(7)/3)/7 + (1 - x**3)**(sympy.S(4)/3)/4 - (1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*log(x**3 + 1)/12 - 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_586():
    f = x**8/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = (1 - x**3)**(sympy.S(4)/3)/4 - 2**(sympy.S(1)/3)*log(x**3 + 1)/12 + 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_587():
    f = x**5/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*log(x**3 + 1)/12 - 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_588():
    f = x**2/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -2**(sympy.S(1)/3)*log(x**3 + 1)/12 + 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_589():
    f = 1/(x*(1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -log(x)/2 + log(1 - (1 - x**3)**(sympy.S(1)/3))/2 + 2**(sympy.S(1)/3)*log(x**3 + 1)/12 - 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_590():
    f = 1/(x**4*(1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = log(x)/6 - log(1 - (1 - x**3)**(sympy.S(1)/3))/6 - 2**(sympy.S(1)/3)*log(x**3 + 1)/12 + 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/9 - (1 - x**3)**(sympy.S(1)/3)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_591():
    f = x**4/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -log(x/(1 - x**3)**(sympy.S(1)/3) + 1)/3 + 2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + log(x**2/(1 - x**3)**(sympy.S(2)/3) - x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_592():
    f = x/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_593():
    f = 1/(x**2*(1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = 2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 - 2**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - (1 - x**3)**(sympy.S(1)/3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_594():
    f = 1/(x**5*(1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/6 + 2**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*x**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/12 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 - (1 - x**3)**(sympy.S(4)/3)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_595():
    f = x**6/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = x**7*appellf1(sympy.S(7)/3, sympy.S(2)/3, 1, sympy.S(10)/3, x**3, -x**3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_596():
    f = x**3/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = x**4*appellf1(sympy.S(4)/3, sympy.S(2)/3, 1, sympy.S(7)/3, x**3, -x**3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_597():
    f = 1/((1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = x*appellf1(sympy.S(1)/3, sympy.S(2)/3, 1, sympy.S(4)/3, x**3, -x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_598():
    f = 1/(x**3*(1 - x**3)**(sympy.S(2)/3)*(x**3 + 1))
    F = -appellf1(sympy.S(-2)/3, sympy.S(2)/3, 1, sympy.S(1)/3, x**3, -x**3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_599():
    f = x**15/((a + b*x**4)*(c + d*x**4))
    F = -a**3*log(a + b*x**4)/(4*b**3*(-a*d + b*c)) + c**3*log(c + d*x**4)/(4*d**3*(-a*d + b*c)) + x**8/(8*b*d) - x**4*(a*d + b*c)/(4*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_600():
    f = x**11/((a + b*x**4)*(c + d*x**4))
    F = a**2*log(a + b*x**4)/(4*b**2*(-a*d + b*c)) - c**2*log(c + d*x**4)/(4*d**2*(-a*d + b*c)) + x**4/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_601():
    f = x**7/((a + b*x**4)*(c + d*x**4))
    F = -a*log(a + b*x**4)/(4*b*(-a*d + b*c)) + c*log(c + d*x**4)/(4*d*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_602():
    f = x**3/((a + b*x**4)*(c + d*x**4))
    F = log(a + b*x**4)/(-4*a*d + 4*b*c) - log(c + d*x**4)/(-4*a*d + 4*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_603():
    f = 1/(x*(a + b*x**4)*(c + d*x**4))
    F = d*log(c + d*x**4)/(4*c*(-a*d + b*c)) - b*log(a + b*x**4)/(4*a*(-a*d + b*c)) + log(x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_604():
    f = 1/(x**5*(a + b*x**4)*(c + d*x**4))
    F = -d**2*log(c + d*x**4)/(4*c**2*(-a*d + b*c)) - 1/(4*a*c*x**4) + b**2*log(a + b*x**4)/(4*a**2*(-a*d + b*c)) - (a*d + b*c)*log(x)/(a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_605():
    f = x**13/((a + b*x**4)*(c + d*x**4))
    F = -a**(sympy.S(5)/2)*atan(sqrt(b)*x**2/sqrt(a))/(2*b**(sympy.S(5)/2)*(-a*d + b*c)) + c**(sympy.S(5)/2)*atan(sqrt(d)*x**2/sqrt(c))/(2*d**(sympy.S(5)/2)*(-a*d + b*c)) + x**6/(6*b*d) - x**2*(a*d + b*c)/(2*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_606():
    f = x**9/((a + b*x**4)*(c + d*x**4))
    F = a**(sympy.S(3)/2)*atan(sqrt(b)*x**2/sqrt(a))/(2*b**(sympy.S(3)/2)*(-a*d + b*c)) - c**(sympy.S(3)/2)*atan(sqrt(d)*x**2/sqrt(c))/(2*d**(sympy.S(3)/2)*(-a*d + b*c)) + x**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_607():
    f = x**5/((a + b*x**4)*(c + d*x**4))
    F = -sqrt(a)*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(b)*(-a*d + b*c)) + sqrt(c)*atan(sqrt(d)*x**2/sqrt(c))/(2*sqrt(d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_608():
    f = x/((a + b*x**4)*(c + d*x**4))
    F = -sqrt(d)*atan(sqrt(d)*x**2/sqrt(c))/(2*sqrt(c)*(-a*d + b*c)) + sqrt(b)*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_609():
    f = 1/(x**3*(a + b*x**4)*(c + d*x**4))
    F = d**(sympy.S(3)/2)*atan(sqrt(d)*x**2/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)) - 1/(2*a*c*x**2) - b**(sympy.S(3)/2)*atan(sqrt(b)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_610():
    f = 1/(x**7*(a + b*x**4)*(c + d*x**4))
    F = -d**(sympy.S(5)/2)*atan(sqrt(d)*x**2/sqrt(c))/(2*c**(sympy.S(5)/2)*(-a*d + b*c)) - 1/(6*a*c*x**6) + (a*d + b*c)/(2*a**2*c**2*x**2) + b**(sympy.S(5)/2)*atan(sqrt(b)*x**2/sqrt(a))/(2*a**(sympy.S(5)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_611():
    f = x**8/((a + b*x**4)*(c + d*x**4))
    F = -sqrt(2)*a**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(5)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(5)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(5)/4)*(-a*d + b*c)) + x/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_612():
    f = x**6/((a + b*x**4)*(c + d*x**4))
    F = -sqrt(2)*a**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_613():
    f = x**4/((a + b*x**4)*(c + d*x**4))
    F = sqrt(2)*a**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*a**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*a**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(1)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(1)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*d**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*c**(sympy.S(1)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*c**(sympy.S(1)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*d**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_614():
    f = x**2/((a + b*x**4)*(c + d*x**4))
    F = -sqrt(2)*d**(sympy.S(1)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(1)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(1)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(1)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(1)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_615():
    f = 1/((a + b*x**4)*(c + d*x**4))
    F = sqrt(2)*d**(sympy.S(3)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(3)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_616():
    f = 1/(x**2*(a + b*x**4)*(c + d*x**4))
    F = sqrt(2)*d**(sympy.S(5)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(5)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(5)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(5)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-a*d + b*c)) - 1/(a*c*x) - sqrt(2)*b**(sympy.S(5)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(5)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(5)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(5)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(5)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(5)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(5)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_617():
    f = 1/(x**4*(a + b*x**4)*(c + d*x**4))
    F = -sqrt(2)*d**(sympy.S(7)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(7)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(7)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(7)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-a*d + b*c)) - 1/(3*a*c*x**3) + sqrt(2)*b**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(7)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(7)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(7)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_618():
    f = 1/(x**6*(a + b*x**4)*(c + d*x**4))
    F = -sqrt(2)*d**(sympy.S(9)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(9)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(9)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(9)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(9)/4)*(-a*d + b*c)) - 1/(5*a*c*x**5) + (a*d + b*c)/(a**2*c**2*x) + sqrt(2)*b**(sympy.S(9)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(9)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(9)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(9)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(9)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(9)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(9)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_619():
    f = x**7*sqrt(c + d*x**4)/(a + b*x**4)
    F = -a*sqrt(c + d*x**4)/(2*b**2) + a*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*b**(sympy.S(5)/2)) + (c + d*x**4)**(sympy.S(3)/2)/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_620():
    f = x**5*sqrt(c + d*x**4)/(a + b*x**4)
    F = -sqrt(a)*sqrt(-a*d + b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*b**2) + x**2*sqrt(c + d*x**4)/(4*b) + (-2*a*d + b*c)*atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(4*b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_621():
    f = x**3*sqrt(c + d*x**4)/(a + b*x**4)
    F = sqrt(c + d*x**4)/(2*b) - sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_622():
    f = x*sqrt(c + d*x**4)/(a + b*x**4)
    F = sqrt(d)*atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(2*b) + sqrt(-a*d + b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*sqrt(a)*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_623():
    f = sqrt(c + d*x**4)/(x*(a + b*x**4))
    F = -sqrt(c)*atanh(sqrt(c + d*x**4)/sqrt(c))/(2*a) + sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*a*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_624():
    f = sqrt(c + d*x**4)/(x**3*(a + b*x**4))
    F = -sqrt(c + d*x**4)/(2*a*x**2) - sqrt(-a*d + b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_625():
    f = sqrt(c + d*x**4)/(x**5*(a + b*x**4))
    F = -sqrt(c + d*x**4)/(4*a*x**4) - sqrt(b)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*a**2) + (-a*d + 2*b*c)*atanh(sqrt(c + d*x**4)/sqrt(c))/(4*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_626():
    f = sqrt(c + d*x**4)/(x**7*(a + b*x**4))
    F = -sqrt(c + d*x**4)/(6*a*x**6) + sqrt(c + d*x**4)*(-a*d + 3*b*c)/(6*a**2*c*x**2) + b*sqrt(-a*d + b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_627():
    f = (e*x)**(sympy.S(3)/2)*sqrt(c + d*x**4)/(a + b*x**4)
    F = 2*(e*x)**(sympy.S(5)/2)*sqrt(c + d*x**4)*appellf1(sympy.S(5)/8, sympy.S(-1)/2, 1, sympy.S(13)/8, -d*x**4/c, -b*x**4/a)/(5*a*e*sqrt(1 + d*x**4/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_628():
    f = sqrt(e*x)*sqrt(c + d*x**4)/(a + b*x**4)
    F = 2*(e*x)**(sympy.S(3)/2)*sqrt(c + d*x**4)*appellf1(sympy.S(3)/8, sympy.S(-1)/2, 1, sympy.S(11)/8, -d*x**4/c, -b*x**4/a)/(3*a*e*sqrt(1 + d*x**4/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_629():
    f = sqrt(c + d*x**4)/(sqrt(e*x)*(a + b*x**4))
    F = 2*sqrt(e*x)*sqrt(c + d*x**4)*appellf1(sympy.S(1)/8, sympy.S(-1)/2, 1, sympy.S(9)/8, -d*x**4/c, -b*x**4/a)/(a*e*sqrt(1 + d*x**4/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_630():
    f = sqrt(c + d*x**4)/((e*x)**(sympy.S(3)/2)*(a + b*x**4))
    F = -2*sqrt(c + d*x**4)*appellf1(sympy.S(-1)/8, sympy.S(-1)/2, 1, sympy.S(7)/8, -d*x**4/c, -b*x**4/a)/(a*e*sqrt(e*x)*sqrt(1 + d*x**4/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_631():
    f = x**11/((a + b*x**4)*sqrt(c + d*x**4))
    F = -a**2*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*b**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + (c + d*x**4)**(sympy.S(3)/2)/(6*b*d**2) - sqrt(c + d*x**4)*(a*d + b*c)/(2*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_632():
    f = x**7/((a + b*x**4)*sqrt(c + d*x**4))
    F = a*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*b**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + sqrt(c + d*x**4)/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_633():
    f = x**3/((a + b*x**4)*sqrt(c + d*x**4))
    F = -atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_634():
    f = 1/(x*(a + b*x**4)*sqrt(c + d*x**4))
    F = sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*a*sqrt(-a*d + b*c)) - atanh(sqrt(c + d*x**4)/sqrt(c))/(2*a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_635():
    f = 1/(x**5*(a + b*x**4)*sqrt(c + d*x**4))
    F = -sqrt(c + d*x**4)/(4*a*c*x**4) - b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(2*a**2*sqrt(-a*d + b*c)) + (a*d + 2*b*c)*atanh(sqrt(c + d*x**4)/sqrt(c))/(4*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_636():
    f = x**9/((a + b*x**4)*sqrt(c + d*x**4))
    F = a**(sympy.S(3)/2)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*b**2*sqrt(-a*d + b*c)) + x**2*sqrt(c + d*x**4)/(4*b*d) - (2*a*d + b*c)*atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(4*b**2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_637():
    f = x**5/((a + b*x**4)*sqrt(c + d*x**4))
    F = -sqrt(a)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*b*sqrt(-a*d + b*c)) + atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(2*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_638():
    f = x/((a + b*x**4)*sqrt(c + d*x**4))
    F = atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*sqrt(a)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_639():
    f = 1/(x**3*(a + b*x**4)*sqrt(c + d*x**4))
    F = -sqrt(c + d*x**4)/(2*a*c*x**2) - b*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*a**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_640():
    f = 1/(x**7*(a + b*x**4)*sqrt(c + d*x**4))
    F = -sqrt(c + d*x**4)/(6*a*c*x**6) + sqrt(c + d*x**4)*(2*a*d + 3*b*c)/(6*a**2*c**2*x**2) + b**2*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(2*a**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_641():
    f = x**8/((a + b*x**4)*sqrt(c + d*x**4))
    F = ((x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((Symbol('a') * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Symbol('a') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((Symbol('a') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_642():
    f = x**15/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = -a**2*(-5*a*d + 6*b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + a*x**8*sqrt(c + d*x**4)/(4*b*(a + b*x**4)*(-a*d + b*c)) - sqrt(c + d*x**4)*(-15*a**2*d**2 + 8*a*b*c*d + 4*b**2*c**2 - b*d*x**4*(-5*a*d + 2*b*c))/(12*b**3*d**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_643():
    f = x**11/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = -a**2*sqrt(c + d*x**4)/(4*b**2*(a + b*x**4)*(-a*d + b*c)) + a*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + sqrt(c + d*x**4)/(2*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_644():
    f = x**7/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = a*sqrt(c + d*x**4)/(4*b*(a + b*x**4)*(-a*d + b*c)) - (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_645():
    f = x**3/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = -sqrt(c + d*x**4)/((a + b*x**4)*(-4*a*d + 4*b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_646():
    f = 1/(x*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = b*sqrt(c + d*x**4)/(4*a*(a + b*x**4)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*a**2*(-a*d + b*c)**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**4)/sqrt(c))/(2*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_647():
    f = 1/(x**5*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = -sqrt(c + d*x**4)/(4*a*c*x**4*(a + b*x**4)) - b*sqrt(c + d*x**4)*(-a*d + 2*b*c)/(4*a**2*c*(a + b*x**4)*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**4)/sqrt(-a*d + b*c))/(4*a**3*(-a*d + b*c)**(sympy.S(3)/2)) + (a*d + 4*b*c)*atanh(sqrt(c + d*x**4)/sqrt(c))/(4*a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_648():
    f = x**13/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = a**(sympy.S(3)/2)*(-4*a*d + 5*b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*b**3*(-a*d + b*c)**(sympy.S(3)/2)) + a*x**6*sqrt(c + d*x**4)/(4*b*(a + b*x**4)*(-a*d + b*c)) + x**2*sqrt(c + d*x**4)*(-2*a*d + b*c)/(4*b**2*d*(-a*d + b*c)) - (4*a*d + b*c)*atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(4*b**3*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_649():
    f = x**9/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = -sqrt(a)*(-2*a*d + 3*b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*b**2*(-a*d + b*c)**(sympy.S(3)/2)) + a*x**2*sqrt(c + d*x**4)/(4*b*(a + b*x**4)*(-a*d + b*c)) + atanh(sqrt(d)*x**2/sqrt(c + d*x**4))/(2*b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_650():
    f = x**5/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = -x**2*sqrt(c + d*x**4)/((a + b*x**4)*(-4*a*d + 4*b*c)) + c*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_651():
    f = x/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = b*x**2*sqrt(c + d*x**4)/(4*a*(a + b*x**4)*(-a*d + b*c)) + (-2*a*d + b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_652():
    f = 1/(x**3*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = b*sqrt(c + d*x**4)/(4*a*x**2*(a + b*x**4)*(-a*d + b*c)) - sqrt(c + d*x**4)*(-2*a*d + 3*b*c)/(4*a**2*c*x**2*(-a*d + b*c)) - b*(-4*a*d + 3*b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_653():
    f = 1/(x**7*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = b*sqrt(c + d*x**4)/(4*a*x**6*(a + b*x**4)*(-a*d + b*c)) - sqrt(c + d*x**4)*(-2*a*d + 5*b*c)/(12*a**2*c*x**6*(-a*d + b*c)) + sqrt(c + d*x**4)*(-4*a**2*d**2 - 8*a*b*c*d + 15*b**2*c**2)/(12*a**3*c**2*x**2*(-a*d + b*c)) + b**2*(-6*a*d + 5*b*c)*atan(x**2*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**4)))/(4*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_654():
    f = x**8/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = ((Symbol('a') * x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_655():
    f = x**4/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = (Integer(-1) * ((x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_656():
    f = 1/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = ((Symbol('b') * x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_657():
    f = 1/(x**4*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = (Integer(-1) * ((((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(11) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(11) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Symbol('b') * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_658():
    f = x**6/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = ((sympy.sqrt(Symbol('d')) * x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1)))) + ((((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * sympy.sqrt((Integer(-1) * Symbol('a'))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * sympy.sqrt((Integer(-1) * Symbol('a'))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_659():
    f = x**2/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('d')) * x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1)))) + ((Symbol('b') * (x)**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_660():
    f = 1/(x**2*(a + b*x**4)**2*sqrt(c + d*x**4))
    F = (Integer(-1) * ((((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * x * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * (Symbol('a') + (Symbol('b') * (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(9) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_661():
    f = (e*x)**m*(a + b*x**4)**2/sqrt(c + d*x**4)
    F = b**2*(e*x)**(m + 5)*sqrt(c + d*x**4)/(d*e**5*(m + 7)) - b*(e*x)**(m + 1)*sqrt(c + d*x**4)*(-2*a*d*(m + 7) + b*c*(m + 5))/(d**2*e*(m + 3)*(m + 7)) + (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*(a**2*d**2*(m + 3)*(m + 7) + b*c*(m + 1)*(-2*a*d*(m + 7) + b*c*(m + 5)))*hyper((sympy.S.Half, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(d**2*e*sqrt(c + d*x**4)*(m + 1)*(m + 3)*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_662():
    f = (e*x)**m*(a + b*x**4)/sqrt(c + d*x**4)
    F = b*(e*x)**(m + 1)*sqrt(c + d*x**4)/(d*e*(m + 3)) - (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*(-a*d*(m + 3) + b*c*(m + 1))*hyper((sympy.S.Half, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(d*e*sqrt(c + d*x**4)*(m + 1)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_663():
    f = (e*x)**m/sqrt(c + d*x**4)
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*hyper((sympy.S.Half, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_664():
    f = (e*x)**m/((a + b*x**4)*sqrt(c + d*x**4))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, sympy.S.Half, 1, m/4 + sympy.S(5)/4, -d*x**4/c, -b*x**4/a)/(a*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_665():
    f = (e*x)**m/((a + b*x**4)**2*sqrt(c + d*x**4))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, sympy.S.Half, 2, m/4 + sympy.S(5)/4, -d*x**4/c, -b*x**4/a)/(a**2*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_666():
    f = (e*x)**m/((a + b*x**4)**3*sqrt(c + d*x**4))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, sympy.S.Half, 3, m/4 + sympy.S(5)/4, -d*x**4/c, -b*x**4/a)/(a**3*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_667():
    f = (e*x)**m*(a + b*x**4)**2/(c + d*x**4)**(sympy.S(3)/2)
    F = b**2*(e*x)**(m + 1)*sqrt(c + d*x**4)/(d**2*e*(m + 3)) - (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*(2*b**2*c**2*(m + 1) - (m + 3)*(2*a**2*d**2 - (m + 1)*(-a*d + b*c)**2))*hyper((sympy.S.Half, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(2*c*d**2*e*sqrt(c + d*x**4)*(m + 1)*(m + 3)) + (e*x)**(m + 1)*(-a*d + b*c)**2/(2*c*d**2*e*sqrt(c + d*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_668():
    f = (e*x)**m*(a + b*x**4)/(c + d*x**4)**(sympy.S(3)/2)
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*(a*d*(1 - m) + b*c*(m + 1))*hyper((sympy.S.Half, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(2*c*d*e*sqrt(c + d*x**4)*(m + 1)) - (e*x)**(m + 1)*(-a*d + b*c)/(2*c*d*e*sqrt(c + d*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_669():
    f = (e*x)**m/(c + d*x**4)**(sympy.S(3)/2)
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*hyper((sympy.S(3)/2, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -d*x**4/c)/(c*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_670():
    f = (e*x)**m/((a + b*x**4)*(c + d*x**4)**(sympy.S(3)/2))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, 1, sympy.S(3)/2, m/4 + sympy.S(5)/4, -b*x**4/a, -d*x**4/c)/(a*c*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_671():
    f = (e*x)**m/((a + b*x**4)**2*(c + d*x**4)**(sympy.S(3)/2))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, sympy.S(3)/2, 2, m/4 + sympy.S(5)/4, -d*x**4/c, -b*x**4/a)/(a**2*c*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_672():
    f = (e*x)**m/((a + b*x**4)**3*(c + d*x**4)**(sympy.S(3)/2))
    F = (e*x)**(m + 1)*sqrt(1 + d*x**4/c)*appellf1(m/4 + sympy.S(1)/4, sympy.S(3)/2, 3, m/4 + sympy.S(5)/4, -d*x**4/c, -b*x**4/a)/(a**3*c*e*sqrt(c + d*x**4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_673():
    f = x**17/((a + b*x**6)*sqrt(c + d*x**6))
    F = -a**2*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(3*b**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + (c + d*x**6)**(sympy.S(3)/2)/(9*b*d**2) - sqrt(c + d*x**6)*(a*d + b*c)/(3*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_674():
    f = x**11/((a + b*x**6)*sqrt(c + d*x**6))
    F = a*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(3*b**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + sqrt(c + d*x**6)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_675():
    f = x**5/((a + b*x**6)*sqrt(c + d*x**6))
    F = -atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(3*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_676():
    f = 1/(x*(a + b*x**6)*sqrt(c + d*x**6))
    F = sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(3*a*sqrt(-a*d + b*c)) - atanh(sqrt(c + d*x**6)/sqrt(c))/(3*a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_677():
    f = 1/(x**7*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(c + d*x**6)/(6*a*c*x**6) - b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(3*a**2*sqrt(-a*d + b*c)) + (a*d + 2*b*c)*atanh(sqrt(c + d*x**6)/sqrt(c))/(6*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_678():
    f = x**14/((a + b*x**6)*sqrt(c + d*x**6))
    F = a**(sympy.S(3)/2)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(3*b**2*sqrt(-a*d + b*c)) + x**3*sqrt(c + d*x**6)/(6*b*d) - (2*a*d + b*c)*atanh(sqrt(d)*x**3/sqrt(c + d*x**6))/(6*b**2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_679():
    f = x**8/((a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(a)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(3*b*sqrt(-a*d + b*c)) + atanh(sqrt(d)*x**3/sqrt(c + d*x**6))/(3*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_680():
    f = x**2/((a + b*x**6)*sqrt(c + d*x**6))
    F = atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(3*sqrt(a)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_681():
    f = 1/(x**4*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(c + d*x**6)/(3*a*c*x**3) - b*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(3*a**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_682():
    f = 1/(x**10*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(c + d*x**6)/(9*a*c*x**9) + sqrt(c + d*x**6)*(2*a*d + 3*b*c)/(9*a**2*c**2*x**3) + b**2*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(3*a**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_683():
    f = x**4/((a + b*x**6)*sqrt(c + d*x**6))
    F = x**5*sqrt(1 + d*x**6/c)*appellf1(sympy.S(5)/6, sympy.S.Half, 1, sympy.S(11)/6, -d*x**6/c, -b*x**6/a)/(5*a*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_684():
    f = x**3/((a + b*x**6)*sqrt(c + d*x**6))
    F = x**4*sqrt(1 + d*x**6/c)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -d*x**6/c, -b*x**6/a)/(4*a*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_685():
    f = x/((a + b*x**6)*sqrt(c + d*x**6))
    F = x**2*sqrt(1 + d*x**6/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -d*x**6/c, -b*x**6/a)/(2*a*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_686():
    f = 1/((a + b*x**6)*sqrt(c + d*x**6))
    F = x*sqrt(1 + d*x**6/c)*appellf1(sympy.S(1)/6, sympy.S.Half, 1, sympy.S(7)/6, -d*x**6/c, -b*x**6/a)/(a*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_687():
    f = 1/(x**2*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-1)/6, sympy.S.Half, 1, sympy.S(5)/6, -d*x**6/c, -b*x**6/a)/(a*x*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_688():
    f = 1/(x**3*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -d*x**6/c, -b*x**6/a)/(2*a*x**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_689():
    f = 1/(x**5*(a + b*x**6)*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, -d*x**6/c, -b*x**6/a)/(4*a*x**4*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_690():
    f = x**17/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = -a**2*sqrt(c + d*x**6)/(6*b**2*(a + b*x**6)*(-a*d + b*c)) + a*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(6*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + sqrt(c + d*x**6)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_691():
    f = x**11/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = a*sqrt(c + d*x**6)/(6*b*(a + b*x**6)*(-a*d + b*c)) - (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(6*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_692():
    f = x**5/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(c + d*x**6)/((a + b*x**6)*(-6*a*d + 6*b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(6*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_693():
    f = 1/(x*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = b*sqrt(c + d*x**6)/(6*a*(a + b*x**6)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(6*a**2*(-a*d + b*c)**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**6)/sqrt(c))/(3*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_694():
    f = 1/(x**7*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(c + d*x**6)/(6*a*c*x**6*(a + b*x**6)) - b*sqrt(c + d*x**6)*(-a*d + 2*b*c)/(6*a**2*c*(a + b*x**6)*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**6)/sqrt(-a*d + b*c))/(6*a**3*(-a*d + b*c)**(sympy.S(3)/2)) + (a*d + 4*b*c)*atanh(sqrt(c + d*x**6)/sqrt(c))/(6*a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_695():
    f = x**14/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(a)*(-2*a*d + 3*b*c)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(6*b**2*(-a*d + b*c)**(sympy.S(3)/2)) + a*x**3*sqrt(c + d*x**6)/(6*b*(a + b*x**6)*(-a*d + b*c)) + atanh(sqrt(d)*x**3/sqrt(c + d*x**6))/(3*b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_696():
    f = x**8/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = -x**3*sqrt(c + d*x**6)/((a + b*x**6)*(-6*a*d + 6*b*c)) + c*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(6*sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_697():
    f = x**2/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = b*x**3*sqrt(c + d*x**6)/(6*a*(a + b*x**6)*(-a*d + b*c)) + (-2*a*d + b*c)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(6*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_698():
    f = 1/(x**4*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = b*sqrt(c + d*x**6)/(6*a*x**3*(a + b*x**6)*(-a*d + b*c)) - sqrt(c + d*x**6)*(-2*a*d + 3*b*c)/(6*a**2*c*x**3*(-a*d + b*c)) - b*(-4*a*d + 3*b*c)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(6*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_699():
    f = 1/(x**10*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = b*sqrt(c + d*x**6)/(6*a*x**9*(a + b*x**6)*(-a*d + b*c)) - sqrt(c + d*x**6)*(-2*a*d + 5*b*c)/(18*a**2*c*x**9*(-a*d + b*c)) + sqrt(c + d*x**6)*(-4*a**2*d**2 - 8*a*b*c*d + 15*b**2*c**2)/(18*a**3*c**2*x**3*(-a*d + b*c)) + b**2*(-6*a*d + 5*b*c)*atan(x**3*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**6)))/(6*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_700():
    f = x**4/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = x**5*sqrt(1 + d*x**6/c)*appellf1(sympy.S(5)/6, sympy.S.Half, 2, sympy.S(11)/6, -d*x**6/c, -b*x**6/a)/(5*a**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_701():
    f = x**3/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = x**4*sqrt(1 + d*x**6/c)*appellf1(sympy.S(2)/3, sympy.S.Half, 2, sympy.S(5)/3, -d*x**6/c, -b*x**6/a)/(4*a**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_702():
    f = x/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = x**2*sqrt(1 + d*x**6/c)*appellf1(sympy.S(1)/3, sympy.S.Half, 2, sympy.S(4)/3, -d*x**6/c, -b*x**6/a)/(2*a**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_703():
    f = 1/((a + b*x**6)**2*sqrt(c + d*x**6))
    F = x*sqrt(1 + d*x**6/c)*appellf1(sympy.S(1)/6, sympy.S.Half, 2, sympy.S(7)/6, -d*x**6/c, -b*x**6/a)/(a**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_704():
    f = 1/(x**2*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-1)/6, sympy.S.Half, 2, sympy.S(5)/6, -d*x**6/c, -b*x**6/a)/(a**2*x*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_705():
    f = 1/(x**3*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-1)/3, sympy.S.Half, 2, sympy.S(2)/3, -d*x**6/c, -b*x**6/a)/(2*a**2*x**2*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_706():
    f = 1/(x**5*(a + b*x**6)**2*sqrt(c + d*x**6))
    F = -sqrt(1 + d*x**6/c)*appellf1(sympy.S(-2)/3, sympy.S.Half, 2, sympy.S(1)/3, -d*x**6/c, -b*x**6/a)/(4*a**2*x**4*sqrt(c + d*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_707():
    f = x**23/((a + b*x**8)*sqrt(c + d*x**8))
    F = -a**2*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(4*b**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + (c + d*x**8)**(sympy.S(3)/2)/(12*b*d**2) - sqrt(c + d*x**8)*(a*d + b*c)/(4*b**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_708():
    f = x**15/((a + b*x**8)*sqrt(c + d*x**8))
    F = a*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + sqrt(c + d*x**8)/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_709():
    f = x**7/((a + b*x**8)*sqrt(c + d*x**8))
    F = -atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(4*sqrt(b)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_710():
    f = 1/(x*(a + b*x**8)*sqrt(c + d*x**8))
    F = sqrt(b)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(4*a*sqrt(-a*d + b*c)) - atanh(sqrt(c + d*x**8)/sqrt(c))/(4*a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_711():
    f = 1/(x**9*(a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(c + d*x**8)/(8*a*c*x**8) - b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(4*a**2*sqrt(-a*d + b*c)) + (a*d + 2*b*c)*atanh(sqrt(c + d*x**8)/sqrt(c))/(8*a**2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_712():
    f = x**19/((a + b*x**8)*sqrt(c + d*x**8))
    F = a**(sympy.S(3)/2)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(4*b**2*sqrt(-a*d + b*c)) + x**4*sqrt(c + d*x**8)/(8*b*d) - (2*a*d + b*c)*atanh(sqrt(d)*x**4/sqrt(c + d*x**8))/(8*b**2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_713():
    f = x**11/((a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(a)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(4*b*sqrt(-a*d + b*c)) + atanh(sqrt(d)*x**4/sqrt(c + d*x**8))/(4*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_714():
    f = x**3/((a + b*x**8)*sqrt(c + d*x**8))
    F = atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(4*sqrt(a)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_715():
    f = 1/(x**5*(a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(c + d*x**8)/(4*a*c*x**4) - b*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(4*a**(sympy.S(3)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_716():
    f = 1/(x**13*(a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(c + d*x**8)/(12*a*c*x**12) + sqrt(c + d*x**8)*(2*a*d + 3*b*c)/(12*a**2*c**2*x**4) + b**2*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(4*a**(sympy.S(5)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_717():
    f = x**9/((a + b*x**8)*sqrt(c + d*x**8))
    F = (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_718():
    f = x/((a + b*x**8)*sqrt(c + d*x**8))
    F = (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_719():
    f = 1/(x**7*(a + b*x**8)*sqrt(c + d*x**8))
    F = (Integer(-1) * (sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8))))) * ((Integer(6) * Symbol('a') * Symbol('c') * (x)**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(12) * Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_720():
    f = x**13/((a + b*x**8)*sqrt(c + d*x**8))
    F = (((x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('d')) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4))))))**(Integer(-1))) + ((((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('a') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('a') * (sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_721():
    f = x**5/((a + b*x**8)*sqrt(c + d*x**8))
    F = (sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1))) + (Integer(-1) * (sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_722():
    f = 1/(x**3*(a + b*x**8)*sqrt(c + d*x**8))
    F = (Integer(-1) * (sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8))))) * ((Integer(2) * Symbol('a') * Symbol('c') * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(2) * Symbol('a') * Symbol('c') * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(8) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_723():
    f = x**4/((a + b*x**8)*sqrt(c + d*x**8))
    F = x**5*sqrt(1 + d*x**8/c)*appellf1(sympy.S(5)/8, sympy.S.Half, 1, sympy.S(13)/8, -d*x**8/c, -b*x**8/a)/(5*a*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_724():
    f = x**2/((a + b*x**8)*sqrt(c + d*x**8))
    F = x**3*sqrt(1 + d*x**8/c)*appellf1(sympy.S(3)/8, sympy.S.Half, 1, sympy.S(11)/8, -d*x**8/c, -b*x**8/a)/(3*a*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_725():
    f = 1/((a + b*x**8)*sqrt(c + d*x**8))
    F = x*sqrt(1 + d*x**8/c)*appellf1(sympy.S(1)/8, sympy.S.Half, 1, sympy.S(9)/8, -d*x**8/c, -b*x**8/a)/(a*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_726():
    f = 1/(x**2*(a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(1 + d*x**8/c)*appellf1(sympy.S(-1)/8, sympy.S.Half, 1, sympy.S(7)/8, -d*x**8/c, -b*x**8/a)/(a*x*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_727():
    f = 1/(x**4*(a + b*x**8)*sqrt(c + d*x**8))
    F = -sqrt(1 + d*x**8/c)*appellf1(sympy.S(-3)/8, sympy.S.Half, 1, sympy.S(5)/8, -d*x**8/c, -b*x**8/a)/(3*a*x**3*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_728():
    f = x**23/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = -a**2*sqrt(c + d*x**8)/(8*b**2*(a + b*x**8)*(-a*d + b*c)) + a*(-3*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(8*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + sqrt(c + d*x**8)/(4*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_729():
    f = x**15/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = a*sqrt(c + d*x**8)/(8*b*(a + b*x**8)*(-a*d + b*c)) - (-a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(8*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_730():
    f = x**7/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = -sqrt(c + d*x**8)/((a + b*x**8)*(-8*a*d + 8*b*c)) + d*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(8*sqrt(b)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_731():
    f = 1/(x*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = b*sqrt(c + d*x**8)/(8*a*(a + b*x**8)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(8*a**2*(-a*d + b*c)**(sympy.S(3)/2)) - atanh(sqrt(c + d*x**8)/sqrt(c))/(4*a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_732():
    f = 1/(x**9*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = -sqrt(c + d*x**8)/(8*a*c*x**8*(a + b*x**8)) - b*sqrt(c + d*x**8)*(-a*d + 2*b*c)/(8*a**2*c*(a + b*x**8)*(-a*d + b*c)) - b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*sqrt(c + d*x**8)/sqrt(-a*d + b*c))/(8*a**3*(-a*d + b*c)**(sympy.S(3)/2)) + (a*d + 4*b*c)*atanh(sqrt(c + d*x**8)/sqrt(c))/(8*a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_733():
    f = x**19/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = -sqrt(a)*(-2*a*d + 3*b*c)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(8*b**2*(-a*d + b*c)**(sympy.S(3)/2)) + a*x**4*sqrt(c + d*x**8)/(8*b*(a + b*x**8)*(-a*d + b*c)) + atanh(sqrt(d)*x**4/sqrt(c + d*x**8))/(4*b**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_734():
    f = x**11/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = -x**4*sqrt(c + d*x**8)/((a + b*x**8)*(-8*a*d + 8*b*c)) + c*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(8*sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_735():
    f = x**3/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = b*x**4*sqrt(c + d*x**8)/(8*a*(a + b*x**8)*(-a*d + b*c)) + (-2*a*d + b*c)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(8*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_736():
    f = 1/(x**5*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = b*sqrt(c + d*x**8)/(8*a*x**4*(a + b*x**8)*(-a*d + b*c)) - sqrt(c + d*x**8)*(-2*a*d + 3*b*c)/(8*a**2*c*x**4*(-a*d + b*c)) - b*(-4*a*d + 3*b*c)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(8*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_737():
    f = 1/(x**13*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = b*sqrt(c + d*x**8)/(8*a*x**12*(a + b*x**8)*(-a*d + b*c)) - sqrt(c + d*x**8)*(-2*a*d + 5*b*c)/(24*a**2*c*x**12*(-a*d + b*c)) + sqrt(c + d*x**8)*(-4*a**2*d**2 - 8*a*b*c*d + 15*b**2*c**2)/(24*a**3*c**2*x**4*(-a*d + b*c)) + b**2*(-6*a*d + 5*b*c)*atan(x**4*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**8)))/(8*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_738():
    f = x**9/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_739():
    f = x/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = ((Symbol('b') * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(7) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + sympy.sqrt(Symbol('d'))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_740():
    f = 1/(x**7*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = (Integer(-1) * ((((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(6))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(6)) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(11) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(11) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(48) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + ((Symbol('b') * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * (Symbol('a'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * (Symbol('a'))**(Integer(3)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_741():
    f = x**13/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = ((sympy.sqrt(Symbol('d')) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(6)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1)))) + ((((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('b') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * sympy.sqrt((Integer(-1) * Symbol('a'))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * sympy.sqrt((Integer(-1) * Symbol('a'))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_742():
    f = x**5/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4))))))**(Integer(-1)))) + ((Symbol('b') * (x)**(Integer(6)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * Symbol('a') * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * ((Integer(-1) * Symbol('a')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_743():
    f = 1/(x**3*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = (Integer(-1) * ((((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))) * ((Integer(8) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * (x)**(Integer(8))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * (x)**(Integer(2))) * ((((Integer(-1) * Symbol('a')))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))))) * ((Integer(32) * ((Integer(-1) * Symbol('a')))**((Integer(9) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + ((Symbol('b') * (sympy.sqrt(Symbol('c')) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))) * (((sympy.sqrt(Symbol('c')) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(4)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(64) * ((Integer(-1) * Symbol('a')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(8)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_744():
    f = x**4/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = x**5*sqrt(1 + d*x**8/c)*appellf1(sympy.S(5)/8, sympy.S.Half, 2, sympy.S(13)/8, -d*x**8/c, -b*x**8/a)/(5*a**2*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_745():
    f = x**2/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = x**3*sqrt(1 + d*x**8/c)*appellf1(sympy.S(3)/8, sympy.S.Half, 2, sympy.S(11)/8, -d*x**8/c, -b*x**8/a)/(3*a**2*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_746():
    f = 1/((a + b*x**8)**2*sqrt(c + d*x**8))
    F = x*sqrt(1 + d*x**8/c)*appellf1(sympy.S(1)/8, sympy.S.Half, 2, sympy.S(9)/8, -d*x**8/c, -b*x**8/a)/(a**2*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_747():
    f = 1/(x**2*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = -sqrt(1 + d*x**8/c)*appellf1(sympy.S(-1)/8, sympy.S.Half, 2, sympy.S(7)/8, -d*x**8/c, -b*x**8/a)/(a**2*x*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_748():
    f = 1/(x**4*(a + b*x**8)**2*sqrt(c + d*x**8))
    F = -sqrt(1 + d*x**8/c)*appellf1(sympy.S(-3)/8, sympy.S.Half, 2, sympy.S(5)/8, -d*x**8/c, -b*x**8/a)/(3*a**2*x**3*sqrt(c + d*x**8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_749():
    f = x**5*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**6*(c + d/x**2)**(sympy.S(3)/2)/(6*c) + x**4*sqrt(c + d/x**2)*(-a*d + 2*b*c)/(8*c) + d*x**2*sqrt(c + d/x**2)*(-a*d + 2*b*c)/(16*c**2) - d**2*(-a*d + 2*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(16*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_750():
    f = x**3*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**4*(c + d/x**2)**(sympy.S(3)/2)/(4*c) + x**2*sqrt(c + d/x**2)*(-a*d + 4*b*c)/(8*c) + d*(-a*d + 4*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(8*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_751():
    f = x*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**2*(c + d/x**2)**(sympy.S(3)/2)/(2*c) - sqrt(c + d/x**2)*(a*d + 2*b*c)/(2*c) + (a*d + 2*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_752():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x
    F = a*sqrt(c)*atanh(sqrt(c + d/x**2)/sqrt(c)) - a*sqrt(c + d/x**2) - b*(c + d/x**2)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_753():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**3
    F = -b*(c + d/x**2)**(sympy.S(5)/2)/(5*d**2) + (c + d/x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_754():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**5
    F = -b*(c + d/x**2)**(sympy.S(7)/2)/(7*d**3) - c*(c + d/x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**3) + (c + d/x**2)**(sympy.S(5)/2)*(-a*d + 2*b*c)/(5*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_755():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**7
    F = -b*(c + d/x**2)**(sympy.S(9)/2)/(9*d**4) + c**2*(c + d/x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**4) - c*(c + d/x**2)**(sympy.S(5)/2)*(-2*a*d + 3*b*c)/(5*d**4) + (c + d/x**2)**(sympy.S(7)/2)*(-a*d + 3*b*c)/(7*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_756():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**9
    F = -b*(c + d/x**2)**(sympy.S(11)/2)/(11*d**5) - c**3*(c + d/x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(3*d**5) + c**2*(c + d/x**2)**(sympy.S(5)/2)*(-3*a*d + 4*b*c)/(5*d**5) - 3*c*(c + d/x**2)**(sympy.S(7)/2)*(-a*d + 2*b*c)/(7*d**5) + (c + d/x**2)**(sympy.S(9)/2)*(-a*d + 4*b*c)/(9*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_757():
    f = x**10*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**11*(c + d/x**2)**(sympy.S(3)/2)/(11*c) + x**9*(c + d/x**2)**(sympy.S(3)/2)*(-8*a*d + 11*b*c)/(99*c**2) - 2*d*x**7*(c + d/x**2)**(sympy.S(3)/2)*(-8*a*d + 11*b*c)/(231*c**3) + 8*d**2*x**5*(c + d/x**2)**(sympy.S(3)/2)*(-8*a*d + 11*b*c)/(1155*c**4) - 16*d**3*x**3*(c + d/x**2)**(sympy.S(3)/2)*(-8*a*d + 11*b*c)/(3465*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_758():
    f = x**8*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**9*(c + d/x**2)**(sympy.S(3)/2)/(9*c) + x**7*(c + d/x**2)**(sympy.S(3)/2)*(-2*a*d + 3*b*c)/(21*c**2) - 4*d*x**5*(c + d/x**2)**(sympy.S(3)/2)*(-2*a*d + 3*b*c)/(105*c**3) + 8*d**2*x**3*(c + d/x**2)**(sympy.S(3)/2)*(-2*a*d + 3*b*c)/(315*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_759():
    f = x**6*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**7*(c + d/x**2)**(sympy.S(3)/2)/(7*c) + x**5*(c + d/x**2)**(sympy.S(3)/2)*(-4*a*d + 7*b*c)/(35*c**2) - 2*d*x**3*(c + d/x**2)**(sympy.S(3)/2)*(-4*a*d + 7*b*c)/(105*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_760():
    f = x**4*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**5*(c + d/x**2)**(sympy.S(3)/2)/(5*c) + x**3*(c + d/x**2)**(sympy.S(3)/2)*(-2*a*d + 5*b*c)/(15*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_761():
    f = x**2*(a + b/x**2)*sqrt(c + d/x**2)
    F = a*x**3*(c + d/x**2)**(sympy.S(3)/2)/(3*c) - b*sqrt(d)*atanh(sqrt(d)/(x*sqrt(c + d/x**2))) + b*x*sqrt(c + d/x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_762():
    f = (a + b/x**2)*sqrt(c + d/x**2)
    F = a*x*(c + d/x**2)**(sympy.S(3)/2)/c - (2*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(2*sqrt(d)) - sqrt(c + d/x**2)*(2*a*d + b*c)/(2*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_763():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**2
    F = -b*(c + d/x**2)**(sympy.S(3)/2)/(4*d*x) + c*(-4*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(8*d**(sympy.S(3)/2)) + sqrt(c + d/x**2)*(-4*a*d + b*c)/(8*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_764():
    f = (a + b/x**2)*sqrt(c + d/x**2)/x**4
    F = -b*(c + d/x**2)**(sympy.S(3)/2)/(6*d*x**3) - c**2*(-2*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(16*d**(sympy.S(5)/2)) + c*sqrt(c + d/x**2)*(-2*a*d + b*c)/(16*d**2*x) + sqrt(c + d/x**2)*(-2*a*d + b*c)/(8*d*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_765():
    f = x**5*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**6*(c + d/x**2)**(sympy.S(5)/2)/(6*c) + d*x**2*sqrt(c + d/x**2)*(-a*d + 6*b*c)/(16*c) + x**4*(c + d/x**2)**(sympy.S(3)/2)*(-a*d + 6*b*c)/(24*c) + d**2*(-a*d + 6*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(16*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_766():
    f = x**3*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**4*(c + d/x**2)**(sympy.S(5)/2)/(4*c) - 3*d*sqrt(c + d/x**2)*(a*d + 4*b*c)/(8*c) + x**2*(c + d/x**2)**(sympy.S(3)/2)*(a*d + 4*b*c)/(8*c) + 3*d*(a*d + 4*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(8*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_767():
    f = x*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**2*(c + d/x**2)**(sympy.S(5)/2)/(2*c) + sqrt(c)*(3*a*d + 2*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/2 + sqrt(c + d/x**2)*(-3*a*d/2 - b*c) - (c + d/x**2)**(sympy.S(3)/2)*(3*a*d + 2*b*c)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_768():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x
    F = a*c**(sympy.S(3)/2)*atanh(sqrt(c + d/x**2)/sqrt(c)) - a*c*sqrt(c + d/x**2) - a*(c + d/x**2)**(sympy.S(3)/2)/3 - b*(c + d/x**2)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_769():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**3
    F = -b*(c + d/x**2)**(sympy.S(7)/2)/(7*d**2) + (c + d/x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_770():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**5
    F = -b*(c + d/x**2)**(sympy.S(9)/2)/(9*d**3) - c*(c + d/x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**3) + (c + d/x**2)**(sympy.S(7)/2)*(-a*d + 2*b*c)/(7*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_771():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**7
    F = -b*(c + d/x**2)**(sympy.S(11)/2)/(11*d**4) + c**2*(c + d/x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**4) - c*(c + d/x**2)**(sympy.S(7)/2)*(-2*a*d + 3*b*c)/(7*d**4) + (c + d/x**2)**(sympy.S(9)/2)*(-a*d + 3*b*c)/(9*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_772():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**9
    F = -b*(c + d/x**2)**(sympy.S(13)/2)/(13*d**5) - c**3*(c + d/x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(5*d**5) + c**2*(c + d/x**2)**(sympy.S(7)/2)*(-3*a*d + 4*b*c)/(7*d**5) - c*(c + d/x**2)**(sympy.S(9)/2)*(-a*d + 2*b*c)/(3*d**5) + (c + d/x**2)**(sympy.S(11)/2)*(-a*d + 4*b*c)/(11*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_773():
    f = x**12*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**13*(c + d/x**2)**(sympy.S(5)/2)/(13*c) + x**11*(c + d/x**2)**(sympy.S(5)/2)*(-8*a*d + 13*b*c)/(143*c**2) - 2*d*x**9*(c + d/x**2)**(sympy.S(5)/2)*(-8*a*d + 13*b*c)/(429*c**3) + 8*d**2*x**7*(c + d/x**2)**(sympy.S(5)/2)*(-8*a*d + 13*b*c)/(3003*c**4) - 16*d**3*x**5*(c + d/x**2)**(sympy.S(5)/2)*(-8*a*d + 13*b*c)/(15015*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_774():
    f = x**10*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**11*(c + d/x**2)**(sympy.S(5)/2)/(11*c) + x**9*(c + d/x**2)**(sympy.S(5)/2)*(-6*a*d + 11*b*c)/(99*c**2) - 4*d*x**7*(c + d/x**2)**(sympy.S(5)/2)*(-6*a*d + 11*b*c)/(693*c**3) + 8*d**2*x**5*(c + d/x**2)**(sympy.S(5)/2)*(-6*a*d + 11*b*c)/(3465*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_775():
    f = x**8*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**9*(c + d/x**2)**(sympy.S(5)/2)/(9*c) + x**7*(c + d/x**2)**(sympy.S(5)/2)*(-4*a*d + 9*b*c)/(63*c**2) - 2*d*x**5*(c + d/x**2)**(sympy.S(5)/2)*(-4*a*d + 9*b*c)/(315*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_776():
    f = x**6*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**7*(c + d/x**2)**(sympy.S(5)/2)/(7*c) + x**5*(c + d/x**2)**(sympy.S(5)/2)*(-2*a*d + 7*b*c)/(35*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_777():
    f = x**4*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**5*(c + d/x**2)**(sympy.S(5)/2)/(5*c) - b*d**(sympy.S(3)/2)*atanh(sqrt(d)/(x*sqrt(c + d/x**2))) + b*d*x*sqrt(c + d/x**2) + b*x**3*(c + d/x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_778():
    f = x**2*(a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**3*(c + d/x**2)**(sympy.S(5)/2)/(3*c) - sqrt(d)*(2*a*d + 3*b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/2 - d*sqrt(c + d/x**2)*(2*a*d + 3*b*c)/(2*c*x) + x*(c + d/x**2)**(sympy.S(3)/2)*(2*a*d + 3*b*c)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_779():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)
    F = a*x*(c + d/x**2)**(sympy.S(5)/2)/c - 3*c*(4*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(8*sqrt(d)) - sqrt(c + d/x**2)*(12*a*d + 3*b*c)/(8*x) - (c + d/x**2)**(sympy.S(3)/2)*(4*a*d + b*c)/(4*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_780():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**2
    F = -b*(c + d/x**2)**(sympy.S(5)/2)/(6*d*x) + c**2*(-6*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(16*d**(sympy.S(3)/2)) + c*sqrt(c + d/x**2)*(-6*a*d + b*c)/(16*d*x) + (c + d/x**2)**(sympy.S(3)/2)*(-6*a*d + b*c)/(24*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_781():
    f = (a + b/x**2)*(c + d/x**2)**(sympy.S(3)/2)/x**4
    F = -b*(c + d/x**2)**(sympy.S(5)/2)/(8*d*x**3) - c**3*(-8*a*d + 3*b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(128*d**(sympy.S(5)/2)) + c**2*sqrt(c + d/x**2)*(-8*a*d + 3*b*c)/(128*d**2*x) + c*sqrt(c + d/x**2)*(-8*a*d + 3*b*c)/(64*d*x**3) + (c + d/x**2)**(sympy.S(3)/2)*(-8*a*d + 3*b*c)/(48*d*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_782():
    f = x**3*(a + b/x**2)/sqrt(c + d/x**2)
    F = a*x**4*sqrt(c + d/x**2)/(4*c) + x**2*sqrt(c + d/x**2)*(-3*a*d + 4*b*c)/(8*c**2) - d*(-3*a*d + 4*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_783():
    f = x*(a + b/x**2)/sqrt(c + d/x**2)
    F = a*x**2*sqrt(c + d/x**2)/(2*c) + (-a*d + 2*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_784():
    f = (a + b/x**2)/(x*sqrt(c + d/x**2))
    F = a*atanh(sqrt(c + d/x**2)/sqrt(c))/sqrt(c) - b*sqrt(c + d/x**2)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_785():
    f = (a + b/x**2)/(x**3*sqrt(c + d/x**2))
    F = -b*(c + d/x**2)**(sympy.S(3)/2)/(3*d**2) + sqrt(c + d/x**2)*(-a*d + b*c)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_786():
    f = (a + b/x**2)/(x**5*sqrt(c + d/x**2))
    F = -b*(c + d/x**2)**(sympy.S(5)/2)/(5*d**3) - c*sqrt(c + d/x**2)*(-a*d + b*c)/d**3 + (c + d/x**2)**(sympy.S(3)/2)*(-a*d + 2*b*c)/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_787():
    f = (a + b/x**2)/(x**7*sqrt(c + d/x**2))
    F = -b*(c + d/x**2)**(sympy.S(7)/2)/(7*d**4) + c**2*sqrt(c + d/x**2)*(-a*d + b*c)/d**4 - c*(c + d/x**2)**(sympy.S(3)/2)*(-2*a*d + 3*b*c)/(3*d**4) + (c + d/x**2)**(sympy.S(5)/2)*(-a*d + 3*b*c)/(5*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_788():
    f = x**4*(a + b/x**2)/sqrt(c + d/x**2)
    F = a*x**5*sqrt(c + d/x**2)/(5*c) + x**3*sqrt(c + d/x**2)*(-4*a*d + 5*b*c)/(15*c**2) - 2*d*x*sqrt(c + d/x**2)*(-4*a*d + 5*b*c)/(15*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_789():
    f = x**2*(a + b/x**2)/sqrt(c + d/x**2)
    F = a*x**3*sqrt(c + d/x**2)/(3*c) + x*sqrt(c + d/x**2)*(-2*a*d + 3*b*c)/(3*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_790():
    f = (a + b/x**2)/sqrt(c + d/x**2)
    F = a*x*sqrt(c + d/x**2)/c - b*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_791():
    f = (a + b/x**2)/(x**2*sqrt(c + d/x**2))
    F = -b*sqrt(c + d/x**2)/(2*d*x) + (-2*a*d + b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_792():
    f = (a + b/x**2)/(x**4*sqrt(c + d/x**2))
    F = -b*sqrt(c + d/x**2)/(4*d*x**3) - c*(-4*a*d + 3*b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(8*d**(sympy.S(5)/2)) + sqrt(c + d/x**2)*(-4*a*d + 3*b*c)/(8*d**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_793():
    f = x**3*(a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**4/(4*c*sqrt(c + d/x**2)) + x**2*(-5*a*d + 4*b*c)/(8*c**2*sqrt(c + d/x**2)) + 3*d*(-5*a*d + 4*b*c)/(8*c**3*sqrt(c + d/x**2)) - 3*d*(-5*a*d + 4*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(8*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_794():
    f = x*(a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**2/(2*c*sqrt(c + d/x**2)) - (-3*a*d + 2*b*c)/(2*c**2*sqrt(c + d/x**2)) + (-3*a*d + 2*b*c)*atanh(sqrt(c + d/x**2)/sqrt(c))/(2*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_795():
    f = (a + b/x**2)/(x*(c + d/x**2)**(sympy.S(3)/2))
    F = a*atanh(sqrt(c + d/x**2)/sqrt(c))/c**(sympy.S(3)/2) + (-a*d + b*c)/(c*d*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_796():
    f = (a + b/x**2)/(x**3*(c + d/x**2)**(sympy.S(3)/2))
    F = -b*sqrt(c + d/x**2)/d**2 - (-a*d + b*c)/(d**2*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_797():
    f = (a + b/x**2)/(x**5*(c + d/x**2)**(sympy.S(3)/2))
    F = -b*(c + d/x**2)**(sympy.S(3)/2)/(3*d**3) + c*(-a*d + b*c)/(d**3*sqrt(c + d/x**2)) + sqrt(c + d/x**2)*(-a*d + 2*b*c)/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_798():
    f = (a + b/x**2)/(x**7*(c + d/x**2)**(sympy.S(3)/2))
    F = -b*(c + d/x**2)**(sympy.S(5)/2)/(5*d**4) - c**2*(-a*d + b*c)/(d**4*sqrt(c + d/x**2)) - c*sqrt(c + d/x**2)*(-2*a*d + 3*b*c)/d**4 + (c + d/x**2)**(sympy.S(3)/2)*(-a*d + 3*b*c)/(3*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_799():
    f = (a + b/x**2)/(x**9*(c + d/x**2)**(sympy.S(3)/2))
    F = -b*(c + d/x**2)**(sympy.S(7)/2)/(7*d**5) + c**3*(-a*d + b*c)/(d**5*sqrt(c + d/x**2)) + c**2*sqrt(c + d/x**2)*(-3*a*d + 4*b*c)/d**5 - c*(c + d/x**2)**(sympy.S(3)/2)*(-a*d + 2*b*c)/d**5 + (c + d/x**2)**(sympy.S(5)/2)*(-a*d + 4*b*c)/(5*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_800():
    f = x**4*(a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**5/(5*c*sqrt(c + d/x**2)) + x**3*(-6*a*d + 5*b*c)/(15*c**2*sqrt(c + d/x**2)) + 4*d*x*(-6*a*d + 5*b*c)/(15*c**3*sqrt(c + d/x**2)) - 8*d*x*sqrt(c + d/x**2)*(-6*a*d + 5*b*c)/(15*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_801():
    f = x**2*(a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = a*x**3/(3*c*sqrt(c + d/x**2)) - x*(-4*a*d + 3*b*c)/(3*c**2*sqrt(c + d/x**2)) + x*sqrt(c + d/x**2)*(-8*a*d + 6*b*c)/(3*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_802():
    f = (a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = a*x/(c*sqrt(c + d/x**2)) - (-2*a*d + b*c)/(c**2*x*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_803():
    f = (a + b/x**2)/(x**2*(c + d/x**2)**(sympy.S(3)/2))
    F = -b*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/d**(sympy.S(3)/2) + (-a*d + b*c)/(c*d*x*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_804():
    f = (a + b/x**2)/(x**4*(c + d/x**2)**(sympy.S(3)/2))
    F = -b/(2*d*x**3*sqrt(c + d/x**2)) - (-2*a*d + 3*b*c)/(2*d**2*x*sqrt(c + d/x**2)) + (-2*a*d + 3*b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(2*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_805():
    f = (a + b/x**2)/(x**6*(c + d/x**2)**(sympy.S(3)/2))
    F = -b/(4*d*x**5*sqrt(c + d/x**2)) - 3*c*(-4*a*d + 5*b*c)*atanh(sqrt(d)/(x*sqrt(c + d/x**2)))/(8*d**(sympy.S(7)/2)) - (-4*a*d + 5*b*c)/(4*d**2*x**3*sqrt(c + d/x**2)) + sqrt(c + d/x**2)*(-12*a*d + 15*b*c)/(8*d**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_806():
    f = x**4*(a + b/x**2)**p*(c + d/x**2)**q
    F = x**5*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-5)/2, -p, -q, sympy.S(-3)/2, -b/(a*x**2), -d/(c*x**2))/(5*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_807():
    f = x**3*(a + b/x**2)**p*(c + d/x**2)**q
    F = b**2*(a + b/x**2)**(p + 1)*(c + d/x**2)**q*appellf1(p + 1, 3, -q, p + 2, (a + b/x**2)/a, -d*(a + b/x**2)/(-a*d + b*c))/(2*a**3*(b*(c + d/x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_808():
    f = x**2*(a + b/x**2)**p*(c + d/x**2)**q
    F = x**3*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-3)/2, -p, -q, sympy.S(-1)/2, -b/(a*x**2), -d/(c*x**2))/(3*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_809():
    f = x*(a + b/x**2)**p*(c + d/x**2)**q
    F = -b*(a + b/x**2)**(p + 1)*(c + d/x**2)**q*appellf1(p + 1, 2, -q, p + 2, (a + b/x**2)/a, -d*(a + b/x**2)/(-a*d + b*c))/(2*a**2*(b*(c + d/x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_810():
    f = (a + b/x**2)**p*(c + d/x**2)**q
    F = x*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-1)/2, -p, -q, sympy.S.Half, -b/(a*x**2), -d/(c*x**2))/((1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_811():
    f = (a + b/x**2)**p*(c + d/x**2)**q/x
    F = (a + b/x**2)**(p + 1)*(c + d/x**2)**q*appellf1(p + 1, 1, -q, p + 2, (a + b/x**2)/a, -d*(a + b/x**2)/(-a*d + b*c))/(2*a*(b*(c + d/x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_812():
    f = (a + b/x**2)**p*(c + d/x**2)**q/x**2
    F = -(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -b/(a*x**2), -d/(c*x**2))/(x*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_813():
    f = (a + b/x**2)**p*(c + d/x**2)**q/x**3
    F = -(a + b/x**2)**(p + 1)*(c + d/x**2)**q*hyper((-q, p + 1), (p + 2,), -d*(a + b/x**2)/(-a*d + b*c))/(2*b*(b*(c + d/x**2)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_814():
    f = (a + b/x**2)**p*(c + d/x**2)**q/x**4
    F = -(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(3)/2, -p, -q, sympy.S(5)/2, -b/(a*x**2), -d/(c*x**2))/(3*x**3*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_815():
    f = (e*x)**(sympy.S(5)/2)*(a + b/x**2)**p*(c + d/x**2)**q
    F = 2*(e*x)**(sympy.S(7)/2)*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-7)/4, -p, -q, sympy.S(-3)/4, -b/(a*x**2), -d/(c*x**2))/(7*e*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_816():
    f = (e*x)**(sympy.S(3)/2)*(a + b/x**2)**p*(c + d/x**2)**q
    F = 2*(e*x)**(sympy.S(5)/2)*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-5)/4, -p, -q, sympy.S(-1)/4, -b/(a*x**2), -d/(c*x**2))/(5*e*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_817():
    f = sqrt(e*x)*(a + b/x**2)**p*(c + d/x**2)**q
    F = 2*(e*x)**(sympy.S(3)/2)*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-3)/4, -p, -q, sympy.S(1)/4, -b/(a*x**2), -d/(c*x**2))/(3*e*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_818():
    f = (a + b/x**2)**p*(c + d/x**2)**q/sqrt(e*x)
    F = 2*sqrt(e*x)*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-1)/4, -p, -q, sympy.S(3)/4, -b/(a*x**2), -d/(c*x**2))/(e*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_819():
    f = (a + b/x**2)**p*(c + d/x**2)**q/(e*x)**(sympy.S(3)/2)
    F = -2*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(1)/4, -p, -q, sympy.S(5)/4, -b/(a*x**2), -d/(c*x**2))/(e*sqrt(e*x)*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_820():
    f = (a + b/x**2)**p*(c + d/x**2)**q/(e*x)**(sympy.S(5)/2)
    F = -2*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(3)/4, -p, -q, sympy.S(7)/4, -b/(a*x**2), -d/(c*x**2))/(3*e*(e*x)**(sympy.S(3)/2)*(1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_821():
    f = x**(sympy.S(5)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)
    F = x**(sympy.S(7)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/4 - x**(sympy.S(5)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/24 - 5*x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/96 - 5*sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/64 - 5*acosh(sqrt(x))/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_822():
    f = x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)
    F = x**(sympy.S(5)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/3 - x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/12 - sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/8 - acosh(sqrt(x))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_823():
    f = sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)
    F = x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/2 - sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/4 - acosh(sqrt(x))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_824():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/sqrt(x)
    F = sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1) - acosh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_825():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/x**(sympy.S(3)/2)
    F = -2*sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1) + 2*acosh(sqrt(x)) + 2*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_826():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/x**(sympy.S(5)/2)
    F = 2*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_827():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/x**(sympy.S(7)/2)
    F = 4*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(15*x**(sympy.S(3)/2)) + 2*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_828():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/x**(sympy.S(9)/2)
    F = 16*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(105*x**(sympy.S(3)/2)) + 8*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(35*x**(sympy.S(5)/2)) + 2*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(7*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_829():
    f = sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/x**(sympy.S(11)/2)
    F = 32*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(315*x**(sympy.S(3)/2)) + 16*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(105*x**(sympy.S(5)/2)) + 4*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(21*x**(sympy.S(7)/2)) + 2*(sqrt(x) - 1)**(sympy.S(3)/2)*(sqrt(x) + 1)**(sympy.S(3)/2)/(9*x**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_830():
    f = x**(sympy.S(5)/2)/(sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = x**(sympy.S(5)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/3 + 5*x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/12 + 5*sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/8 + 5*acosh(sqrt(x))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_831():
    f = x**(sympy.S(3)/2)/(sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/2 + 3*sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/4 + 3*acosh(sqrt(x))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_832():
    f = sqrt(x)/(sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1) + acosh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_833():
    f = 1/(sqrt(x)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = 2*acosh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_834():
    f = 1/(x**(sympy.S(3)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = 2*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_835():
    f = 1/(x**(sympy.S(5)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = 4*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/(3*sqrt(x)) + 2*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_836():
    f = 1/(x**(sympy.S(7)/2)*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1))
    F = 16*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/(15*sqrt(x)) + 8*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/(15*x**(sympy.S(3)/2)) + 2*sqrt(sqrt(x) - 1)*sqrt(sqrt(x) + 1)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_837():
    f = x**2*(-a + b*x**n)**p*(a + b*x**n)**p
    F = x**3*(-a + b*x**n)**p*(a + b*x**n)**p*hyper((-p, 3/(2*n)), (1 + 3/(2*n),), b**2*x**(2*n)/a**2)/(3*(1 - b**2*x**(2*n)/a**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_838():
    f = x*(-a + b*x**n)**p*(a + b*x**n)**p
    F = x**2*(-a + b*x**n)**p*(a + b*x**n)**p*hyper((1/n, -p), (1 + 1/n,), b**2*x**(2*n)/a**2)/(2*(1 - b**2*x**(2*n)/a**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_839():
    f = (-a + b*x**n)**p*(a + b*x**n)**p
    F = x*(-a + b*x**n)**p*(a + b*x**n)**p*hyper((-p, 1/(2*n)), (1 + 1/(2*n),), b**2*x**(2*n)/a**2)/(1 - b**2*x**(2*n)/a**2)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_840():
    f = (-a + b*x**n)**p*(a + b*x**n)**p/x
    F = -(-a + b*x**n)**p*(a + b*x**n)**p*(a**2 - b**2*x**(2*n))*hyper((1, p + 1), (p + 2,), 1 - b**2*x**(2*n)/a**2)/(2*a**2*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_841():
    f = (-a + b*x**n)**p*(a + b*x**n)**p/x**2
    F = -(-a + b*x**n)**p*(a + b*x**n)**p*hyper((-p, -1/(2*n)), (1 - 1/(2*n),), b**2*x**(2*n)/a**2)/(x*(1 - b**2*x**(2*n)/a**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_842():
    f = (x**6 + 1)/(x*(1 - x**6))
    F = log(x) - log(1 - x**6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_843():
    f = (e*x)**m*(a + b*x**n)**p*(a*(m + 1) + b*x**n*(m + n*p + n + 1))
    F = (e*x)**(m + 1)*(a + b*x**n)**(p + 1)/e
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_844():
    f = (e*x)**m/((a + b*x**n)*(c + d*x**n))
    F = -d*(e*x)**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*e*(m + 1)*(-a*d + b*c)) + b*(e*x)**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*e*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_845():
    f = x**2/((a + b*x**n)*(c + d*x**n))
    F = -d*x**3*hyper((1, 3/n), ((n + 3)/n,), -d*x**n/c)/(3*c*(-a*d + b*c)) + b*x**3*hyper((1, 3/n), ((n + 3)/n,), -b*x**n/a)/(3*a*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_846():
    f = x/((a + b*x**n)*(c + d*x**n))
    F = -d*x**2*hyper((1, 2/n), ((n + 2)/n,), -d*x**n/c)/(2*c*(-a*d + b*c)) + b*x**2*hyper((1, 2/n), ((n + 2)/n,), -b*x**n/a)/(2*a*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_847():
    f = 1/((a + b*x**n)*(c + d*x**n))
    F = -d*x*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*(-a*d + b*c)) + b*x*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_848():
    f = 1/(x*(a + b*x**n)*(c + d*x**n))
    F = d*log(c + d*x**n)/(c*n*(-a*d + b*c)) - b*log(a + b*x**n)/(a*n*(-a*d + b*c)) + log(x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_849():
    f = 1/(x**2*(a + b*x**n)*(c + d*x**n))
    F = d*hyper((1, -1/n), (-(1 - n)/n,), -d*x**n/c)/(c*x*(-a*d + b*c)) - b*hyper((1, -1/n), (-(1 - n)/n,), -b*x**n/a)/(a*x*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_850():
    f = 1/(x**3*(a + b*x**n)*(c + d*x**n))
    F = d*hyper((1, -2/n), (-(2 - n)/n,), -d*x**n/c)/(2*c*x**2*(-a*d + b*c)) - b*hyper((1, -2/n), (-(2 - n)/n,), -b*x**n/a)/(2*a*x**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_851():
    f = (e*x)**m/((a + b*x**n)**2*(c + d*x**n))
    F = d**2*(e*x)**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*e*(m + 1)*(-a*d + b*c)**2) + b*(e*x)**(m + 1)/(a*e*n*(a + b*x**n)*(-a*d + b*c)) + b*(e*x)**(m + 1)*(a*d*(m - 2*n + 1) - b*c*(m - n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*e*n*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_852():
    f = x**2/((a + b*x**n)**2*(c + d*x**n))
    F = d**2*x**3*hyper((1, 3/n), ((n + 3)/n,), -d*x**n/c)/(3*c*(-a*d + b*c)**2) + b*x**3/(a*n*(a + b*x**n)*(-a*d + b*c)) + b*x**3*(a*d*(3 - 2*n) - b*c*(3 - n))*hyper((1, 3/n), ((n + 3)/n,), -b*x**n/a)/(3*a**2*n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_853():
    f = x/((a + b*x**n)**2*(c + d*x**n))
    F = d**2*x**2*hyper((1, 2/n), ((n + 2)/n,), -d*x**n/c)/(2*c*(-a*d + b*c)**2) + b*x**2/(a*n*(a + b*x**n)*(-a*d + b*c)) + b*x**2*(2*a*d*(1 - n) - b*c*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -b*x**n/a)/(2*a**2*n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_854():
    f = 1/((a + b*x**n)**2*(c + d*x**n))
    F = d**2*x*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*(-a*d + b*c)**2) + b*x/(a*n*(a + b*x**n)*(-a*d + b*c)) + b*x*(a*d*(1 - 2*n) - b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_855():
    f = 1/(x*(a + b*x**n)**2*(c + d*x**n))
    F = -d**2*log(c + d*x**n)/(c*n*(-a*d + b*c)**2) + b/(a*n*(a + b*x**n)*(-a*d + b*c)) - b*(-2*a*d + b*c)*log(a + b*x**n)/(a**2*n*(-a*d + b*c)**2) + log(x)/(a**2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_856():
    f = 1/(x**2*(a + b*x**n)**2*(c + d*x**n))
    F = -d**2*hyper((1, -1/n), (-(1 - n)/n,), -d*x**n/c)/(c*x*(-a*d + b*c)**2) + b/(a*n*x*(a + b*x**n)*(-a*d + b*c)) - b*(-a*d*(2*n + 1) + b*c*(n + 1))*hyper((1, -1/n), (-(1 - n)/n,), -b*x**n/a)/(a**2*n*x*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_857():
    f = 1/(x**3*(a + b*x**n)**2*(c + d*x**n))
    F = -d**2*hyper((1, -2/n), (-(2 - n)/n,), -d*x**n/c)/(2*c*x**2*(-a*d + b*c)**2) + b/(a*n*x**2*(a + b*x**n)*(-a*d + b*c)) + b*(2*a*d*(n + 1) - b*c*(n + 2))*hyper((1, -2/n), (-(2 - n)/n,), -b*x**n/a)/(2*a**2*n*x**2*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_858():
    f = x**(2*n - 1)*(a + b*x**n)**3/(c + d*x**n)
    F = b**3*x**(4*n)/(4*d*n) - b**2*x**(3*n)*(-3*a*d + b*c)/(3*d**2*n) + b*x**(2*n)*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/(2*d**3*n) + c*(-a*d + b*c)**3*log(c + d*x**n)/(d**5*n) - x**n*(-a*d + b*c)**3/(d**4*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_859():
    f = x**(2*n - 1)*(a + b*x**n)**2/(c + d*x**n)
    F = b**2*x**(3*n)/(3*d*n) - b*x**(2*n)*(-2*a*d + b*c)/(2*d**2*n) - c*(-a*d + b*c)**2*log(c + d*x**n)/(d**4*n) + x**n*(-a*d + b*c)**2/(d**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_860():
    f = x**(2*n - 1)*(a + b*x**n)/(c + d*x**n)
    F = b*x**(2*n)/(2*d*n) + c*(-a*d + b*c)*log(c + d*x**n)/(d**3*n) - x**n*(-a*d + b*c)/(d**2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_861():
    f = x**(2*n - 1)/((a + b*x**n)*(c + d*x**n))
    F = -a*log(a + b*x**n)/(b*n*(-a*d + b*c)) + c*log(c + d*x**n)/(d*n*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_862():
    f = x**(2*n - 1)/((a + b*x**n)**2*(c + d*x**n))
    F = a/(b*n*(a + b*x**n)*(-a*d + b*c)) + c*log(a + b*x**n)/(n*(-a*d + b*c)**2) - c*log(c + d*x**n)/(n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_863():
    f = x**(2*n - 1)/((a + b*x**n)**3*(c + d*x**n))
    F = a/(2*b*n*(a + b*x**n)**2*(-a*d + b*c)) - c*d*log(a + b*x**n)/(n*(-a*d + b*c)**3) + c*d*log(c + d*x**n)/(n*(-a*d + b*c)**3) - c/(n*(a + b*x**n)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_864():
    f = x**(3*n - 1)*(a + b*x**n)**3/(c + d*x**n)
    F = b**3*x**(5*n)/(5*d*n) - b**2*x**(4*n)*(-3*a*d + b*c)/(4*d**2*n) + b*x**(3*n)*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/(3*d**3*n) - c**2*(-a*d + b*c)**3*log(c + d*x**n)/(d**6*n) + c*x**n*(-a*d + b*c)**3/(d**5*n) - x**(2*n)*(-a*d + b*c)**3/(2*d**4*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_865():
    f = x**(3*n - 1)*(a + b*x**n)**2/(c + d*x**n)
    F = b**2*x**(4*n)/(4*d*n) - b*x**(3*n)*(-2*a*d + b*c)/(3*d**2*n) + c**2*(-a*d + b*c)**2*log(c + d*x**n)/(d**5*n) - c*x**n*(-a*d + b*c)**2/(d**4*n) + x**(2*n)*(-a*d + b*c)**2/(2*d**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_866():
    f = x**(3*n - 1)*(a + b*x**n)/(c + d*x**n)
    F = b*x**(3*n)/(3*d*n) - c**2*(-a*d + b*c)*log(c + d*x**n)/(d**4*n) + c*x**n*(-a*d + b*c)/(d**3*n) - x**(2*n)*(-a*d + b*c)/(2*d**2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_867():
    f = x**(3*n - 1)/((a + b*x**n)*(c + d*x**n))
    F = a**2*log(a + b*x**n)/(b**2*n*(-a*d + b*c)) - c**2*log(c + d*x**n)/(d**2*n*(-a*d + b*c)) + x**n/(b*d*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_868():
    f = x**(3*n - 1)/((a + b*x**n)**2*(c + d*x**n))
    F = -a**2/(b**2*n*(a + b*x**n)*(-a*d + b*c)) - a*(-a*d + 2*b*c)*log(a + b*x**n)/(b**2*n*(-a*d + b*c)**2) + c**2*log(c + d*x**n)/(d*n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_869():
    f = x**(3*n - 1)/((a + b*x**n)**3*(c + d*x**n))
    F = -a**2/(2*b**2*n*(a + b*x**n)**2*(-a*d + b*c)) + a*(-a*d + 2*b*c)/(b**2*n*(a + b*x**n)*(-a*d + b*c)**2) + c**2*log(a + b*x**n)/(n*(-a*d + b*c)**3) - c**2*log(c + d*x**n)/(n*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_870():
    f = x**13*(b + c*x)**13*(b + 2*c*x)
    F = x**14*(b + c*x)**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_871():
    f = x**27*(b + c*x**2)**13*(b + 2*c*x**2)
    F = x**28*(b + c*x**2)**14/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_872():
    f = x**41*(b + c*x**3)**13*(b + 2*c*x**3)
    F = x**42*(b + c*x**3)**14/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_873():
    f = x**(14*n - 1)*(b + c*x**n)**13*(b + 2*c*x**n)
    F = x**(14*n)*(b + c*x**n)**14/(14*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_874():
    f = x**(m - 1)*(a + b*x**n)**(p - 1)*(a*m + b*x**n*(m + n*p))
    F = x**m*(a + b*x**n)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_875():
    f = (b + 2*c*x**2)/(x*(b + c*x**2))
    F = log(x) + log(b + c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_876():
    f = (b + 2*c*x**3)/(x*(b + c*x**3))
    F = log(x) + log(b + c*x**3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_877():
    f = (b + 2*c*x**n)/(x*(b + c*x**n))
    F = log(x) + log(b + c*x**n)/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_878():
    f = (b + 2*c*x)/(x**8*(b + c*x)**8)
    F = -1/(7*x**7*(b + c*x)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_879():
    f = (b + 2*c*x**2)/(x**15*(b + c*x**2)**8)
    F = -1/(14*x**14*(b + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_880():
    f = (b + 2*c*x**3)/(x**22*(b + c*x**3)**8)
    F = -1/(21*x**21*(b + c*x**3)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_881():
    f = x**(-7*n - 1)*(b + 2*c*x**n)/(b + c*x**n)**8
    F = -1/(7*n*x**(7*n)*(b + c*x**n)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_882():
    f = x**31*sqrt(x**16 + 1)/(1 - x**16)
    F = -(x**16 + 1)**(sympy.S(3)/2)/24 - sqrt(x**16 + 1)/8 + sqrt(2)*atanh(sqrt(2)*sqrt(x**16 + 1)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_883():
    f = sqrt(c + d/x)/(x*sqrt(a + b/x))
    F = -2*sqrt(d)*atanh(sqrt(d)*sqrt(a + b/x)/(sqrt(b)*sqrt(c + d/x)))/sqrt(b) + 2*sqrt(c)*atanh(sqrt(c)*sqrt(a + b/x)/(sqrt(a)*sqrt(c + d/x)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_884():
    f = x**(2*n - 1)*(a + b*x**n)**(sympy.S(5)/2)/sqrt(c + d*x**n)
    F = (a + b*x**n)**(sympy.S(7)/2)*sqrt(c + d*x**n)/(4*b*d*n) - (a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n)*(a*d + 7*b*c)/(24*b*d**2*n) + (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)*(-5*a*d + 5*b*c)*(a*d + 7*b*c)/(96*b*d**3*n) - 5*sqrt(a + b*x**n)*sqrt(c + d*x**n)*(-a*d + b*c)**2*(a*d + 7*b*c)/(64*b*d**4*n) + 5*(-a*d + b*c)**3*(a*d + 7*b*c)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(64*b**(sympy.S(3)/2)*d**(sympy.S(9)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_885():
    f = x**(2*n - 1)*(a + b*x**n)**(sympy.S(3)/2)/sqrt(c + d*x**n)
    F = (a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n)/(3*b*d*n) - (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)*(a*d + 5*b*c)/(12*b*d**2*n) + sqrt(a + b*x**n)*sqrt(c + d*x**n)*(-a*d + b*c)*(a*d + 5*b*c)/(8*b*d**3*n) - (-a*d + b*c)**2*(a*d + 5*b*c)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(8*b**(sympy.S(3)/2)*d**(sympy.S(7)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_886():
    f = x**(2*n - 1)*sqrt(a + b*x**n)/sqrt(c + d*x**n)
    F = (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)/(2*b*d*n) - sqrt(a + b*x**n)*sqrt(c + d*x**n)*(a*d + 3*b*c)/(4*b*d**2*n) + (-a*d + b*c)*(a*d + 3*b*c)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(4*b**(sympy.S(3)/2)*d**(sympy.S(5)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_887():
    f = x**(2*n - 1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n))
    F = sqrt(a + b*x**n)*sqrt(c + d*x**n)/(b*d*n) - (a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(b**(sympy.S(3)/2)*d**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_888():
    f = x**(2*n - 1)/((a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n))
    F = 2*a*sqrt(c + d*x**n)/(b*n*sqrt(a + b*x**n)*(-a*d + b*c)) + 2*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(b**(sympy.S(3)/2)*sqrt(d)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_889():
    f = x**(2*n - 1)/((a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n))
    F = 2*a*sqrt(c + d*x**n)/(3*b*n*(a + b*x**n)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt(c + d*x**n)*(-2*a*d + 6*b*c)/(3*b*n*sqrt(a + b*x**n)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_890():
    f = x**(3*n - 1)*(a + b*x**n)**(sympy.S(5)/2)/sqrt(c + d*x**n)
    F = x**n*(a + b*x**n)**(sympy.S(7)/2)*sqrt(c + d*x**n)/(5*b*d*n) - (a + b*x**n)**(sympy.S(7)/2)*sqrt(c + d*x**n)*(3*a*d + 9*b*c)/(40*b**2*d**2*n) + (a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n)*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(240*b**2*d**3*n) - (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)*(-a*d + b*c)*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(192*b**2*d**4*n) + sqrt(a + b*x**n)*sqrt(c + d*x**n)*(-a*d + b*c)**2*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)/(128*b**2*d**5*n) - (-a*d + b*c)**3*(3*a**2*d**2 + 14*a*b*c*d + 63*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(128*b**(sympy.S(5)/2)*d**(sympy.S(11)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_891():
    f = x**(3*n - 1)*(a + b*x**n)**(sympy.S(3)/2)/sqrt(c + d*x**n)
    F = x**n*(a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n)/(4*b*d*n) - (a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n)*(3*a*d + 7*b*c)/(24*b**2*d**2*n) + (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)/(96*b**2*d**3*n) - sqrt(a + b*x**n)*sqrt(c + d*x**n)*(-a*d + b*c)*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)/(64*b**2*d**4*n) + (-a*d + b*c)**2*(3*a**2*d**2 + 10*a*b*c*d + 35*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(64*b**(sympy.S(5)/2)*d**(sympy.S(9)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_892():
    f = x**(3*n - 1)*sqrt(a + b*x**n)/sqrt(c + d*x**n)
    F = x**n*(a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)/(3*b*d*n) - (a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n)*(3*a*d + 5*b*c)/(12*b**2*d**2*n) + sqrt(a + b*x**n)*sqrt(c + d*x**n)*(a**2*d**2 + 2*a*b*c*d + 5*b**2*c**2)/(8*b**2*d**3*n) - (-a*d + b*c)*(a**2*d**2 + 2*a*b*c*d + 5*b**2*c**2)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(8*b**(sympy.S(5)/2)*d**(sympy.S(7)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_893():
    f = x**(3*n - 1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n))
    F = x**n*sqrt(a + b*x**n)*sqrt(c + d*x**n)/(2*b*d*n) - sqrt(a + b*x**n)*sqrt(c + d*x**n)*(3*a*d + 3*b*c)/(4*b**2*d**2*n) - (4*a*b*c*d - 3*(a*d + b*c)**2)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(4*b**(sympy.S(5)/2)*d**(sympy.S(5)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_894():
    f = x**(3*n - 1)/((a + b*x**n)**(sympy.S(3)/2)*sqrt(c + d*x**n))
    F = -2*a**2*sqrt(c + d*x**n)/(b**2*n*sqrt(a + b*x**n)*(-a*d + b*c)) + sqrt(a + b*x**n)*sqrt(c + d*x**n)/(b**2*d*n) - (3*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(b**(sympy.S(5)/2)*d**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_895():
    f = x**(3*n - 1)/((a + b*x**n)**(sympy.S(5)/2)*sqrt(c + d*x**n))
    F = -2*a**2*sqrt(c + d*x**n)/(3*b**2*n*(a + b*x**n)**(sympy.S(3)/2)*(-a*d + b*c)) + 4*a*sqrt(c + d*x**n)*(-2*a*d + 3*b*c)/(3*b**2*n*sqrt(a + b*x**n)*(-a*d + b*c)**2) + 2*atanh(sqrt(d)*sqrt(a + b*x**n)/(sqrt(b)*sqrt(c + d*x**n)))/(b**(sympy.S(5)/2)*sqrt(d)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_896():
    f = x**p*(b + c*x)**p*(b + 2*c*x)
    F = x**(p + 1)*(b + c*x)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_897():
    f = x**(2*p + 1)*(b + c*x**2)**p*(b + 2*c*x**2)
    F = x**(2*p + 2)*(b + c*x**2)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_898():
    f = x**(3*p + 2)*(b + c*x**3)**p*(b + 2*c*x**3)
    F = x**(3*p + 3)*(b + c*x**3)**(p + 1)/(3*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_4_e_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_899():
    f = x**(n*(p + 1) - 1)*(b + c*x**n)**p*(b + 2*c*x**n)
    F = x**(n*(p + 1))*(b + c*x**n)**(p + 1)/(n*(p + 1))
    assert integrate(f, x) == F

