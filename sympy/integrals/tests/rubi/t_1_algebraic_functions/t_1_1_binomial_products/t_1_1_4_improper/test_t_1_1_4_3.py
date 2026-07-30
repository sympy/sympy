"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.4 Improper/1.1.4.3 (e x)^m (a x^j+b x^k)^p (c+d x^n)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, j, m, n, p, q = symbols('A B a b c d e j m n p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_1():
    f = x**2*(A + B*x**2)*(b*x**2 + c*x**4)
    F = A*b*x**5/5 + B*c*x**9/9 + x**7*(A*c/7 + B*b/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_2():
    f = x*(A + B*x**2)*(b*x**2 + c*x**4)
    F = A*b*x**4/4 + B*c*x**8/8 + x**6*(A*c/6 + B*b/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_3():
    f = (A + B*x**2)*(b*x**2 + c*x**4)
    F = A*b*x**3/3 + B*c*x**7/7 + x**5*(A*c/5 + B*b/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_4():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x
    F = A*b*x**2/2 + B*c*x**6/6 + x**4*(A*c/4 + B*b/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_5():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**2
    F = A*b*x + B*c*x**5/5 + x**3*(A*c/3 + B*b/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_6():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**3
    F = A*b*log(x) + B*c*x**4/4 + x**2*(A*c/2 + B*b/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_7():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**4
    F = -A*b/x + B*c*x**3/3 + x*(A*c + B*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_8():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**5
    F = -A*b/(2*x**2) + B*c*x**2/2 + (A*c + B*b)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_9():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**6
    F = -A*b/(3*x**3) + B*c*x - (A*c + B*b)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_10():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**7
    F = -A*b/(4*x**4) + B*c*log(x) - (A*c + B*b)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_11():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**8
    F = -A*b/(5*x**5) - B*c/x - (A*c + B*b)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_12():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2
    F = A*b**2*x**5/5 + B*c**2*x**11/11 + b*x**7*(2*A*c + B*b)/7 + c*x**9*(A*c + 2*B*b)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_13():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x
    F = A*b**2*x**4/4 + B*c**2*x**10/10 + b*x**6*(2*A*c + B*b)/6 + c*x**8*(A*c + 2*B*b)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_14():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**2
    F = A*b**2*x**3/3 + B*c**2*x**9/9 + b*x**5*(2*A*c + B*b)/5 + c*x**7*(A*c + 2*B*b)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_15():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**3
    F = B*(b + c*x**2)**4/(8*c**2) - (b + c*x**2)**3*(-A*c + B*b)/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_16():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**4
    F = A*b**2*x + B*c**2*x**7/7 + b*x**3*(2*A*c + B*b)/3 + c*x**5*(A*c + 2*B*b)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_17():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**5
    F = A*b**2*log(x) + A*b*c*x**2 + A*c**2*x**4/4 + B*(b + c*x**2)**3/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_18():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**6
    F = -A*b**2/x + B*c**2*x**5/5 + b*x*(2*A*c + B*b) + c*x**3*(A*c + 2*B*b)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_19():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**7
    F = -A*b**2/(2*x**2) + B*c**2*x**4/4 + b*(2*A*c + B*b)*log(x) + c*x**2*(A*c + 2*B*b)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_20():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**8
    F = -A*b**2/(3*x**3) + B*c**2*x**3/3 - b*(2*A*c + B*b)/x + c*x*(A*c + 2*B*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_21():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**9
    F = -A*b**2/(4*x**4) + B*c**2*x**2/2 - b*(2*A*c + B*b)/(2*x**2) + c*(A*c + 2*B*b)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_22():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**10
    F = -A*b**2/(5*x**5) + B*c**2*x - b*(2*A*c + B*b)/(3*x**3) - c*(A*c + 2*B*b)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_23():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**11
    F = -A*b**2/(6*x**6) + B*c**2*log(x) - b*(2*A*c + B*b)/(4*x**4) - c*(A*c + 2*B*b)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_24():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**12
    F = -A*b**2/(7*x**7) - B*c**2/x - b*(2*A*c + B*b)/(5*x**5) - c*(A*c + 2*B*b)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_25():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**2
    F = A*b**3*x**5/5 + B*c**3*x**13/13 + b**2*x**7*(3*A*c + B*b)/7 + b*c*x**9*(A*c + B*b)/3 + c**2*x**11*(A*c + 3*B*b)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_26():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**3
    F = B*(b + c*x**2)**6/(12*c**3) + b*(b + c*x**2)**4*(-A*c + B*b)/(8*c**3) - (b + c*x**2)**5*(-A*c + 2*B*b)/(10*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_27():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**4
    F = A*b**3*x**3/3 + B*c**3*x**11/11 + b**2*x**5*(3*A*c + B*b)/5 + 3*b*c*x**7*(A*c + B*b)/7 + c**2*x**9*(A*c + 3*B*b)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_28():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**5
    F = B*(b + c*x**2)**5/(10*c**2) - (b + c*x**2)**4*(-A*c + B*b)/(8*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_29():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**6
    F = A*b**3*x + B*c**3*x**9/9 + b**2*x**3*(3*A*c + B*b)/3 + 3*b*c*x**5*(A*c + B*b)/5 + c**2*x**7*(A*c + 3*B*b)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_30():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**7
    F = A*b**3*log(x) + 3*A*b**2*c*x**2/2 + 3*A*b*c**2*x**4/4 + A*c**3*x**6/6 + B*(b + c*x**2)**4/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_31():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**8
    F = -A*b**3/x + B*c**3*x**7/7 + b**2*x*(3*A*c + B*b) + b*c*x**3*(A*c + B*b) + c**2*x**5*(A*c + 3*B*b)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_32():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**9
    F = -A*b**3/(2*x**2) + B*c**3*x**6/6 + b**2*(3*A*c + B*b)*log(x) + 3*b*c*x**2*(A*c + B*b)/2 + c**2*x**4*(A*c + 3*B*b)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_33():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**10
    F = -A*b**3/(3*x**3) + B*c**3*x**5/5 - b**2*(3*A*c + B*b)/x + 3*b*c*x*(A*c + B*b) + c**2*x**3*(A*c + 3*B*b)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_34():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**11
    F = -A*b**3/(4*x**4) + B*c**3*x**4/4 - b**2*(3*A*c + B*b)/(2*x**2) + 3*b*c*(A*c + B*b)*log(x) + c**2*x**2*(A*c + 3*B*b)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_35():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**12
    F = -A*b**3/(5*x**5) + B*c**3*x**3/3 - b**2*(3*A*c + B*b)/(3*x**3) - 3*b*c*(A*c + B*b)/x + c**2*x*(A*c + 3*B*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_36():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**13
    F = -A*b**3/(6*x**6) + B*c**3*x**2/2 - b**2*(3*A*c + B*b)/(4*x**4) - 3*b*c*(A*c + B*b)/(2*x**2) + c**2*(A*c + 3*B*b)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_37():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**14
    F = -A*b**3/(7*x**7) + B*c**3*x - b**2*(3*A*c + B*b)/(5*x**5) - b*c*(A*c + B*b)/x**3 - c**2*(A*c + 3*B*b)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_38():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**15
    F = -A*(b + c*x**2)**4/(8*b*x**8) - B*b**3/(6*x**6) - 3*B*b**2*c/(4*x**4) - 3*B*b*c**2/(2*x**2) + B*c**3*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_39():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**16
    F = -A*b**3/(9*x**9) - B*c**3/x - b**2*(3*A*c + B*b)/(7*x**7) - 3*b*c*(A*c + B*b)/(5*x**5) - c**2*(A*c + 3*B*b)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_40():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**17
    F = -A*(b + c*x**2)**4/(10*b*x**10) - (b + c*x**2)**4*(-A*c + 5*B*b)/(40*b**2*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_41():
    f = x**10*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**9/(9*c) - b**(sympy.S(7)/2)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(11)/2) + b**3*x*(-A*c + B*b)/c**5 - b**2*x**3*(-A*c + B*b)/(3*c**4) + b*x**5*(-A*c + B*b)/(5*c**3) - x**7*(-A*c + B*b)/(7*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_42():
    f = x**9*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**8/(8*c) + b**3*(-A*c + B*b)*log(b + c*x**2)/(2*c**5) - b**2*x**2*(-A*c + B*b)/(2*c**4) + b*x**4*(-A*c + B*b)/(4*c**3) - x**6*(-A*c + B*b)/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_43():
    f = x**8*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**7/(7*c) + b**(sympy.S(5)/2)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(9)/2) - b**2*x*(-A*c + B*b)/c**4 + b*x**3*(-A*c + B*b)/(3*c**3) - x**5*(-A*c + B*b)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_44():
    f = x**7*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**6/(6*c) - b**2*(-A*c + B*b)*log(b + c*x**2)/(2*c**4) + b*x**2*(-A*c + B*b)/(2*c**3) - x**4*(-A*c + B*b)/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_45():
    f = x**6*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**5/(5*c) - b**(sympy.S(3)/2)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(7)/2) + b*x*(-A*c + B*b)/c**3 - x**3*(-A*c + B*b)/(3*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_46():
    f = x**5*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**4/(4*c) + b*(-A*c + B*b)*log(b + c*x**2)/(2*c**3) - x**2*(-A*c + B*b)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_47():
    f = x**4*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**3/(3*c) + sqrt(b)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/c**(sympy.S(5)/2) - x*(-A*c + B*b)/c**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_48():
    f = x**3*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x**2/(2*c) - (-A*c + B*b)*log(b + c*x**2)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_49():
    f = x**2*(A + B*x**2)/(b*x**2 + c*x**4)
    F = B*x/c - (-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/(sqrt(b)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_50():
    f = x*(A + B*x**2)/(b*x**2 + c*x**4)
    F = A*log(x)/b + (-A*c + B*b)*log(b + c*x**2)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_51():
    f = (A + B*x**2)/(b*x**2 + c*x**4)
    F = -A/(b*x) + (-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/(b**(sympy.S(3)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_52():
    f = (A + B*x**2)/(b*x**2 - c*x**4)
    F = -A/(b*x) + (A*c + B*b)*atanh(sqrt(c)*x/sqrt(b))/(b**(sympy.S(3)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_53():
    f = (A + B*x**2)/(x*(b*x**2 + c*x**4))
    F = -A/(2*b*x**2) + (-A*c + B*b)*log(x)/b**2 - (-A*c + B*b)*log(b + c*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_54():
    f = (A + B*x**2)/(x**2*(b*x**2 + c*x**4))
    F = -A/(3*b*x**3) - (-A*c + B*b)/(b**2*x) - sqrt(c)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_55():
    f = (A + B*x**2)/(x**3*(b*x**2 + c*x**4))
    F = -A/(4*b*x**4) - (-A*c + B*b)/(2*b**2*x**2) - c*(-A*c + B*b)*log(x)/b**3 + c*(-A*c + B*b)*log(b + c*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_56():
    f = (A + B*x**2)/(x**4*(b*x**2 + c*x**4))
    F = -A/(5*b*x**5) - (-A*c + B*b)/(3*b**2*x**3) + c*(-A*c + B*b)/(b**3*x) + c**(sympy.S(3)/2)*(-A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_57():
    f = (A + B*x**2)/(x**5*(b*x**2 + c*x**4))
    F = -A/(6*b*x**6) - (-A*c + B*b)/(4*b**2*x**4) + c*(-A*c + B*b)/(2*b**3*x**2) + c**2*(-A*c + B*b)*log(x)/b**4 - c**2*(-A*c + B*b)*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_58():
    f = x**12*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**7/(7*c**2) + b**(sympy.S(5)/2)*(-7*A*c + 9*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(11)/2)) - b**3*x*(-A*c + B*b)/(2*c**5*(b + c*x**2)) - b**2*x*(-3*A*c + 4*B*b)/c**5 + b*x**3*(-2*A*c + 3*B*b)/(3*c**4) - x**5*(-A*c + 2*B*b)/(5*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_59():
    f = x**11*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**6/(6*c**2) - b**3*(-A*c + B*b)/(2*c**5*(b + c*x**2)) - b**2*(-3*A*c + 4*B*b)*log(b + c*x**2)/(2*c**5) + b*x**2*(-2*A*c + 3*B*b)/(2*c**4) - x**4*(-A*c + 2*B*b)/(4*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_60():
    f = x**10*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**5/(5*c**2) - b**(sympy.S(3)/2)*(-5*A*c + 7*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(9)/2)) + b**2*x*(-A*c + B*b)/(2*c**4*(b + c*x**2)) + b*x*(-2*A*c + 3*B*b)/c**4 - x**3*(-A*c + 2*B*b)/(3*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_61():
    f = x**9*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**4/(4*c**2) + b**2*(-A*c + B*b)/(2*c**4*(b + c*x**2)) + b*(-2*A*c + 3*B*b)*log(b + c*x**2)/(2*c**4) - x**2*(-A*c + 2*B*b)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_62():
    f = x**8*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**3/(3*c**2) + sqrt(b)*(-3*A*c + 5*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*c**(sympy.S(7)/2)) - b*x*(-A*c + B*b)/(2*c**3*(b + c*x**2)) - x*(-A*c + 2*B*b)/c**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_63():
    f = x**7*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x**2/(2*c**2) - b*(-A*c + B*b)/(2*c**3*(b + c*x**2)) - (-A*c + 2*B*b)*log(b + c*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_64():
    f = x**6*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*x/c**2 + x*(-A*c + B*b)/(2*c**2*(b + c*x**2)) - (-A*c + 3*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*sqrt(b)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_65():
    f = x**5*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = B*log(b + c*x**2)/(2*c**2) + (-A*c + B*b)/(2*c**2*(b + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_66():
    f = x**4*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -x*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + (A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_67():
    f = x**3*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = A*log(x)/b**2 - A*log(b + c*x**2)/(2*b**2) - (-A*c + B*b)/(2*b*c*(b + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_68():
    f = x**2*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -A/(b**2*x) + x*(-A*c + B*b)/(2*b**2*(b + c*x**2)) + (-3*A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(5)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_69():
    f = x*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -A/(2*b**2*x**2) + (-A*c + B*b)/(2*b**2*(b + c*x**2)) + (-2*A*c + B*b)*log(x)/b**3 - (-2*A*c + B*b)*log(b + c*x**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_70():
    f = (A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -A/(3*b**2*x**3) - c*x*(-A*c + B*b)/(2*b**3*(b + c*x**2)) - (-2*A*c + B*b)/(b**3*x) - sqrt(c)*(-5*A*c + 3*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_71():
    f = (A + B*x**2)/(x*(b*x**2 + c*x**4)**2)
    F = -A/(4*b**2*x**4) - c*(-A*c + B*b)/(2*b**3*(b + c*x**2)) - (-2*A*c + B*b)/(2*b**3*x**2) - c*(-3*A*c + 2*B*b)*log(x)/b**4 + c*(-3*A*c + 2*B*b)*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_72():
    f = (A + B*x**2)/(x**2*(b*x**2 + c*x**4)**2)
    F = -A/(5*b**2*x**5) - (-2*A*c + B*b)/(3*b**3*x**3) + c**2*x*(-A*c + B*b)/(2*b**4*(b + c*x**2)) + c*(-3*A*c + 2*B*b)/(b**4*x) + c**(sympy.S(3)/2)*(-7*A*c + 5*B*b)*atan(sqrt(c)*x/sqrt(b))/(2*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_73():
    f = x**14*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*x**5/(5*c**3) - 7*b**(sympy.S(3)/2)*(-5*A*c + 9*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*c**(sympy.S(11)/2)) - b**3*x*(-A*c + B*b)/(4*c**5*(b + c*x**2)**2) + b**2*x*(-13*A*c + 17*B*b)/(8*c**5*(b + c*x**2)) + 3*b*x*(-A*c + 2*B*b)/c**5 - x**3*(-A*c + 3*B*b)/(3*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_74():
    f = x**13*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*x**4/(4*c**3) - b**3*(-A*c + B*b)/(4*c**5*(b + c*x**2)**2) + b**2*(-3*A*c + 4*B*b)/(2*c**5*(b + c*x**2)) + 3*b*(-A*c + 2*B*b)*log(b + c*x**2)/(2*c**5) - x**2*(-A*c + 3*B*b)/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_75():
    f = x**12*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*x**3/(3*c**3) + 5*sqrt(b)*(-3*A*c + 7*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*c**(sympy.S(9)/2)) + b**2*x*(-A*c + B*b)/(4*c**4*(b + c*x**2)**2) - b*x*(-9*A*c + 13*B*b)/(8*c**4*(b + c*x**2)) - x*(-A*c + 3*B*b)/c**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_76():
    f = x**11*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*x**2/(2*c**3) + b**2*(-A*c + B*b)/(4*c**4*(b + c*x**2)**2) - b*(-2*A*c + 3*B*b)/(2*c**4*(b + c*x**2)) - (-A*c + 3*B*b)*log(b + c*x**2)/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_77():
    f = x**10*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*x/c**3 - b*x*(-A*c + B*b)/(4*c**3*(b + c*x**2)**2) + x*(-5*A*c + 9*B*b)/(8*c**3*(b + c*x**2)) - (-3*A*c + 15*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*sqrt(b)*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_78():
    f = x**9*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = B*log(b + c*x**2)/(2*c**3) - b*(-A*c + B*b)/(4*c**3*(b + c*x**2)**2) + (-A*c + 2*B*b)/(2*c**3*(b + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_79():
    f = x**8*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = x*(-A*c + B*b)/(4*c**2*(b + c*x**2)**2) - x*(-A*c + 5*B*b)/(8*b*c**2*(b + c*x**2)) + (A*c + 3*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(3)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_80():
    f = x**7*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = (A + B*x**2)**2/((b + c*x**2)**2*(-4*A*c + 4*B*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_81():
    f = x**6*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) + x*(3*A*c + B*b)/(8*b**2*c*(b + c*x**2)) + (3*A*c + B*b)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(5)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_82():
    f = x**5*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = A/(2*b**2*(b + c*x**2)) + A*log(x)/b**3 - A*log(b + c*x**2)/(2*b**3) - (-A*c + B*b)/(4*b*c*(b + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_83():
    f = x**4*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -A/(b**3*x) + x*(-A*c + B*b)/(4*b**2*(b + c*x**2)**2) + x*(-7*A*c + 3*B*b)/(8*b**3*(b + c*x**2)) + (-15*A*c + 3*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(7)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_84():
    f = x**3*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -A/(2*b**3*x**2) + (-A*c + B*b)/(4*b**2*(b + c*x**2)**2) + (-2*A*c + B*b)/(2*b**3*(b + c*x**2)) + (-3*A*c + B*b)*log(x)/b**4 - (-3*A*c + B*b)*log(b + c*x**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_85():
    f = x**2*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -A/(3*b**3*x**3) - c*x*(-A*c + B*b)/(4*b**3*(b + c*x**2)**2) - c*x*(-11*A*c + 7*B*b)/(8*b**4*(b + c*x**2)) - (-3*A*c + B*b)/(b**4*x) - 5*sqrt(c)*(-7*A*c + 3*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_86():
    f = x*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -A/(4*b**3*x**4) - c*(-A*c + B*b)/(4*b**3*(b + c*x**2)**2) - c*(-3*A*c + 2*B*b)/(2*b**4*(b + c*x**2)) - (-3*A*c + B*b)/(2*b**4*x**2) - 3*c*(-2*A*c + B*b)*log(x)/b**5 + 3*c*(-2*A*c + B*b)*log(b + c*x**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_87():
    f = (A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -A/(5*b**3*x**5) + c**2*x*(-A*c + B*b)/(4*b**4*(b + c*x**2)**2) - (-3*A*c + B*b)/(3*b**4*x**3) + c**2*x*(-15*A*c + 11*B*b)/(8*b**5*(b + c*x**2)) + 3*c*(-2*A*c + B*b)/(b**5*x) + 7*c**(sympy.S(3)/2)*(-9*A*c + 5*B*b)*atan(sqrt(c)*x/sqrt(b))/(8*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_88():
    f = (A + B*x**2)/(x*(b*x**2 + c*x**4)**3)
    F = -A/(6*b**3*x**6) + c**2*(-A*c + B*b)/(4*b**4*(b + c*x**2)**2) - (-3*A*c + B*b)/(4*b**4*x**4) + c**2*(-4*A*c + 3*B*b)/(2*b**5*(b + c*x**2)) + 3*c*(-2*A*c + B*b)/(2*b**5*x**2) + 2*c**2*(-5*A*c + 3*B*b)*log(x)/b**6 - c**2*(-5*A*c + 3*B*b)*log(b + c*x**2)/b**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_89():
    f = x**7*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*x**6*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(12*c) - 7*b**5*(-4*A*c + 3*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(1024*c**(sympy.S(11)/2)) + 7*b**3*(b + 2*c*x**2)*(-4*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(1024*c**5) - 7*b**2*(-4*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(384*c**4) + 7*b*x**2*(-4*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(320*c**3) - x**4*(-4*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(40*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_90():
    f = x**5*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*x**4*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(10*c) + b**4*(-10*A*c + 7*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(256*c**(sympy.S(9)/2)) - b**2*(b + 2*c*x**2)*(-10*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(256*c**4) + b*(-10*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(96*c**3) - x**2*(-10*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(80*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_91():
    f = x**3*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = -b**3*(-8*A*c + 5*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(128*c**(sympy.S(7)/2)) + b*(b + 2*c*x**2)*(-8*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(128*c**3) - (b*x**2 + c*x**4)**(sympy.S(3)/2)*(-8*A*c + 5*B*b - 6*B*c*x**2)/(48*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_92():
    f = x*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*c) + b**2*(-2*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(5)/2)) - (b + 2*c*x**2)*(-2*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(16*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_93():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x
    F = B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*c*x**2) - b*(-4*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(3)/2)) - (-4*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_94():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**3
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(b*x**4) + (2*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*sqrt(c)) + (2*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_95():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**5
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*b*x**6) + B*sqrt(c)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4)) - B*sqrt(b*x**2 + c*x**4)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_96():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**7
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*b*x**8) - (-2*A*c + 5*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*b**2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_97():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**9
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*b*x**10) - (-4*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(35*b**2*x**8) + 2*c*(-4*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*b**3*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_98():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**11
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*b*x**12) - (-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(21*b**2*x**10) + 4*c*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*b**3*x**8) - 8*c**2*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(315*b**4*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_99():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**13
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*b*x**14) - (-8*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(99*b**2*x**12) + 2*c*(-8*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(231*b**3*x**10) - 8*c**2*(-8*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(1155*b**4*x**8) + 16*c**3*(-8*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3465*b**5*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_100():
    f = x**4*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*x**3*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*c) - 8*b**2*(-3*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(315*c**4*x**3) + 4*b*(-3*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*c**3*x) - x*(-3*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(21*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_101():
    f = x**2*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*x*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*c) + 2*b*(-7*A*c + 4*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*c**3*x**3) - (-7*A*c + 4*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(35*c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_102():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*c*x) - (-5*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*c**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_103():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**2
    F = -A*sqrt(b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4)) + A*sqrt(b*x**2 + c*x**4)/x + B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_104():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**4
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(2*b*x**5) + (A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(2*b*x) - (A*c + 2*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_105():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**6
    F = -A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*b*x**7) - (-A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(8*b*x**3) - c*(-A*c + 4*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_106():
    f = x**5*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*x**4*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(14*c) - b**6*(-14*A*c + 9*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2048*c**(sympy.S(11)/2)) + b**4*(b + 2*c*x**2)*(-14*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(2048*c**5) - b**2*(b + 2*c*x**2)*(-14*A*c + 9*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(768*c**4) + b*(-14*A*c + 9*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(240*c**3) - x**2*(-14*A*c + 9*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(168*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_107():
    f = x**3*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = b**5*(-12*A*c + 7*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(1024*c**(sympy.S(9)/2)) - b**3*(b + 2*c*x**2)*(-12*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(1024*c**4) + b*(b + 2*c*x**2)*(-12*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(384*c**3) - (b*x**2 + c*x**4)**(sympy.S(5)/2)*(-12*A*c + 7*B*b - 10*B*c*x**2)/(120*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_108():
    f = x*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*c) - 3*b**4*(-2*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(256*c**(sympy.S(7)/2)) + 3*b**2*(b + 2*c*x**2)*(-2*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(256*c**3) - (b + 2*c*x**2)*(-2*A*c + B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(32*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_109():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x
    F = B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(8*c*x**2) + b**3*(-8*A*c + 3*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(128*c**(sympy.S(5)/2)) - b*(b + 2*c*x**2)*(-8*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(128*c**2) - (-8*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(48*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_110():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**3
    F = A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(b*x**4) - b**2*(-6*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(3)/2)) + (b + 2*c*x**2)*(-6*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(16*c) + (-6*A*c + B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_111():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**5
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(b*x**6) + 3*b*(4*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*sqrt(c)) + (3*A*c/2 + 3*B*b/8)*sqrt(b*x**2 + c*x**4) + (4*A*c + B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(4*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_112():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**7
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3*b*x**8) + sqrt(c)*(2*A*c + 3*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/2 + c*(2*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(2*b) - (2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_113():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**9
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(5*b*x**10) + B*c**(sympy.S(3)/2)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4)) - B*c*sqrt(b*x**2 + c*x**4)/x**2 - B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_114():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**11
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(7*b*x**12) - (-2*A*c + 7*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(35*b**2*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_115():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**13
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(9*b*x**14) - (-4*A*c + 9*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(63*b**2*x**12) + 2*c*(-4*A*c + 9*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(315*b**3*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_116():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**15
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(11*b*x**16) - (-6*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(99*b**2*x**14) + 4*c*(-6*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(693*b**3*x**12) - 8*c**2*(-6*A*c + 11*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3465*b**4*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_117():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**17
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(13*b*x**18) - (-8*A*c + 13*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(143*b**2*x**16) + 2*c*(-8*A*c + 13*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(429*b**3*x**14) - 8*c**2*(-8*A*c + 13*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3003*b**4*x**12) + 16*c**3*(-8*A*c + 13*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15015*b**5*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_118():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**19
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15*b*x**20) - (-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(39*b**2*x**18) + 8*c*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(429*b**3*x**16) - 16*c**2*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(1287*b**4*x**14) + 64*c**3*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(9009*b**5*x**12) - 128*c**4*(-2*A*c + 3*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(45045*b**6*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_119():
    f = x**4*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*x**3*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(13*c) + 16*b**3*(-13*A*c + 8*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15015*c**5*x**5) - 8*b**2*(-13*A*c + 8*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3003*c**4*x**3) + 2*b*(-13*A*c + 8*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(429*c**3*x) - x*(-13*A*c + 8*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(143*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_120():
    f = x**2*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*x*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(11*c) - 8*b**2*(-11*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3465*c**4*x**5) + 4*b*(-11*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(693*c**3*x**3) - (-11*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(99*c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_121():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(9*c*x) + 2*b*(-9*A*c + 4*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(315*c**3*x**5) - (-9*A*c + 4*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(63*c**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_122():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**2
    F = B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(7*c*x**3) - (-7*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(35*c**2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_123():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**4
    F = -A*b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4)) + A*b*sqrt(b*x**2 + c*x**4)/x + A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*x**3) + B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(5*c*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_124():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**6
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(2*b*x**7) - sqrt(b)*(3*A*c + 2*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/2 + (3*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(2*x) + (3*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_125():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**8
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(4*b*x**9) + 3*c*(A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(8*b*x) - (A*c + 4*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(8*b*x**5) - 3*c*(A*c + 4*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_126():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**10
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(6*b*x**11) - c*(-A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(16*b*x**3) - (-A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(24*b*x**7) - c**2*(-A*c + 6*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(16*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_127():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**12
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(8*b*x**13) - c*(-3*A*c + 8*B*b)*sqrt(b*x**2 + c*x**4)/(64*b*x**5) - (-3*A*c + 8*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(48*b*x**9) - c**2*(-3*A*c + 8*B*b)*sqrt(b*x**2 + c*x**4)/(128*b**2*x**3) + c**3*(-3*A*c + 8*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(128*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_128():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**14
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*b*x**15) - c*(-A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(32*b*x**7) - (-A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(16*b*x**11) - c**2*(-A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(128*b**2*x**5) + 3*c**3*(-A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(256*b**3*x**3) - 3*c**4*(-A*c + 2*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(256*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_129():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**16
    F = -A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(12*b*x**17) - c*(-7*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(320*b*x**9) - (-7*A*c + 12*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(120*b*x**13) - c**2*(-7*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(1920*b**2*x**7) + c**3*(-7*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(1536*b**3*x**5) - c**4*(-7*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(1024*b**4*x**3) + c**5*(-7*A*c + 12*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(1024*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_130():
    f = x**7*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*x**6*sqrt(b*x**2 + c*x**4)/(8*c) + 5*b**3*(-8*A*c + 7*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(128*c**(sympy.S(9)/2)) - 5*b**2*(-8*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(128*c**4) + 5*b*x**2*(-8*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(192*c**3) - x**4*(-8*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(48*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_131():
    f = x**5*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*x**4*sqrt(b*x**2 + c*x**4)/(6*c) - b**2*(-6*A*c + 5*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(7)/2)) + b*(-6*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(16*c**3) - x**2*(-6*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(24*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_132():
    f = x**3*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = b*(-4*A*c + 3*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(5)/2)) - sqrt(b*x**2 + c*x**4)*(-4*A*c + 3*B*b - 2*B*c*x**2)/(8*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_133():
    f = x*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*sqrt(b*x**2 + c*x**4)/(2*c) - (-2*A*c + B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_134():
    f = (A + B*x**2)/(x*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(b*x**2) + B*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_135():
    f = (A + B*x**2)/(x**3*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(3*b*x**4) - (-2*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(3*b**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_136():
    f = (A + B*x**2)/(x**5*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(5*b*x**6) - (-4*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(15*b**2*x**4) + 2*c*(-4*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(15*b**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_137():
    f = (A + B*x**2)/(x**7*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(7*b*x**8) - (-6*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(35*b**2*x**6) + 4*c*(-6*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(105*b**3*x**4) - 8*c**2*(-6*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(105*b**4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_138():
    f = (A + B*x**2)/(x**9*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(9*b*x**10) - (-8*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(63*b**2*x**8) + 2*c*(-8*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(105*b**3*x**6) - 8*c**2*(-8*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(315*b**4*x**4) + 16*c**3*(-8*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(315*b**5*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_139():
    f = x**6*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*x**5*sqrt(b*x**2 + c*x**4)/(7*c) - 8*b**2*(-7*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(105*c**4*x) + 4*b*x*(-7*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(105*c**3) - x**3*(-7*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(35*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_140():
    f = x**4*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*x**3*sqrt(b*x**2 + c*x**4)/(5*c) + 2*b*(-5*A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(15*c**3*x) - x*(-5*A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(15*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_141():
    f = x**2*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = B*x*sqrt(b*x**2 + c*x**4)/(3*c) - (-3*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(3*c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_142():
    f = (A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = -A*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/sqrt(b) + B*sqrt(b*x**2 + c*x**4)/(c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_143():
    f = (A + B*x**2)/(x**2*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(2*b*x**3) - (-A*c + 2*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_144():
    f = (A + B*x**2)/(x**4*sqrt(b*x**2 + c*x**4))
    F = -A*sqrt(b*x**2 + c*x**4)/(4*b*x**5) - (-3*A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(8*b**2*x**3) + c*(-3*A*c + 4*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_145():
    f = x**9*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -5*b**2*(-6*A*c + 7*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(16*c**(sympy.S(9)/2)) + 5*b*(-6*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(16*c**4) - x**2*(-30*A*c + 35*B*b)*sqrt(b*x**2 + c*x**4)/(24*c**3) - x**8*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**4*(-6*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(6*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_146():
    f = x**7*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 3*b*(-4*A*c + 5*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(8*c**(sympy.S(7)/2)) - (-12*A*c + 15*B*b)*sqrt(b*x**2 + c*x**4)/(8*c**3) - x**6*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**2*(-4*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(4*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_147():
    f = x**5*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(-2*A*c + 3*B*b)*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/(2*c**(sympy.S(5)/2)) - x**4*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + (-2*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(2*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_148():
    f = x**3*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = B*atanh(sqrt(c)*x**2/sqrt(b*x**2 + c*x**4))/c**(sympy.S(3)/2) - x**2*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_149():
    f = x*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(A*b - x**2*(-2*A*c + B*b))/(b**2*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_150():
    f = (A + B*x**2)/(x*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -A/(3*b*x**2*sqrt(b*x**2 + c*x**4)) - (b + 2*c*x**2)*(-4*A*c + 3*B*b)/(3*b**3*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_151():
    f = (A + B*x**2)/(x**3*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -A/(5*b*x**4*sqrt(b*x**2 + c*x**4)) - (-6*A*c + 5*B*b)/(15*b**2*x**2*sqrt(b*x**2 + c*x**4)) + 4*c*(b + 2*c*x**2)*(-6*A*c + 5*B*b)/(15*b**4*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_152():
    f = (A + B*x**2)/(x**5*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -A/(7*b*x**6*sqrt(b*x**2 + c*x**4)) - (-8*A*c + 7*B*b)/(35*b**2*x**4*sqrt(b*x**2 + c*x**4)) + 2*c*(-8*A*c + 7*B*b)/(35*b**3*x**2*sqrt(b*x**2 + c*x**4)) - 8*c**2*(b + 2*c*x**2)*(-8*A*c + 7*B*b)/(35*b**5*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_153():
    f = x**8*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 8*b*(-5*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(15*c**4*x) - x*(-20*A*c + 24*B*b)*sqrt(b*x**2 + c*x**4)/(15*c**3) - x**7*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**3*(-5*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(5*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_154():
    f = x**6*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(-6*A*c + 8*B*b)*sqrt(b*x**2 + c*x**4)/(3*c**3*x) - x**5*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x*(-3*A*c + 4*B*b)*sqrt(b*x**2 + c*x**4)/(3*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_155():
    f = x**4*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**3*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + (-A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(b*c**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_156():
    f = x**2*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -A*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/b**(sympy.S(3)/2) - x*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_157():
    f = (A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -B/(3*c*x*sqrt(b*x**2 + c*x**4)) - (-3*A*c + 2*B*b)/(3*b*c*x*sqrt(b*x**2 + c*x**4)) + (-3*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(2*b**2*c*x**3) - (-3*A*c + 2*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(2*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_158():
    f = (A + B*x**2)/(x**2*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -A/(4*b*x**3*sqrt(b*x**2 + c*x**4)) + (-5*A*c + 4*B*b)/(4*b**2*x*sqrt(b*x**2 + c*x**4)) - (-15*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(8*b**3*x**3) + 3*c*(-5*A*c + 4*B*b)*atanh(sqrt(b)*x/sqrt(b*x**2 + c*x**4))/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_159():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(b*x**2 + c*x**4)
    F = 2*A*b*x**(sympy.S(13)/2)/13 + 2*B*c*x**(sympy.S(21)/2)/21 + x**(sympy.S(17)/2)*(2*A*c + 2*B*b)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_160():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(b*x**2 + c*x**4)
    F = 2*A*b*x**(sympy.S(11)/2)/11 + 2*B*c*x**(sympy.S(19)/2)/19 + x**(sympy.S(15)/2)*(2*A*c + 2*B*b)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_161():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(b*x**2 + c*x**4)
    F = 2*A*b*x**(sympy.S(9)/2)/9 + 2*B*c*x**(sympy.S(17)/2)/17 + x**(sympy.S(13)/2)*(2*A*c + 2*B*b)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_162():
    f = sqrt(x)*(A + B*x**2)*(b*x**2 + c*x**4)
    F = 2*A*b*x**(sympy.S(7)/2)/7 + 2*B*c*x**(sympy.S(15)/2)/15 + x**(sympy.S(11)/2)*(2*A*c + 2*B*b)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_163():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/sqrt(x)
    F = 2*A*b*x**(sympy.S(5)/2)/5 + 2*B*c*x**(sympy.S(13)/2)/13 + x**(sympy.S(9)/2)*(2*A*c + 2*B*b)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_164():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    F = 2*A*b*x**(sympy.S(3)/2)/3 + 2*B*c*x**(sympy.S(11)/2)/11 + x**(sympy.S(7)/2)*(2*A*c + 2*B*b)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_165():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**(sympy.S(5)/2)
    F = 2*A*b*sqrt(x) + 2*B*c*x**(sympy.S(9)/2)/9 + x**(sympy.S(5)/2)*(2*A*c + 2*B*b)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_166():
    f = (A + B*x**2)*(b*x**2 + c*x**4)/x**(sympy.S(7)/2)
    F = -2*A*b/sqrt(x) + 2*B*c*x**(sympy.S(7)/2)/7 + x**(sympy.S(3)/2)*(2*A*c + 2*B*b)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_167():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**2
    F = 2*A*b**2*x**(sympy.S(17)/2)/17 + 2*B*c**2*x**(sympy.S(29)/2)/29 + 2*b*x**(sympy.S(21)/2)*(2*A*c + B*b)/21 + 2*c*x**(sympy.S(25)/2)*(A*c + 2*B*b)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_168():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**2
    F = 2*A*b**2*x**(sympy.S(15)/2)/15 + 2*B*c**2*x**(sympy.S(27)/2)/27 + 2*b*x**(sympy.S(19)/2)*(2*A*c + B*b)/19 + 2*c*x**(sympy.S(23)/2)*(A*c + 2*B*b)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_169():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**2
    F = 2*A*b**2*x**(sympy.S(13)/2)/13 + 2*B*c**2*x**(sympy.S(25)/2)/25 + 2*b*x**(sympy.S(17)/2)*(2*A*c + B*b)/17 + 2*c*x**(sympy.S(21)/2)*(A*c + 2*B*b)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_170():
    f = sqrt(x)*(A + B*x**2)*(b*x**2 + c*x**4)**2
    F = 2*A*b**2*x**(sympy.S(11)/2)/11 + 2*B*c**2*x**(sympy.S(23)/2)/23 + 2*b*x**(sympy.S(15)/2)*(2*A*c + B*b)/15 + 2*c*x**(sympy.S(19)/2)*(A*c + 2*B*b)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_171():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/sqrt(x)
    F = 2*A*b**2*x**(sympy.S(9)/2)/9 + 2*B*c**2*x**(sympy.S(21)/2)/21 + 2*b*x**(sympy.S(13)/2)*(2*A*c + B*b)/13 + 2*c*x**(sympy.S(17)/2)*(A*c + 2*B*b)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_172():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**(sympy.S(3)/2)
    F = 2*A*b**2*x**(sympy.S(7)/2)/7 + 2*B*c**2*x**(sympy.S(19)/2)/19 + 2*b*x**(sympy.S(11)/2)*(2*A*c + B*b)/11 + 2*c*x**(sympy.S(15)/2)*(A*c + 2*B*b)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_173():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**(sympy.S(5)/2)
    F = 2*A*b**2*x**(sympy.S(5)/2)/5 + 2*B*c**2*x**(sympy.S(17)/2)/17 + 2*b*x**(sympy.S(9)/2)*(2*A*c + B*b)/9 + 2*c*x**(sympy.S(13)/2)*(A*c + 2*B*b)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_174():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**2/x**(sympy.S(7)/2)
    F = 2*A*b**2*x**(sympy.S(3)/2)/3 + 2*B*c**2*x**(sympy.S(15)/2)/15 + 2*b*x**(sympy.S(7)/2)*(2*A*c + B*b)/7 + 2*c*x**(sympy.S(11)/2)*(A*c + 2*B*b)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_175():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**3
    F = 2*A*b**3*x**(sympy.S(21)/2)/21 + 2*B*c**3*x**(sympy.S(37)/2)/37 + 2*b**2*x**(sympy.S(25)/2)*(3*A*c + B*b)/25 + 6*b*c*x**(sympy.S(29)/2)*(A*c + B*b)/29 + 2*c**2*x**(sympy.S(33)/2)*(A*c + 3*B*b)/33
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_176():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**3
    F = 2*A*b**3*x**(sympy.S(19)/2)/19 + 2*B*c**3*x**(sympy.S(35)/2)/35 + 2*b**2*x**(sympy.S(23)/2)*(3*A*c + B*b)/23 + 2*b*c*x**(sympy.S(27)/2)*(A*c + B*b)/9 + 2*c**2*x**(sympy.S(31)/2)*(A*c + 3*B*b)/31
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_177():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**3
    F = 2*A*b**3*x**(sympy.S(17)/2)/17 + 2*B*c**3*x**(sympy.S(33)/2)/33 + 2*b**2*x**(sympy.S(21)/2)*(3*A*c + B*b)/21 + 6*b*c*x**(sympy.S(25)/2)*(A*c + B*b)/25 + 2*c**2*x**(sympy.S(29)/2)*(A*c + 3*B*b)/29
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_178():
    f = sqrt(x)*(A + B*x**2)*(b*x**2 + c*x**4)**3
    F = 2*A*b**3*x**(sympy.S(15)/2)/15 + 2*B*c**3*x**(sympy.S(31)/2)/31 + 2*b**2*x**(sympy.S(19)/2)*(3*A*c + B*b)/19 + 6*b*c*x**(sympy.S(23)/2)*(A*c + B*b)/23 + 2*c**2*x**(sympy.S(27)/2)*(A*c + 3*B*b)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_179():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/sqrt(x)
    F = 2*A*b**3*x**(sympy.S(13)/2)/13 + 2*B*c**3*x**(sympy.S(29)/2)/29 + 2*b**2*x**(sympy.S(17)/2)*(3*A*c + B*b)/17 + 2*b*c*x**(sympy.S(21)/2)*(A*c + B*b)/7 + 2*c**2*x**(sympy.S(25)/2)*(A*c + 3*B*b)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_180():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**(sympy.S(3)/2)
    F = 2*A*b**3*x**(sympy.S(11)/2)/11 + 2*B*c**3*x**(sympy.S(27)/2)/27 + 2*b**2*x**(sympy.S(15)/2)*(3*A*c + B*b)/15 + 6*b*c*x**(sympy.S(19)/2)*(A*c + B*b)/19 + 2*c**2*x**(sympy.S(23)/2)*(A*c + 3*B*b)/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_181():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**(sympy.S(5)/2)
    F = 2*A*b**3*x**(sympy.S(9)/2)/9 + 2*B*c**3*x**(sympy.S(25)/2)/25 + 2*b**2*x**(sympy.S(13)/2)*(3*A*c + B*b)/13 + 6*b*c*x**(sympy.S(17)/2)*(A*c + B*b)/17 + 2*c**2*x**(sympy.S(21)/2)*(A*c + 3*B*b)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_182():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**3/x**(sympy.S(7)/2)
    F = 2*A*b**3*x**(sympy.S(7)/2)/7 + 2*B*c**3*x**(sympy.S(23)/2)/23 + 2*b**2*x**(sympy.S(11)/2)*(3*A*c + B*b)/11 + 2*b*c*x**(sympy.S(15)/2)*(A*c + B*b)/5 + 2*c**2*x**(sympy.S(19)/2)*(A*c + 3*B*b)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_183():
    f = x**(sympy.S(13)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(11)/2)/(11*c) - sqrt(2)*b**(sympy.S(7)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(15)/4)) + sqrt(2)*b**(sympy.S(7)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(15)/4)) + sqrt(2)*b**(sympy.S(7)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(15)/4)) - sqrt(2)*b**(sympy.S(7)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(15)/4)) + 2*b*x**(sympy.S(3)/2)*(-A*c + B*b)/(3*c**3) - x**(sympy.S(7)/2)*(-2*A*c + 2*B*b)/(7*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_184():
    f = x**(sympy.S(11)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(9)/2)/(9*c) + sqrt(2)*b**(sympy.S(5)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(5)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(5)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(5)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(13)/4)) + 2*b*sqrt(x)*(-A*c + B*b)/c**3 - x**(sympy.S(5)/2)*(-2*A*c + 2*B*b)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_185():
    f = x**(sympy.S(9)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(7)/2)/(7*c) + sqrt(2)*b**(sympy.S(3)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(11)/4)) - sqrt(2)*b**(sympy.S(3)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(11)/4)) - sqrt(2)*b**(sympy.S(3)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)) + sqrt(2)*b**(sympy.S(3)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(11)/4)) + x**(sympy.S(3)/2)*(2*A*c - 2*B*b)/(3*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_186():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(5)/2)/(5*c) - sqrt(2)*b**(sympy.S(1)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*c**(sympy.S(9)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*c**(sympy.S(9)/4)) + sqrt(x)*(2*A*c - 2*B*b)/c**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_187():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(3)/2)/(3*c) - sqrt(2)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(1)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_188():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = 2*B*sqrt(x)/c + sqrt(2)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(3)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_189():
    f = sqrt(x)*(A + B*x**2)/(b*x**2 + c*x**4)
    F = -2*A/(b*sqrt(x)) + sqrt(2)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(5)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_190():
    f = (A + B*x**2)/(sqrt(x)*(b*x**2 + c*x**4))
    F = -2*A/(3*b*x**(sympy.S(3)/2)) - sqrt(2)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) - sqrt(2)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_191():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4))
    F = -2*A/(5*b*x**(sympy.S(5)/2)) - (-2*A*c + 2*B*b)/(b**2*sqrt(x)) - sqrt(2)*c**(sympy.S(1)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(9)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(9)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_192():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*(b*x**2 + c*x**4))
    F = -2*A/(7*b*x**(sympy.S(7)/2)) - (-2*A*c + 2*B*b)/(3*b**2*x**(sympy.S(3)/2)) + sqrt(2)*c**(sympy.S(3)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(11)/4)) - sqrt(2)*c**(sympy.S(3)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(11)/4)) + sqrt(2)*c**(sympy.S(3)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4)) - sqrt(2)*c**(sympy.S(3)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_193():
    f = (A + B*x**2)/(x**(sympy.S(7)/2)*(b*x**2 + c*x**4))
    F = -2*A/(9*b*x**(sympy.S(9)/2)) - (-2*A*c + 2*B*b)/(5*b**2*x**(sympy.S(5)/2)) + 2*c*(-A*c + B*b)/(b**3*sqrt(x)) + sqrt(2)*c**(sympy.S(5)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(5)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(5)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(5)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_194():
    f = (A + B*x**2)/(x**(sympy.S(9)/2)*(b*x**2 + c*x**4))
    F = -2*A/(11*b*x**(sympy.S(11)/2)) - (-2*A*c + 2*B*b)/(7*b**2*x**(sympy.S(7)/2)) + 2*c*(-A*c + B*b)/(3*b**3*x**(sympy.S(3)/2)) - sqrt(2)*c**(sympy.S(7)/4)*(-A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(15)/4)) + sqrt(2)*c**(sympy.S(7)/4)*(-A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(4*b**(sympy.S(15)/4)) - sqrt(2)*c**(sympy.S(7)/4)*(-A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(15)/4)) + sqrt(2)*c**(sympy.S(7)/4)*(-A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(2*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_195():
    f = x**(sympy.S(19)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = sqrt(2)*b**(sympy.S(5)/4)*(-9*A*c + 13*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(17)/4)) - sqrt(2)*b**(sympy.S(5)/4)*(-9*A*c + 13*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(17)/4)) + sqrt(2)*b**(sympy.S(5)/4)*(-9*A*c + 13*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(17)/4)) - sqrt(2)*b**(sympy.S(5)/4)*(-9*A*c + 13*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(17)/4)) + b*sqrt(x)*(-9*A*c + 13*B*b)/(2*c**4) - x**(sympy.S(5)/2)*(-9*A*c + 13*B*b)/(10*c**3) - x**(sympy.S(13)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + x**(sympy.S(9)/2)*(-9*A*c + 13*B*b)/(18*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_196():
    f = x**(sympy.S(17)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = sqrt(2)*b**(sympy.S(3)/4)*(-7*A*c + 11*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(15)/4)) - sqrt(2)*b**(sympy.S(3)/4)*(-7*A*c + 11*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(15)/4)) - sqrt(2)*b**(sympy.S(3)/4)*(-7*A*c + 11*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(15)/4)) + sqrt(2)*b**(sympy.S(3)/4)*(-7*A*c + 11*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(15)/4)) - x**(sympy.S(3)/2)*(-7*A*c + 11*B*b)/(6*c**3) - x**(sympy.S(11)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + x**(sympy.S(7)/2)*(-7*A*c + 11*B*b)/(14*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_197():
    f = x**(sympy.S(15)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 9*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 9*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*c**(sympy.S(13)/4)) - sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 9*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)) + sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 9*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*c**(sympy.S(13)/4)) - sqrt(x)*(-5*A*c + 9*B*b)/(2*c**3) - x**(sympy.S(9)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + x**(sympy.S(5)/2)*(-5*A*c + 9*B*b)/(10*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_198():
    f = x**(sympy.S(13)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -x**(sympy.S(7)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + x**(sympy.S(3)/2)*(-3*A*c + 7*B*b)/(6*b*c**2) - sqrt(2)*(-3*A*c + 7*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) + sqrt(2)*(-3*A*c + 7*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) + sqrt(2)*(-3*A*c + 7*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*c**(sympy.S(11)/4)) - sqrt(2)*(-3*A*c + 7*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(1)/4)*c**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_199():
    f = x**(sympy.S(11)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -x**(sympy.S(5)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + sqrt(x)*(-A*c + 5*B*b)/(2*b*c**2) + sqrt(2)*(-A*c + 5*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) - sqrt(2)*(-A*c + 5*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) + sqrt(2)*(-A*c + 5*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(3)/4)*c**(sympy.S(9)/4)) - sqrt(2)*(-A*c + 5*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(3)/4)*c**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_200():
    f = x**(sympy.S(9)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -x**(sympy.S(3)/2)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + sqrt(2)*(A*c + 3*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(A*c + 3*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(A*c + 3*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(5)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(A*c + 3*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(5)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_201():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -sqrt(x)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) - sqrt(2)*(3*A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(3*A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(3*A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(3*A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(7)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_202():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -(-A*c + B*b)/(2*b*c*sqrt(x)*(b + c*x**2)) + (-5*A*c + B*b)/(2*b**2*c*sqrt(x)) + sqrt(2)*(-5*A*c + B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-5*A*c + B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-5*A*c + B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-5*A*c + B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(9)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_203():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -(-A*c + B*b)/(2*b*c*x**(sympy.S(3)/2)*(b + c*x**2)) + (-7*A*c + 3*B*b)/(6*b**2*c*x**(sympy.S(3)/2)) - sqrt(2)*(-7*A*c + 3*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-7*A*c + 3*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) - sqrt(2)*(-7*A*c + 3*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-7*A*c + 3*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(11)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_204():
    f = sqrt(x)*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -(-A*c + B*b)/(2*b*c*x**(sympy.S(5)/2)*(b + c*x**2)) + (-9*A*c + 5*B*b)/(10*b**2*c*x**(sympy.S(5)/2)) - (-9*A*c + 5*B*b)/(2*b**3*sqrt(x)) - sqrt(2)*c**(sympy.S(1)/4)*(-9*A*c + 5*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-9*A*c + 5*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(13)/4)) + sqrt(2)*c**(sympy.S(1)/4)*(-9*A*c + 5*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4)) - sqrt(2)*c**(sympy.S(1)/4)*(-9*A*c + 5*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_205():
    f = (A + B*x**2)/(sqrt(x)*(b*x**2 + c*x**4)**2)
    F = -(-A*c + B*b)/(2*b*c*x**(sympy.S(7)/2)*(b + c*x**2)) + (-11*A*c + 7*B*b)/(14*b**2*c*x**(sympy.S(7)/2)) - (-11*A*c + 7*B*b)/(6*b**3*x**(sympy.S(3)/2)) + sqrt(2)*c**(sympy.S(3)/4)*(-11*A*c + 7*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(15)/4)) - sqrt(2)*c**(sympy.S(3)/4)*(-11*A*c + 7*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(15)/4)) + sqrt(2)*c**(sympy.S(3)/4)*(-11*A*c + 7*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(15)/4)) - sqrt(2)*c**(sympy.S(3)/4)*(-11*A*c + 7*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_206():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**2)
    F = -(-A*c + B*b)/(2*b*c*x**(sympy.S(9)/2)*(b + c*x**2)) + (-13*A*c + 9*B*b)/(18*b**2*c*x**(sympy.S(9)/2)) - (-13*A*c + 9*B*b)/(10*b**3*x**(sympy.S(5)/2)) + c*(-13*A*c + 9*B*b)/(2*b**4*sqrt(x)) + sqrt(2)*c**(sympy.S(5)/4)*(-13*A*c + 9*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(17)/4)) - sqrt(2)*c**(sympy.S(5)/4)*(-13*A*c + 9*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(16*b**(sympy.S(17)/4)) - sqrt(2)*c**(sympy.S(5)/4)*(-13*A*c + 9*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(17)/4)) + sqrt(2)*c**(sympy.S(5)/4)*(-13*A*c + 9*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(8*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_207():
    f = x**(sympy.S(23)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -9*sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 13*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*c**(sympy.S(17)/4)) + 9*sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 13*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*c**(sympy.S(17)/4)) - 9*sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 13*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)) + 9*sqrt(2)*b**(sympy.S(1)/4)*(-5*A*c + 13*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*c**(sympy.S(17)/4)) + sqrt(x)*(45*A*c - 117*B*b)/(16*c**4) - x**(sympy.S(13)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) - x**(sympy.S(9)/2)*(-5*A*c + 13*B*b)/(16*b*c**2*(b + c*x**2)) + x**(sympy.S(5)/2)*(-45*A*c + 117*B*b)/(80*b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_208():
    f = x**(sympy.S(21)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(11)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) - x**(sympy.S(7)/2)*(-3*A*c + 11*B*b)/(16*b*c**2*(b + c*x**2)) + x**(sympy.S(3)/2)*(-21*A*c + 77*B*b)/(48*b*c**3) - sqrt(2)*(-21*A*c + 77*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(1)/4)*c**(sympy.S(15)/4)) + sqrt(2)*(-21*A*c + 77*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(1)/4)*c**(sympy.S(15)/4)) + sqrt(2)*(-21*A*c + 77*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(1)/4)*c**(sympy.S(15)/4)) - sqrt(2)*(-21*A*c + 77*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(1)/4)*c**(sympy.S(15)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_209():
    f = x**(sympy.S(19)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(9)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) - x**(sympy.S(5)/2)*(-A*c + 9*B*b)/(16*b*c**2*(b + c*x**2)) + sqrt(x)*(-5*A*c + 45*B*b)/(16*b*c**3) + sqrt(2)*(-5*A*c + 45*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(3)/4)*c**(sympy.S(13)/4)) - sqrt(2)*(-5*A*c + 45*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(3)/4)*c**(sympy.S(13)/4)) + sqrt(2)*(-5*A*c + 45*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(3)/4)*c**(sympy.S(13)/4)) - sqrt(2)*(-5*A*c + 45*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(3)/4)*c**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_210():
    f = x**(sympy.S(17)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(7)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) - x**(sympy.S(3)/2)*(A*c + 7*B*b)/(16*b*c**2*(b + c*x**2)) + sqrt(2)*(3*A*c + 21*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(5)/4)*c**(sympy.S(11)/4)) - sqrt(2)*(3*A*c + 21*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(5)/4)*c**(sympy.S(11)/4)) - sqrt(2)*(3*A*c + 21*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(5)/4)*c**(sympy.S(11)/4)) + sqrt(2)*(3*A*c + 21*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(5)/4)*c**(sympy.S(11)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_211():
    f = x**(sympy.S(15)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(5)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) - sqrt(x)*(3*A*c + 5*B*b)/(16*b*c**2*(b + c*x**2)) - sqrt(2)*(3*A*c + 5*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(7)/4)*c**(sympy.S(9)/4)) + sqrt(2)*(3*A*c + 5*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(7)/4)*c**(sympy.S(9)/4)) - sqrt(2)*(3*A*c + 5*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(7)/4)*c**(sympy.S(9)/4)) + sqrt(2)*(3*A*c + 5*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(7)/4)*c**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_212():
    f = x**(sympy.S(13)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -x**(sympy.S(3)/2)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) + x**(sympy.S(3)/2)*(5*A*c + 3*B*b)/(16*b**2*c*(b + c*x**2)) + sqrt(2)*(5*A*c + 3*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(9)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(5*A*c + 3*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(9)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(5*A*c + 3*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(9)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(5*A*c + 3*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(9)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_213():
    f = x**(sympy.S(11)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -sqrt(x)*(-A*c + B*b)/(4*b*c*(b + c*x**2)**2) + sqrt(x)*(7*A*c + B*b)/(16*b**2*c*(b + c*x**2)) - sqrt(2)*(21*A*c + 3*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(11)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(21*A*c + 3*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(11)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(21*A*c + 3*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(11)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(21*A*c + 3*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(11)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_214():
    f = x**(sympy.S(9)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -(-A*c + B*b)/(4*b*c*sqrt(x)*(b + c*x**2)**2) - (-9*A*c + B*b)/(16*b**2*c*sqrt(x)*(b + c*x**2)) + (-45*A*c + 5*B*b)/(16*b**3*c*sqrt(x)) + sqrt(2)*(-45*A*c + 5*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(13)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-45*A*c + 5*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(13)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(-45*A*c + 5*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(13)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-45*A*c + 5*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(13)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_215():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -(-A*c + B*b)/(4*b*c*x**(sympy.S(3)/2)*(b + c*x**2)**2) - (-11*A*c + 3*B*b)/(16*b**2*c*x**(sympy.S(3)/2)*(b + c*x**2)) + (-77*A*c + 21*B*b)/(48*b**3*c*x**(sympy.S(3)/2)) - sqrt(2)*(-77*A*c + 21*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(15)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-77*A*c + 21*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(15)/4)*c**(sympy.S(1)/4)) - sqrt(2)*(-77*A*c + 21*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(15)/4)*c**(sympy.S(1)/4)) + sqrt(2)*(-77*A*c + 21*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(15)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_216():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -(-A*c + B*b)/(4*b*c*x**(sympy.S(5)/2)*(b + c*x**2)**2) - (-13*A*c + 5*B*b)/(16*b**2*c*x**(sympy.S(5)/2)*(b + c*x**2)) + (-117*A*c + 45*B*b)/(80*b**3*c*x**(sympy.S(5)/2)) - (-117*A*c + 45*B*b)/(16*b**4*sqrt(x)) - 9*sqrt(2)*c**(sympy.S(1)/4)*(-13*A*c + 5*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(17)/4)) + 9*sqrt(2)*c**(sympy.S(1)/4)*(-13*A*c + 5*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(17)/4)) + 9*sqrt(2)*c**(sympy.S(1)/4)*(-13*A*c + 5*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(17)/4)) - 9*sqrt(2)*c**(sympy.S(1)/4)*(-13*A*c + 5*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_217():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -(-A*c + B*b)/(4*b*c*x**(sympy.S(7)/2)*(b + c*x**2)**2) - (-15*A*c + 7*B*b)/(16*b**2*c*x**(sympy.S(7)/2)*(b + c*x**2)) + (-165*A*c + 77*B*b)/(112*b**3*c*x**(sympy.S(7)/2)) - (-165*A*c + 77*B*b)/(48*b**4*x**(sympy.S(3)/2)) + 11*sqrt(2)*c**(sympy.S(3)/4)*(-15*A*c + 7*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(19)/4)) - 11*sqrt(2)*c**(sympy.S(3)/4)*(-15*A*c + 7*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(19)/4)) + 11*sqrt(2)*c**(sympy.S(3)/4)*(-15*A*c + 7*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(19)/4)) - 11*sqrt(2)*c**(sympy.S(3)/4)*(-15*A*c + 7*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(19)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_218():
    f = sqrt(x)*(A + B*x**2)/(b*x**2 + c*x**4)**3
    F = -(-A*c + B*b)/(4*b*c*x**(sympy.S(9)/2)*(b + c*x**2)**2) - (-17*A*c + 9*B*b)/(16*b**2*c*x**(sympy.S(9)/2)*(b + c*x**2)) + (-221*A*c + 117*B*b)/(144*b**3*c*x**(sympy.S(9)/2)) - (-221*A*c + 117*B*b)/(80*b**4*x**(sympy.S(5)/2)) + 13*c*(-17*A*c + 9*B*b)/(16*b**5*sqrt(x)) + 13*sqrt(2)*c**(sympy.S(5)/4)*(-17*A*c + 9*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(21)/4)) - 13*sqrt(2)*c**(sympy.S(5)/4)*(-17*A*c + 9*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(21)/4)) - 13*sqrt(2)*c**(sympy.S(5)/4)*(-17*A*c + 9*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(21)/4)) + 13*sqrt(2)*c**(sympy.S(5)/4)*(-17*A*c + 9*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(21)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_219():
    f = (A + B*x**2)/(sqrt(x)*(b*x**2 + c*x**4)**3)
    F = -(-A*c + B*b)/(4*b*c*x**(sympy.S(11)/2)*(b + c*x**2)**2) - (-19*A*c + 11*B*b)/(16*b**2*c*x**(sympy.S(11)/2)*(b + c*x**2)) + (-285*A*c + 165*B*b)/(176*b**3*c*x**(sympy.S(11)/2)) - (-285*A*c + 165*B*b)/(112*b**4*x**(sympy.S(7)/2)) + 5*c*(-19*A*c + 11*B*b)/(16*b**5*x**(sympy.S(3)/2)) - 15*sqrt(2)*c**(sympy.S(7)/4)*(-19*A*c + 11*B*b)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(23)/4)) + 15*sqrt(2)*c**(sympy.S(7)/4)*(-19*A*c + 11*B*b)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(x) + sqrt(b) + sqrt(c)*x)/(128*b**(sympy.S(23)/4)) - 15*sqrt(2)*c**(sympy.S(7)/4)*(-19*A*c + 11*B*b)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(23)/4)) + 15*sqrt(2)*c**(sympy.S(7)/4)*(-19*A*c + 11*B*b)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4))/(64*b**(sympy.S(23)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_220():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(15*c) - 2*b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) + 4*b**2*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(231*c**3*sqrt(x)) - 4*b*x**(sympy.S(3)/2)*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(385*c**2) - x**(sympy.S(7)/2)*(-10*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(55*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_221():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = 2*B*sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(13*c) - 4*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 7*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 2*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 4*b**2*x**(sympy.S(3)/2)*(b + c*x**2)*(-13*A*c + 7*B*b)/(195*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 4*b*sqrt(x)*(-13*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(585*c**2) - x**(sympy.S(5)/2)*(-26*A*c + 14*B*b)*sqrt(b*x**2 + c*x**4)/(117*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_222():
    f = sqrt(x)*(A + B*x**2)*sqrt(b*x**2 + c*x**4)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*c*sqrt(x)) + 2*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - 4*b*(-11*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(231*c**2*sqrt(x)) - x**(sympy.S(3)/2)*(-22*A*c + 10*B*b)*sqrt(b*x**2 + c*x**4)/(77*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_223():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/sqrt(x)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*c*x**(sympy.S(3)/2)) + 4*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 2*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 4*b*x**(sympy.S(3)/2)*(b + c*x**2)*(-3*A*c + B*b)/(15*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - sqrt(x)*(-6*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(15*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_224():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(3)/2)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*c*x**(sympy.S(5)/2)) - 2*b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4)) - (-14*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(21*c*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_225():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(5)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(b*x**(sympy.S(7)/2)) - 4*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(5*A*c + B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 2*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(5*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(3)/2)*(b + c*x**2)*(20*A*c + 4*B*b)/(5*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + sqrt(x)*(10*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_226():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(7)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(3*b*x**(sympy.S(9)/2)) + (2*A*c + 2*B*b)*sqrt(b*x**2 + c*x**4)/(3*b*sqrt(x)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(2*A*c + 2*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_227():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(9)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(5*b*x**(sympy.S(11)/2)) + 4*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)*(A*c + 5*B*b)/(5*b*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - (2*A*c + 10*B*b)*sqrt(b*x**2 + c*x**4)/(5*b*x**(sympy.S(3)/2)) - 4*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + 5*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 2*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_228():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(11)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(7*b*x**(sympy.S(13)/2)) - (-2*A*c + 14*B*b)*sqrt(b*x**2 + c*x**4)/(21*b*x**(sympy.S(5)/2)) + 2*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_229():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(13)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*b*x**(sympy.S(15)/2)) - (-2*A*c + 6*B*b)*sqrt(b*x**2 + c*x**4)/(15*b*x**(sympy.S(7)/2)) + 4*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)*(-A*c + 3*B*b)/(15*b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 4*c*(-A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(15*b**2*x**(sympy.S(3)/2)) - 4*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 3*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) + 2*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_230():
    f = (A + B*x**2)*sqrt(b*x**2 + c*x**4)/x**(sympy.S(15)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(11*b*x**(sympy.S(17)/2)) - (-10*A*c + 22*B*b)*sqrt(b*x**2 + c*x**4)/(77*b*x**(sympy.S(9)/2)) - 4*c*(-5*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(231*b**2*x**(sympy.S(5)/2)) - 2*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 11*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_231():
    f = x**(sympy.S(7)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*B*x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(25*c) - 88*b**(sympy.S(21)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(16575*c**(sympy.S(19)/4)*sqrt(b*x**2 + c*x**4)) + 44*b**(sympy.S(21)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(16575*c**(sympy.S(19)/4)*sqrt(b*x**2 + c*x**4)) + 88*b**5*x**(sympy.S(3)/2)*(b + c*x**2)*(-5*A*c + 3*B*b)/(16575*c**(sympy.S(9)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 88*b**4*sqrt(x)*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(49725*c**4) + 88*b**3*x**(sympy.S(5)/2)*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(69615*c**3) - 8*b**2*x**(sympy.S(9)/2)*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(7735*c**2) - 4*b*x**(sympy.S(13)/2)*(-5*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(595*c) - x**(sympy.S(9)/2)*(-10*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(105*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_232():
    f = x**(sympy.S(5)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*B*x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(23*c) + 12*b**(sympy.S(19)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-23*A*c + 13*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(33649*c**(sympy.S(17)/4)*sqrt(b*x**2 + c*x**4)) - 24*b**4*(-23*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(33649*c**4*sqrt(x)) + 72*b**3*x**(sympy.S(3)/2)*(-23*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(168245*c**3) - 8*b**2*x**(sympy.S(7)/2)*(-23*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(24035*c**2) - 4*b*x**(sympy.S(11)/2)*(-23*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(2185*c) - x**(sympy.S(7)/2)*(-46*A*c + 26*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(437*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_233():
    f = x**(sympy.S(3)/2)*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*B*sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(21*c) + 8*b**(sympy.S(17)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-21*A*c + 11*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3315*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) - 4*b**(sympy.S(17)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-21*A*c + 11*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3315*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) - 8*b**4*x**(sympy.S(3)/2)*(b + c*x**2)*(-21*A*c + 11*B*b)/(3315*c**(sympy.S(7)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 8*b**3*sqrt(x)*(-21*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(9945*c**3) - 8*b**2*x**(sympy.S(5)/2)*(-21*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(13923*c**2) - 4*b*x**(sympy.S(9)/2)*(-21*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(1547*c) - x**(sympy.S(5)/2)*(-42*A*c + 22*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(357*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_234():
    f = sqrt(x)*(A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(19*c*sqrt(x)) - 4*b**(sympy.S(15)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-19*A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(4389*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) + 8*b**3*(-19*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(4389*c**3*sqrt(x)) - 8*b**2*x**(sympy.S(3)/2)*(-19*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(7315*c**2) - 4*b*x**(sympy.S(7)/2)*(-19*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(1045*c) - x**(sympy.S(3)/2)*(-38*A*c + 18*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(285*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_235():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/sqrt(x)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(17*c*x**(sympy.S(3)/2)) - 8*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-17*A*c + 7*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 4*b**(sympy.S(13)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-17*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 8*b**3*x**(sympy.S(3)/2)*(b + c*x**2)*(-17*A*c + 7*B*b)/(1105*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 8*b**2*sqrt(x)*(-17*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(3315*c**2) - 4*b*x**(sympy.S(5)/2)*(-17*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(663*c) - sqrt(x)*(-34*A*c + 14*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(221*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_236():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(15*c*x**(sympy.S(5)/2)) + 4*b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - 8*b**2*(-3*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(231*c**2*sqrt(x)) - 4*b*x**(sympy.S(3)/2)*(-3*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(77*c) - (-6*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(33*c*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_237():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(13*c*x**(sympy.S(7)/2)) + 8*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 3*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 4*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - 8*b**2*x**(sympy.S(3)/2)*(b + c*x**2)*(-13*A*c + 3*B*b)/(195*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 4*b*sqrt(x)*(-13*A*c + 3*B*b)*sqrt(b*x**2 + c*x**4)/(195*c) - (-26*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(117*c*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_238():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(7)/2)
    F = 2*B*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(11*c*x**(sympy.S(9)/2)) - 4*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4)) - 4*b*(-11*A*c + B*b)*sqrt(b*x**2 + c*x**4)/(77*c*sqrt(x)) - (-22*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(77*c*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_239():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(9)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(b*x**(sympy.S(11)/2)) - 8*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(9*A*c + B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 4*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(9*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 8*b*x**(sympy.S(3)/2)*(b + c*x**2)*(9*A*c + B*b)/(15*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + sqrt(x)*(12*A*c/5 + 4*B*b/15)*sqrt(b*x**2 + c*x**4) + (18*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(9*b*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_240():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(11)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(3*b*x**(sympy.S(13)/2)) + 4*b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(7*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4)) + (28*A*c + 12*B*b)*sqrt(b*x**2 + c*x**4)/(21*sqrt(x)) + (14*A*c + 6*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(21*b*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_241():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(13)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(5*b*x**(sympy.S(15)/2)) - 24*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 12*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(b*x**2 + c*x**4)) + 24*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)*(A*c + B*b)/((5*sqrt(b) + 5*sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 12*c*sqrt(x)*(A*c + B*b)*sqrt(b*x**2 + c*x**4)/(5*b) - (2*A*c + 2*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(b*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_242():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(15)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(7*b*x**(sympy.S(17)/2)) + 4*c*(3*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(21*b*sqrt(x)) - (6*A*c + 14*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(21*b*x**(sympy.S(9)/2)) + 4*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(3*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_243():
    f = (A + B*x**2)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/x**(sympy.S(17)/2)
    F = -2*A*(b*x**2 + c*x**4)**(sympy.S(5)/2)/(9*b*x**(sympy.S(19)/2)) + 8*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)*(A*c + 9*B*b)/(15*b*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - 4*c*(A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(15*b*x**(sympy.S(3)/2)) - (2*A*c + 18*B*b)*(b*x**2 + c*x**4)**(sympy.S(3)/2)/(45*b*x**(sympy.S(11)/2)) - 8*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + 9*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) + 4*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_244():
    f = x**(sympy.S(13)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(11)/2)*sqrt(b*x**2 + c*x**4)/(15*c) + b**(sympy.S(11)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-15*A*c + 13*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*c**(sympy.S(17)/4)*sqrt(b*x**2 + c*x**4)) - 2*b**2*(-15*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(77*c**4*sqrt(x)) + 6*b*x**(sympy.S(3)/2)*(-15*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(385*c**3) - x**(sympy.S(7)/2)*(-30*A*c + 26*B*b)*sqrt(b*x**2 + c*x**4)/(165*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_245():
    f = x**(sympy.S(11)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(9)/2)*sqrt(b*x**2 + c*x**4)/(13*c) + 14*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 11*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) - 7*b**(sympy.S(9)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-13*A*c + 11*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(195*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) - 14*b**2*x**(sympy.S(3)/2)*(b + c*x**2)*(-13*A*c + 11*B*b)/(195*c**(sympy.S(7)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 14*b*sqrt(x)*(-13*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(585*c**3) - x**(sympy.S(5)/2)*(-26*A*c + 22*B*b)*sqrt(b*x**2 + c*x**4)/(117*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_246():
    f = x**(sympy.S(9)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)/(11*c) - 5*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) + 10*b*(-11*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(231*c**3*sqrt(x)) - x**(sympy.S(3)/2)*(-22*A*c + 18*B*b)*sqrt(b*x**2 + c*x**4)/(77*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_247():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)/(9*c) - 2*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 7*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 2*b*x**(sympy.S(3)/2)*(b + c*x**2)*(-9*A*c + 7*B*b)/(15*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - sqrt(x)*(-18*A*c + 14*B*b)*sqrt(b*x**2 + c*x**4)/(45*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_248():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)/(7*c) + b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4)) - (-14*A*c + 10*B*b)*sqrt(b*x**2 + c*x**4)/(21*c**2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_249():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*sqrt(x)*sqrt(b*x**2 + c*x**4)/(5*c) + 2*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(3)/2)*(b + c*x**2)*(-10*A*c + 6*B*b)/(5*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_250():
    f = sqrt(x)*(A + B*x**2)/sqrt(b*x**2 + c*x**4)
    F = 2*B*sqrt(b*x**2 + c*x**4)/(3*c*sqrt(x)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_251():
    f = (A + B*x**2)/(sqrt(x)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(b*x**(sympy.S(3)/2)) + x**(sympy.S(3)/2)*(b + c*x**2)*(2*A*c + 2*B*b)/(b*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(2*A*c + 2*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_252():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(3*b*x**(sympy.S(5)/2)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(5)/4)*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_253():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(5*b*x**(sympy.S(7)/2)) + 2*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)*(-3*A*c + 5*B*b)/(5*b**2*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - (-6*A*c + 10*B*b)*sqrt(b*x**2 + c*x**4)/(5*b**2*x**(sympy.S(3)/2)) - 2*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + 5*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_254():
    f = (A + B*x**2)/(x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(7*b*x**(sympy.S(9)/2)) - (-10*A*c + 14*B*b)*sqrt(b*x**2 + c*x**4)/(21*b**2*x**(sympy.S(5)/2)) - c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_255():
    f = (A + B*x**2)/(x**(sympy.S(9)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(9*b*x**(sympy.S(11)/2)) - (-14*A*c + 18*B*b)*sqrt(b*x**2 + c*x**4)/(45*b**2*x**(sympy.S(7)/2)) - 2*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)*(-7*A*c + 9*B*b)/(15*b**3*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 2*c*(-7*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(15*b**3*x**(sympy.S(3)/2)) + 2*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 9*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_256():
    f = (A + B*x**2)/(x**(sympy.S(11)/2)*sqrt(b*x**2 + c*x**4))
    F = -2*A*sqrt(b*x**2 + c*x**4)/(11*b*x**(sympy.S(13)/2)) - (-18*A*c + 22*B*b)*sqrt(b*x**2 + c*x**4)/(77*b**2*x**(sympy.S(9)/2)) + 10*c*(-9*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(231*b**3*x**(sympy.S(5)/2)) + 5*c**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 11*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_257():
    f = x**(sympy.S(17)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -15*b**(sympy.S(7)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + 13*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(154*c**(sympy.S(17)/4)*sqrt(b*x**2 + c*x**4)) + 15*b*(-11*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(77*c**4*sqrt(x)) - x**(sympy.S(3)/2)*(-99*A*c + 117*B*b)*sqrt(b*x**2 + c*x**4)/(77*c**3) - x**(sympy.S(15)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(7)/2)*(-11*A*c + 13*B*b)*sqrt(b*x**2 + c*x**4)/(11*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_258():
    f = x**(sympy.S(15)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -7*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 11*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) + 7*b**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 11*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(30*c**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) + 7*b*x**(sympy.S(3)/2)*(b + c*x**2)*(-9*A*c + 11*B*b)/(15*c**(sympy.S(7)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - sqrt(x)*(-63*A*c + 77*B*b)*sqrt(b*x**2 + c*x**4)/(45*c**3) - x**(sympy.S(13)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(5)/2)*(-9*A*c + 11*B*b)*sqrt(b*x**2 + c*x**4)/(9*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_259():
    f = x**(sympy.S(13)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 5*b**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(42*c**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4)) - (-35*A*c + 45*B*b)*sqrt(b*x**2 + c*x**4)/(21*c**3*sqrt(x)) - x**(sympy.S(11)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(3)/2)*(-7*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(7*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_260():
    f = x**(sympy.S(11)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 3*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 7*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - 3*b**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*c**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(3)/2)*(b + c*x**2)*(-15*A*c + 21*B*b)/(5*c**(sympy.S(5)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(9)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + sqrt(x)*(-5*A*c + 7*B*b)*sqrt(b*x**2 + c*x**4)/(5*b*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_261():
    f = x**(sympy.S(9)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**(sympy.S(7)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + (-3*A*c + 5*B*b)*sqrt(b*x**2 + c*x**4)/(3*b*c**2*sqrt(x)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(1)/4)*c**(sympy.S(9)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_262():
    f = x**(sympy.S(7)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**(sympy.S(5)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(3)/2)*(b + c*x**2)*(-A*c + 3*B*b)/(b*c**(sympy.S(3)/2)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 3*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*c**(sympy.S(7)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_263():
    f = x**(sympy.S(5)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -x**(sympy.S(3)/2)*(-A*c + B*b)/(b*c*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(5)/4)*c**(sympy.S(5)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_264():
    f = x**(sympy.S(3)/2)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -2*A*sqrt(x)/(b*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(5)/2)*(-3*A*c + B*b)/(b**2*sqrt(b*x**2 + c*x**4)) - x**(sympy.S(3)/2)*(b + c*x**2)*(-3*A*c + B*b)/(b**2*sqrt(c)*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(7)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4)) - x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-3*A*c + B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*c**(sympy.S(3)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_265():
    f = sqrt(x)*(A + B*x**2)/(b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -2*A/(3*b*sqrt(x)*sqrt(b*x**2 + c*x**4)) + x**(sympy.S(3)/2)*(-5*A*c + 3*B*b)/(3*b**2*sqrt(b*x**2 + c*x**4)) + x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-5*A*c + 3*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(9)/4)*c**(sympy.S(1)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_266():
    f = (A + B*x**2)/(sqrt(x)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -2*A/(5*b*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)) + sqrt(x)*(-7*A*c + 5*B*b)/(5*b**2*sqrt(b*x**2 + c*x**4)) + 3*sqrt(c)*x**(sympy.S(3)/2)*(b + c*x**2)*(-7*A*c + 5*B*b)/(5*b**3*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) - (-21*A*c + 15*B*b)*sqrt(b*x**2 + c*x**4)/(5*b**3*x**(sympy.S(3)/2)) - 3*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 5*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4)) + 3*c**(sympy.S(1)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-7*A*c + 5*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*b**(sympy.S(11)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_267():
    f = (A + B*x**2)/(x**(sympy.S(3)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -2*A/(7*b*x**(sympy.S(5)/2)*sqrt(b*x**2 + c*x**4)) + (-9*A*c + 7*B*b)/(7*b**2*sqrt(x)*sqrt(b*x**2 + c*x**4)) - (-45*A*c + 35*B*b)*sqrt(b*x**2 + c*x**4)/(21*b**3*x**(sympy.S(5)/2)) - 5*c**(sympy.S(3)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-9*A*c + 7*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(42*b**(sympy.S(13)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_268():
    f = (A + B*x**2)/(x**(sympy.S(5)/2)*(b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -2*A/(9*b*x**(sympy.S(7)/2)*sqrt(b*x**2 + c*x**4)) + (-11*A*c + 9*B*b)/(9*b**2*x**(sympy.S(3)/2)*sqrt(b*x**2 + c*x**4)) - (-77*A*c + 63*B*b)*sqrt(b*x**2 + c*x**4)/(45*b**3*x**(sympy.S(7)/2)) - 7*c**(sympy.S(3)/2)*x**(sympy.S(3)/2)*(b + c*x**2)*(-11*A*c + 9*B*b)/(15*b**4*(sqrt(b) + sqrt(c)*x)*sqrt(b*x**2 + c*x**4)) + 7*c*(-11*A*c + 9*B*b)*sqrt(b*x**2 + c*x**4)/(15*b**4*x**(sympy.S(3)/2)) + 7*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + 9*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4)) - 7*c**(sympy.S(5)/4)*x*sqrt((b + c*x**2)/(sqrt(b) + sqrt(c)*x)**2)*(sqrt(b) + sqrt(c)*x)*(-11*A*c + 9*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*sqrt(x)/b**(sympy.S(1)/4)), sympy.S.Half)/(30*b**(sympy.S(15)/4)*sqrt(b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_269():
    f = x**m*(A + B*x**2)*(b*x**2 + c*x**4)**3
    F = A*b**3*x**(m + 7)/(m + 7) + B*c**3*x**(m + 15)/(m + 15) + b**2*x**(m + 9)*(3*A*c + B*b)/(m + 9) + 3*b*c*x**(m + 11)*(A*c + B*b)/(m + 11) + c**2*x**(m + 13)*(A*c + 3*B*b)/(m + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_270():
    f = x**m*(A + B*x**2)*(b*x**2 + c*x**4)**2
    F = A*b**2*x**(m + 5)/(m + 5) + B*c**2*x**(m + 11)/(m + 11) + b*x**(m + 7)*(2*A*c + B*b)/(m + 7) + c*x**(m + 9)*(A*c + 2*B*b)/(m + 9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_271():
    f = x**m*(A + B*x**2)*(b*x**2 + c*x**4)
    F = A*b*x**(m + 3)/(m + 3) + B*c*x**(m + 7)/(m + 7) + x**(m + 5)*(A*c + B*b)/(m + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_272():
    f = x**m*(A + B*x**2)/(b*x**2 + c*x**4)
    F = -B*x**(m - 1)/(c*(1 - m)) + x**(m - 1)*(-A*c + B*b)*hyper((1, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), -c*x**2/b)/(b*c*(1 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_273():
    f = x**m*(A + B*x**2)/(b*x**2 + c*x**4)**2
    F = -x**(m - 3)*(-A*c + B*b)/(2*b*c*(b + c*x**2)) + x**(m - 3)*(-A*c*(5 - m) + B*b*(3 - m))*hyper((1, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), -c*x**2/b)/(2*b**2*c*(3 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_274():
    f = x**m*(A + B*x**2)*(b*x**2 + c*x**4)**p
    F = B*x**(m - 1)*(b*x**2 + c*x**4)**(p + 1)/(c*(m + 4*p + 3)) - x**(m + 1)*(b*x**2 + c*x**4)**p*(-A*c*(m + 4*p + 3) + B*b*(m + 2*p + 1))*hyper((-p, m/2 + p + sympy.S.Half), (m/2 + p + sympy.S(3)/2,), -c*x**2/b)/(c*(1 + c*x**2/b)**p*(m + 2*p + 1)*(m + 4*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_275():
    f = x**(-j*p + n - 1)*(c + d*x**n)*(a*x**j + b*x**(j + n))**p
    F = d*x**(-j*(p + 1) + n)*(a*x**j + b*x**(j + n))**(p + 1)/(b*n*(p + 2)) - (a*d - b*c*(p + 2))*(a*x**j + b*x**(j + n))**(p + 1)/(b**2*n*x**(j*(p + 1))*(p + 1)*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_276():
    f = (e*x)**m*(c + d*x**n)**q*(a*x**j + b*x**(j + n))**p
    F = x*(e*x)**m*(c + d*x**n)**q*(a*x**j + b*x**(j + n))**p*appellf1((j*p + m + 1)/n, -p, -q, (j*p + m + n + 1)/n, -b*x**n/a, -d*x**n/c)/((1 + b*x**n/a)**p*(1 + d*x**n/c)**q*(j*p + m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_277():
    f = (e*x)**(sympy.S(7)/4)*(c + d*x**n)**q*(a*x**j + b*x**(j + n))**(sympy.S(5)/3)
    F = 12*a*e*x**(j + 2)*(e*x)**(sympy.S(3)/4)*(c + d*x**n)**q*(a*x**j + b*x**(j + n))**(sympy.S(2)/3)*appellf1((20*j + 33)/(12*n), sympy.S(-5)/3, -q, 1 + (5*j/3 + sympy.S(11)/4)/n, -b*x**n/a, -d*x**n/c)/((1 + b*x**n/a)**(sympy.S(2)/3)*(1 + d*x**n/c)**q*(20*j + 33))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_278():
    f = (3*x**4 + 4)/(2*x**5 + 5*x)
    F = 4*log(x)/5 + 7*log(2*x**4 + 5)/40
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_279():
    f = (x**6 + 1)/(-x**7 + x)
    F = log(x) - log(1 - x**6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_280():
    f = (5*x**10 + 8)/(-x**11 + 2*x)
    F = 4*log(x) - 9*log(2 - x**10)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_281():
    f = (2*x - 3)/(x**3 - x**2)
    F = log(x) - log(1 - x) - 3/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_282():
    f = (a*x**m + b*x**n)/(c*x**m + d*x**n)
    F = a*x/c + x*(-a*d + b*c)*hyper((1, 1/(m - n)), (1 + 1/(m - n),), -c*x**(m - n)/d)/(c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_283():
    f = x**m*(a + b*x**n)**p*(a*x**q*(m + q + 1) + b*x**(n + q)*(m + n*(p + 1) + q + 1))
    F = x**(m + q + 1)*(a + b*x**n)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_284():
    f = x**m*(a + b/x)**n/(c + d*x)
    F = x**m*(a + b/x)**n*appellf1(-m, 1, -n, 1 - m, -c/(d*x), -b/(a*x))/(d*m*(1 + b/(a*x))**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_285():
    f = x**2*(a + b/x)**n/(c + d*x)
    F = -c**3*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(d**3*(n + 1)*(a*c - b*d)) + x**2*(a + b/x)**(n + 1)/(2*a*d) - x*(a + b/x)**(n + 1)*(2*a*c + b*d*(1 - n))/(2*a**2*d**2) + (a + b/x)**(n + 1)*(2*a**2*c**2 - 2*a*b*c*d*n - b**2*d**2*n*(1 - n))*hyper((1, n + 1), (n + 2,), 1 + b/(a*x))/(2*a**3*d**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_286():
    f = x*(a + b/x)**n/(c + d*x)
    F = c**2*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(d**2*(n + 1)*(a*c - b*d)) + x*(a + b/x)**(n + 1)/(a*d) - (a + b/x)**(n + 1)*(a*c - b*d*n)*hyper((1, n + 1), (n + 2,), 1 + b/(a*x))/(a**2*d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_287():
    f = (a + b/x)**n/(c + d*x)
    F = -c*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(d*(n + 1)*(a*c - b*d)) + (a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b/(a*x))/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_288():
    f = (a + b/x)**n/(x*(c + d*x))
    F = (a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/((n + 1)*(a*c - b*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_289():
    f = (a + b/x)**n/(x**2*(c + d*x))
    F = -d*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c*(n + 1)*(a*c - b*d)) - (a + b/x)**(n + 1)/(b*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_290():
    f = (a + b/x)**n/(x**3*(c + d*x))
    F = d**2*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c**2*(n + 1)*(a*c - b*d)) - (a + b/x)**(n + 2)/(b**2*c*(n + 2)) + (a + b/x)**(n + 1)*(a*c + b*d)/(b**2*c**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_291():
    f = (a + b/x)**n/(x**5*(c + d*x))
    F = d**4*(a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c**4*(n + 1)*(a*c - b*d)) - (a + b/x)**(n + 4)/(b**4*c*(n + 4)) + (a + b/x)**(n + 3)*(3*a*c + b*d)/(b**4*c**2*(n + 3)) - (a + b/x)**(n + 2)*(3*a**2*c**2 + 2*a*b*c*d + b**2*d**2)/(b**4*c**3*(n + 2)) + (a + b/x)**(n + 1)*(a*c + b*d)*(a**2*c**2 + b**2*d**2)/(b**4*c**4*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_292():
    f = x**m*(a + b/x)**n/(c + d*x)**2
    F = -x**(m - 1)*(a + b/x)**n*appellf1(1 - m, 2, -n, 2 - m, -c/(d*x), -b/(a*x))/(d**2*(1 - m)*(1 + b/(a*x))**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_293():
    f = x**2*(a + b/x)**n/(c + d*x)**2
    F = c**2*(a + b/x)**(n + 1)*(2*a*c - b*d*(2 - n))*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(d**3*(n + 1)*(a*c - b*d)**2) + c*(a + b/x)**(n + 1)*(2*a*c - b*d)/(a*d**2*(a*c - b*d)*(c/x + d)) + x*(a + b/x)**(n + 1)/(a*d*(c/x + d)) - (a + b/x)**(n + 1)*(2*a*c - b*d*n)*hyper((1, n + 1), (n + 2,), 1 + b/(a*x))/(a**2*d**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_294():
    f = x*(a + b/x)**n/(c + d*x)**2
    F = -c*(a + b/x)**(n + 1)/(d*(a*c - b*d)*(c/x + d)) - c*(a + b/x)**(n + 1)*(a*c - b*d*(1 - n))*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(d**2*(n + 1)*(a*c - b*d)**2) + (a + b/x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b/(a*x))/(a*d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_295():
    f = (a + b/x)**n/(c + d*x)**2
    F = -b*(a + b/x)**(n + 1)*hyper((2, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/((n + 1)*(a*c - b*d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_296():
    f = (a + b/x)**n/(x*(c + d*x)**2)
    F = -d*(a + b/x)**(n + 1)/(c*(a*c - b*d)*(c/x + d)) + (a + b/x)**(n + 1)*(a*c - b*d*(n + 1))*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c*(n + 1)*(a*c - b*d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_297():
    f = (a + b/x)**n/(x**2*(c + d*x)**2)
    F = d**2*(a + b/x)**(n + 1)/(c**2*(a*c - b*d)*(c/x + d)) - d*(a + b/x)**(n + 1)*(2*a*c - b*d*(n + 2))*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c**2*(n + 1)*(a*c - b*d)**2) - (a + b/x)**(n + 1)/(b*c**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_3_e_x_pow_m_a_x_pow_j_plus_b_x_pow_k_pow_p_c_plus_d_x_pow_n_pow_q_298():
    f = (a + b/x)**n/(x**3*(c + d*x)**2)
    F = d**2*(a + b/x)**(n + 1)*(3*a*c - b*d*(n + 3))*hyper((1, n + 1), (n + 2,), c*(a + b/x)/(a*c - b*d))/(c**3*(n + 1)*(a*c - b*d)**2) - (a + b/x)**(n + 1)/(b*c*x**2*(n + 2)*(c/x + d)) - (a + b/x)**(n + 1)*(-c*(a*c - b*d)*(a*c + b*d*(n + 3))/x + d*(-a*c*(a*c + b*d*(3*n + 5)) + b*d*(n + 2)*(a*c + b*d*(n + 3))))/(b**2*c**3*(n + 1)*(n + 2)*(a*c - b*d)*(c/x + d))
    assert integrate(f, x) == F

