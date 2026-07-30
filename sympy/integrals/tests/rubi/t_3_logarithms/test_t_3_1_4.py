"""Generated from MathematicaSyntaxTestSuite.

Source: 3 Logarithms/3.1.4 (f x)^m (d+e x^r)^q (a+b log(c x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, n, p, q, r = symbols('a b c d e f g m n p q r')

def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_1():
    f = x**3*(a + b*log(c*x**n))*(d + e*x)
    F = -b*d*n*x**4/16 - b*e*n*x**5/25 + (a + b*log(c*x**n))*(d*x**4/4 + e*x**5/5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_2():
    f = x**2*(a + b*log(c*x**n))*(d + e*x)
    F = -b*d*n*x**3/9 - b*e*n*x**4/16 + (a + b*log(c*x**n))*(d*x**3/3 + e*x**4/4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_3():
    f = x*(a + b*log(c*x**n))*(d + e*x)
    F = -b*d*n*x**2/4 - b*e*n*x**3/9 + (a + b*log(c*x**n))*(d*x**2/2 + e*x**3/3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_4():
    f = (a + b*log(c*x**n))*(d + e*x)
    F = -b*d*n*x - b*e*n*x**2/4 + d*x*(a + b*log(c*x**n)) + e*x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_5():
    f = (a + b*log(c*x**n))*(d + e*x)/x
    F = a*e*x - b*e*n*x + b*e*x*log(c*x**n) + d*(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_6():
    f = (a + b*log(c*x**n))*(d + e*x)/x**3
    F = -b*d*n/(4*x**2) - b*e*n/x + b*e**2*n*log(x)/(2*d) - (a + b*log(c*x**n))*(d + e*x)**2/(2*d*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_7():
    f = (a + b*log(c*x**n))*(d + e*x)/x**4
    F = -b*d*n/(9*x**3) - b*e*n/(4*x**2) - d*(a + b*log(c*x**n))/(3*x**3) - e*(a + b*log(c*x**n))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_8():
    f = x**3*(a + b*log(c*x**n))*(d + e*x)**2
    F = -b*d**2*n*x**4/16 - 2*b*d*e*n*x**5/25 - b*e**2*n*x**6/36 + (a + b*log(c*x**n))*(d**2*x**4/4 + 2*d*e*x**5/5 + e**2*x**6/6)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_9():
    f = x**2*(a + b*log(c*x**n))*(d + e*x)**2
    F = -b*d**2*n*x**3/9 - b*d*e*n*x**4/8 - b*e**2*n*x**5/25 + (a + b*log(c*x**n))*(d**2*x**3/3 + d*e*x**4/2 + e**2*x**5/5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_10():
    f = x*(a + b*log(c*x**n))*(d + e*x)**2
    F = -b*d**2*n*x**2/4 - 2*b*d*e*n*x**3/9 - b*e**2*n*x**4/16 + (a + b*log(c*x**n))*(d**2*x**2/2 + 2*d*e*x**3/3 + e**2*x**4/4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_11():
    f = (a + b*log(c*x**n))*(d + e*x)**2
    F = -b*d**3*n*log(x)/(3*e) - b*d**2*n*x - b*d*e*n*x**2/2 - b*e**2*n*x**3/9 + (a + b*log(c*x**n))*(d + e*x)**3/(3*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_12():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x
    F = -b*d**2*n*log(x)**2/2 - b*n*(4*d + e*x)**2/4 + d**2*(a + b*log(c*x**n))*log(x) + 2*d*e*x*(a + b*log(c*x**n)) + e**2*x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_13():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x**2
    F = -b*d**2*n/x - b*d*e*n*log(x)**2 - b*e**2*n*x - d**2*(a + b*log(c*x**n))/x + 2*d*e*(a + b*log(c*x**n))*log(x) + e**2*x*(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_14():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x**3
    F = -b*e**2*n*log(x)**2/2 - b*n*(d + 4*e*x)**2/(4*x**2) - d**2*(a + b*log(c*x**n))/(2*x**2) - 2*d*e*(a + b*log(c*x**n))/x + e**2*(a + b*log(c*x**n))*log(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_15():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x**4
    F = -b*d**2*n/(9*x**3) - b*d*e*n/(2*x**2) - b*e**2*n/x + b*e**3*n*log(x)/(3*d) - (a + b*log(c*x**n))*(d + e*x)**3/(3*d*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_16():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x**5
    F = -b*d**2*n/(16*x**4) - 2*b*d*e*n/(9*x**3) - b*e**2*n/(4*x**2) - d**2*(a + b*log(c*x**n))/(4*x**4) - 2*d*e*(a + b*log(c*x**n))/(3*x**3) - e**2*(a + b*log(c*x**n))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_17():
    f = (a + b*log(c*x**n))*(d + e*x)**2/x**6
    F = -b*d**2*n/(25*x**5) - b*d*e*n/(8*x**4) - b*e**2*n/(9*x**3) - d**2*(a + b*log(c*x**n))/(5*x**5) - d*e*(a + b*log(c*x**n))/(2*x**4) - e**2*(a + b*log(c*x**n))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_18():
    f = x**3*(a + b*log(c*x**n))*(d + e*x)**3
    F = -b*d**3*n*x**4/16 - 3*b*d**2*e*n*x**5/25 - b*d*e**2*n*x**6/12 - b*e**3*n*x**7/49 + (a + b*log(c*x**n))*(d**3*x**4/4 + 3*d**2*e*x**5/5 + d*e**2*x**6/2 + e**3*x**7/7)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_19():
    f = x**2*(a + b*log(c*x**n))*(d + e*x)**3
    F = -b*d**3*n*x**3/9 - 3*b*d**2*e*n*x**4/16 - 3*b*d*e**2*n*x**5/25 - b*e**3*n*x**6/36 + (a + b*log(c*x**n))*(d**3*x**3/3 + 3*d**2*e*x**4/4 + 3*d*e**2*x**5/5 + e**3*x**6/6)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_20():
    f = x*(a + b*log(c*x**n))*(d + e*x)**3
    F = b*d**5*n*log(x)/(20*e**2) + b*d**4*n*x/(5*e) + 3*b*d**3*n*x**2/20 + b*d**2*e*n*x**3/15 + b*d*e**2*n*x**4/80 - b*n*(d + e*x)**5/(25*e**2) - (a + b*log(c*x**n))*(d*(d + e*x)**4/(4*e**2) - (d + e*x)**5/(5*e**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_21():
    f = (a + b*log(c*x**n))*(d + e*x)**3
    F = -b*d**4*n*log(x)/(4*e) - b*d**3*n*x - 3*b*d**2*e*n*x**2/4 - b*d*e**2*n*x**3/3 - b*e**3*n*x**4/16 + (a + b*log(c*x**n))*(d + e*x)**4/(4*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_22():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x
    F = -b*d**3*n*log(x)**2/2 - 3*b*d**2*e*n*x - 3*b*d*e**2*n*x**2/4 - b*e**3*n*x**3/9 + d**3*(a + b*log(c*x**n))*log(x) + 3*d**2*e*x*(a + b*log(c*x**n)) + 3*d*e**2*x**2*(a + b*log(c*x**n))/2 + e**3*x**3*(a + b*log(c*x**n))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_23():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**2
    F = -b*d**3*n/x - 3*b*d**2*e*n*log(x)**2/2 - 3*b*d*e**2*n*x - b*e**3*n*x**2/4 - d**3*(a + b*log(c*x**n))/x + 3*d**2*e*(a + b*log(c*x**n))*log(x) + 3*d*e**2*x*(a + b*log(c*x**n)) + e**3*x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_24():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**3
    F = -b*d**3*n/(4*x**2) - 3*b*d**2*e*n/x - 3*b*d*e**2*n*log(x)**2/2 - b*e**3*n*x - d**3*(a + b*log(c*x**n))/(2*x**2) - 3*d**2*e*(a + b*log(c*x**n))/x + 3*d*e**2*(a + b*log(c*x**n))*log(x) + e**3*x*(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_25():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**4
    F = -b*d**3*n/(9*x**3) - 3*b*d**2*e*n/(4*x**2) - 3*b*d*e**2*n/x - b*e**3*n*log(x)**2/2 - d**3*(a + b*log(c*x**n))/(3*x**3) - 3*d**2*e*(a + b*log(c*x**n))/(2*x**2) - 3*d*e**2*(a + b*log(c*x**n))/x + e**3*(a + b*log(c*x**n))*log(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_26():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**5
    F = -b*d**3*n/(16*x**4) - b*d**2*e*n/(3*x**3) - 3*b*d*e**2*n/(4*x**2) - b*e**3*n/x + b*e**4*n*log(x)/(4*d) - (a + b*log(c*x**n))*(d + e*x)**4/(4*d*x**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_27():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**6
    F = b*d**2*e*n/(80*x**4) + b*d*e**2*n/(15*x**3) + 3*b*e**3*n/(20*x**2) + b*e**4*n/(5*d*x) - b*e**5*n*log(x)/(20*d**2) - b*n*(d + e*x)**5/(25*d**2*x**5) - (a + b*log(c*x**n))*(d + e*x)**4/(5*d*x**5) + e*(a + b*log(c*x**n))*(d + e*x)**4/(20*d**2*x**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_28():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**7
    F = -b*d**3*n/(36*x**6) - 3*b*d**2*e*n/(25*x**5) - 3*b*d*e**2*n/(16*x**4) - b*e**3*n/(9*x**3) - d**3*(a + b*log(c*x**n))/(6*x**6) - 3*d**2*e*(a + b*log(c*x**n))/(5*x**5) - 3*d*e**2*(a + b*log(c*x**n))/(4*x**4) - e**3*(a + b*log(c*x**n))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_29():
    f = (a + b*log(c*x**n))*(d + e*x)**3/x**8
    F = -b*d**3*n/(49*x**7) - b*d**2*e*n/(12*x**6) - 3*b*d*e**2*n/(25*x**5) - b*e**3*n/(16*x**4) - d**3*(a + b*log(c*x**n))/(7*x**7) - d**2*e*(a + b*log(c*x**n))/(2*x**6) - 3*d*e**2*(a + b*log(c*x**n))/(5*x**5) - e**3*(a + b*log(c*x**n))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_30():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)
    F = ((Symbol('a') * (Symbol('d'))**(Integer(2)) * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(3))) * ((Integer(9) * Symbol('e')))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_31():
    f = x**2*(a + b*log(c*x**n))/(d + e*x)
    F = (Integer(-1) * ((Symbol('a') * Symbol('d') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_32():
    f = x*(a + b*log(c*x**n))/(d + e*x)
    F = ((Symbol('a') * x) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * (Symbol('e'))**(Integer(-1)))) + ((Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_33():
    f = (a + b*log(c*x**n))/(d + e*x)
    F = (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_34():
    f = (a + b*log(c*x**n))/(x*(d + e*x))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_35():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('d') * x))**(Integer(-1)))) + ((Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_36():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_37():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(9) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(3) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(3)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(3)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_38():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)**2
    F = ((Integer(3) * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Integer(3) * Symbol('a')) + (Symbol('b') * Symbol('n'))) * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('d') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Integer(3) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Integer(3) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_39():
    f = x*(a + b*log(c*x**n))/(d + e*x)**2
    F = (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * Symbol('n')) + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_40():
    f = (a + b*log(c*x**n))/(d + e*x)**2
    F = -b*n*log(d + e*x)/(d*e) + x*(a + b*log(c*x**n))/(d*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_41():
    f = (a + b*log(c*x**n))/(x*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('e') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_42():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_43():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(3)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_44():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)**3
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n') * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((((Integer(6) * Symbol('a')) + (Integer(5) * Symbol('b') * Symbol('n'))) * x) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * ((Integer(3) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Integer(6) * Symbol('a')) + (Integer(5) * Symbol('b') * Symbol('n')) + (Integer(6) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_45():
    f = x**2*(a + b*log(c*x**n))/(d + e*x)**3
    F = (Integer(-1) * (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Integer(2) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(2) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((((Integer(2) * Symbol('a')) + (Integer(3) * Symbol('b') * Symbol('n')) + (Integer(2) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_46():
    f = x*(a + b*log(c*x**n))/(d + e*x)**3
    F = -b*n/(2*e**2*(d + e*x)) - b*n*log(d + e*x)/(2*d*e**2) + x**2*(a + b*log(c*x**n))/(2*d*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_47():
    f = (a + b*log(c*x**n))/(d + e*x)**3
    F = b*n/(2*d*e*(d + e*x)) + b*n*log(x)/(2*d**2*e) - b*n*log(d + e*x)/(2*d**2*e) - (a + b*log(c*x**n))/(2*e*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_48():
    f = (a + b*log(c*x**n))/(x*(d + e*x)**3)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.log(x)) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_49():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x)**3)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.log(x)) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('e'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(3) * Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_50():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x)**3)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(4)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log(x)) * ((Integer(2) * (Symbol('d'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(4)) * x))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(3)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + ((Integer(7) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(5))))**(Integer(-1))) + ((Integer(6) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_51():
    f = x**5*(a + b*log(c*x**n))/(d + e*x)**4
    F = ((Integer(10) * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Integer(60) * Symbol('a')) + (Integer(47) * Symbol('b') * Symbol('n'))) * x) * ((Integer(6) * (Symbol('e'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(10) * Symbol('b') * Symbol('d') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * ((Integer(5) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(5) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * ((Integer(20) * Symbol('a')) + (Integer(9) * Symbol('b') * Symbol('n')) + (Integer(20) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Integer(60) * Symbol('a')) + (Integer(47) * Symbol('b') * Symbol('n')) + (Integer(60) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(12) * (Symbol('e'))**(Integer(4))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Integer(60) * Symbol('a')) + (Integer(47) * Symbol('b') * Symbol('n')) + (Integer(60) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('e'))**(Integer(6))))**(Integer(-1))) + ((Integer(10) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(6)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_52():
    f = x**4*(a + b*log(c*x**n))/(d + e*x)**4
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * x) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))) + ((((Integer(12) * Symbol('a')) + (Integer(13) * Symbol('b') * Symbol('n'))) * x) * ((Integer(3) * (Symbol('e'))**(Integer(4))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * ((Integer(4) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(4) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * ((Integer(12) * Symbol('a')) + (Integer(7) * Symbol('b') * Symbol('n')) + (Integer(12) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Integer(12) * Symbol('a')) + (Integer(13) * Symbol('b') * Symbol('n')) + (Integer(12) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * (Symbol('e'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_53():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)**4
    F = (Integer(-1) * (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * ((Integer(3) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Integer(6) * Symbol('a')) + (Integer(5) * Symbol('b') * Symbol('n')) + (Integer(6) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((((Integer(6) * Symbol('a')) + (Integer(11) * Symbol('b') * Symbol('n')) + (Integer(6) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('e'))**(Integer(4))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_54():
    f = x**2*(a + b*log(c*x**n))/(d + e*x)**4
    F = b*d*n/(6*e**3*(d + e*x)**2) - 2*b*n/(3*e**3*(d + e*x)) - b*n*log(d + e*x)/(3*d*e**3) + x**3*(a + b*log(c*x**n))/(3*d*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_55():
    f = x*(a + b*log(c*x**n))/(d + e*x)**4
    F = -b*n/(6*e**2*(d + e*x)**2) + b*n/(6*d*e**2*(d + e*x)) + b*n*log(x)/(6*d**2*e**2) - b*n*log(d + e*x)/(6*d**2*e**2) + d*(a + b*log(c*x**n))/(3*e**2*(d + e*x)**3) - (a + b*log(c*x**n))/(2*e**2*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_56():
    f = (a + b*log(c*x**n))/(d + e*x)**4
    F = b*n/(6*d*e*(d + e*x)**2) + b*n/(3*d**2*e*(d + e*x)) + b*n*log(x)/(3*d**3*e) - b*n*log(d + e*x)/(3*d**3*e) - (a + b*log(c*x**n))/(3*e*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_57():
    f = (a + b*log(c*x**n))/(x*(d + e*x)**4)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('n')) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('n') * sympy.log(x)) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(3) * Symbol('d') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(11) * Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_58():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x)**4)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * (((Symbol('d'))**(Integer(4)) * x))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(3) * (Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.log(x)) * ((Integer(3) * (Symbol('d'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((Symbol('d'))**(Integer(4)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(4) * Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(13) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_59():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x)**4)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(4)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(5)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(6) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(6) * (Symbol('d'))**(Integer(5)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log(x)) * ((Integer(6) * (Symbol('d'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(4)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(5)) * x))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(6)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(10) * (Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1)))) + ((Integer(47) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(6))))**(Integer(-1))) + ((Integer(10) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_60():
    f = x**8*(a + b*log(c*x**n))/(d + e*x)**7
    F = ((Integer(28) * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(8)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Integer(280) * Symbol('a')) + (Integer(341) * Symbol('b') * Symbol('n'))) * x) * ((Integer(10) * (Symbol('e'))**(Integer(8))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Symbol('e'))**(Integer(7)))**(Integer(-1)))) + (Integer(-1) * ((Integer(28) * Symbol('b') * Symbol('d') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(8)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(8)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(7)) * ((Integer(8) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(8) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(30) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(6)) * ((Integer(56) * Symbol('a')) + (Integer(15) * Symbol('b') * Symbol('n')) + (Integer(56) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(120) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * ((Integer(168) * Symbol('a')) + (Integer(73) * Symbol('b') * Symbol('n')) + (Integer(168) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(180) * (Symbol('e'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Integer(280) * Symbol('a')) + (Integer(341) * Symbol('b') * Symbol('n')) + (Integer(280) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(20) * (Symbol('e'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(4)) * ((Integer(840) * Symbol('a')) + (Integer(533) * Symbol('b') * Symbol('n')) + (Integer(840) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(360) * (Symbol('e'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * ((Integer(840) * Symbol('a')) + (Integer(743) * Symbol('b') * Symbol('n')) + (Integer(840) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(90) * (Symbol('e'))**(Integer(6)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Integer(280) * Symbol('a')) + (Integer(341) * Symbol('b') * Symbol('n')) + (Integer(280) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(10) * (Symbol('e'))**(Integer(9))))**(Integer(-1))) + ((Integer(28) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(9)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_61():
    f = x**7*(a + b*log(c*x**n))/(d + e*x)**7
    F = (Integer(-1) * ((Integer(7) * Symbol('b') * Symbol('n') * x) * ((Symbol('e'))**(Integer(7)))**(Integer(-1)))) + ((((Integer(140) * Symbol('a')) + (Integer(223) * Symbol('b') * Symbol('n'))) * x) * ((Integer(20) * (Symbol('e'))**(Integer(7))))**(Integer(-1))) + ((Integer(7) * Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(7)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(7)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(6)) * ((Integer(7) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(7) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(30) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * ((Integer(42) * Symbol('a')) + (Integer(13) * Symbol('b') * Symbol('n')) + (Integer(42) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(120) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * ((Integer(140) * Symbol('a')) + (Integer(153) * Symbol('b') * Symbol('n')) + (Integer(140) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(40) * (Symbol('e'))**(Integer(6)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * ((Integer(210) * Symbol('a')) + (Integer(107) * Symbol('b') * Symbol('n')) + (Integer(210) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(360) * (Symbol('e'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * ((Integer(420) * Symbol('a')) + (Integer(319) * Symbol('b') * Symbol('n')) + (Integer(420) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(360) * (Symbol('e'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Integer(140) * Symbol('a')) + (Integer(223) * Symbol('b') * Symbol('n')) + (Integer(140) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(20) * (Symbol('e'))**(Integer(8))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(8)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_62():
    f = x**6*(a + b*log(c*x**n))/(d + e*x)**7
    F = (Integer(-1) * (((x)**(Integer(6)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * ((Integer(6) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(6) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(30) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * ((Integer(20) * Symbol('a')) + (Integer(19) * Symbol('b') * Symbol('n')) + (Integer(20) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(40) * (Symbol('e'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Integer(20) * Symbol('a')) + (Integer(29) * Symbol('b') * Symbol('n')) + (Integer(20) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(20) * (Symbol('e'))**(Integer(6)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * ((Integer(30) * Symbol('a')) + (Integer(11) * Symbol('b') * Symbol('n')) + (Integer(30) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(120) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * ((Integer(60) * Symbol('a')) + (Integer(37) * Symbol('b') * Symbol('n')) + (Integer(60) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(180) * (Symbol('e'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + ((((Integer(20) * Symbol('a')) + (Integer(49) * Symbol('b') * Symbol('n')) + (Integer(20) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(20) * (Symbol('e'))**(Integer(7))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(7)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_63():
    f = x**5*(a + b*log(c*x**n))/(d + e*x)**7
    F = -b*d**4*n/(30*e**6*(d + e*x)**5) + 5*b*d**3*n/(24*e**6*(d + e*x)**4) - 5*b*d**2*n/(9*e**6*(d + e*x)**3) + 5*b*d*n/(6*e**6*(d + e*x)**2) - 5*b*n/(6*e**6*(d + e*x)) - b*n*log(d + e*x)/(6*d*e**6) + x**6*(a + b*log(c*x**n))/(6*d*(d + e*x)**6)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_64():
    f = x**4*(a + b*log(c*x**n))/(d + e*x)**7
    F = b*d**2*n/(120*e**5*(d + e*x)**4) - 2*b*d*n/(45*e**5*(d + e*x)**3) + b*n/(10*e**5*(d + e*x)**2) - 2*b*n/(15*d*e**5*(d + e*x)) - b*n*x**5/(30*d**2*(d + e*x)**5) - b*n*log(d + e*x)/(30*d**2*e**5) + x**5*(a + b*log(c*x**n))/(6*d*(d + e*x)**6) + x**5*(a + b*log(c*x**n))/(30*d**2*(d + e*x)**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_65():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)**7
    F = -b*d**2*n/(30*e**4*(d + e*x)**5) + 13*b*d*n/(120*e**4*(d + e*x)**4) - 19*b*n/(180*e**4*(d + e*x)**3) + b*n/(120*d*e**4*(d + e*x)**2) + b*n/(60*d**2*e**4*(d + e*x)) + b*n*log(x)/(60*d**3*e**4) - b*n*log(d + e*x)/(60*d**3*e**4) + d**3*(a + b*log(c*x**n))/(6*e**4*(d + e*x)**6) - 3*d**2*(a + b*log(c*x**n))/(5*e**4*(d + e*x)**5) + 3*d*(a + b*log(c*x**n))/(4*e**4*(d + e*x)**4) - (a + b*log(c*x**n))/(3*e**4*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_66():
    f = x**2*(a + b*log(c*x**n))/(d + e*x)**7
    F = b*d*n/(30*e**3*(d + e*x)**5) - 7*b*n/(120*e**3*(d + e*x)**4) + b*n/(180*d*e**3*(d + e*x)**3) + b*n/(120*d**2*e**3*(d + e*x)**2) + b*n/(60*d**3*e**3*(d + e*x)) + b*n*log(x)/(60*d**4*e**3) - b*n*log(d + e*x)/(60*d**4*e**3) - d**2*(a + b*log(c*x**n))/(6*e**3*(d + e*x)**6) + 2*d*(a + b*log(c*x**n))/(5*e**3*(d + e*x)**5) - (a + b*log(c*x**n))/(4*e**3*(d + e*x)**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_67():
    f = x*(a + b*log(c*x**n))/(d + e*x)**7
    F = -b*n/(30*e**2*(d + e*x)**5) + b*n/(120*d*e**2*(d + e*x)**4) + b*n/(90*d**2*e**2*(d + e*x)**3) + b*n/(60*d**3*e**2*(d + e*x)**2) + b*n/(30*d**4*e**2*(d + e*x)) + b*n*log(x)/(30*d**5*e**2) - b*n*log(d + e*x)/(30*d**5*e**2) + d*(a + b*log(c*x**n))/(6*e**2*(d + e*x)**6) - (a + b*log(c*x**n))/(5*e**2*(d + e*x)**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_68():
    f = (a + b*log(c*x**n))/(d + e*x)**7
    F = b*n/(30*d*e*(d + e*x)**5) + b*n/(24*d**2*e*(d + e*x)**4) + b*n/(18*d**3*e*(d + e*x)**3) + b*n/(12*d**4*e*(d + e*x)**2) + b*n/(6*d**5*e*(d + e*x)) + b*n*log(x)/(6*d**6*e) - b*n*log(d + e*x)/(6*d**6*e) - (a + b*log(c*x**n))/(6*e*(d + e*x)**6)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_69():
    f = (a + b*log(c*x**n))/(x*(d + e*x)**7)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(30) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('b') * Symbol('n')) * ((Integer(120) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(37) * Symbol('b') * Symbol('n')) * ((Integer(180) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(19) * Symbol('b') * Symbol('n')) * ((Integer(40) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(29) * Symbol('b') * Symbol('n')) * ((Integer(20) * (Symbol('d'))**(Integer(6)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(29) * Symbol('b') * Symbol('n') * sympy.log(x)) * ((Integer(20) * (Symbol('d'))**(Integer(7))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(6) * Symbol('d') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(3) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(7)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(7)))**(Integer(-1)))) + ((Integer(49) * Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(20) * (Symbol('d'))**(Integer(7))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(7)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_70():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x)**7)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * (((Symbol('d'))**(Integer(7)) * x))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(30) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1))) + ((Integer(17) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(120) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1))) + ((Integer(79) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(180) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(53) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(40) * (Symbol('d'))**(Integer(6)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(103) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(20) * (Symbol('d'))**(Integer(7)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(103) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.log(x)) * ((Integer(20) * (Symbol('d'))**(Integer(8))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((Symbol('d'))**(Integer(7)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(5) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(6)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * (Symbol('e'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(8)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(7) * Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(8)))**(Integer(-1))) + (Integer(-1) * ((Integer(223) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(20) * (Symbol('d'))**(Integer(8))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(8)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_71():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x)**7)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(7)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(7) * Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(8)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(30) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(23) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(120) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(34) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(45) * (Symbol('d'))**(Integer(6)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(14) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(5) * (Symbol('d'))**(Integer(7)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(131) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')) * ((Integer(10) * (Symbol('d'))**(Integer(8)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(131) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log(x)) * ((Integer(10) * (Symbol('d'))**(Integer(9))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(7)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(7) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(8)) * x))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(6))))**(Integer(-1))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(5) * (Symbol('d'))**(Integer(4)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(5)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1))) + ((Integer(10) * (Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(6)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(15) * (Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(7)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(21) * (Symbol('e'))**(Integer(3)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(9)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(28) * (Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**(Integer(9)))**(Integer(-1)))) + ((Integer(341) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * x)))) * ((Integer(10) * (Symbol('d'))**(Integer(9))))**(Integer(-1))) + ((Integer(28) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(9)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_72():
    f = log(c*x)/(-c*x + 1)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (Symbol('c'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_73():
    f = log(x/c)/(c - x)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (x * (Symbol('c'))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_74():
    f = x**2*(a + b*log(c*x**n))**2*(d + e*x)
    F = 2*b**2*d*n**2*x**3/27 + b**2*e*n**2*x**4/32 - 2*b*d*n*x**3*(a + b*log(c*x**n))/9 - b*e*n*x**4*(a + b*log(c*x**n))/8 + d*x**3*(a + b*log(c*x**n))**2/3 + e*x**4*(a + b*log(c*x**n))**2/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_75():
    f = x*(a + b*log(c*x**n))**2*(d + e*x)
    F = b**2*d*n**2*x**2/4 + 2*b**2*e*n**2*x**3/27 - b*d*n*x**2*(a + b*log(c*x**n))/2 - 2*b*e*n*x**3*(a + b*log(c*x**n))/9 + d*x**2*(a + b*log(c*x**n))**2/2 + e*x**3*(a + b*log(c*x**n))**2/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_76():
    f = (a + b*log(c*x**n))**2*(d + e*x)
    F = -2*a*b*d*n*x + 2*b**2*d*n**2*x - 2*b**2*d*n*x*log(c*x**n) + b**2*e*n**2*x**2/4 - b*e*n*x**2*(a + b*log(c*x**n))/2 + d*x*(a + b*log(c*x**n))**2 + e*x**2*(a + b*log(c*x**n))**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_77():
    f = (a + b*log(c*x**n))**2*(d + e*x)/x
    F = -2*a*b*e*n*x + 2*b**2*e*n**2*x - 2*b**2*e*n*x*log(c*x**n) + e*x*(a + b*log(c*x**n))**2 + d*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_78():
    f = (a + b*log(c*x**n))**2*(d + e*x)/x**2
    F = -2*b**2*d*n**2/x - 2*b*d*n*(a + b*log(c*x**n))/x - d*(a + b*log(c*x**n))**2/x + e*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_79():
    f = (a + b*log(c*x**n))**2*(d + e*x)/x**3
    F = -b**2*d*n**2/(4*x**2) - 2*b**2*e*n**2/x - b*d*n*(a + b*log(c*x**n))/(2*x**2) - 2*b*e*n*(a + b*log(c*x**n))/x - d*(a + b*log(c*x**n))**2/(2*x**2) - e*(a + b*log(c*x**n))**2/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_80():
    f = (a + b*log(c*x**n))**2*(d + e*x)/x**4
    F = -2*b**2*d*n**2/(27*x**3) - b**2*e*n**2/(4*x**2) - 2*b*d*n*(a + b*log(c*x**n))/(9*x**3) - b*e*n*(a + b*log(c*x**n))/(2*x**2) - d*(a + b*log(c*x**n))**2/(3*x**3) - e*(a + b*log(c*x**n))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_81():
    f = (a + b*log(c*x**n))**2*(d + e*x)/x**5
    F = -b**2*d*n**2/(32*x**4) - 2*b**2*e*n**2/(27*x**3) - b*d*n*(a + b*log(c*x**n))/(8*x**4) - 2*b*e*n*(a + b*log(c*x**n))/(9*x**3) - d*(a + b*log(c*x**n))**2/(4*x**4) - e*(a + b*log(c*x**n))**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_82():
    f = x**2*(a + b*log(c*x**n))**2*(d + e*x)**2
    F = 2*b**2*d**2*n**2*x**3/27 + b**2*d*e*n**2*x**4/16 + 2*b**2*e**2*n**2*x**5/125 - 2*b*d**2*n*x**3*(a + b*log(c*x**n))/9 - b*d*e*n*x**4*(a + b*log(c*x**n))/4 - 2*b*e**2*n*x**5*(a + b*log(c*x**n))/25 + d**2*x**3*(a + b*log(c*x**n))**2/3 + d*e*x**4*(a + b*log(c*x**n))**2/2 + e**2*x**5*(a + b*log(c*x**n))**2/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_83():
    f = x*(a + b*log(c*x**n))**2*(d + e*x)**2
    F = b**2*d**2*n**2*x**2/4 + 4*b**2*d*e*n**2*x**3/27 + b**2*e**2*n**2*x**4/32 - b*d**2*n*x**2*(a + b*log(c*x**n))/2 - 4*b*d*e*n*x**3*(a + b*log(c*x**n))/9 - b*e**2*n*x**4*(a + b*log(c*x**n))/8 + d**2*x**2*(a + b*log(c*x**n))**2/2 + 2*d*e*x**3*(a + b*log(c*x**n))**2/3 + e**2*x**4*(a + b*log(c*x**n))**2/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_84():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2
    F = b**2*d**3*n**2*log(x)**2/(3*e) + 2*b**2*d**2*n**2*x + b**2*d*e*n**2*x**2/2 + 2*b**2*e**2*n**2*x**3/27 - 2*b*d**3*n*(a + b*log(c*x**n))*log(x)/(3*e) - 2*b*d**2*n*x*(a + b*log(c*x**n)) - b*d*e*n*x**2*(a + b*log(c*x**n)) - 2*b*e**2*n*x**3*(a + b*log(c*x**n))/9 + (a + b*log(c*x**n))**2*(d + e*x)**3/(3*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_85():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2/x
    F = -4*a*b*d*e*n*x + 4*b**2*d*e*n**2*x - 4*b**2*d*e*n*x*log(c*x**n) + b**2*e**2*n**2*x**2/4 - b*e**2*n*x**2*(a + b*log(c*x**n))/2 + 2*d*e*x*(a + b*log(c*x**n))**2 + e**2*x**2*(a + b*log(c*x**n))**2/2 + d**2*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_86():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2/x**2
    F = -2*a*b*e**2*n*x - 2*b**2*d**2*n**2/x + 2*b**2*e**2*n**2*x - 2*b**2*e**2*n*x*log(c*x**n) - 2*b*d**2*n*(a + b*log(c*x**n))/x - d**2*(a + b*log(c*x**n))**2/x + e**2*x*(a + b*log(c*x**n))**2 + 2*d*e*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_87():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2/x**3
    F = -b**2*d**2*n**2/(4*x**2) - 4*b**2*d*e*n**2/x - b*d**2*n*(a + b*log(c*x**n))/(2*x**2) - 4*b*d*e*n*(a + b*log(c*x**n))/x - d**2*(a + b*log(c*x**n))**2/(2*x**2) - 2*d*e*(a + b*log(c*x**n))**2/x + e**2*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_88():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2/x**4
    F = -2*b**2*d**2*n**2/(27*x**3) - b**2*d*e*n**2/(2*x**2) - 2*b**2*e**2*n**2/x - 2*b*d**2*n*(a + b*log(c*x**n))/(9*x**3) - b*d*e*n*(a + b*log(c*x**n))/x**2 - 2*b*e**2*n*(a + b*log(c*x**n))/x - d**2*(a + b*log(c*x**n))**2/(3*x**3) - d*e*(a + b*log(c*x**n))**2/x**2 - e**2*(a + b*log(c*x**n))**2/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_89():
    f = (a + b*log(c*x**n))**2*(d + e*x)**2/x**5
    F = -b**2*d**2*n**2/(32*x**4) - 4*b**2*d*e*n**2/(27*x**3) - b**2*e**2*n**2/(4*x**2) - b*d**2*n*(a + b*log(c*x**n))/(8*x**4) - 4*b*d*e*n*(a + b*log(c*x**n))/(9*x**3) - b*e**2*n*(a + b*log(c*x**n))/(2*x**2) - d**2*(a + b*log(c*x**n))**2/(4*x**4) - 2*d*e*(a + b*log(c*x**n))**2/(3*x**3) - e**2*(a + b*log(c*x**n))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_90():
    f = x**3*(a + b*log(c*x**n))**2/(d + e*x)
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(27) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * Symbol('n') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(9) * Symbol('e')))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(3) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_91():
    f = x**2*(a + b*log(c*x**n))**2/(d + e*x)
    F = ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(4) * Symbol('e')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * Symbol('n') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_92():
    f = x*(a + b*log(c*x**n))**2/(d + e*x)
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('n') * x) * (Symbol('e'))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * x) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('n') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('e'))**(Integer(-1)))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_93():
    f = (a + b*log(c*x**n))**2/(d + e*x)
    F = ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_94():
    f = (a + b*log(c*x**n))**2/(x*(d + e*x))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_95():
    f = (a + b*log(c*x**n))**2/(x**2*(d + e*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Symbol('d') * x))**(Integer(-1)))) + ((Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_96():
    f = (a + b*log(c*x**n))**2/(x**3*(d + e*x))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * ((Integer(4) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_97():
    f = (a + b*log(c*x**n))**2/(x**4*(d + e*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * ((Integer(27) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(9) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(3) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(3)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(3)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_98():
    f = x**3*(a + b*log(c*x**n))**2/(d + e*x)**2
    F = ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * x) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * Symbol('n') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + ((Integer(6) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_99():
    f = x**2*(a + b*log(c*x**n))**2/(d + e*x)**2
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('n') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('n') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_100():
    f = x*(a + b*log(c*x**n))**2/(d + e*x)**2
    F = (Integer(-1) * ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_101():
    f = (a + b*log(c*x**n))**2/(d + e*x)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_102():
    f = (a + b*log(c*x**n))**2/(x*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('e') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_103():
    f = (a + b*log(c*x**n))**2/(x**2*(d + e*x)**2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_104():
    f = (a + b*log(c*x**n))**2/(x**3*(d + e*x)**2)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (Symbol('n'))**(Integer(2))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(3)) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * (((Symbol('d'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Integer(6) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_105():
    f = x*(a + b*log(c*x**n))**2/(d + e*x)**3
    F = ((Symbol('b') * Symbol('n') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') * Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * Symbol('n')) + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('e'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_106():
    f = (a + b*log(c*x**n))**2/(d + e*x)**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('e') * x)))) * (((Symbol('d'))**(Integer(2)) * Symbol('e')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_107():
    f = x**2*(a + b*log(c*x**n))**2/(d + e*x)**4
    F = ((Symbol('b') * Symbol('n') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('d') * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(3) * Symbol('d') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * x * ((Integer(2) * Symbol('a')) + (Symbol('b') * Symbol('n')) + (Integer(2) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('d') * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * ((Integer(2) * Symbol('a')) + (Integer(3) * Symbol('b') * Symbol('n')) + (Integer(2) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * (Symbol('e'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * (Symbol('e'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_108():
    f = (a + b*log(c*x**n))**2/(d + e*x)**4
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.log(x)) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('d') * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('e') * x)))) * (((Symbol('d'))**(Integer(3)) * Symbol('e')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_109():
    f = (a + b*log(c*x**n))**3/(x*(d + e*x))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_110():
    f = (a + b*log(c*x**n))**3/(x*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('e') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * (((Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * x))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (Symbol('d') * ((Symbol('e') * x))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_111():
    f = sqrt(a + b*log(c*x**n))*(d + e*x)
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * Symbol('d') * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Integer(4))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + ((Integer(2))**(Integer(-1)) * Symbol('e') * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_112():
    f = sqrt(a + b*log(c*x**n))*(d + e*x)**2
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * Symbol('d') * Symbol('e') * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(6))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('d'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Symbol('d') * Symbol('e') * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + ((Integer(3))**(Integer(-1)) * (Symbol('e'))**(Integer(2)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_113():
    f = sqrt(a + b*log(c*x**n))*(d + e*x)**3
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3)) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2)) * Symbol('e') * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('d'))**(Integer(3)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * Symbol('e') * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Symbol('d') * (Symbol('e'))**(Integer(2)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + ((Integer(4))**(Integer(-1)) * (Symbol('e'))**(Integer(3)) * (x)**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_114():
    f = sqrt(a + b*log(c*x**n))/(d + e*x)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_115():
    f = sqrt(a + b*log(c*x**n))/(d + e*x)**2
    F = ((x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Symbol('d') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * x)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))**(Integer(-1)), x)) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_116():
    f = sqrt(a + b*log(c*x**n))/(d + e*x)**3
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('Unintegrable')(((x * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))**(Integer(-1)), x)) * ((Integer(4) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_117():
    f = (a + b*log(c*x**n))*sqrt(d + e*x)
    F = 4*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(3*e) - 4*b*d*n*sqrt(d + e*x)/(3*e) - 4*b*n*(d + e*x)**(sympy.S(3)/2)/(9*e) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/(3*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_118():
    f = (a + b*log(c*x**n))*sqrt(d + e*x)/x
    F = (Integer(-4) * Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) + (Integer(4) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + (Integer(2) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + (Integer(2) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('d')) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Integer(-1) * (Integer(4) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_119():
    f = (a + b*log(c*x**n))*sqrt(d + e*x)/x**2
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_120():
    f = (a + b*log(c*x**n))*sqrt(d + e*x)/x**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(8) * Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('d') * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_121():
    f = x**3*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)
    F = -64*b*d**(sympy.S(11)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(1155*e**4) + 64*b*d**5*n*sqrt(d + e*x)/(1155*e**4) + 64*b*d**4*n*(d + e*x)**(sympy.S(3)/2)/(3465*e**4) + 64*b*d**3*n*(d + e*x)**(sympy.S(5)/2)/(5775*e**4) - 172*b*d**2*n*(d + e*x)**(sympy.S(7)/2)/(1617*e**4) + 32*b*d*n*(d + e*x)**(sympy.S(9)/2)/(297*e**4) - 4*b*n*(d + e*x)**(sympy.S(11)/2)/(121*e**4) - 2*d**3*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**4) + 6*d**2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(7)/2)/(7*e**4) - 2*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(9)/2)/(3*e**4) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(11)/2)/(11*e**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_122():
    f = x**2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)
    F = 32*b*d**(sympy.S(9)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(315*e**3) - 32*b*d**4*n*sqrt(d + e*x)/(315*e**3) - 32*b*d**3*n*(d + e*x)**(sympy.S(3)/2)/(945*e**3) - 32*b*d**2*n*(d + e*x)**(sympy.S(5)/2)/(1575*e**3) + 44*b*d*n*(d + e*x)**(sympy.S(7)/2)/(441*e**3) - 4*b*n*(d + e*x)**(sympy.S(9)/2)/(81*e**3) + 2*d**2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**3) - 4*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(7)/2)/(7*e**3) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(9)/2)/(9*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_123():
    f = x*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)
    F = -8*b*d**(sympy.S(7)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(35*e**2) + 8*b*d**3*n*sqrt(d + e*x)/(35*e**2) + 8*b*d**2*n*(d + e*x)**(sympy.S(3)/2)/(105*e**2) + 8*b*d*n*(d + e*x)**(sympy.S(5)/2)/(175*e**2) - 4*b*n*(d + e*x)**(sympy.S(7)/2)/(49*e**2) - 2*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**2) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(7)/2)/(7*e**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_124():
    f = (a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)
    F = 4*b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(5*e) - 4*b*d**2*n*sqrt(d + e*x)/(5*e) - 4*b*d*n*(d + e*x)**(sympy.S(3)/2)/(15*e) - 4*b*n*(d + e*x)**(sympy.S(5)/2)/(25*e) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_125():
    f = (a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/x
    F = ((Integer(-1) * (Integer(16) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) + (Integer(-1) * ((Integer(4) * (Integer(9))**(Integer(-1))) * Symbol('b') * Symbol('n') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(16) * (Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + (Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + (Integer(2) * Symbol('d') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Integer(-1) * (Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_126():
    f = (a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/x**2
    F = (Integer(-4) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (x)**(Integer(-1)))) + (Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + (Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + (Integer(3) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Symbol('d')) * Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Integer(-1) * (Integer(6) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_127():
    f = (a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/x**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_128():
    f = x**3*(a + b*log(c*x**n))/sqrt(d + e*x)
    F = -64*b*d**(sympy.S(7)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(35*e**4) + 64*b*d**3*n*sqrt(d + e*x)/(35*e**4) - 76*b*d**2*n*(d + e*x)**(sympy.S(3)/2)/(105*e**4) + 64*b*d*n*(d + e*x)**(sympy.S(5)/2)/(175*e**4) - 4*b*n*(d + e*x)**(sympy.S(7)/2)/(49*e**4) - 2*d**3*(a + b*log(c*x**n))*sqrt(d + e*x)/e**4 + 2*d**2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/e**4 - 6*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**4) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(7)/2)/(7*e**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_129():
    f = x**2*(a + b*log(c*x**n))/sqrt(d + e*x)
    F = 32*b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(15*e**3) - 32*b*d**2*n*sqrt(d + e*x)/(15*e**3) + 28*b*d*n*(d + e*x)**(sympy.S(3)/2)/(45*e**3) - 4*b*n*(d + e*x)**(sympy.S(5)/2)/(25*e**3) + 2*d**2*(a + b*log(c*x**n))*sqrt(d + e*x)/e**3 - 4*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/(3*e**3) + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_130():
    f = x*(a + b*log(c*x**n))/sqrt(d + e*x)
    F = -8*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(3*e**2) + 8*b*d*n*sqrt(d + e*x)/(3*e**2) - 4*b*n*(d + e*x)**(sympy.S(3)/2)/(9*e**2) - 2*d*(a + b*log(c*x**n))*sqrt(d + e*x)/e**2 + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_131():
    f = (a + b*log(c*x**n))/sqrt(d + e*x)
    F = 4*b*sqrt(d)*n*atanh(sqrt(d + e*x)/sqrt(d))/e - 4*b*n*sqrt(d + e*x)/e + 2*(a + b*log(c*x**n))*sqrt(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_132():
    f = (a + b*log(c*x**n))/(x*sqrt(d + e*x))
    F = ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_133():
    f = (a + b*log(c*x**n))/(x**2*sqrt(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') * x))**(Integer(-1)))) + ((Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_134():
    f = (a + b*log(c*x**n))/(x**3*sqrt(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(4) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + ((Integer(7) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_135():
    f = x**3*(a + b*log(c*x**n))/(d + e*x)**(sympy.S(3)/2)
    F = 64*b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(5*e**4) - 44*b*d**2*n*sqrt(d + e*x)/(5*e**4) + 16*b*d*n*(d + e*x)**(sympy.S(3)/2)/(15*e**4) - 4*b*n*(d + e*x)**(sympy.S(5)/2)/(25*e**4) + 2*d**3*(a + b*log(c*x**n))/(e**4*sqrt(d + e*x)) + 6*d**2*(a + b*log(c*x**n))*sqrt(d + e*x)/e**4 - 2*d*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/e**4 + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(5)/2)/(5*e**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_136():
    f = x**2*(a + b*log(c*x**n))/(d + e*x)**(sympy.S(3)/2)
    F = -32*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x)/sqrt(d))/(3*e**3) + 20*b*d*n*sqrt(d + e*x)/(3*e**3) - 4*b*n*(d + e*x)**(sympy.S(3)/2)/(9*e**3) - 2*d**2*(a + b*log(c*x**n))/(e**3*sqrt(d + e*x)) - 4*d*(a + b*log(c*x**n))*sqrt(d + e*x)/e**3 + 2*(a + b*log(c*x**n))*(d + e*x)**(sympy.S(3)/2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_137():
    f = x*(a + b*log(c*x**n))/(d + e*x)**(sympy.S(3)/2)
    F = 8*b*sqrt(d)*n*atanh(sqrt(d + e*x)/sqrt(d))/e**2 - 4*b*n*sqrt(d + e*x)/e**2 + 2*d*(a + b*log(c*x**n))/(e**2*sqrt(d + e*x)) + 2*(a + b*log(c*x**n))*sqrt(d + e*x)/e**2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_138():
    f = (a + b*log(c*x**n))/(d + e*x)**(sympy.S(3)/2)
    F = -4*b*n*atanh(sqrt(d + e*x)/sqrt(d))/(sqrt(d)*e) - (2*a + 2*b*log(c*x**n))/(e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_139():
    f = (a + b*log(c*x**n))/(x*(d + e*x)**(sympy.S(3)/2))
    F = ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))))))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_140():
    f = x**2/((a + b*log(c*x**n))*(d + e*x))
    F = sympy.Function('Unintegrable')(((x)**(Integer(2)) * (((Symbol('d') + (Symbol('e') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_141():
    f = x/((a + b*log(c*x**n))*(d + e*x))
    F = sympy.Function('Unintegrable')((x * (((Symbol('d') + (Symbol('e') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_142():
    f = 1/((a + b*log(c*x**n))*(d + e*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_143():
    f = 1/(x*(a + b*log(c*x**n))*(d + e*x))
    F = sympy.Function('Unintegrable')(((x * (Symbol('d') + (Symbol('e') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_144():
    f = 1/(x**2*(a + b*log(c*x**n))*(d + e*x))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('d') + (Symbol('e') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_145():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x)**3
    F = -b*d**3*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 3*b*d**2*e*n*(f*x)**(m + 2)/(f**2*(m + 2)**2) - 3*b*d*e**2*n*(f*x)**(m + 3)/(f**3*(m + 3)**2) - b*e**3*n*(f*x)**(m + 4)/(f**4*(m + 4)**2) + d**3*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 3*d**2*e*(f*x)**(m + 2)*(a + b*log(c*x**n))/(f**2*(m + 2)) + 3*d*e**2*(f*x)**(m + 3)*(a + b*log(c*x**n))/(f**3*(m + 3)) + e**3*(f*x)**(m + 4)*(a + b*log(c*x**n))/(f**4*(m + 4))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_146():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x)**2
    F = -b*d**2*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 2*b*d*e*n*(f*x)**(m + 2)/(f**2*(m + 2)**2) - b*e**2*n*(f*x)**(m + 3)/(f**3*(m + 3)**2) + d**2*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 2*d*e*(f*x)**(m + 2)*(a + b*log(c*x**n))/(f**2*(m + 2)) + e**2*(f*x)**(m + 3)*(a + b*log(c*x**n))/(f**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_147():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x)
    F = -b*d*n*(f*x)**(m + 1)/(f*(m + 1)**2) - b*e*n*(f*x)**(m + 2)/(f**2*(m + 2)**2) + d*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + e*(f*x)**(m + 2)*(a + b*log(c*x**n))/(f**2*(m + 2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_148():
    f = (f*x)**m*(a + b*log(c*x**n))
    F = -b*n*(f*x)**(m + 1)/(f*(m + 1)**2) + (f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_149():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_150():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_151():
    f = x*(a + b*x)**m*log(c*x**n)
    F = sympy.Function('Unintegrable')((x * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('m')) * sympy.log((Symbol('c') * (x)**(Symbol('n'))))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_152():
    f = (a + b*x)**m*log(c*x**n)
    F = (a + b*x)**(m + 1)*log(c*x**n)/(b*(m + 1)) + n*(a + b*x)**(m + 2)*hyper((1, m + 2), (m + 3,), 1 + b*x/a)/(a*b*(m**2 + 3*m + 2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_153():
    f = (a + b*x)**m*log(c*x**n)/x
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * x)))**(Symbol('m')) * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_154():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x**6/36 - b*e*n*x**8/64 + (a + b*log(c*x**n))*(d*x**6/6 + e*x**8/8)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_155():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x**4/16 - b*e*n*x**6/36 + (a + b*log(c*x**n))*(d*x**4/4 + e*x**6/6)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_156():
    f = x*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x**2/4 - b*e*n*x**4/16 + (a + b*log(c*x**n))*(d*x**2/2 + e*x**4/4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_157():
    f = (a + b*log(c*x**n))*(d + e*x**2)/x
    F = -b*e*n*x**2/4 + e*x**2*(a + b*log(c*x**n))/2 + d*(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_158():
    f = (a + b*log(c*x**n))*(d + e*x**2)/x**5
    F = -b*d*n/(16*x**4) - b*e*n/(4*x**2) - d*(a + b*log(c*x**n))/(4*x**4) - e*(a + b*log(c*x**n))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_159():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x**5/25 - b*e*n*x**7/49 + (a + b*log(c*x**n))*(d*x**5/5 + e*x**7/7)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_160():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x**3/9 - b*e*n*x**5/25 + (a + b*log(c*x**n))*(d*x**3/3 + e*x**5/5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_161():
    f = (a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*x - b*e*n*x**3/9 + d*x*(a + b*log(c*x**n)) + e*x**3*(a + b*log(c*x**n))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_162():
    f = (a + b*log(c*x**n))*(d + e*x**2)/x**2
    F = -b*d*n/x - b*e*n*x - d*(a + b*log(c*x**n))/x + e*x*(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_163():
    f = (a + b*log(c*x**n))*(d + e*x**2)/x**4
    F = -b*d*n/(9*x**3) - b*e*n/x - d*(a + b*log(c*x**n))/(3*x**3) - e*(a + b*log(c*x**n))/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_164():
    f = (a + b*log(c*x**n))*(d + e*x**2)/x**6
    F = -b*d*n/(25*x**5) - b*e*n/(9*x**3) - d*(a + b*log(c*x**n))/(5*x**5) - e*(a + b*log(c*x**n))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_165():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*x**6/36 - b*d*e*n*x**8/32 - b*e**2*n*x**10/100 + (a + b*log(c*x**n))*(d**2*x**6/6 + d*e*x**8/4 + e**2*x**10/10)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_166():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*x**4/16 - b*d*e*n*x**6/18 - b*e**2*n*x**8/64 + (a + b*log(c*x**n))*(d**2*x**4/4 + d*e*x**6/3 + e**2*x**8/8)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_167():
    f = x*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**3*n*log(x)/(6*e) - b*d**2*n*x**2/4 - b*d*e*n*x**4/8 - b*e**2*n*x**6/36 + (a + b*log(c*x**n))*(d + e*x**2)**3/(6*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_168():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x
    F = -b*d**2*n*log(x)**2/2 - b*d*e*n*x**2/2 - b*e**2*n*x**4/16 + d**2*(a + b*log(c*x**n))*log(x) + d*e*x**2*(a + b*log(c*x**n)) + e**2*x**4*(a + b*log(c*x**n))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_169():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**3
    F = -b*d**2*n/(4*x**2) - b*d*e*n*log(x)**2 - b*e**2*n*x**2/4 - d**2*(a + b*log(c*x**n))/(2*x**2) + 2*d*e*(a + b*log(c*x**n))*log(x) + e**2*x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_170():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**5
    F = -b*d**2*n/(16*x**4) - b*d*e*n/(2*x**2) - b*e**2*n*log(x)**2/2 - d**2*(a + b*log(c*x**n))/(4*x**4) - d*e*(a + b*log(c*x**n))/x**2 + e**2*(a + b*log(c*x**n))*log(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_171():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*x**5/25 - 2*b*d*e*n*x**7/49 - b*e**2*n*x**9/81 + (a + b*log(c*x**n))*(d**2*x**5/5 + 2*d*e*x**7/7 + e**2*x**9/9)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_172():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*x**3/9 - 2*b*d*e*n*x**5/25 - b*e**2*n*x**7/49 + (a + b*log(c*x**n))*(d**2*x**3/3 + 2*d*e*x**5/5 + e**2*x**7/7)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_173():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*x - 2*b*d*e*n*x**3/9 - b*e**2*n*x**5/25 + d**2*x*(a + b*log(c*x**n)) + 2*d*e*x**3*(a + b*log(c*x**n))/3 + e**2*x**5*(a + b*log(c*x**n))/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_174():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**2
    F = -b*d**2*n/x - 2*b*d*e*n*x - b*e**2*n*x**3/9 - d**2*(a + b*log(c*x**n))/x + 2*d*e*x*(a + b*log(c*x**n)) + e**2*x**3*(a + b*log(c*x**n))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_175():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**4
    F = -b*d**2*n/(9*x**3) - 2*b*d*e*n/x - b*e**2*n*x - d**2*(a + b*log(c*x**n))/(3*x**3) - 2*d*e*(a + b*log(c*x**n))/x + e**2*x*(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_176():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**6
    F = -b*d**2*n/(25*x**5) - 2*b*d*e*n/(9*x**3) - b*e**2*n/x - d**2*(a + b*log(c*x**n))/(5*x**5) - 2*d*e*(a + b*log(c*x**n))/(3*x**3) - e**2*(a + b*log(c*x**n))/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_177():
    f = (a + b*log(c*x**n))*(d + e*x**2)**2/x**8
    F = -b*d**2*n/(49*x**7) - 2*b*d*e*n/(25*x**5) - b*e**2*n/(9*x**3) - d**2*(a + b*log(c*x**n))/(7*x**7) - 2*d*e*(a + b*log(c*x**n))/(5*x**5) - e**2*(a + b*log(c*x**n))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_178():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**3*n*x**6/36 - 3*b*d**2*e*n*x**8/64 - 3*b*d*e**2*n*x**10/100 - b*e**3*n*x**12/144 + (a + b*log(c*x**n))*(d**3*x**6/6 + 3*d**2*e*x**8/8 + 3*d*e**2*x**10/10 + e**3*x**12/12)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_179():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = b*d**5*n*log(x)/(40*e**2) + b*d**4*n*x**2/(20*e) + 3*b*d**3*n*x**4/80 + b*d**2*e*n*x**6/60 + b*d*e**2*n*x**8/320 - b*n*(d + e*x**2)**5/(100*e**2) - (a + b*log(c*x**n))*(d*(d + e*x**2)**4/(8*e**2) - (d + e*x**2)**5/(10*e**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_180():
    f = x*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**4*n*log(x)/(8*e) - b*d**3*n*x**2/4 - 3*b*d**2*e*n*x**4/16 - b*d*e**2*n*x**6/12 - b*e**3*n*x**8/64 + (a + b*log(c*x**n))*(d + e*x**2)**4/(8*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_181():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x
    F = -b*d**3*n*log(x)**2/2 - 3*b*d**2*e*n*x**2/4 - 3*b*d*e**2*n*x**4/16 - b*e**3*n*x**6/36 + d**3*(a + b*log(c*x**n))*log(x) + 3*d**2*e*x**2*(a + b*log(c*x**n))/2 + 3*d*e**2*x**4*(a + b*log(c*x**n))/4 + e**3*x**6*(a + b*log(c*x**n))/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_182():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**3
    F = -b*d**3*n/(4*x**2) - 3*b*d**2*e*n*log(x)**2/2 - 3*b*d*e**2*n*x**2/4 - b*e**3*n*x**4/16 - d**3*(a + b*log(c*x**n))/(2*x**2) + 3*d**2*e*(a + b*log(c*x**n))*log(x) + 3*d*e**2*x**2*(a + b*log(c*x**n))/2 + e**3*x**4*(a + b*log(c*x**n))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_183():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**5
    F = -b*d**3*n/(16*x**4) - 3*b*d**2*e*n/(4*x**2) - 3*b*d*e**2*n*log(x)**2/2 - b*e**3*n*x**2/4 - d**3*(a + b*log(c*x**n))/(4*x**4) - 3*d**2*e*(a + b*log(c*x**n))/(2*x**2) + 3*d*e**2*(a + b*log(c*x**n))*log(x) + e**3*x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_184():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**3*n*x**5/25 - 3*b*d**2*e*n*x**7/49 - b*d*e**2*n*x**9/27 - b*e**3*n*x**11/121 + (a + b*log(c*x**n))*(231*d**3*x**5 + 495*d**2*e*x**7 + 385*d*e**2*x**9 + 105*e**3*x**11)/1155
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_185():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**3*n*x**3/9 - 3*b*d**2*e*n*x**5/25 - 3*b*d*e**2*n*x**7/49 - b*e**3*n*x**9/81 + (a + b*log(c*x**n))*(d**3*x**3/3 + 3*d**2*e*x**5/5 + 3*d*e**2*x**7/7 + e**3*x**9/9)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_186():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**3*n*x - b*d**2*e*n*x**3/3 - 3*b*d*e**2*n*x**5/25 - b*e**3*n*x**7/49 + d**3*x*(a + b*log(c*x**n)) + d**2*e*x**3*(a + b*log(c*x**n)) + 3*d*e**2*x**5*(a + b*log(c*x**n))/5 + e**3*x**7*(a + b*log(c*x**n))/7
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_187():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**2
    F = -b*d**3*n/x - 3*b*d**2*e*n*x - b*d*e**2*n*x**3/3 - b*e**3*n*x**5/25 - d**3*(a + b*log(c*x**n))/x + 3*d**2*e*x*(a + b*log(c*x**n)) + d*e**2*x**3*(a + b*log(c*x**n)) + e**3*x**5*(a + b*log(c*x**n))/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_188():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**4
    F = -b*d**3*n/(9*x**3) - 3*b*d**2*e*n/x - 3*b*d*e**2*n*x - b*e**3*n*x**3/9 - d**3*(a + b*log(c*x**n))/(3*x**3) - 3*d**2*e*(a + b*log(c*x**n))/x + 3*d*e**2*x*(a + b*log(c*x**n)) + e**3*x**3*(a + b*log(c*x**n))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_189():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**6
    F = -b*d**3*n/(25*x**5) - b*d**2*e*n/(3*x**3) - 3*b*d*e**2*n/x - b*e**3*n*x - d**3*(a + b*log(c*x**n))/(5*x**5) - d**2*e*(a + b*log(c*x**n))/x**3 - 3*d*e**2*(a + b*log(c*x**n))/x + e**3*x*(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_190():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**8
    F = -b*d**3*n/(49*x**7) - 3*b*d**2*e*n/(25*x**5) - b*d*e**2*n/(3*x**3) - b*e**3*n/x - d**3*(a + b*log(c*x**n))/(7*x**7) - 3*d**2*e*(a + b*log(c*x**n))/(5*x**5) - d*e**2*(a + b*log(c*x**n))/x**3 - e**3*(a + b*log(c*x**n))/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_191():
    f = (a + b*log(c*x**n))*(d + e*x**2)**3/x**10
    F = -b*d**3*n/(81*x**9) - 3*b*d**2*e*n/(49*x**7) - 3*b*d*e**2*n/(25*x**5) - b*e**3*n/(9*x**3) - d**3*(a + b*log(c*x**n))/(9*x**9) - 3*d**2*e*(a + b*log(c*x**n))/(7*x**7) - 3*d*e**2*(a + b*log(c*x**n))/(5*x**5) - e**3*(a + b*log(c*x**n))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_192():
    f = x**5*(a + b*log(c*x**n))/(d + e*x**2)
    F = ((Symbol('b') * Symbol('d') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(4))) * ((Integer(16) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('e')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_193():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * Symbol('e')))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_194():
    f = x*(a + b*log(c*x**n))/(d + e*x**2)
    F = (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_195():
    f = (a + b*log(c*x**n))/(x*(d + e*x**2))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_196():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_197():
    f = (a + b*log(c*x**n))/(x**5*(d + e*x**2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(16) * Symbol('d') * (x)**(Integer(4))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * Symbol('d') * (x)**(Integer(4))))**(Integer(-1)))) + ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_198():
    f = x**4*(a + b*log(c*x**n))/(d + e*x**2)
    F = (Integer(-1) * ((Symbol('a') * Symbol('d') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(3))) * ((Integer(9) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e')))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_199():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**2)
    F = ((Symbol('a') * x) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * (Symbol('e'))**(Integer(-1)))) + ((Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_200():
    f = (a + b*log(c*x**n))/(d + e*x**2)
    F = ((sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_201():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_202():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x**2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(9) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n')) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(3) * Symbol('d') * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_203():
    f = x**5*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('d') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_204():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = (Integer(-1) * (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_205():
    f = x*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = -b*n*log(d + e*x**2)/(4*d*e) + x**2*(a + b*log(c*x**n))/(2*d*(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_206():
    f = (a + b*log(c*x**n))/(x*(d + e*x**2)**2)
    F = ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * ((Integer(2) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(2) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_207():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**2)**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(4) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(4) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(4) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_208():
    f = x**4*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = ((Symbol('a') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('d')) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_209():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = ((Symbol('b') * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_210():
    f = (a + b*log(c*x**n))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_211():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**2)**2)
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * x * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * ((Integer(3) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(3) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_212():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x**2)**2)
    F = (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('n')) * ((Integer(18) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(5) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(5) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(5) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(5) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * x))**(Integer(-1))) + (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * ((Integer(5) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(5) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_213():
    f = x**5*(a + b*log(c*x**n))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('d') * Symbol('n')) * ((Integer(8) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.log(x)) * ((Integer(4) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(8) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_214():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**2)**3
    F = -b*n/(8*e**2*(d + e*x**2)) - b*n*log(d + e*x**2)/(8*d*e**2) + x**4*(a + b*log(c*x**n))/(4*d*(d + e*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_215():
    f = x*(a + b*log(c*x**n))/(d + e*x**2)**3
    F = b*n/(8*d*e*(d + e*x**2)) + b*n*log(x)/(4*d**2*e) - b*n*log(d + e*x**2)/(8*d**2*e) - (a + b*log(c*x**n))/(4*e*(d + e*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_216():
    f = (a + b*log(c*x**n))/(x*(d + e*x**2)**3)
    F = ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * Symbol('d') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')) + (Integer(-1) * (Integer(3) * Symbol('b') * Symbol('n'))) + (Integer(4) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (((Integer(4) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(4) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_217():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**2)**3)
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n')) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * Symbol('d') * (x)**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (((Integer(6) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(6) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(12) * Symbol('a')) + (Integer(-1) * (Integer(5) * Symbol('b') * Symbol('n'))) + (Integer(12) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1))))) * ((Integer(12) * Symbol('a')) + (Integer(-1) * (Integer(5) * Symbol('b') * Symbol('n'))) + (Integer(12) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_218():
    f = x**4*(a + b*log(c*x**n))/(d + e*x**2)**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * ((Integer(8) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(3) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(16) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('d')) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_219():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('n') * x) * ((Integer(8) * Symbol('d') * Symbol('e') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('e') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * Symbol('d') * Symbol('e') * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_220():
    f = (a + b*log(c*x**n))/(d + e*x**2)**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('d') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_221():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**2)**3)
    F = (Integer(-1) * ((Integer(15) * Symbol('b') * Symbol('n')) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * Symbol('d') * x * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (((Integer(5) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(5) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * x * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(15) * Symbol('a')) + (Integer(-1) * (Integer(8) * Symbol('b') * Symbol('n'))) + (Integer(15) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * ((Integer(15) * Symbol('a')) + (Integer(-1) * (Integer(8) * Symbol('b') * Symbol('n'))) + (Integer(15) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * sympy.I * Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_222():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x**2)**3)
    F = (Integer(-1) * ((Integer(35) * Symbol('b') * Symbol('n')) * ((Integer(72) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(35) * Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(8) * (Symbol('d'))**(Integer(4)) * x))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(4) * Symbol('d') * (x)**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (((Integer(7) * Symbol('a')) + (Integer(-1) * (Symbol('b') * Symbol('n'))) + (Integer(7) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Integer(35) * Symbol('a')) + (Integer(-1) * (Integer(12) * Symbol('b') * Symbol('n'))) + (Integer(35) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(24) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(35) * Symbol('a')) + (Integer(-1) * (Integer(12) * Symbol('b') * Symbol('n'))) + (Integer(35) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**(Integer(4)) * x))**(Integer(-1))) + (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * ((Integer(35) * Symbol('a')) + (Integer(-1) * (Integer(12) * Symbol('b') * Symbol('n'))) + (Integer(35) * Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(35) * sympy.I * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_223():
    f = x*log(c*x**2)/(-c*x**2 + 1)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_224():
    f = x*log(x**2/c)/(c - x**2)
    F = (Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((x)**(Integer(2)) * (Symbol('c'))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_225():
    f = log(x)/(1 - x**2)
    F = (sympy.atanh(x) * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_226():
    f = log(x)/(x**2 + 1)
    F = (sympy.atan(x) * sympy.log(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_227():
    f = (a + b*log(c*x))/(-e*x**2 + 1)
    F = ((sympy.atanh((sympy.sqrt(Symbol('e')) * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.sqrt(Symbol('e'))) * x))) * ((Integer(2) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.sqrt(Symbol('e')) * x))) * ((Integer(2) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_228():
    f = (a + b*log(c*x**n))/(-e*x**2 + 1)
    F = ((sympy.atanh((sympy.sqrt(Symbol('e')) * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.sqrt(Symbol('e'))) * x))) * ((Integer(2) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.sqrt(Symbol('e')) * x))) * ((Integer(2) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_229():
    f = (a + b*log(c*x**n))**2/(d + e*x**2)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))))**(Integer(-1))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_230():
    f = (a + b*log(c*x**n))**3/(d + e*x**2)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))))**(Integer(-1))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_231():
    f = 1/((a + b*log(c*x**n))*(d + e*x**2)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_232():
    f = 1/((a + b*log(c*x**n))**2*(d + e*x**2)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_233():
    f = x*(a + b*log(c*x**n))*sqrt(d + e*x**2)
    F = b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*e) - b*d*n*sqrt(d + e*x**2)/(3*e) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*e) + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_234():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**2)/x
    F = ((Integer(-1) * Symbol('b')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + (Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + ((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_235():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**2)/x**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_236():
    f = x**4*(a + b*log(c*x**n))*sqrt(d + e*x**2)
    F = ((Integer(7) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(192) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('d') * Symbol('n') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(288) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(36))**(Integer(-1)) * Symbol('b') * Symbol('n') * (x)**(Integer(5)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))) + ((Integer(5) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(192) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(32) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(16) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(16) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(24) * Symbol('e')))**(Integer(-1))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(5)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(16) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(32) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_237():
    f = x**2*(a + b*log(c*x**n))*sqrt(d + e*x**2)
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('d') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(32) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * Symbol('b') * Symbol('n') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(32) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(8) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('d') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * Symbol('e')))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_238():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**2)
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * Symbol('b') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.atanh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(4) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_239():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**2)/x**2
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * (x)**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (x)**(Integer(-1)))) + ((sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_240():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**2)/x**4
    F = b*e**(sympy.S(3)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d) - b*e*n*sqrt(d + e*x**2)/(3*d*x) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*d*x**3) - (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*d*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_241():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)
    F = 8*b*d**(sympy.S(9)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(315*e**3) - 8*b*d**4*n*sqrt(d + e*x**2)/(315*e**3) - 8*b*d**3*n*(d + e*x**2)**(sympy.S(3)/2)/(945*e**3) - 8*b*d**2*n*(d + e*x**2)**(sympy.S(5)/2)/(1575*e**3) + 11*b*d*n*(d + e*x**2)**(sympy.S(7)/2)/(441*e**3) - b*n*(d + e*x**2)**(sympy.S(9)/2)/(81*e**3) + d**2*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**3) - 2*d*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(7)/2)/(7*e**3) + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(9)/2)/(9*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_242():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)
    F = -2*b*d**(sympy.S(7)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(35*e**2) + 2*b*d**3*n*sqrt(d + e*x**2)/(35*e**2) + 2*b*d**2*n*(d + e*x**2)**(sympy.S(3)/2)/(105*e**2) + 2*b*d*n*(d + e*x**2)**(sympy.S(5)/2)/(175*e**2) - b*n*(d + e*x**2)**(sympy.S(7)/2)/(49*e**2) - d*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**2) + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(7)/2)/(7*e**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_243():
    f = x*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)
    F = b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(5*e) - b*d**2*n*sqrt(d + e*x**2)/(5*e) - b*d*n*(d + e*x**2)**(sympy.S(3)/2)/(15*e) - b*n*(d + e*x**2)**(sympy.S(5)/2)/(25*e) + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_244():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x
    F = ((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(9))**(Integer(-1)) * Symbol('b') * Symbol('n') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + ((Integer(3))**(Integer(-1)) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_245():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**3
    F = ((Integer(-1) * Symbol('b')) * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('d')) * Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_246():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(11) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(192) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(23) * (Integer(288))**(Integer(-1))) * Symbol('b') * Symbol('d') * Symbol('n') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))) + (Integer(-1) * ((Integer(36))**(Integer(-1)) * Symbol('b') * Symbol('e') * Symbol('n') * (x)**(Integer(5)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(192) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(32) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(16) * Symbol('e')))**(Integer(-1))) + ((Integer(8))**(Integer(-1)) * Symbol('d') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(32) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_247():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)
    F = ((Integer(-1) * (Integer(9) * (Integer(32))**(Integer(-1)))) * Symbol('b') * Symbol('d') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * Symbol('b') * Symbol('n') * x * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(3) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(16) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.atanh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))))) * ((Integer(32) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(8) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('d') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + ((Integer(4))**(Integer(-1)) * x * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(8) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(16) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_248():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**2
    F = (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * Symbol('e') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))) + ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(2) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('e') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (x)**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(4) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_249():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**4
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_250():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**6
    F = b*e**(sympy.S(5)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(5*d) - b*e**2*n*sqrt(d + e*x**2)/(5*d*x) - b*e*n*(d + e*x**2)**(sympy.S(3)/2)/(15*d*x**3) - b*n*(d + e*x**2)**(sympy.S(5)/2)/(25*d*x**5) - (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*d*x**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_251():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**8
    F = -2*b*e**(sympy.S(7)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(35*d**2) + 2*b*e**3*n*sqrt(d + e*x**2)/(35*d**2*x) + 2*b*e**2*n*(d + e*x**2)**(sympy.S(3)/2)/(105*d**2*x**3) + 2*b*e*n*(d + e*x**2)**(sympy.S(5)/2)/(175*d**2*x**5) - b*n*(d + e*x**2)**(sympy.S(7)/2)/(49*d**2*x**7) - (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(7*d*x**7) + 2*e*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(35*d**2*x**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_252():
    f = (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/x**10
    F = -b*n*(d + e*x**2)**(sympy.S(7)/2)/(81*d**2*x**9) + 8*b*e**(sympy.S(9)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(315*d**3) - 8*b*e**4*n*sqrt(d + e*x**2)/(315*d**3*x) - 8*b*e**3*n*(d + e*x**2)**(sympy.S(3)/2)/(945*d**3*x**3) - 8*b*e**2*n*(d + e*x**2)**(sympy.S(5)/2)/(1575*d**3*x**5) + 50*b*e*n*(d + e*x**2)**(sympy.S(7)/2)/(3969*d**3*x**7) - (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(9*d*x**9) + 4*e*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(63*d**2*x**7) - 8*e**2*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(315*d**3*x**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_253():
    f = x*sqrt(x**2 + 4)*log(x)
    F = (x**2 + 4)**(sympy.S(3)/2)*log(x)/3 - (x**2 + 4)**(sympy.S(3)/2)/9 - 4*sqrt(x**2 + 4)/3 + 8*atanh(sqrt(x**2 + 4)/2)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_254():
    f = x**5*(a + b*log(c*x**n))/sqrt(d + e*x**2)
    F = 8*b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(15*e**3) - 8*b*d**2*n*sqrt(d + e*x**2)/(15*e**3) + 7*b*d*n*(d + e*x**2)**(sympy.S(3)/2)/(45*e**3) - b*n*(d + e*x**2)**(sympy.S(5)/2)/(25*e**3) + d**2*(a + b*log(c*x**n))*sqrt(d + e*x**2)/e**3 - 2*d*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**3) + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_255():
    f = x**3*(a + b*log(c*x**n))/sqrt(d + e*x**2)
    F = -2*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*e**2) + 2*b*d*n*sqrt(d + e*x**2)/(3*e**2) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*e**2) - d*(a + b*log(c*x**n))*sqrt(d + e*x**2)/e**2 + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_256():
    f = x*(a + b*log(c*x**n))/sqrt(d + e*x**2)
    F = b*sqrt(d)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/e - b*n*sqrt(d + e*x**2)/e + (a + b*log(c*x**n))*sqrt(d + e*x**2)/e
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_257():
    f = (a + b*log(c*x**n))/(x*sqrt(d + e*x**2))
    F = ((Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_258():
    f = (a + b*log(c*x**n))/(x**3*sqrt(d + e*x**2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_259():
    f = x**2*(a + b*log(c*x**n))/sqrt(d + e*x**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + ((x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_260():
    f = (a + b*log(c*x**n))/sqrt(d + e*x**2)
    F = ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_261():
    f = (a + b*log(c*x**n))/(x**2*sqrt(d + e*x**2))
    F = b*sqrt(e)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/d - b*n*sqrt(d + e*x**2)/(d*x) - (a + b*log(c*x**n))*sqrt(d + e*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_262():
    f = (a + b*log(c*x**n))/(x**4*sqrt(d + e*x**2))
    F = -2*b*e**(sympy.S(3)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d**2) + 2*b*e*n*sqrt(d + e*x**2)/(3*d**2*x) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*d**2*x**3) - (a + b*log(c*x**n))*sqrt(d + e*x**2)/(3*d*x**3) + 2*e*(a + b*log(c*x**n))*sqrt(d + e*x**2)/(3*d**2*x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_263():
    f = (a + b*log(c*x**n))/(x**6*sqrt(d + e*x**2))
    F = -b*n*(d + e*x**2)**(sympy.S(3)/2)/(25*d**2*x**5) + 8*b*e**(sympy.S(5)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(15*d**3) - 8*b*e**2*n*sqrt(d + e*x**2)/(15*d**3*x) + 26*b*e*n*(d + e*x**2)**(sympy.S(3)/2)/(225*d**3*x**3) - (a + b*log(c*x**n))*sqrt(d + e*x**2)/(5*d*x**5) + 4*e*(a + b*log(c*x**n))*sqrt(d + e*x**2)/(15*d**2*x**3) - 8*e**2*(a + b*log(c*x**n))*sqrt(d + e*x**2)/(15*d**3*x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_264():
    f = x**7*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = 16*b*d**(sympy.S(5)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(5*e**4) - 11*b*d**2*n*sqrt(d + e*x**2)/(5*e**4) + 4*b*d*n*(d + e*x**2)**(sympy.S(3)/2)/(15*e**4) - b*n*(d + e*x**2)**(sympy.S(5)/2)/(25*e**4) + d**3*(a + b*log(c*x**n))/(e**4*sqrt(d + e*x**2)) + 3*d**2*(a + b*log(c*x**n))*sqrt(d + e*x**2)/e**4 - d*(a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/e**4 + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_265():
    f = x**5*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = -8*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*e**3) + 5*b*d*n*sqrt(d + e*x**2)/(3*e**3) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*e**3) - d**2*(a + b*log(c*x**n))/(e**3*sqrt(d + e*x**2)) - 2*d*(a + b*log(c*x**n))*sqrt(d + e*x**2)/e**3 + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_266():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = 2*b*sqrt(d)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/e**2 - b*n*sqrt(d + e*x**2)/e**2 + d*(a + b*log(c*x**n))/(e**2*sqrt(d + e*x**2)) + (a + b*log(c*x**n))*sqrt(d + e*x**2)/e**2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_267():
    f = x*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = -b*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(sqrt(d)*e) - (a + b*log(c*x**n))/(e*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_268():
    f = (a + b*log(c*x**n))/(x*(d + e*x**2)**(sympy.S(3)/2))
    F = ((Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('d') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)) + (Integer(-1) * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_269():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_270():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_271():
    f = (a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(3)/2)
    F = -b*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(d*sqrt(e)) + x*(a + b*log(c*x**n))/(d*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_272():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**2)**(sympy.S(3)/2))
    F = 2*b*sqrt(e)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/d**2 - b*n*sqrt(d + e*x**2)/(d**2*x) - (a + b*log(c*x**n))/(d*x*sqrt(d + e*x**2)) - 2*e*x*(a + b*log(c*x**n))/(d**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_273():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x**2)**(sympy.S(3)/2))
    F = -b*n*sqrt(d + e*x**2)/(9*d**2*x**3) - 8*b*e**(sympy.S(3)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d**3) + 14*b*e*n*sqrt(d + e*x**2)/(9*d**3*x) - (a + b*log(c*x**n))/(3*d*x**3*sqrt(d + e*x**2)) + 4*e*(a + b*log(c*x**n))/(3*d**2*x*sqrt(d + e*x**2)) + 8*e**2*x*(a + b*log(c*x**n))/(3*d**3*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_274():
    f = (a + b*log(c*x**n))/(x**6*(d + e*x**2)**(sympy.S(3)/2))
    F = -b*n*sqrt(d + e*x**2)/(25*d**2*x**5) + 14*b*e*n*sqrt(d + e*x**2)/(75*d**3*x**3) + 16*b*e**(sympy.S(5)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(5*d**4) - 148*b*e**2*n*sqrt(d + e*x**2)/(75*d**4*x) - (a + b*log(c*x**n))/(5*d*x**5*sqrt(d + e*x**2)) + 2*e*(a + b*log(c*x**n))/(5*d**2*x**3*sqrt(d + e*x**2)) - 8*e**2*(a + b*log(c*x**n))/(5*d**3*x*sqrt(d + e*x**2)) - 16*e**3*x*(a + b*log(c*x**n))/(5*d**4*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_275():
    f = x**7*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = -16*b*d**(sympy.S(3)/2)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*e**4) - b*d**2*n/(3*e**4*sqrt(d + e*x**2)) + 8*b*d*n*sqrt(d + e*x**2)/(3*e**4) - b*n*(d + e*x**2)**(sympy.S(3)/2)/(9*e**4) + d**3*(a + b*log(c*x**n))/(3*e**4*(d + e*x**2)**(sympy.S(3)/2)) - 3*d**2*(a + b*log(c*x**n))/(e**4*sqrt(d + e*x**2)) - 3*d*(a + b*log(c*x**n))*sqrt(d + e*x**2)/e**4 + (a + b*log(c*x**n))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_276():
    f = x**5*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = 8*b*sqrt(d)*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*e**3) + b*d*n/(3*e**3*sqrt(d + e*x**2)) - b*n*sqrt(d + e*x**2)/e**3 - d**2*(a + b*log(c*x**n))/(3*e**3*(d + e*x**2)**(sympy.S(3)/2)) + 2*d*(a + b*log(c*x**n))/(e**3*sqrt(d + e*x**2)) + (a + b*log(c*x**n))*sqrt(d + e*x**2)/e**3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_277():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = -b*n/(3*e**2*sqrt(d + e*x**2)) - 2*b*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*sqrt(d)*e**2) + d*(a + b*log(c*x**n))/(3*e**2*(d + e*x**2)**(sympy.S(3)/2)) - (a + b*log(c*x**n))/(e**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_278():
    f = x*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = b*n/(3*d*e*sqrt(d + e*x**2)) - b*n*atanh(sqrt(d + e*x**2)/sqrt(d))/(3*d**(sympy.S(3)/2)*e) - (a + b*log(c*x**n))/(3*e*(d + e*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_279():
    f = (a + b*log(c*x**n))/(x*(d + e*x**2)**(sympy.S(5)/2))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (((Symbol('d') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)) + (Integer(3) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_280():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**2)**(sympy.S(5)/2))
    F = ((Symbol('b') * Symbol('e') * Symbol('n')) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(31) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(12) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(5) * Symbol('e') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_281():
    f = x**6*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = ((Integer(5) * Symbol('b') * Symbol('d') * Symbol('n') * x) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * (x)**(Integer(3))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('n') * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(31) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(12) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.atanh((sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * (((Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))))))) * ((Integer(2) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(5)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(5) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_282():
    f = x**4*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * (((Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_283():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = b*n*x/(3*d*e*sqrt(d + e*x**2)) - b*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d*e**(sympy.S(3)/2)) + x**3*(a + b*log(c*x**n))/(3*d*(d + e*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_284():
    f = (a + b*log(c*x**n))/(d + e*x**2)**(sympy.S(5)/2)
    F = -b*n*x/(3*d**2*sqrt(d + e*x**2)) - 2*b*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d**2*sqrt(e)) + x*(a + b*log(c*x**n))/(3*d*(d + e*x**2)**(sympy.S(3)/2)) + 2*x*(a + b*log(c*x**n))/(3*d**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_285():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**2)**(sympy.S(5)/2))
    F = -b*n/(d**2*x*sqrt(d + e*x**2)) + 8*b*sqrt(e)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d**3) - 2*b*e*n*x/(3*d**3*sqrt(d + e*x**2)) - (a + b*log(c*x**n))/(d*x*(d + e*x**2)**(sympy.S(3)/2)) - 4*e*x*(a + b*log(c*x**n))/(3*d**2*(d + e*x**2)**(sympy.S(3)/2)) - 8*e*x*(a + b*log(c*x**n))/(3*d**3*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_286():
    f = (a + b*log(c*x**n))/(x**4*(d + e*x**2)**(sympy.S(5)/2))
    F = -b*n*sqrt(d + e*x**2)/(9*d**3*x**3) - 16*b*e**(sympy.S(3)/2)*n*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*d**4) - b*e**2*n*x/(3*d**4*sqrt(d + e*x**2)) + 23*b*e*n*sqrt(d + e*x**2)/(9*d**4*x) - (a + b*log(c*x**n))/(3*d*x**3*(d + e*x**2)**(sympy.S(3)/2)) + 2*e*(a + b*log(c*x**n))/(d**2*x*(d + e*x**2)**(sympy.S(3)/2)) + 8*e**2*x*(a + b*log(c*x**n))/(3*d**3*(d + e*x**2)**(sympy.S(3)/2)) + 16*e**2*x*(a + b*log(c*x**n))/(3*d**4*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_287():
    f = x**3*(a + b*log(c*x**n))/(sqrt(d - e*x)*sqrt(d + e*x))
    F = -2*b*d**4*n*sqrt(1 - e**2*x**2/d**2)*atanh(sqrt(1 - e**2*x**2/d**2))/(3*e**4*sqrt(d - e*x)*sqrt(d + e*x)) + 2*b*d**2*n*(d**2 - e**2*x**2)/(3*e**4*sqrt(d - e*x)*sqrt(d + e*x)) - b*n*(d**2 - e**2*x**2)**2/(9*e**4*sqrt(d - e*x)*sqrt(d + e*x)) - d**2*(a + b*log(c*x**n))*(d**2 - e**2*x**2)/(e**4*sqrt(d - e*x)*sqrt(d + e*x)) + (a + b*log(c*x**n))*(d**2 - e**2*x**2)**2/(3*e**4*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_288():
    f = x*(a + b*log(c*x**n))/(sqrt(d - e*x)*sqrt(d + e*x))
    F = -b*d**2*n*sqrt(1 - e**2*x**2/d**2)*atanh(sqrt(1 - e**2*x**2/d**2))/(e**2*sqrt(d - e*x)*sqrt(d + e*x)) + b*n*(d**2 - e**2*x**2)/(e**2*sqrt(d - e*x)*sqrt(d + e*x)) - (a + b*log(c*x**n))*(d**2 - e**2*x**2)/(e**2*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_289():
    f = (a + b*log(c*x**n))/(x*sqrt(d - e*x)*sqrt(d + e*x))
    F = ((Symbol('b') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * (sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_290():
    f = (a + b*log(c*x**n))/(x**3*sqrt(d - e*x)*sqrt(d + e*x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * (sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_291():
    f = x**2*(a + b*log(c*x**n))/(sqrt(d - e*x)*sqrt(d + e*x))
    F = ((Symbol('b') * Symbol('n') * x * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * (sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))**(Integer(2))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(3)) * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_292():
    f = (a + b*log(c*x**n))/(sqrt(d - e*x)*sqrt(d + e*x))
    F = ((sympy.I * Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * (sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * Symbol('e') * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin(((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))))) * ((Integer(2) * Symbol('e') * sympy.sqrt((Symbol('d') + (Integer(-1) * (Symbol('e') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_293():
    f = (a + b*log(c*x**n))/(x**2*sqrt(d - e*x)*sqrt(d + e*x))
    F = -b*e*n*sqrt(1 - e**2*x**2/d**2)*asin(e*x/d)/(d*sqrt(d - e*x)*sqrt(d + e*x)) - b*n*(d**2 - e**2*x**2)/(d**2*x*sqrt(d - e*x)*sqrt(d + e*x)) - (a + b*log(c*x**n))*(d**2 - e**2*x**2)/(d**2*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_294():
    f = (a + b*log(c*x**n))/(x**4*sqrt(d - e*x)*sqrt(d + e*x))
    F = -2*b*e**3*n*sqrt(1 - e**2*x**2/d**2)*asin(e*x/d)/(3*d**3*sqrt(d - e*x)*sqrt(d + e*x)) - 2*b*e**2*n*(d**2 - e**2*x**2)/(3*d**4*x*sqrt(d - e*x)*sqrt(d + e*x)) - b*n*(d**2 - e**2*x**2)**2/(9*d**4*x**3*sqrt(d - e*x)*sqrt(d + e*x)) - (a + b*log(c*x**n))*(d**2 - e**2*x**2)/(3*d**2*x**3*sqrt(d - e*x)*sqrt(d + e*x)) - 2*e**2*(a + b*log(c*x**n))*(d**2 - e**2*x**2)/(3*d**4*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_295():
    f = x*log(x)/sqrt(x**2 - 1)
    F = sqrt(x**2 - 1)*log(x) - sqrt(x**2 - 1) + atan(sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_296():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**2)**3
    F = -b*d**3*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 3*b*d**2*e*n*(f*x)**(m + 3)/(f**3*(m + 3)**2) - 3*b*d*e**2*n*(f*x)**(m + 5)/(f**5*(m + 5)**2) - b*e**3*n*(f*x)**(m + 7)/(f**7*(m + 7)**2) + d**3*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 3*d**2*e*(f*x)**(m + 3)*(a + b*log(c*x**n))/(f**3*(m + 3)) + 3*d*e**2*(f*x)**(m + 5)*(a + b*log(c*x**n))/(f**5*(m + 5)) + e**3*(f*x)**(m + 7)*(a + b*log(c*x**n))/(f**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_297():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**2)**2
    F = -b*d**2*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 2*b*d*e*n*(f*x)**(m + 3)/(f**3*(m + 3)**2) - b*e**2*n*(f*x)**(m + 5)/(f**5*(m + 5)**2) + d**2*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 2*d*e*(f*x)**(m + 3)*(a + b*log(c*x**n))/(f**3*(m + 3)) + e**2*(f*x)**(m + 5)*(a + b*log(c*x**n))/(f**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_298():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**2)
    F = -b*d*n*(f*x)**(m + 1)/(f*(m + 1)**2) - b*e*n*(f*x)**(m + 3)/(f**3*(m + 3)**2) + d*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + e*(f*x)**(m + 3)*(a + b*log(c*x**n))/(f**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_299():
    f = (f*x)**m*(a + b*log(c*x**n))
    F = -b*n*(f*x)**(m + 1)/(f*(m + 1)**2) + (f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_300():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x**2)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_301():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x**2)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_302():
    f = (a + b*log(c*x**n))**3/(d + e*x**3)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1)))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3)) * sympy.log((Integer(1) + (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3) * (Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3)) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(6) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.sqrt(Integer(3)) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(6) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(12) * sympy.I * sympy.sqrt(Integer(3)) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * sympy.I * sympy.sqrt(Integer(3)) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(12) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_303():
    f = (a + b*log(c*x**n))**2/(d + e*x**3)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1)))) + ((x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * sympy.sqrt(Integer(3)) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(4) * sympy.I * sympy.sqrt(Integer(3)) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_304():
    f = (a + b*log(c*x**n))/(d + e*x**3)**2
    F = ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1))) + (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * sympy.log((((Integer(-1) * (Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1))))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.log(((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Symbol('e'))**((Integer(3))**(Integer(-1))) * x)))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('n') * sympy.log(((Symbol('d'))**((Integer(3))**(Integer(-1))) + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x)))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((Integer(9) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Integer(3)) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(5)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('d'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))))) * ((((Integer(1) + (Integer(-1))**((Integer(3))**(Integer(-1)))))**(Integer(4)) * (Symbol('d'))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (Symbol('e'))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_305():
    f = 1/((a + b*log(c*x**n))*(d + e*x**3)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(3)))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_306():
    f = 1/((a + b*log(c*x**n))**2*(d + e*x**3)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(3)))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_307():
    f = x**3*(a + b*log(c*x**n))/(d + e/x)
    F = (Integer(-1) * ((Symbol('a') * (Symbol('e'))**(Integer(3)) * x) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(3)) * Symbol('n') * x) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * (x)**(Integer(3))) * ((Integer(9) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(4))) * ((Integer(16) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(3)) * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('e'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(4)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_308():
    f = x**2*(a + b*log(c*x**n))/(d + e/x)
    F = ((Symbol('a') * (Symbol('e'))**(Integer(2)) * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(3))) * ((Integer(9) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(3)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_309():
    f = x*(a + b*log(c*x**n))/(d + e/x)
    F = (Integer(-1) * ((Symbol('a') * Symbol('e') * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('n') * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_310():
    f = (a + b*log(c*x**n))/(d + e/x)
    F = ((Symbol('a') * x) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * x) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('b') * x * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_311():
    f = (a + b*log(c*x**n))/(x*(d + e/x))
    F = (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_312():
    f = (a + b*log(c*x**n))/(x**2*(d + e/x))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('e') * ((Symbol('d') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (Symbol('e'))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('e') * ((Symbol('d') * x))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_313():
    f = (a + b*log(c*x**n))/(x**3*(d + e/x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Symbol('e') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('e') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * Symbol('n')))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_314():
    f = (a + b*log(c*x**n))/(x**4*(d + e/x))
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(4) * Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * Symbol('n')) * (((Symbol('e'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('e'))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(3)) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_315():
    f = x**3*(a + b*log(c*x))/(d + e/x)
    F = (Integer(-1) * ((Symbol('a') * (Symbol('e'))**(Integer(3)) * x) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(3)) * x) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * (x)**(Integer(3))) * ((Integer(9) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * ((Integer(16) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(3)) * x * sympy.log((Symbol('c') * x))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('e'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_316():
    f = x**2*(a + b*log(c*x))/(d + e/x)
    F = ((Symbol('a') * (Symbol('e'))**(Integer(2)) * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * (x)**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(3))) * ((Integer(9) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * x * sympy.log((Symbol('c') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_317():
    f = x*(a + b*log(c*x))/(d + e/x)
    F = (Integer(-1) * ((Symbol('a') * Symbol('e') * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * x * sympy.log((Symbol('c') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_318():
    f = (a + b*log(c*x))/(d + e/x)
    F = ((Symbol('a') * x) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * x) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('b') * x * sympy.log((Symbol('c') * x))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_319():
    f = (a + b*log(c*x))/(x*(d + e/x))
    F = (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_320():
    f = (a + b*log(c*x))/(x**2*(d + e/x))
    F = (Integer(-1) * ((sympy.log((Integer(1) + (Symbol('e') * ((Symbol('d') * x))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('e'))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('e') * ((Symbol('d') * x))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_321():
    f = (a + b*log(c*x))/(x**3*(d + e/x))
    F = (Integer(-1) * (Symbol('b') * ((Symbol('e') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * ((Symbol('e') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_322():
    f = (a + b*log(c*x))/(x**4*(d + e/x))
    F = (Integer(-1) * (Symbol('b') * ((Integer(4) * Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d')) * (((Symbol('e'))**(Integer(2)) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * ((Integer(2) * Symbol('e') * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (((Symbol('e'))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('d') * x) * (Symbol('e'))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_323():
    f = x**(n - 1)*log(e*x**n)/(-e*x**n + 1)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('e') * (x)**(Symbol('n')))))) * ((Symbol('e') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_324():
    f = x**(n - 1)*log(x**n/d)/(d - x**n)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((x)**(Symbol('n')) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('n'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_325():
    f = x**(n - 1)*log(-e*x**n/d)/(d + e*x**n)
    F = Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Symbol('e') * (x)**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('e') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_326():
    f = log(a/x)/(a*x - x**2)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('a') * (x)**(Integer(-1)))))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_327():
    f = log(a/x**2)/(a*x - x**3)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_328():
    f = log(a*x**(1 - n))/(a*x - x**n)
    F = Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('a') * (x)**((Integer(1) + (Integer(-1) * Symbol('n')))))))) * ((Symbol('a') * (Integer(1) + (Integer(-1) * Symbol('n')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_329():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))*(d + e*x**m)**3
    F = -b*d**4*n*x**(1 - m)*(f*x)**(m - 1)*log(x)/(4*e*m) - b*d**3*n*x*(f*x)**(m - 1)/m**2 - 3*b*d**2*e*n*x**(m + 1)*(f*x)**(m - 1)/(4*m**2) - b*d*e**2*n*x**(2*m + 1)*(f*x)**(m - 1)/(3*m**2) - b*e**3*n*x**(3*m + 1)*(f*x)**(m - 1)/(16*m**2) + x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))*(d + e*x**m)**4/(4*e*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_330():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))*(d + e*x**m)**2
    F = -b*d**3*n*x**(1 - m)*(f*x)**(m - 1)*log(x)/(3*e*m) - b*d**2*n*x*(f*x)**(m - 1)/m**2 - b*d*e*n*x**(m + 1)*(f*x)**(m - 1)/(2*m**2) - b*e**2*n*x**(2*m + 1)*(f*x)**(m - 1)/(9*m**2) + x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))*(d + e*x**m)**3/(3*e*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_331():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))
    F = -b*n*(f*x)**m/(f*m**2) + (f*x)**m*(a + b*log(c*x**n))/(f*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_332():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))/(d + e*x**m)
    F = (((x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Symbol('m'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e') * Symbol('m')))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Symbol('m'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_333():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))/(d + e*x**m)**2
    F = -b*n*(f*x)**m*log(d + e*x**m)/(d*e*f*m**2*x**m) + (f*x)**m*(a + b*log(c*x**n))/(d*f*m*(d + e*x**m))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_334():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))/(d + e*x**m)**3
    F = b*n*x**(1 - m)*(f*x)**(m - 1)/(2*d*e*m**2*(d + e*x**m)) + b*n*x**(1 - m)*(f*x)**(m - 1)*log(x)/(2*d**2*e*m) - b*n*x**(1 - m)*(f*x)**(m - 1)*log(d + e*x**m)/(2*d**2*e*m**2) - x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(2*e*m*(d + e*x**m)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_335():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))/(d + e*x**m)**4
    F = b*n*x**(1 - m)*(f*x)**(m - 1)/(6*d*e*m**2*(d + e*x**m)**2) + b*n*x**(1 - m)*(f*x)**(m - 1)/(3*d**2*e*m**2*(d + e*x**m)) + b*n*x**(1 - m)*(f*x)**(m - 1)*log(x)/(3*d**3*e*m) - b*n*x**(1 - m)*(f*x)**(m - 1)*log(d + e*x**m)/(3*d**3*e*m**2) - x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(3*e*m*(d + e*x**m)**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_336():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)**3
    F = b**2*d**4*n**2*x**(1 - m)*(f*x)**(m - 1)*log(x)**2/(4*e*m) + 2*b**2*d**3*n**2*x*(f*x)**(m - 1)/m**3 + 3*b**2*d**2*e*n**2*x**(m + 1)*(f*x)**(m - 1)/(4*m**3) + 2*b**2*d*e**2*n**2*x**(2*m + 1)*(f*x)**(m - 1)/(9*m**3) + b**2*e**3*n**2*x**(3*m + 1)*(f*x)**(m - 1)/(32*m**3) - b*d**4*n*x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))*log(x)/(2*e*m) - 2*b*d**3*n*x*(f*x)**(m - 1)*(a + b*log(c*x**n))/m**2 - 3*b*d**2*e*n*x**(m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(2*m**2) - 2*b*d*e**2*n*x**(2*m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(3*m**2) - b*e**3*n*x**(3*m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(8*m**2) + x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)**4/(4*e*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_337():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)**2
    F = b**2*d**3*n**2*x**(1 - m)*(f*x)**(m - 1)*log(x)**2/(3*e*m) + 2*b**2*d**2*n**2*x*(f*x)**(m - 1)/m**3 + b**2*d*e*n**2*x**(m + 1)*(f*x)**(m - 1)/(2*m**3) + 2*b**2*e**2*n**2*x**(2*m + 1)*(f*x)**(m - 1)/(27*m**3) - 2*b*d**3*n*x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))*log(x)/(3*e*m) - 2*b*d**2*n*x*(f*x)**(m - 1)*(a + b*log(c*x**n))/m**2 - b*d*e*n*x**(m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/m**2 - 2*b*e**2*n*x**(2*m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(9*m**2) + x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)**3/(3*e*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_338():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)
    F = b**2*d**2*n**2*x**(1 - m)*(f*x)**(m - 1)*log(x)**2/(2*e*m) + 2*b**2*d*n**2*x*(f*x)**(m - 1)/m**3 + b**2*e*n**2*x**(m + 1)*(f*x)**(m - 1)/(4*m**3) - b*d**2*n*x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))*log(x)/(e*m) - 2*b*d*n*x*(f*x)**(m - 1)*(a + b*log(c*x**n))/m**2 - b*e*n*x**(m + 1)*(f*x)**(m - 1)*(a + b*log(c*x**n))/(2*m**2) + x**(1 - m)*(f*x)**(m - 1)*(a + b*log(c*x**n))**2*(d + e*x**m)**2/(2*e*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_339():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*(f*x)**m/(f*m**3) - 2*b*n*(f*x)**m*(a + b*log(c*x**n))/(f*m**2) + (f*x)**m*(a + b*log(c*x**n))**2/(f*m)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_340():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2/(d + e*x**m)
    F = (((x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * (x)**(Symbol('m'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e') * Symbol('m')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * (x)**(Symbol('m'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * (x)**(Symbol('m'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_341():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2/(d + e*x**m)**2
    F = (Integer(-1) * (((x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('e') * Symbol('m') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('m'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_342():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2/(d + e*x**m)**3
    F = (Integer(-1) * ((Symbol('b') * Symbol('n') * x * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('m'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('m'))))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * Symbol('e') * Symbol('m') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('m')))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('m')))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_343():
    f = (f*x)**(m - 1)*(a + b*log(c*x**n))**2/(d + e*x**m)**4
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m')))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('m'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('m'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.log(x)) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * Symbol('d') * Symbol('e') * (Symbol('m'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('m')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * x * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('m'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('m'))))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(3) * Symbol('e') * Symbol('m') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('m')))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e') * (Symbol('m'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('m')))))) * (((Symbol('d'))**(Integer(3)) * Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**((Integer(1) + (Integer(-1) * Symbol('m')))) * ((Symbol('f') * x))**((Integer(-1) + Symbol('m'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('m')) * Symbol('e')))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e') * (Symbol('m'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_344():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x**6/36 - b*e*n*x**(r + 6)/(r + 6)**2 + (a + b*log(c*x**n))*(d*x**6 + 6*e*x**(r + 6)/(r + 6))/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_345():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x**4/16 - b*e*n*x**(r + 4)/(r + 4)**2 + (a + b*log(c*x**n))*(d*x**4 + 4*e*x**(r + 4)/(r + 4))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_346():
    f = x*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x**2/4 - b*e*n*x**(r + 2)/(r + 2)**2 + (a + b*log(c*x**n))*(d*x**2 + 2*e*x**(r + 2)/(r + 2))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_347():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x
    F = -b*e*n*x**r/r**2 + e*x**r*(a + b*log(c*x**n))/r + d*(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_348():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x**3
    F = -b*d*n/(4*x**2) - b*e*n*x**(r - 2)/(2 - r)**2 - d*(a + b*log(c*x**n))/(2*x**2) - e*x**(r - 2)*(a + b*log(c*x**n))/(2 - r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_349():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x**5
    F = -b*d*n/(16*x**4) - b*e*n*x**(r - 4)/(4 - r)**2 - d*(a + b*log(c*x**n))/(4*x**4) - e*x**(r - 4)*(a + b*log(c*x**n))/(4 - r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_350():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x**5/25 - b*e*n*x**(r + 5)/(r + 5)**2 + (a + b*log(c*x**n))*(d*x**5 + 5*e*x**(r + 5)/(r + 5))/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_351():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x**3/9 - b*e*n*x**(r + 3)/(r + 3)**2 + (a + b*log(c*x**n))*(d*x**3 + 3*e*x**(r + 3)/(r + 3))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_352():
    f = (a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*x - b*e*n*x**(r + 1)/(r + 1)**2 + d*x*(a + b*log(c*x**n)) + e*x**(r + 1)*(a + b*log(c*x**n))/(r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_353():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x**2
    F = -b*d*n/x - b*e*n*x**(r - 1)/(1 - r)**2 - d*(a + b*log(c*x**n))/x - e*x**(r - 1)*(a + b*log(c*x**n))/(1 - r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_354():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x**4
    F = -b*d*n/(9*x**3) - b*e*n*x**(r - 3)/(3 - r)**2 - d*(a + b*log(c*x**n))/(3*x**3) - e*x**(r - 3)*(a + b*log(c*x**n))/(3 - r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_355():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x**6
    F = -b*d*n/(25*x**5) - b*e*n*x**(r - 5)/(5 - r)**2 - d*(a + b*log(c*x**n))/(5*x**5) - e*x**(r - 5)*(a + b*log(c*x**n))/(5 - r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_356():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x**6/36 - 2*b*d*e*n*x**(r + 6)/(r + 6)**2 - b*e**2*n*x**(2*r + 6)/(4*(r + 3)**2) + (a + b*log(c*x**n))*(d**2*x**6 + 12*d*e*x**(r + 6)/(r + 6) + 3*e**2*x**(2*r + 6)/(r + 3))/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_357():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x**4/16 - 2*b*d*e*n*x**(r + 4)/(r + 4)**2 - b*e**2*n*x**(2*r + 4)/(4*(r + 2)**2) + (a + b*log(c*x**n))*(d**2*x**4 + 8*d*e*x**(r + 4)/(r + 4) + 2*e**2*x**(2*r + 4)/(r + 2))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_358():
    f = x*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x**2/4 - 2*b*d*e*n*x**(r + 2)/(r + 2)**2 - b*e**2*n*x**(2*r + 2)/(4*(r + 1)**2) + (a + b*log(c*x**n))*(d**2*x**2 + 4*d*e*x**(r + 2)/(r + 2) + e**2*x**(2*r + 2)/(r + 1))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_359():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x
    F = -b*d**2*n*log(x)**2/2 - 2*b*d*e*n*x**r/r**2 - b*e**2*n*x**(2*r)/(4*r**2) + d**2*(a + b*log(c*x**n))*log(x) + 2*d*e*x**r*(a + b*log(c*x**n))/r + e**2*x**(2*r)*(a + b*log(c*x**n))/(2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_360():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**3
    F = -b*d**2*n/(4*x**2) - 2*b*d*e*n*x**(r - 2)/(2 - r)**2 - b*e**2*n*x**(2*r - 2)/(4*(1 - r)**2) - d**2*(a + b*log(c*x**n))/(2*x**2) - 2*d*e*x**(r - 2)*(a + b*log(c*x**n))/(2 - r) - e**2*x**(2*r - 2)*(a + b*log(c*x**n))/(2 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_361():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**5
    F = -b*d**2*n/(16*x**4) - 2*b*d*e*n*x**(r - 4)/(4 - r)**2 - b*e**2*n*x**(2*r - 4)/(4*(2 - r)**2) - d**2*(a + b*log(c*x**n))/(4*x**4) - 2*d*e*x**(r - 4)*(a + b*log(c*x**n))/(4 - r) - e**2*x**(2*r - 4)*(a + b*log(c*x**n))/(4 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_362():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x**5/25 - 2*b*d*e*n*x**(r + 5)/(r + 5)**2 - b*e**2*n*x**(2*r + 5)/(2*r + 5)**2 + (a + b*log(c*x**n))*(d**2*x**5 + 10*d*e*x**(r + 5)/(r + 5) + 5*e**2*x**(2*r + 5)/(2*r + 5))/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_363():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x**3/9 - 2*b*d*e*n*x**(r + 3)/(r + 3)**2 - b*e**2*n*x**(2*r + 3)/(2*r + 3)**2 + (a + b*log(c*x**n))*(d**2*x**3 + 6*d*e*x**(r + 3)/(r + 3) + 3*e**2*x**(2*r + 3)/(2*r + 3))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_364():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*x - 2*b*d*e*n*x**(r + 1)/(r + 1)**2 - b*e**2*n*x**(2*r + 1)/(2*r + 1)**2 + d**2*x*(a + b*log(c*x**n)) + 2*d*e*x**(r + 1)*(a + b*log(c*x**n))/(r + 1) + e**2*x**(2*r + 1)*(a + b*log(c*x**n))/(2*r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_365():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**2
    F = -b*d**2*n/x - 2*b*d*e*n*x**(r - 1)/(1 - r)**2 - b*e**2*n*x**(2*r - 1)/(1 - 2*r)**2 - d**2*(a + b*log(c*x**n))/x - 2*d*e*x**(r - 1)*(a + b*log(c*x**n))/(1 - r) - e**2*x**(2*r - 1)*(a + b*log(c*x**n))/(1 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_366():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**4
    F = -b*d**2*n/(9*x**3) - 2*b*d*e*n*x**(r - 3)/(3 - r)**2 - b*e**2*n*x**(2*r - 3)/(3 - 2*r)**2 - d**2*(a + b*log(c*x**n))/(3*x**3) - 2*d*e*x**(r - 3)*(a + b*log(c*x**n))/(3 - r) - e**2*x**(2*r - 3)*(a + b*log(c*x**n))/(3 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_367():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**6
    F = -b*d**2*n/(25*x**5) - 2*b*d*e*n*x**(r - 5)/(5 - r)**2 - b*e**2*n*x**(2*r - 5)/(5 - 2*r)**2 - d**2*(a + b*log(c*x**n))/(5*x**5) - 2*d*e*x**(r - 5)*(a + b*log(c*x**n))/(5 - r) - e**2*x**(2*r - 5)*(a + b*log(c*x**n))/(5 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_368():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x**8
    F = -b*d**2*n/(49*x**7) - 2*b*d*e*n*x**(r - 7)/(7 - r)**2 - b*e**2*n*x**(2*r - 7)/(7 - 2*r)**2 - d**2*(a + b*log(c*x**n))/(7*x**7) - 2*d*e*x**(r - 7)*(a + b*log(c*x**n))/(7 - r) - e**2*x**(2*r - 7)*(a + b*log(c*x**n))/(7 - 2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_369():
    f = x**5*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x**6/36 - 3*b*d**2*e*n*x**(r + 6)/(r + 6)**2 - 3*b*d*e**2*n*x**(2*r + 6)/(4*(r + 3)**2) - b*e**3*n*x**(3*r + 6)/(9*(r + 2)**2) + (a + b*log(c*x**n))*(d**3*x**6 + 18*d**2*e*x**(r + 6)/(r + 6) + 9*d*e**2*x**(2*r + 6)/(r + 3) + 2*e**3*x**(3*r + 6)/(r + 2))/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_370():
    f = x**3*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x**4/16 - 3*b*d**2*e*n*x**(r + 4)/(r + 4)**2 - 3*b*d*e**2*n*x**(2*r + 4)/(4*(r + 2)**2) - b*e**3*n*x**(3*r + 4)/(3*r + 4)**2 + (a + b*log(c*x**n))*(d**3*x**4 + 12*d**2*e*x**(r + 4)/(r + 4) + 6*d*e**2*x**(2*r + 4)/(r + 2) + 4*e**3*x**(3*r + 4)/(3*r + 4))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_371():
    f = x*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x**2/4 - 3*b*d**2*e*n*x**(r + 2)/(r + 2)**2 - 3*b*d*e**2*n*x**(2*r + 2)/(4*(r + 1)**2) - b*e**3*n*x**(3*r + 2)/(3*r + 2)**2 + (a + b*log(c*x**n))*(d**3*x**2 + 6*d**2*e*x**(r + 2)/(r + 2) + 3*d*e**2*x**(2*r + 2)/(r + 1) + 2*e**3*x**(3*r + 2)/(3*r + 2))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_372():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x
    F = -b*d**3*n*log(x)**2/2 - 3*b*d**2*e*n*x**r/r**2 - 3*b*d*e**2*n*x**(2*r)/(4*r**2) - b*e**3*n*x**(3*r)/(9*r**2) + d**3*(a + b*log(c*x**n))*log(x) + 3*d**2*e*x**r*(a + b*log(c*x**n))/r + 3*d*e**2*x**(2*r)*(a + b*log(c*x**n))/(2*r) + e**3*x**(3*r)*(a + b*log(c*x**n))/(3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_373():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**3
    F = -b*d**3*n/(4*x**2) - 3*b*d**2*e*n*x**(r - 2)/(2 - r)**2 - 3*b*d*e**2*n*x**(2*r - 2)/(4*(1 - r)**2) - b*e**3*n*x**(3*r - 2)/(2 - 3*r)**2 - d**3*(a + b*log(c*x**n))/(2*x**2) - 3*d**2*e*x**(r - 2)*(a + b*log(c*x**n))/(2 - r) - 3*d*e**2*x**(2*r - 2)*(a + b*log(c*x**n))/(2 - 2*r) - e**3*x**(3*r - 2)*(a + b*log(c*x**n))/(2 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_374():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**5
    F = -b*d**3*n/(16*x**4) - 3*b*d**2*e*n*x**(r - 4)/(4 - r)**2 - 3*b*d*e**2*n*x**(2*r - 4)/(4*(2 - r)**2) - b*e**3*n*x**(3*r - 4)/(4 - 3*r)**2 - d**3*(a + b*log(c*x**n))/(4*x**4) - 3*d**2*e*x**(r - 4)*(a + b*log(c*x**n))/(4 - r) - 3*d*e**2*x**(2*r - 4)*(a + b*log(c*x**n))/(4 - 2*r) - e**3*x**(3*r - 4)*(a + b*log(c*x**n))/(4 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_375():
    f = x**4*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x**5/25 - 3*b*d**2*e*n*x**(r + 5)/(r + 5)**2 - 3*b*d*e**2*n*x**(2*r + 5)/(2*r + 5)**2 - b*e**3*n*x**(3*r + 5)/(3*r + 5)**2 + (a + b*log(c*x**n))*(d**3*x**5 + 15*d**2*e*x**(r + 5)/(r + 5) + 15*d*e**2*x**(2*r + 5)/(2*r + 5) + 5*e**3*x**(3*r + 5)/(3*r + 5))/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_376():
    f = x**2*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x**3/9 - 3*b*d**2*e*n*x**(r + 3)/(r + 3)**2 - 3*b*d*e**2*n*x**(2*r + 3)/(2*r + 3)**2 - b*e**3*n*x**(3*r + 3)/(9*(r + 1)**2) + (a + b*log(c*x**n))*(d**3*x**3 + 9*d**2*e*x**(r + 3)/(r + 3) + 9*d*e**2*x**(2*r + 3)/(2*r + 3) + e**3*x**(3*r + 3)/(r + 1))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_377():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*x - 3*b*d**2*e*n*x**(r + 1)/(r + 1)**2 - 3*b*d*e**2*n*x**(2*r + 1)/(2*r + 1)**2 - b*e**3*n*x**(3*r + 1)/(3*r + 1)**2 + d**3*x*(a + b*log(c*x**n)) + 3*d**2*e*x**(r + 1)*(a + b*log(c*x**n))/(r + 1) + 3*d*e**2*x**(2*r + 1)*(a + b*log(c*x**n))/(2*r + 1) + e**3*x**(3*r + 1)*(a + b*log(c*x**n))/(3*r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_378():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**2
    F = -b*d**3*n/x - 3*b*d**2*e*n*x**(r - 1)/(1 - r)**2 - 3*b*d*e**2*n*x**(2*r - 1)/(1 - 2*r)**2 - b*e**3*n*x**(3*r - 1)/(1 - 3*r)**2 - d**3*(a + b*log(c*x**n))/x - 3*d**2*e*x**(r - 1)*(a + b*log(c*x**n))/(1 - r) - 3*d*e**2*x**(2*r - 1)*(a + b*log(c*x**n))/(1 - 2*r) - e**3*x**(3*r - 1)*(a + b*log(c*x**n))/(1 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_379():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**4
    F = -b*d**3*n/(9*x**3) - 3*b*d**2*e*n*x**(r - 3)/(3 - r)**2 - 3*b*d*e**2*n*x**(2*r - 3)/(3 - 2*r)**2 - b*e**3*n*x**(3*r - 3)/(9*(1 - r)**2) - d**3*(a + b*log(c*x**n))/(3*x**3) - 3*d**2*e*x**(r - 3)*(a + b*log(c*x**n))/(3 - r) - 3*d*e**2*x**(2*r - 3)*(a + b*log(c*x**n))/(3 - 2*r) - e**3*x**(3*r - 3)*(a + b*log(c*x**n))/(3 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_380():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**6
    F = -b*d**3*n/(25*x**5) - 3*b*d**2*e*n*x**(r - 5)/(5 - r)**2 - 3*b*d*e**2*n*x**(2*r - 5)/(5 - 2*r)**2 - b*e**3*n*x**(3*r - 5)/(5 - 3*r)**2 - d**3*(a + b*log(c*x**n))/(5*x**5) - 3*d**2*e*x**(r - 5)*(a + b*log(c*x**n))/(5 - r) - 3*d*e**2*x**(2*r - 5)*(a + b*log(c*x**n))/(5 - 2*r) - e**3*x**(3*r - 5)*(a + b*log(c*x**n))/(5 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_381():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**8
    F = -b*d**3*n/(49*x**7) - 3*b*d**2*e*n*x**(r - 7)/(7 - r)**2 - 3*b*d*e**2*n*x**(2*r - 7)/(7 - 2*r)**2 - b*e**3*n*x**(3*r - 7)/(7 - 3*r)**2 - d**3*(a + b*log(c*x**n))/(7*x**7) - 3*d**2*e*x**(r - 7)*(a + b*log(c*x**n))/(7 - r) - 3*d*e**2*x**(2*r - 7)*(a + b*log(c*x**n))/(7 - 2*r) - e**3*x**(3*r - 7)*(a + b*log(c*x**n))/(7 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_382():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x**10
    F = -b*d**3*n/(81*x**9) - 3*b*d**2*e*n*x**(r - 9)/(9 - r)**2 - 3*b*d*e**2*n*x**(2*r - 9)/(9 - 2*r)**2 - b*e**3*n*x**(3*r - 9)/(9*(3 - r)**2) - d**3*(a + b*log(c*x**n))/(9*x**9) - 3*d**2*e*x**(r - 9)*(a + b*log(c*x**n))/(9 - r) - 3*d*e**2*x**(2*r - 9)*(a + b*log(c*x**n))/(9 - 2*r) - e**3*x**(3*r - 9)*(a + b*log(c*x**n))/(9 - 3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_383():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**r)
    F = sympy.Function('Unintegrable')((((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_384():
    f = x*(a + b*log(c*x**n))/(d + e*x**r)
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_385():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r))
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * ((Symbol('e') * (x)**(Symbol('r'))))**(Integer(-1)))))) * ((Symbol('d') * Symbol('r')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * ((Symbol('e') * (x)**(Symbol('r'))))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('r'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_386():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**r))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((x)**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_387():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**r)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_388():
    f = (a + b*log(c*x**n))/(d + e*x**r)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_389():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**r))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((x)**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_390():
    f = x**3*(a + b*log(c*x**n))/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')((((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_391():
    f = x*(a + b*log(c*x**n))/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_392():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**2)
    F = (Integer(-1) * ((Symbol('e') * (x)**(Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * Symbol('r') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('r')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_393():
    f = (a + b*log(c*x**n))/(x**3*(d + e*x**r)**2)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((x)**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_394():
    f = x**2*(a + b*log(c*x**n))/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_395():
    f = (a + b*log(c*x**n))/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_396():
    f = (a + b*log(c*x**n))/(x**2*(d + e*x**r)**2)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (((x)**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_397():
    f = (a + b*log(c*x**n))/(x*(c - 1/x**n))
    F = ((Symbol('a') * sympy.log((Integer(1) + (Integer(-1) * (Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('c') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('c') * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_398():
    f = (a + b*log(c*x**n))*(d + e*x**r)**3/x
    F = -b*d**3*n*log(x)**2/2 - 3*b*d**2*e*n*x**r/r**2 - 3*b*d*e**2*n*x**(2*r)/(4*r**2) - b*e**3*n*x**(3*r)/(9*r**2) + d**3*(a + b*log(c*x**n))*log(x) + 3*d**2*e*x**r*(a + b*log(c*x**n))/r + 3*d*e**2*x**(2*r)*(a + b*log(c*x**n))/(2*r) + e**3*x**(3*r)*(a + b*log(c*x**n))/(3*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_399():
    f = (a + b*log(c*x**n))*(d + e*x**r)**2/x
    F = -b*d**2*n*log(x)**2/2 - 2*b*d*e*n*x**r/r**2 - b*e**2*n*x**(2*r)/(4*r**2) + d**2*(a + b*log(c*x**n))*log(x) + 2*d*e*x**r*(a + b*log(c*x**n))/r + e**2*x**(2*r)*(a + b*log(c*x**n))/(2*r)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_400():
    f = (a + b*log(c*x**n))*(d + e*x**r)/x
    F = -b*e*n*x**r/r**2 + e*x**r*(a + b*log(c*x**n))/r + d*(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_401():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r))
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * Symbol('r')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('r'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_402():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**2)
    F = (Integer(-1) * ((Symbol('e') * (x)**(Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * Symbol('r') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('r')))**(Integer(-1)))) + ((Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_403():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**3)
    F = (Integer(-1) * ((Symbol('b') * Symbol('n')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * sympy.log(x)) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * Symbol('r')))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('d') * Symbol('r') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (x)**(Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * Symbol('r') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * Symbol('r')))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('n') * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * ((Integer(2) * (Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_404():
    f = (a + b*log(c*x**n))**2*(d + e*x**r)**3/x
    F = 6*b**2*d**2*e*n**2*x**r/r**3 + 3*b**2*d*e**2*n**2*x**(2*r)/(4*r**3) + 2*b**2*e**3*n**2*x**(3*r)/(27*r**3) - 6*b*d**2*e*n*x**r*(a + b*log(c*x**n))/r**2 - 3*b*d*e**2*n*x**(2*r)*(a + b*log(c*x**n))/(2*r**2) - 2*b*e**3*n*x**(3*r)*(a + b*log(c*x**n))/(9*r**2) + 3*d**2*e*x**r*(a + b*log(c*x**n))**2/r + 3*d*e**2*x**(2*r)*(a + b*log(c*x**n))**2/(2*r) + e**3*x**(3*r)*(a + b*log(c*x**n))**2/(3*r) + d**3*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_405():
    f = (a + b*log(c*x**n))**2*(d + e*x**r)**2/x
    F = 4*b**2*d*e*n**2*x**r/r**3 + b**2*e**2*n**2*x**(2*r)/(4*r**3) - 4*b*d*e*n*x**r*(a + b*log(c*x**n))/r**2 - b*e**2*n*x**(2*r)*(a + b*log(c*x**n))/(2*r**2) + 2*d*e*x**r*(a + b*log(c*x**n))**2/r + e**2*x**(2*r)*(a + b*log(c*x**n))**2/(2*r) + d**2*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_406():
    f = (a + b*log(c*x**n))**2*(d + e*x**r)/x
    F = 2*b**2*e*n**2*x**r/r**3 - 2*b*e*n*x**r*(a + b*log(c*x**n))/r**2 + e*x**r*(a + b*log(c*x**n))**2/r + d*(a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_407():
    f = (a + b*log(c*x**n))**2/(x*(d + e*x**r))
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * Symbol('r')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * ((Symbol('d') * (Symbol('r'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_408():
    f = (a + b*log(c*x**n))**2/(x*(d + e*x**r)**2)
    F = (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Symbol('d') * Symbol('r') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * Symbol('r')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_409():
    f = (a + b*log(c*x**n))**2/(x*(d + e*x**r)**3)
    F = ((Symbol('b') * Symbol('e') * Symbol('n') * (x)**(Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * ((Integer(2) * Symbol('d') * Symbol('r') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * (((Symbol('d'))**(Integer(2)) * Symbol('r') * (Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * Symbol('r')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('d') * (((x)**(Symbol('r')) * Symbol('e')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_410():
    f = (a + b*log(c*x**n))*(d + e*x**r)**(sympy.S(5)/2)/x
    F = (Integer(-1) * ((Integer(92) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * ((Integer(15) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * Symbol('b') * Symbol('d') * Symbol('n') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(45) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(25) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + ((Integer(92) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Integer(15))**(Integer(-1))) * (((Integer(15) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * (Symbol('r'))**(Integer(-1))) + ((Integer(5) * Symbol('d') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Symbol('r'))**(Integer(-1))) + ((Integer(3) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Symbol('r'))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (Symbol('r'))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_411():
    f = (a + b*log(c*x**n))*(d + e*x**r)**(sympy.S(3)/2)/x
    F = (Integer(-1) * ((Integer(16) * Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * ((Integer(3) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * (Symbol('r'))**(Integer(-1))) + (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (Symbol('r'))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_412():
    f = (a + b*log(c*x**n))*sqrt(d + e*x**r)/x
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))) + (Integer(2) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (Symbol('r'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (Symbol('r'))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('d')) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * ((Symbol('r'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_413():
    f = (a + b*log(c*x**n))/(x*sqrt(d + e*x**r))
    F = ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((sympy.sqrt(Symbol('d')) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('d')) * Symbol('r')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * ((sympy.sqrt(Symbol('d')) * (Symbol('r'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_414():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**(sympy.S(3)/2))
    F = ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + (Integer(2) * (((Symbol('d') * Symbol('r') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))))**(Integer(-1)) + (Integer(-1) * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('r')))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_415():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**(sympy.S(5)/2))
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n')) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))))**(Integer(-1)))) + ((Integer(16) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (((Symbol('d') * Symbol('r') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)) + (Integer(3) * (((Symbol('d'))**(Integer(2)) * Symbol('r') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('r')))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_416():
    f = (a + b*log(c*x**n))/(x*(d + e*x**r)**(sympy.S(7)/2))
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n')) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * (Symbol('r'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * Symbol('b') * Symbol('n')) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * (Symbol('r'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))))**(Integer(-1)))) + ((Integer(92) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('n') * (sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * (((Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Integer(15))**(Integer(-1))) * ((Integer(3) * ((Symbol('d') * Symbol('r') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(5) * (((Symbol('d'))**(Integer(2)) * Symbol('r') * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(15) * (((Symbol('d'))**(Integer(3)) * Symbol('r') * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * (((Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('r')))**(Integer(-1))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * sympy.atanh((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))) * (((Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('d')) + (Integer(-1) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Symbol('r'))))))))**(Integer(-1))))))) * (((Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('r'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_417():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**r)**3
    F = -b*d**3*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 3*b*d**2*e*n*x**(r + 1)*(f*x)**m/(m + r + 1)**2 - 3*b*d*e**2*n*x**(2*r + 1)*(f*x)**m/(m + 2*r + 1)**2 - b*e**3*n*x**(3*r + 1)*(f*x)**m/(m + 3*r + 1)**2 + d**3*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 3*d**2*e*x**(r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + r + 1) + 3*d*e**2*x**(2*r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + 2*r + 1) + e**3*x**(3*r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + 3*r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_418():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**r)**2
    F = -b*d**2*n*(f*x)**(m + 1)/(f*(m + 1)**2) - 2*b*d*e*n*x**(r + 1)*(f*x)**m/(m + r + 1)**2 - b*e**2*n*x**(2*r + 1)*(f*x)**m/(m + 2*r + 1)**2 + d**2*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + 2*d*e*x**(r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + r + 1) + e**2*x**(2*r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + 2*r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_419():
    f = (f*x)**m*(a + b*log(c*x**n))*(d + e*x**r)
    F = -b*d*n*(f*x)**(m + 1)/(f*(m + 1)**2) - b*e*n*x**(r + 1)*(f*x)**m/(m + r + 1)**2 + d*(f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1)) + e*x**(r + 1)*(f*x)**m*(a + b*log(c*x**n))/(m + r + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_420():
    f = (f*x)**m*(a + b*log(c*x**n))
    F = -b*n*(f*x)**(m + 1)/(f*(m + 1)**2) + (f*x)**(m + 1)*(a + b*log(c*x**n))/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_421():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x**r)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_422():
    f = (f*x)**m*(a + b*log(c*x**n))/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_423():
    f = (a + b*log(c*x**n))*(d + e/x**(1/(q + 1)))**q
    F = -b*n*x*(d + e/x**(1/(q + 1)))**q*hyper((-q - 1, -q - 1), (-q,), -e/(d*x**(1/(q + 1))))/(1 + e/(d*x**(1/(q + 1))))**q + x*(a + b*log(c*x**n))*(d + e/x**(1/(q + 1)))**(q + 1)/d
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_424():
    f = (f*x)**(-r*(q + 1) - 1)*(a + b*log(c*x**n))*(d + e*x**r)**q
    F = -b*n*(d + e*x**r)**q*hyper((-q - 1, -q - 1), (-q,), -e*x**r/d)/(f*r**2*(f*x)**(r*(q + 1))*(1 + e*x**r/d)**q*(q + 1)**2) - (a + b*log(c*x**n))*(d + e*x**r)**(q + 1)/(d*f*r*(f*x)**(r*(q + 1))*(q + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_425():
    f = (f*x)**m*(a + b*log(c*x**n))**p*(d + e*x**r)**3
    F = (((Symbol('d'))**(Integer(3)) * ((Symbol('f') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('f') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * Symbol('e') * (x)**((Integer(1) + Symbol('r'))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + Symbol('r'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + Symbol('r'))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * (Symbol('e'))**(Integer(2)) * (x)**((Integer(1) + (Integer(2) * Symbol('r')))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + (Integer(2) * Symbol('r')))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + (Integer(2) * Symbol('r')))))**(Integer(-1))) + (((Symbol('e'))**(Integer(3)) * (x)**((Integer(1) + (Integer(3) * Symbol('r')))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(3) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + (Integer(3) * Symbol('r')))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + (Integer(3) * Symbol('r'))) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(3) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + (Integer(3) * Symbol('r')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_426():
    f = (f*x)**m*(a + b*log(c*x**n))**p*(d + e*x**r)**2
    F = (((Symbol('d'))**(Integer(2)) * ((Symbol('f') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('f') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * Symbol('e') * (x)**((Integer(1) + Symbol('r'))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + Symbol('r'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + Symbol('r'))))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (x)**((Integer(1) + (Integer(2) * Symbol('r')))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + (Integer(2) * Symbol('r')))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(2) * Symbol('r'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + (Integer(2) * Symbol('r')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_427():
    f = (f*x)**m*(a + b*log(c*x**n))**p*(d + e*x**r)
    F = ((Symbol('d') * ((Symbol('f') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('f') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + ((Symbol('e') * (x)**((Integer(1) + Symbol('r'))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m') + Symbol('r'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + Symbol('r')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + Symbol('m') + Symbol('r'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_428():
    f = (f*x)**m*(a + b*log(c*x**n))**p
    F = (((Symbol('f') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('f') * (Integer(1) + Symbol('m')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_429():
    f = (f*x)**m*(a + b*log(c*x**n))**p/(d + e*x**r)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_430():
    f = (f*x)**m*(a + b*log(c*x**n))**p/(d + e*x**r)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('r')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_431():
    f = (a + b*log(c*x**n))*(f + g*x)/(d + e*x)**3
    F = b*n*(-d*g + e*f)/(2*d*e**2*(d + e*x)) + b*f**2*n*log(x)/(2*d**2*(-d*g + e*f)) - b*n*(d*g + e*f)*log(d + e*x)/(2*d**2*e**2) - (a + b*log(c*x**n))*(f + g*x)**2/((d + e*x)**2*(-2*d*g + 2*e*f))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_432():
    f = (a + b*log(c*x**n))**2*(f + g*x)/(d + e*x)**3
    F = (Integer(-1) * ((Symbol('b') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * Symbol('n') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (((Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('f'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('e') * x)))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_4_f_x_pow_m_d_plus_e_x_pow_r_pow_q_a_plus_b_log_c_x_pow_n_pow_p_433():
    f = (a + b*log(c*x**n))**3*(f + g*x)/(d + e*x)**3
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * Symbol('n') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('f'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(3))) * ((Integer(2) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * (Symbol('n'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * x) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F

