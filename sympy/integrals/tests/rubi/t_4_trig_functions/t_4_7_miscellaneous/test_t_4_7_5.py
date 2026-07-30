"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.7 Miscellaneous/4.7.5 x^m trig(a+b log(c x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, m, n, p = symbols('a b c d e m n p')

def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_1():
    f = x**2*sin(a + b*log(c*x**n))
    F = -b*n*x**3*cos(a + b*log(c*x**n))/(b**2*n**2 + 9) + 3*x**3*sin(a + b*log(c*x**n))/(b**2*n**2 + 9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_2():
    f = x*sin(a + b*log(c*x**n))
    F = -b*n*x**2*cos(a + b*log(c*x**n))/(b**2*n**2 + 4) + 2*x**2*sin(a + b*log(c*x**n))/(b**2*n**2 + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_3():
    f = sin(a + b*log(c*x**n))
    F = -b*n*x*cos(a + b*log(c*x**n))/(b**2*n**2 + 1) + x*sin(a + b*log(c*x**n))/(b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_4():
    f = sin(a + b*log(c*x**n))/x
    F = -cos(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_5():
    f = sin(a + b*log(c*x**n))/x**2
    F = -b*n*cos(a + b*log(c*x**n))/(x*(b**2*n**2 + 1)) - sin(a + b*log(c*x**n))/(x*(b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_6():
    f = sin(a + b*log(c*x**n))/x**3
    F = -b*n*cos(a + b*log(c*x**n))/(x**2*(b**2*n**2 + 4)) - 2*sin(a + b*log(c*x**n))/(x**2*(b**2*n**2 + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_7():
    f = x**2*sin(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x**3/(12*b**2*n**2 + 27) - 2*b*n*x**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(4*b**2*n**2 + 9) + 3*x**3*sin(a + b*log(c*x**n))**2/(4*b**2*n**2 + 9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_8():
    f = x*sin(a + b*log(c*x**n))**2
    F = b**2*n**2*x**2/(4*b**2*n**2 + 4) - b*n*x**2*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(2*b**2*n**2 + 2) + x**2*sin(a + b*log(c*x**n))**2/(2*b**2*n**2 + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_9():
    f = sin(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x/(4*b**2*n**2 + 1) - 2*b*n*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(4*b**2*n**2 + 1) + x*sin(a + b*log(c*x**n))**2/(4*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_10():
    f = sin(a + b*log(c*x**n))**2/x
    F = log(x)/2 - sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_11():
    f = sin(a + b*log(c*x**n))**2/x**2
    F = -2*b**2*n**2/(x*(4*b**2*n**2 + 1)) - 2*b*n*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(x*(4*b**2*n**2 + 1)) - sin(a + b*log(c*x**n))**2/(x*(4*b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_12():
    f = sin(a + b*log(c*x**n))**2/x**3
    F = -b**2*n**2/(x**2*(4*b**2*n**2 + 4)) - b*n*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(x**2*(2*b**2*n**2 + 2)) - sin(a + b*log(c*x**n))**2/(x**2*(2*b**2*n**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_13():
    f = x**2*sin(a + b*log(c*x**n))**3
    F = -2*b**3*n**3*x**3*cos(a + b*log(c*x**n))/(3*b**4*n**4 + 30*b**2*n**2 + 27) + 2*b**2*n**2*x**3*sin(a + b*log(c*x**n))/(b**4*n**4 + 10*b**2*n**2 + 9) - b*n*x**3*sin(a + b*log(c*x**n))**2*cos(a + b*log(c*x**n))/(3*b**2*n**2 + 3) + x**3*sin(a + b*log(c*x**n))**3/(3*b**2*n**2 + 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_14():
    f = x*sin(a + b*log(c*x**n))**3
    F = -6*b**3*n**3*x**2*cos(a + b*log(c*x**n))/(9*b**4*n**4 + 40*b**2*n**2 + 16) + 12*b**2*n**2*x**2*sin(a + b*log(c*x**n))/(9*b**4*n**4 + 40*b**2*n**2 + 16) - 3*b*n*x**2*sin(a + b*log(c*x**n))**2*cos(a + b*log(c*x**n))/(9*b**2*n**2 + 4) + 2*x**2*sin(a + b*log(c*x**n))**3/(9*b**2*n**2 + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_15():
    f = sin(a + b*log(c*x**n))**3
    F = -6*b**3*n**3*x*cos(a + b*log(c*x**n))/(9*b**4*n**4 + 10*b**2*n**2 + 1) + 6*b**2*n**2*x*sin(a + b*log(c*x**n))/(9*b**4*n**4 + 10*b**2*n**2 + 1) - 3*b*n*x*sin(a + b*log(c*x**n))**2*cos(a + b*log(c*x**n))/(9*b**2*n**2 + 1) + x*sin(a + b*log(c*x**n))**3/(9*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_16():
    f = sin(a + b*log(c*x**n))**3/x
    F = cos(a + b*log(c*x**n))**3/(3*b*n) - cos(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_17():
    f = sin(a + b*log(c*x**n))**3/x**2
    F = -6*b**3*n**3*cos(a + b*log(c*x**n))/(x*(9*b**4*n**4 + 10*b**2*n**2 + 1)) - 6*b**2*n**2*sin(a + b*log(c*x**n))/(x*(9*b**4*n**4 + 10*b**2*n**2 + 1)) - 3*b*n*sin(a + b*log(c*x**n))**2*cos(a + b*log(c*x**n))/(x*(9*b**2*n**2 + 1)) - sin(a + b*log(c*x**n))**3/(x*(9*b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_18():
    f = sin(a + b*log(c*x**n))**3/x**3
    F = -6*b**3*n**3*cos(a + b*log(c*x**n))/(x**2*(9*b**4*n**4 + 40*b**2*n**2 + 16)) - 12*b**2*n**2*sin(a + b*log(c*x**n))/(x**2*(9*b**4*n**4 + 40*b**2*n**2 + 16)) - 3*b*n*sin(a + b*log(c*x**n))**2*cos(a + b*log(c*x**n))/(x**2*(9*b**2*n**2 + 4)) - 2*sin(a + b*log(c*x**n))**3/(x**2*(9*b**2*n**2 + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_19():
    f = x**2*sin(a + b*log(c*x**n))**4
    F = 8*b**4*n**4*x**3/(64*b**4*n**4 + 180*b**2*n**2 + 81) - 24*b**3*n**3*x**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(64*b**4*n**4 + 180*b**2*n**2 + 81) + 36*b**2*n**2*x**3*sin(a + b*log(c*x**n))**2/(64*b**4*n**4 + 180*b**2*n**2 + 81) - 4*b*n*x**3*sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(16*b**2*n**2 + 9) + 3*x**3*sin(a + b*log(c*x**n))**4/(16*b**2*n**2 + 9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_20():
    f = x*sin(a + b*log(c*x**n))**4
    F = 3*b**4*n**4*x**2/(16*b**4*n**4 + 20*b**2*n**2 + 4) - 3*b**3*n**3*x**2*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(8*b**4*n**4 + 10*b**2*n**2 + 2) + 3*b**2*n**2*x**2*sin(a + b*log(c*x**n))**2/(8*b**4*n**4 + 10*b**2*n**2 + 2) - b*n*x**2*sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(4*b**2*n**2 + 1) + x**2*sin(a + b*log(c*x**n))**4/(8*b**2*n**2 + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_21():
    f = sin(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x/(64*b**4*n**4 + 20*b**2*n**2 + 1) - 24*b**3*n**3*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(64*b**4*n**4 + 20*b**2*n**2 + 1) + 12*b**2*n**2*x*sin(a + b*log(c*x**n))**2/(64*b**4*n**4 + 20*b**2*n**2 + 1) - 4*b*n*x*sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(16*b**2*n**2 + 1) + x*sin(a + b*log(c*x**n))**4/(16*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_22():
    f = sin(a + b*log(c*x**n))**4/x
    F = 3*log(x)/8 - sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(4*b*n) - 3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_23():
    f = sin(a + b*log(c*x**n))**4/x**2
    F = -24*b**4*n**4/(x*(64*b**4*n**4 + 20*b**2*n**2 + 1)) - 24*b**3*n**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(x*(64*b**4*n**4 + 20*b**2*n**2 + 1)) - 12*b**2*n**2*sin(a + b*log(c*x**n))**2/(x*(64*b**4*n**4 + 20*b**2*n**2 + 1)) - 4*b*n*sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(x*(16*b**2*n**2 + 1)) - sin(a + b*log(c*x**n))**4/(x*(16*b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_24():
    f = sin(a + b*log(c*x**n))**4/x**3
    F = -3*b**4*n**4/(x**2*(16*b**4*n**4 + 20*b**2*n**2 + 4)) - 3*b**3*n**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(x**2*(8*b**4*n**4 + 10*b**2*n**2 + 2)) - 3*b**2*n**2*sin(a + b*log(c*x**n))**2/(x**2*(8*b**4*n**4 + 10*b**2*n**2 + 2)) - b*n*sin(a + b*log(c*x**n))**3*cos(a + b*log(c*x**n))/(x**2*(4*b**2*n**2 + 1)) - sin(a + b*log(c*x**n))**4/(x**2*(8*b**2*n**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_25():
    f = sin(log(a + b*x))
    F = (a + b*x)*sin(log(a + b*x))/(2*b) - (a + b*x)*cos(log(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_26():
    f = x**m*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n))
    F = -x**(m + 1)*(c*x**n)**((m + 1)/n)*exp(a*(m + 1)/(n*sqrt(-(m + 1)**2/n**2)))/(4*n*sqrt(-(m + 1)**2/n**2)) + x**(m + 1)*(m + 1)*exp(a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))*log(x)/(2*n*(c*x**n)**((m + 1)/n)*sqrt(-(m + 1)**2/n**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_27():
    f = x**2*sin(a + 3*sqrt(-1/n**2)*log(c*x**n))
    F = n*x**3*sqrt(-1/n**2)*(c*x**n)**(3/n)*exp(-a*n*sqrt(-1/n**2))/12 - n*x**3*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))*log(x)/(2*(c*x**n)**(3/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_28():
    f = x*sin(a + 2*sqrt(-1/n**2)*log(c*x**n))
    F = n*x**2*sqrt(-1/n**2)*(c*x**n)**(2/n)*exp(-a*n*sqrt(-1/n**2))/8 - n*x**2*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))*log(x)/(2*(c*x**n)**(2/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_29():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n))
    F = n*x*sqrt(-1/n**2)*(c*x**n)**(1/n)*exp(-a*n*sqrt(-1/n**2))/4 - n*x*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))*log(x)/(2*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_30():
    f = sin(a)/x
    F = log(x)*sin(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_31():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n))/x**2
    F = n*sqrt(-1/n**2)*(c*x**n)**(1/n)*exp(-a*n*sqrt(-1/n**2))*log(x)/(2*x) + n*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(4*x*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_32():
    f = sin(a + 2*sqrt(-1/n**2)*log(c*x**n))/x**3
    F = n*sqrt(-1/n**2)*(c*x**n)**(2/n)*exp(-a*n*sqrt(-1/n**2))*log(x)/(2*x**2) + n*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(8*x**2*(c*x**n)**(2/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_33():
    f = x**m*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**2
    F = -x**(m + 1)*(c*x**n)**((m + 1)/n)*exp(-2*a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))/(8*m + 8) + x**(m + 1)/(2*m + 2) - x**(m + 1)*exp(2*a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))*log(x)/(4*(c*x**n)**((m + 1)/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_34():
    f = x**2*sin(a + 3*sqrt(-1/n**2)*log(c*x**n)/2)**2
    F = -x**3*(c*x**n)**(3/n)*exp(-2*a*n*sqrt(-1/n**2))/24 + x**3/6 - x**3*exp(2*a*n*sqrt(-1/n**2))*log(x)/(4*(c*x**n)**(3/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_35():
    f = x*sin(a + sqrt(-1/n**2)*log(c*x**n))**2
    F = -x**2*(c*x**n)**(2/n)*exp(-2*a*n*sqrt(-1/n**2))/16 + x**2/4 - x**2*exp(2*a*n*sqrt(-1/n**2))*log(x)/(4*(c*x**n)**(2/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_36():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n)/2)**2
    F = -x*(c*x**n)**(1/n)*exp(-2*a*n*sqrt(-1/n**2))/8 + x/2 - x*exp(2*a*n*sqrt(-1/n**2))*log(x)/(4*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_37():
    f = sin(a)**2/x
    F = log(x)*sin(a)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_38():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n)/2)**2/x**2
    F = -(c*x**n)**(1/n)*exp(-2*a*n*sqrt(-1/n**2))*log(x)/(4*x) - 1/(2*x) + exp(2*a*n*sqrt(-1/n**2))/(8*x*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_39():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n))**2/x**3
    F = -(c*x**n)**(2/n)*exp(-2*a*n*sqrt(-1/n**2))*log(x)/(4*x**2) - 1/(4*x**2) + exp(2*a*n*sqrt(-1/n**2))/(16*x**2*(c*x**n)**(2/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_40():
    f = x**m*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**3
    F = 6*n*x**(m + 1)*sqrt(-(m + 1)**2/n**2)*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**2*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)/(5*(m + 1)**2) - 4*n*x**(m + 1)*sqrt(-(m + 1)**2/n**2)*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)/(5*(m + 1)**2) - 4*x**(m + 1)*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**3/(5*m + 5) + 8*x**(m + 1)*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)/(5*m + 5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_41():
    f = x**2*sin(a + sqrt(-1/n**2)*log(c*x**n))**3
    F = -n*x**3*sqrt(-1/n**2)*(c*x**n)**(3/n)*exp(-3*a*n*sqrt(-1/n**2))/48 + 3*n*x**3*sqrt(-1/n**2)*(c*x**n)**(1/n)*exp(-a*n*sqrt(-1/n**2))/32 - 3*n*x**3*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(16*(c*x**n)**(1/n)) + n*x**3*sqrt(-1/n**2)*exp(3*a*n*sqrt(-1/n**2))*log(x)/(8*(c*x**n)**(3/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_42():
    f = x*sin(a + 2*sqrt(-1/n**2)*log(c*x**n)/3)**3
    F = 9*n*x**2*sqrt(-1/n**2)*(c*x**n)**(2/(3*n))*exp(-a*n*sqrt(-1/n**2))/64 - n*x**2*sqrt(-1/n**2)*(c*x**n)**(2/n)*exp(-3*a*n*sqrt(-1/n**2))/32 + n*x**2*sqrt(-1/n**2)*exp(3*a*n*sqrt(-1/n**2))*log(x)/(8*(c*x**n)**(2/n)) - 9*n*x**2*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(32*(c*x**n)**(2/(3*n)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_43():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n)/3)**3
    F = 9*n*x*sqrt(-1/n**2)*(c*x**n)**(1/(3*n))*exp(-a*n*sqrt(-1/n**2))/32 - n*x*sqrt(-1/n**2)*(c*x**n)**(1/n)*exp(-3*a*n*sqrt(-1/n**2))/16 + n*x*sqrt(-1/n**2)*exp(3*a*n*sqrt(-1/n**2))*log(x)/(8*(c*x**n)**(1/n)) - 9*n*x*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(16*(c*x**n)**(1/(3*n)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_44():
    f = sin(a)**3/x
    F = log(x)*sin(a)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_45():
    f = sin(a + sqrt(-1/n**2)*log(c*x**n)/3)**3/x**2
    F = -9*n*sqrt(-1/n**2)*(c*x**n)**(1/(3*n))*exp(-a*n*sqrt(-1/n**2))/(16*x) - n*sqrt(-1/n**2)*(c*x**n)**(1/n)*exp(-3*a*n*sqrt(-1/n**2))*log(x)/(8*x) - n*sqrt(-1/n**2)*exp(3*a*n*sqrt(-1/n**2))/(16*x*(c*x**n)**(1/n)) + 9*n*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(32*x*(c*x**n)**(1/(3*n)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_46():
    f = sin(a + 2*sqrt(-1/n**2)*log(c*x**n)/3)**3/x**3
    F = -9*n*sqrt(-1/n**2)*(c*x**n)**(2/(3*n))*exp(-a*n*sqrt(-1/n**2))/(32*x**2) - n*sqrt(-1/n**2)*(c*x**n)**(2/n)*exp(-3*a*n*sqrt(-1/n**2))*log(x)/(8*x**2) - n*sqrt(-1/n**2)*exp(3*a*n*sqrt(-1/n**2))/(32*x**2*(c*x**n)**(2/n)) + 9*n*sqrt(-1/n**2)*exp(a*n*sqrt(-1/n**2))/(64*x**2*(c*x**n)**(2/(3*n)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_47():
    f = x**m*sin(a + sqrt(-(m + 1)**2)*log(c*x**2)/2)
    F = x**(m + 1)*(c*x**2)**(-m/2 + sympy.S(-1)/2)*(m + 1)*exp(a*sqrt(-(m + 1)**2)/(m + 1))*log(x)/(2*sqrt(-(m + 1)**2)) - x**(m + 1)*(c*x**2)**(m/2 + sympy.S.Half)*exp(a*(m + 1)/sqrt(-(m + 1)**2))/(4*sqrt(-(m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_48():
    f = sin(a + I*log(c*x**2)/2)
    F = I*c*x**3*exp(-I*a)/(4*sqrt(c*x**2)) - I*x*exp(I*a)*log(x)/(2*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_49():
    f = x**m*sin(a + sqrt(-(m + 1)**2)*log(c*x**2)/4)**2
    F = -x**(m + 1)*(c*x**2)**(-m/2 + sympy.S(-1)/2)*exp(-2*a*(m + 1)/sqrt(-(m + 1)**2))*log(x)/4 - x**(m + 1)*(c*x**2)**(m/2 + sympy.S.Half)*exp(2*a*(m + 1)/sqrt(-(m + 1)**2))/(8*m + 8) + x**(m + 1)/(2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_50():
    f = sin(a + I*log(c*x**2)/4)**2
    F = -c*x**3*exp(-2*I*a)/(8*sqrt(c*x**2)) + x/2 - x*exp(2*I*a)*log(x)/(4*sqrt(c*x**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_51():
    f = x**m*sin(a + sqrt(-(m + 1)**2)*log(c*x**2)/6)**3
    F = -x**(m + 1)*(c*x**2)**(-m/2 + sympy.S(-1)/2)*(m + 1)*exp(-3*a*(m + 1)/sqrt(-(m + 1)**2))*log(x)/(8*sqrt(-(m + 1)**2)) + 9*x**(m + 1)*(c*x**2)**(-m/6 + sympy.S(-1)/6)*exp(a*sqrt(-(m + 1)**2)/(m + 1))/(16*sqrt(-(m + 1)**2)) - 9*x**(m + 1)*(c*x**2)**(m/6 + sympy.S(1)/6)*exp(a*(m + 1)/sqrt(-(m + 1)**2))/(32*sqrt(-(m + 1)**2)) + x**(m + 1)*(c*x**2)**(m/2 + sympy.S.Half)*exp(3*a*(m + 1)/sqrt(-(m + 1)**2))/(16*sqrt(-(m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_52():
    f = sin(a + I*log(c*x**2)/6)**3
    F = -I*c*x**3*exp(-3*I*a)/(16*sqrt(c*x**2)) + 9*I*x*(c*x**2)**(sympy.S(1)/6)*exp(-I*a)/32 + I*x*exp(3*I*a)*log(x)/(8*sqrt(c*x**2)) - 9*I*x*exp(I*a)/(16*(c*x**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_53():
    f = x*sqrt(sin(a + b*log(c*x**n)))
    F = 2*x**2*sqrt(sin(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, sympy.S(-1)/4 - I/(b*n)), (sympy.S(3)/4 - I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_54():
    f = sqrt(sin(a + b*log(c*x**n)))
    F = 2*x*sqrt(sin(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, -(b*n + 2*I)/(4*b*n)), (sympy.S(3)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_55():
    f = sqrt(sin(a + b*log(c*x**n)))/x
    F = 2*elliptic_e(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_56():
    f = sqrt(sin(a + b*log(c*x**n)))/x**2
    F = -2*sqrt(sin(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, sympy.S(-1)/4 + I/(2*b*n)), (sympy.S(3)/4 + I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x*sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_57():
    f = sqrt(sin(a + b*log(c*x**n)))/x**3
    F = -2*sqrt(sin(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, sympy.S(-1)/4 + I/(b*n)), (sympy.S(3)/4 + I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(I*b*n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_58():
    f = x*sin(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**2*sin(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, sympy.S(-3)/4 - I/(b*n)), (sympy.S(1)/4 - I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_59():
    f = sin(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x*sin(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, sympy.S(-3)/4 - I/(2*b*n)), (sympy.S(1)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_60():
    f = sin(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -2*sqrt(sin(a + b*log(c*x**n)))*cos(a + b*log(c*x**n))/(3*b*n) + 2*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_61():
    f = sin(a + b*log(c*x**n))**(sympy.S(3)/2)/x**2
    F = -2*sin(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, sympy.S(-3)/4 + I/(2*b*n)), (sympy.S(1)/4 + I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(3*I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_62():
    f = sin(a + b*log(c*x**n))**(sympy.S(3)/2)/x**3
    F = -2*sin(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, sympy.S(-3)/4 + I/(b*n)), (sympy.S(1)/4 + I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(3*I*b*n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_63():
    f = 1/sqrt(sin(a + b*log(c*x**n)))
    F = 2*x*sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*hyper((sympy.S.Half, sympy.S(1)/4 - I/(2*b*n)), (sympy.S(5)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((I*b*n + 2)*sqrt(sin(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_64():
    f = 1/(x*sqrt(sin(a + b*log(c*x**n))))
    F = 2*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_65():
    f = sin(a + b*log(c*x**n))**(sympy.S(-3)/2)
    F = 2*x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, sympy.S(3)/4 - I/(2*b*n)), (sympy.S(7)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((3*I*b*n + 2)*sin(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_66():
    f = 1/(x*sin(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = -2*elliptic_e(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n) - 2*cos(a + b*log(c*x**n))/(b*n*sqrt(sin(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_67():
    f = sin(a + b*log(c*x**n))**(sympy.S(-5)/2)
    F = 2*x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, sympy.S(5)/4 - I/(2*b*n)), (sympy.S(9)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((5*I*b*n + 2)*sin(a + b*log(c*x**n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_68():
    f = 1/(x*sin(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(3*b*n) - 2*cos(a + b*log(c*x**n))/(3*b*n*sin(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_69():
    f = sin(a - 2*I*log(c*x))**(sympy.S(-3)/2)
    F = (-c**4*x**4*exp(2*I*a) + 1)*exp(-2*I*a)/(2*c**4*x**3*sin(a - 2*I*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_70():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))**4
    F = 24*b**4*d**4*n**4*(e*x)**(m + 1)/(e*(m + 1)*(4*b**2*d**2*n**2 + (m + 1)**2)*(16*b**2*d**2*n**2 + (m + 1)**2)) - 24*b**3*d**3*n**3*(e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))*cos(d*(a + b*log(c*x**n)))/(e*(4*b**2*d**2*n**2 + (m + 1)**2)*(16*b**2*d**2*n**2 + (m + 1)**2)) + 12*b**2*d**2*n**2*(e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))**2/(e*(4*b**2*d**2*n**2 + (m + 1)**2)*(16*b**2*d**2*n**2 + (m + 1)**2)) - 4*b*d*n*(e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))**3*cos(d*(a + b*log(c*x**n)))/(e*(16*b**2*d**2*n**2 + (m + 1)**2)) + (e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))**4/(e*(16*b**2*d**2*n**2 + (m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_71():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))**3
    F = -6*b**3*d**3*n**3*(e*x)**(m + 1)*cos(d*(a + b*log(c*x**n)))/(e*(b**2*d**2*n**2 + (m + 1)**2)*(9*b**2*d**2*n**2 + (m + 1)**2)) + 6*b**2*d**2*n**2*(e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))/(e*(b**2*d**2*n**2 + (m + 1)**2)*(9*b**2*d**2*n**2 + (m + 1)**2)) - 3*b*d*n*(e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))**2*cos(d*(a + b*log(c*x**n)))/(e*(9*b**2*d**2*n**2 + (m + 1)**2)) + (e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))**3/(e*(9*b**2*d**2*n**2 + (m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_72():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))**2
    F = 2*b**2*d**2*n**2*(e*x)**(m + 1)/(e*(m + 1)*(4*b**2*d**2*n**2 + (m + 1)**2)) - 2*b*d*n*(e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))*cos(d*(a + b*log(c*x**n)))/(e*(4*b**2*d**2*n**2 + (m + 1)**2)) + (e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))**2/(e*(4*b**2*d**2*n**2 + (m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_73():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))
    F = -b*d*n*(e*x)**(m + 1)*cos(d*(a + b*log(c*x**n)))/(e*(b**2*d**2*n**2 + (m + 1)**2)) + (e*x)**(m + 1)*(m + 1)*sin(d*(a + b*log(c*x**n)))/(e*(b**2*d**2*n**2 + (m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_74():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))**(sympy.S(3)/2)
    F = 2*(e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, -(3*b*d*n + 2*I*m + 2*I)/(4*b*d*n)), (-(-b*d*n + 2*I*m + 2*I)/(4*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**(sympy.S(3)/2)*(-3*I*b*d*n + 2*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_75():
    f = (e*x)**m*sqrt(sin(d*(a + b*log(c*x**n))))
    F = 2*(e*x)**(m + 1)*sqrt(sin(d*(a + b*log(c*x**n))))*hyper((sympy.S(-1)/2, -(b*d*n + 2*I*m + 2*I)/(4*b*d*n)), (-(-3*b*d*n + 2*I*m + 2*I)/(4*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*sqrt(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)*(-I*b*d*n + 2*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_76():
    f = (e*x)**m/sqrt(sin(d*(a + b*log(c*x**n))))
    F = 2*(e*x)**(m + 1)*sqrt(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)*hyper((sympy.S.Half, -(-b*d*n + 2*I*m + 2*I)/(4*b*d*n)), (-(-5*b*d*n + 2*I*m + 2*I)/(4*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(I*b*d*n + 2*m + 2)*sqrt(sin(d*(a + b*log(c*x**n)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_77():
    f = (e*x)**m/sin(d*(a + b*log(c*x**n)))**(sympy.S(3)/2)
    F = 2*(e*x)**(m + 1)*(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, -(-3*b*d*n + 2*I*m + 2*I)/(4*b*d*n)), (-(-7*b*d*n + 2*I*m + 2*I)/(4*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(3*I*b*d*n + 2*m + 2)*sin(d*(a + b*log(c*x**n)))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_78():
    f = (e*x)**m/sin(d*(a + b*log(c*x**n)))**(sympy.S(5)/2)
    F = 2*(e*x)**(m + 1)*(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, -(-5*b*d*n + 2*I*m + 2*I)/(4*b*d*n)), (-(-9*b*d*n + 2*I*m + 2*I)/(4*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(5*I*b*d*n + 2*m + 2)*sin(d*(a + b*log(c*x**n)))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_79():
    f = (e*x)**m*sin(d*(a + b*log(c*x**n)))**p
    F = (e*x)**(m + 1)*sin(d*(a + b*log(c*x**n)))**p*hyper((-p, -(b*d*n*p + I*m + I)/(2*b*d*n)), (-p/2 + 1 - I*(m + 1)/(2*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**p*(-I*b*d*n*p + m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_80():
    f = x**2*sin(a + b*log(c*x**n))**p
    F = x**3*sin(a + b*log(c*x**n))**p*hyper((-p, -(b*n*p + 3*I)/(2*b*n)), (-p/2 + 1 - 3*I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(-I*b*n*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_81():
    f = x*sin(a + b*log(c*x**n))**p
    F = x**2*sin(a + b*log(c*x**n))**p*hyper((-p, -p/2 - I/(b*n)), (-p/2 + 1 - I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(-I*b*n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_82():
    f = sin(a + b*log(c*x**n))**p
    F = x*sin(a + b*log(c*x**n))**p*hyper((-p, -(b*n*p + I)/(2*b*n)), (-p/2 + 1 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(-I*b*n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_83():
    f = sin(a + b*log(c*x**n))**p/x
    F = sin(a + b*log(c*x**n))**(p + 1)*cos(a + b*log(c*x**n))*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sin(a + b*log(c*x**n))**2)/(b*n*(p + 1)*sqrt(cos(a + b*log(c*x**n))**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_84():
    f = sin(a + b*log(c*x**n))**p/x**2
    F = -sin(a + b*log(c*x**n))**p*hyper((-p, -p/2 + I/(2*b*n)), (-p/2 + 1 + I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(I*b*n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_85():
    f = sin(a + b*log(c*x**n))**p/x**3
    F = -sin(a + b*log(c*x**n))**p*hyper((-p, -p/2 + I/(b*n)), (-p/2 + 1 + I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(I*b*n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_86():
    f = x**2*cos(a + b*log(c*x**n))
    F = b*n*x**3*sin(a + b*log(c*x**n))/(b**2*n**2 + 9) + 3*x**3*cos(a + b*log(c*x**n))/(b**2*n**2 + 9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_87():
    f = x*cos(a + b*log(c*x**n))
    F = b*n*x**2*sin(a + b*log(c*x**n))/(b**2*n**2 + 4) + 2*x**2*cos(a + b*log(c*x**n))/(b**2*n**2 + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_88():
    f = cos(a + b*log(c*x**n))
    F = b*n*x*sin(a + b*log(c*x**n))/(b**2*n**2 + 1) + x*cos(a + b*log(c*x**n))/(b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_89():
    f = cos(a + b*log(c*x**n))/x
    F = sin(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_90():
    f = cos(a + b*log(c*x**n))/x**2
    F = b*n*sin(a + b*log(c*x**n))/(x*(b**2*n**2 + 1)) - cos(a + b*log(c*x**n))/(x*(b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_91():
    f = x**2*cos(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x**3/(12*b**2*n**2 + 27) + 2*b*n*x**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(4*b**2*n**2 + 9) + 3*x**3*cos(a + b*log(c*x**n))**2/(4*b**2*n**2 + 9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_92():
    f = x*cos(a + b*log(c*x**n))**2
    F = b**2*n**2*x**2/(4*b**2*n**2 + 4) + b*n*x**2*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(2*b**2*n**2 + 2) + x**2*cos(a + b*log(c*x**n))**2/(2*b**2*n**2 + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_93():
    f = cos(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x/(4*b**2*n**2 + 1) + 2*b*n*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(4*b**2*n**2 + 1) + x*cos(a + b*log(c*x**n))**2/(4*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_94():
    f = cos(a + b*log(c*x**n))**2/x
    F = log(x)/2 + sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_95():
    f = cos(a + b*log(c*x**n))**2/x**2
    F = -2*b**2*n**2/(x*(4*b**2*n**2 + 1)) + 2*b*n*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(x*(4*b**2*n**2 + 1)) - cos(a + b*log(c*x**n))**2/(x*(4*b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_96():
    f = x**2*cos(a + b*log(c*x**n))**3
    F = 2*b**3*n**3*x**3*sin(a + b*log(c*x**n))/(3*b**4*n**4 + 30*b**2*n**2 + 27) + 2*b**2*n**2*x**3*cos(a + b*log(c*x**n))/(b**4*n**4 + 10*b**2*n**2 + 9) + b*n*x**3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**2/(3*b**2*n**2 + 3) + x**3*cos(a + b*log(c*x**n))**3/(3*b**2*n**2 + 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_97():
    f = x*cos(a + b*log(c*x**n))**3
    F = 6*b**3*n**3*x**2*sin(a + b*log(c*x**n))/(9*b**4*n**4 + 40*b**2*n**2 + 16) + 12*b**2*n**2*x**2*cos(a + b*log(c*x**n))/(9*b**4*n**4 + 40*b**2*n**2 + 16) + 3*b*n*x**2*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**2/(9*b**2*n**2 + 4) + 2*x**2*cos(a + b*log(c*x**n))**3/(9*b**2*n**2 + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_98():
    f = cos(a + b*log(c*x**n))**3
    F = 6*b**3*n**3*x*sin(a + b*log(c*x**n))/(9*b**4*n**4 + 10*b**2*n**2 + 1) + 6*b**2*n**2*x*cos(a + b*log(c*x**n))/(9*b**4*n**4 + 10*b**2*n**2 + 1) + 3*b*n*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**2/(9*b**2*n**2 + 1) + x*cos(a + b*log(c*x**n))**3/(9*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_99():
    f = cos(a + b*log(c*x**n))**3/x
    F = -sin(a + b*log(c*x**n))**3/(3*b*n) + sin(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_100():
    f = cos(a + b*log(c*x**n))**3/x**2
    F = 6*b**3*n**3*sin(a + b*log(c*x**n))/(x*(9*b**4*n**4 + 10*b**2*n**2 + 1)) - 6*b**2*n**2*cos(a + b*log(c*x**n))/(x*(9*b**4*n**4 + 10*b**2*n**2 + 1)) + 3*b*n*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**2/(x*(9*b**2*n**2 + 1)) - cos(a + b*log(c*x**n))**3/(x*(9*b**2*n**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_101():
    f = cos(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x/(64*b**4*n**4 + 20*b**2*n**2 + 1) + 24*b**3*n**3*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(64*b**4*n**4 + 20*b**2*n**2 + 1) + 12*b**2*n**2*x*cos(a + b*log(c*x**n))**2/(64*b**4*n**4 + 20*b**2*n**2 + 1) + 4*b*n*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**3/(16*b**2*n**2 + 1) + x*cos(a + b*log(c*x**n))**4/(16*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_102():
    f = cos(a + b*log(c*x**n))**4/x
    F = 3*log(x)/8 + sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**3/(4*b*n) + 3*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_103():
    f = cos(log(3*x + 6))
    F = (x/2 + 1)*sin(log(3*x + 6)) + (x/2 + 1)*cos(log(3*x + 6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_104():
    f = x**m*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n))
    F = x**(m + 1)*(c*x**n)**((m + 1)/n)*exp(a*(m + 1)/(n*sqrt(-(m + 1)**2/n**2)))/(4*m + 4) + x**(m + 1)*exp(a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))*log(x)/(2*(c*x**n)**((m + 1)/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_105():
    f = cos(a + sqrt(-1/n**2)*log(c*x**n))
    F = x*(c*x**n)**(1/n)*exp(-a*n*sqrt(-1/n**2))/4 + x*exp(a*n*sqrt(-1/n**2))*log(x)/(2*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_106():
    f = x**m*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**2
    F = x**(m + 1)*(c*x**n)**((m + 1)/n)*exp(-2*a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))/(8*m + 8) + x**(m + 1)/(2*m + 2) + x**(m + 1)*exp(2*a*n*sqrt(-(m + 1)**2/n**2)/(m + 1))*log(x)/(4*(c*x**n)**((m + 1)/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_107():
    f = cos(a + sqrt(-1/n**2)*log(c*x**n)/2)**2
    F = x*(c*x**n)**(1/n)*exp(-2*a*n*sqrt(-1/n**2))/8 + x/2 + x*exp(2*a*n*sqrt(-1/n**2))*log(x)/(4*(c*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_108():
    f = x**m*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**3
    F = -6*n*x**(m + 1)*sqrt(-(m + 1)**2/n**2)*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**2/(5*(m + 1)**2) + 4*n*x**(m + 1)*sqrt(-(m + 1)**2/n**2)*sin(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)/(5*(m + 1)**2) - 4*x**(m + 1)*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)**3/(5*m + 5) + 8*x**(m + 1)*cos(a + sqrt(-(m + 1)**2/n**2)*log(c*x**n)/2)/(5*m + 5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_109():
    f = cos(a + sqrt(-1/n**2)*log(c*x**n)/3)**3
    F = 9*x*(c*x**n)**(1/(3*n))*exp(-a*n*sqrt(-1/n**2))/32 + x*(c*x**n)**(1/n)*exp(-3*a*n*sqrt(-1/n**2))/16 + x*exp(3*a*n*sqrt(-1/n**2))*log(x)/(8*(c*x**n)**(1/n)) + 9*x*exp(a*n*sqrt(-1/n**2))/(16*(c*x**n)**(1/(3*n)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_110():
    f = sqrt(cos(a + b*log(c*x**n)))
    F = 2*x*sqrt(cos(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, -(b*n + 2*I)/(4*b*n)), (sympy.S(3)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_111():
    f = sqrt(cos(a + b*log(c*x**n)))/x
    F = 2*elliptic_e(a/2 + b*log(c*x**n)/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_112():
    f = cos(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x*cos(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, sympy.S(-3)/4 - I/(2*b*n)), (sympy.S(1)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_113():
    f = cos(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = 2*sin(a + b*log(c*x**n))*sqrt(cos(a + b*log(c*x**n)))/(3*b*n) + 2*elliptic_f(a/2 + b*log(c*x**n)/2, 2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_114():
    f = cos(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x*cos(a + b*log(c*x**n))**(sympy.S(5)/2)*hyper((sympy.S(-5)/2, sympy.S(-5)/4 - I/(2*b*n)), (-(b*n + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*(-5*I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_115():
    f = cos(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**(sympy.S(3)/2)/(5*b*n) + 6*elliptic_e(a/2 + b*log(c*x**n)/2, 2)/(5*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_116():
    f = 1/sqrt(cos(a + b*log(c*x**n)))
    F = 2*x*sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*hyper((sympy.S.Half, sympy.S(1)/4 - I/(2*b*n)), (sympy.S(5)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((I*b*n + 2)*sqrt(cos(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_117():
    f = 1/(x*sqrt(cos(a + b*log(c*x**n))))
    F = 2*elliptic_f(a/2 + b*log(c*x**n)/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_118():
    f = cos(a + b*log(c*x**n))**(sympy.S(-3)/2)
    F = 2*x*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, sympy.S(3)/4 - I/(2*b*n)), (sympy.S(7)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((3*I*b*n + 2)*cos(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_119():
    f = 1/(x*cos(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*sin(a + b*log(c*x**n))/(b*n*sqrt(cos(a + b*log(c*x**n)))) - 2*elliptic_e(a/2 + b*log(c*x**n)/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_120():
    f = cos(a + b*log(c*x**n))**(sympy.S(-5)/2)
    F = 2*x*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, sympy.S(5)/4 - I/(2*b*n)), (sympy.S(9)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((5*I*b*n + 2)*cos(a + b*log(c*x**n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_121():
    f = 1/(x*cos(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*sin(a + b*log(c*x**n))/(3*b*n*cos(a + b*log(c*x**n))**(sympy.S(3)/2)) + 2*elliptic_f(a/2 + b*log(c*x**n)/2, 2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_122():
    f = cos(a - 2*I*log(c*x))**(sympy.S(-3)/2)
    F = -(c**4*x**4*exp(2*I*a) + 1)*exp(-2*I*a)/(2*c**4*x**3*cos(a - 2*I*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_123():
    f = x**m*cos(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x**(m + 1)/((m + 1)*(4*b**2*n**2 + (m + 1)**2)*(16*b**2*n**2 + (m + 1)**2)) + 24*b**3*n**3*x**(m + 1)*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/((4*b**2*n**2 + (m + 1)**2)*(16*b**2*n**2 + (m + 1)**2)) + 12*b**2*n**2*x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))**2/((4*b**2*n**2 + (m + 1)**2)*(16*b**2*n**2 + (m + 1)**2)) + 4*b*n*x**(m + 1)*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**3/(16*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))**4/(16*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_124():
    f = x**m*cos(a + b*log(c*x**n))**3
    F = 6*b**3*n**3*x**(m + 1)*sin(a + b*log(c*x**n))/((b**2*n**2 + (m + 1)**2)*(9*b**2*n**2 + (m + 1)**2)) + 6*b**2*n**2*x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))/((b**2*n**2 + (m + 1)**2)*(9*b**2*n**2 + (m + 1)**2)) + 3*b*n*x**(m + 1)*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**2/(9*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))**3/(9*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_125():
    f = x**m*cos(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x**(m + 1)/((m + 1)*(4*b**2*n**2 + (m + 1)**2)) + 2*b*n*x**(m + 1)*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))/(4*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))**2/(4*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_126():
    f = x**m*cos(a + b*log(c*x**n))
    F = b*n*x**(m + 1)*sin(a + b*log(c*x**n))/(b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cos(a + b*log(c*x**n))/(b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_127():
    f = x**m*cos(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*cos(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(-3)/2, -(3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_128():
    f = x**m*sqrt(cos(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*sqrt(cos(a + b*log(c*x**n)))*hyper((sympy.S(-1)/2, -(b*n + 2*I*m + 2*I)/(4*b*n)), (-(-3*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_129():
    f = x**m/sqrt(cos(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*hyper((sympy.S.Half, -(-b*n + 2*I*m + 2*I)/(4*b*n)), (-(-5*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((I*b*n + 2*m + 2)*sqrt(cos(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_130():
    f = x**m/cos(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, -(-3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-7*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((3*I*b*n + 2*m + 2)*cos(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_131():
    f = x**m/cos(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x**(m + 1)*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, -(-5*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-9*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/((5*I*b*n + 2*m + 2)*cos(a + b*log(c*x**n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_132():
    f = (e*x)**m*cos(d*(a + b*log(c*x**n)))**p
    F = (e*x)**(m + 1)*cos(d*(a + b*log(c*x**n)))**p*hyper((-p, -(b*d*n*p + I*m + I)/(2*b*d*n)), (-p/2 + 1 - I*(m + 1)/(2*b*d*n),), -(c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*((c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**p*(-I*b*d*n*p + m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_133():
    f = x*cos(a + b*log(c*x**n))**p
    F = x**2*cos(a + b*log(c*x**n))**p*hyper((-p, -p/2 - I/(b*n)), (-p/2 + 1 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(-I*b*n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_134():
    f = cos(a + b*log(c*x**n))**p
    F = x*cos(a + b*log(c*x**n))**p*hyper((-p, -(b*n*p + I)/(2*b*n)), (-p/2 + 1 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*(-I*b*n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_135():
    f = tan(a + b*log(c*x**n))/x
    F = -log(cos(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_136():
    f = tan(a + b*log(c*x**n))**2/x
    F = -log(x) + tan(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_137():
    f = tan(a + b*log(c*x**n))**3/x
    F = log(cos(a + b*log(c*x**n)))/(b*n) + tan(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_138():
    f = tan(a + b*log(c*x**n))**4/x
    F = log(x) + tan(a + b*log(c*x**n))**3/(3*b*n) - tan(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_139():
    f = tan(a + b*log(c*x**n))**5/x
    F = -log(cos(a + b*log(c*x**n)))/(b*n) + tan(a + b*log(c*x**n))**4/(4*b*n) - tan(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_140():
    f = tan(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + 2*tan(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_141():
    f = tan(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + 2*sqrt(tan(a + b*log(c*x**n)))/(b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_142():
    f = sqrt(tan(a + b*log(c*x**n)))/x
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_143():
    f = 1/(x*sqrt(tan(a + b*log(c*x**n))))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_144():
    f = 1/(x*tan(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n) - 2/(b*n*sqrt(tan(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_145():
    f = 1/(x*tan(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + tan(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(tan(a + b*log(c*x**n))) + 1)/(2*b*n) - 2/(3*b*n*tan(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_146():
    f = cot(a + b*log(c*x**n))/x
    F = log(sin(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_147():
    f = cot(a + b*log(c*x**n))**2/x
    F = -log(x) - cot(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_148():
    f = cot(a + b*log(c*x**n))**3/x
    F = -log(sin(a + b*log(c*x**n)))/(b*n) - cot(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_149():
    f = cot(a + b*log(c*x**n))**4/x
    F = log(x) - cot(a + b*log(c*x**n))**3/(3*b*n) + cot(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_150():
    f = cot(a + b*log(c*x**n))**5/x
    F = log(sin(a + b*log(c*x**n)))/(b*n) - cot(a + b*log(c*x**n))**4/(4*b*n) + cot(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_151():
    f = cot(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - 2*cot(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_152():
    f = cot(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - 2*sqrt(cot(a + b*log(c*x**n)))/(b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_153():
    f = sqrt(cot(a + b*log(c*x**n)))/x
    F = -sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_154():
    f = 1/(x*sqrt(cot(a + b*log(c*x**n))))
    F = sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) - sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_155():
    f = 1/(x*cot(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) - sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n) + 2/(b*n*sqrt(cot(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_156():
    f = 1/(x*cot(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*log(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + cot(a + b*log(c*x**n)) + 1)/(4*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) - 1)/(2*b*n) + sqrt(2)*atan(sqrt(2)*sqrt(cot(a + b*log(c*x**n))) + 1)/(2*b*n) + 2/(3*b*n*cot(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_157():
    f = x**2*sec(a + b*log(c*x**n))
    F = 2*x**3*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - 3*I/(2*b*n)), (sympy.S(3)/2 - 3*I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_158():
    f = x*sec(a + b*log(c*x**n))
    F = 2*x**2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - I/(b*n)), (sympy.S(3)/2 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_159():
    f = sec(a + b*log(c*x**n))
    F = 2*x*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - I/(2*b*n)), (sympy.S(3)/2 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_160():
    f = sec(a + b*log(c*x**n))/x
    F = atanh(sin(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_161():
    f = sec(a + b*log(c*x**n))/x**2
    F = -2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half + I/(2*b*n)), (sympy.S(3)/2 + I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_162():
    f = sec(a + b*log(c*x**n))/x**3
    F = -2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half + I/(b*n)), (sympy.S(3)/2 + I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_163():
    f = x**2*sec(a + b*log(c*x**n))**2
    F = 4*x**3*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 - 3*I/(2*b*n)), (2 - 3*I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(2*I*b*n + 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_164():
    f = x*sec(a + b*log(c*x**n))**2
    F = 2*x**2*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 - I/(b*n)), (2 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_165():
    f = sec(a + b*log(c*x**n))**2
    F = 4*x*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 - I/(2*b*n)), (2 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(2*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_166():
    f = sec(a + b*log(c*x**n))**2/x
    F = tan(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_167():
    f = sec(a + b*log(c*x**n))**2/x**2
    F = -4*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 + I/(2*b*n)), (2 + I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-2*I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_168():
    f = sec(a + b*log(c*x**n))**2/x**3
    F = -2*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 + I/(b*n)), (2 + I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_169():
    f = x*sec(a + b*log(c*x**n))**3
    F = 8*x**2*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, sympy.S(3)/2 - I/(b*n)), (sympy.S(5)/2 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(3*I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_170():
    f = sec(a + b*log(c*x**n))**3
    F = 8*x*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, sympy.S(3)/2 - I/(2*b*n)), (sympy.S(5)/2 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(3*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_171():
    f = sec(a + b*log(c*x**n))**3/x
    F = tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n))/(2*b*n) + atanh(sin(a + b*log(c*x**n)))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_172():
    f = sec(a + b*log(c*x**n))**3/x**2
    F = -8*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, sympy.S(3)/2 + I/(2*b*n)), (sympy.S(5)/2 + I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-3*I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_173():
    f = sec(a + b*log(c*x**n))**3/x**3
    F = -8*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, sympy.S(3)/2 + I/(b*n)), (sympy.S(5)/2 + I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-3*I*b*n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_174():
    f = x*sec(a + b*log(c*x**n))**4
    F = 8*x**2*(c*x**n)**(4*I*b)*exp(4*I*a)*hyper((4, 2 - I/(b*n)), (3 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(2*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_175():
    f = sec(a + b*log(c*x**n))**4
    F = 16*x*(c*x**n)**(4*I*b)*exp(4*I*a)*hyper((4, 2 - I/(2*b*n)), (3 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(4*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_176():
    f = sec(a + b*log(c*x**n))**4/x
    F = tan(a + b*log(c*x**n))**3/(3*b*n) + tan(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_177():
    f = sec(a + b*log(c*x**n))**4/x**2
    F = -16*(c*x**n)**(4*I*b)*exp(4*I*a)*hyper((4, 2 + I/(2*b*n)), (3 + I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x*(-4*I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_178():
    f = sec(a + b*log(c*x**n))**4/x**3
    F = -8*(c*x**n)**(4*I*b)*exp(4*I*a)*hyper((4, 2 + I/(b*n)), (3 + I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(-2*I*b*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_179():
    f = 2*b**2*n**2*sec(a + b*log(c*x**n))**3 - (b**2*n**2 + 1)*sec(a + b*log(c*x**n))
    F = b*n*x*tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n)) - x*sec(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_180():
    f = x**m*sec(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))**3
    F = x**(m + 1)*sec(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))/(2*m + 2) + x**(m + 1)*tan(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))*sec(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))/(2*sqrt(-(m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_181():
    f = x*sec(a + 2*log(c*x**I))**3
    F = x**2*(c*x**I)**(2*I)*exp(I*a)/((c*x**I)**(4*I)*exp(2*I*a) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_182():
    f = sec(a + 2*log(c/x**(I/2)))**3
    F = 2*x*(c/x**(I/2))**(6*I)*exp(3*I*a)/((c/x**(I/2))**(4*I)*exp(2*I*a) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_183():
    f = sec(a + I*log(c*x**n)/(n*(p - 2)))**p
    F = x*(2 - p)*((c*x**n)**(2/(n*(2 - p)))*exp(2*I*a) + 1)*exp(-2*I*a)*sec(a - I*log(c*x**n)/(n*(2 - p)))**p/((c*x**n)**(2/(n*(2 - p)))*(2 - 2*p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_184():
    f = sec(a - I*log(c*x**n)/(n*(p - 2)))**p
    F = x*(1 + exp(2*I*a)/(c*x**n)**(2/(n*(2 - p))))*(2 - p)*sec(a + I*log(c*x**n)/(n*(2 - p)))**p/(2 - 2*p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_185():
    f = sqrt(sec(a + b*log(c*x**n)))
    F = 2*x*sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*hyper((sympy.S.Half, sympy.S(1)/4 - I/(2*b*n)), (sympy.S(5)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sqrt(sec(a + b*log(c*x**n)))/(I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_186():
    f = sqrt(sec(a + b*log(c*x**n)))/x
    F = 2*sqrt(cos(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_187():
    f = sec(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, sympy.S(3)/4 - I/(2*b*n)), (sympy.S(7)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_188():
    f = sec(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = 2*sin(a + b*log(c*x**n))*sqrt(sec(a + b*log(c*x**n)))/(b*n) - 2*sqrt(cos(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_189():
    f = sec(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, sympy.S(5)/4 - I/(2*b*n)), (sympy.S(9)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**(sympy.S(5)/2)/(5*I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_190():
    f = sec(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sin(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) + 2*sqrt(cos(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_191():
    f = 1/sqrt(sec(a + b*log(c*x**n)))
    F = 2*x*hyper((sympy.S(-1)/2, -(b*n + 2*I)/(4*b*n)), (sympy.S(3)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2)*sqrt(sec(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_192():
    f = 1/(x*sqrt(sec(a + b*log(c*x**n))))
    F = 2*sqrt(cos(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_193():
    f = sec(a + b*log(c*x**n))**(sympy.S(-3)/2)
    F = 2*x*hyper((sympy.S(-3)/2, sympy.S(-3)/4 - I/(2*b*n)), (sympy.S(1)/4 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2)*sec(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_194():
    f = 1/(x*sec(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*sin(a + b*log(c*x**n))/(3*b*n*sqrt(sec(a + b*log(c*x**n)))) + 2*sqrt(cos(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_195():
    f = sec(a + b*log(c*x**n))**(sympy.S(-5)/2)
    F = 2*x*hyper((sympy.S(-5)/2, sympy.S(-5)/4 - I/(2*b*n)), (-(b*n + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*(-5*I*b*n + 2)*sec(a + b*log(c*x**n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_196():
    f = 1/(x*sec(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*sin(a + b*log(c*x**n))/(5*b*n*sec(a + b*log(c*x**n))**(sympy.S(3)/2)) + 6*sqrt(cos(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2, 2)*sqrt(sec(a + b*log(c*x**n)))/(5*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_197():
    f = x**m*sec(a + b*log(c*x**n))**3
    F = 8*x**(m + 1)*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, -(-3*b*n + I*(m + 1))/(2*b*n)), (-(-5*b*n + I*(m + 1))/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(3*I*b*n + m + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_198():
    f = x**m*sec(a + b*log(c*x**n))**2
    F = 4*x**(m + 1)*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, -(-2*b*n + I*(m + 1))/(2*b*n)), (-(-4*b*n + I*(m + 1))/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(2*I*b*n + m + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_199():
    f = x**m*sec(a + b*log(c*x**n))
    F = 2*x**(m + 1)*(c*x**n)**(I*b)*exp(I*a)*hyper((1, -(-b*n + I*m + I)/(2*b*n)), (-(-3*b*n + I*(m + 1))/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + m + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_200():
    f = x**m*sec(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x**(m + 1)*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, -(-5*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-9*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**(sympy.S(5)/2)/(5*I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_201():
    f = x**m*sec(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, -(-3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-7*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_202():
    f = x**m*sqrt(sec(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*hyper((sympy.S.Half, -(-b*n + 2*I*m + 2*I)/(4*b*n)), (-(-5*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sqrt(sec(a + b*log(c*x**n)))/(I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_203():
    f = x**m/sqrt(sec(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*hyper((sympy.S(-1)/2, -(b*n + 2*I*m + 2*I)/(4*b*n)), (-(-3*b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt((c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2*m + 2)*sqrt(sec(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_204():
    f = x**m/sec(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*hyper((sympy.S(-3)/2, -(3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-b*n + 2*I*m + 2*I)/(4*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))/(((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2*m + 2)*sec(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_205():
    f = (e*x)**m*sec(d*(a + b*log(c*x**n)))**p
    F = (e*x)**(m + 1)*((c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**p*hyper((p, -(-b*d*n*p + I*m + I)/(2*b*d*n)), (p/2 + 1 - I*(m + 1)/(2*b*d*n),), -(c*x**n)**(2*I*b*d)*exp(2*I*a*d))*sec(d*(a + b*log(c*x**n)))**p/(e*(I*b*d*n*p + m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_206():
    f = x*sec(a + b*log(c*x**n))**p
    F = x**2*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*hyper((p, p/2 - I/(b*n)), (p/2 + 1 - I/(b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**p/(I*b*n*p + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_207():
    f = sec(a + b*log(c*x**n))**p
    F = x*((c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*hyper((p, -(-b*n*p + I)/(2*b*n)), (p/2 + 1 - I/(2*b*n),), -(c*x**n)**(2*I*b)*exp(2*I*a))*sec(a + b*log(c*x**n))**p/(I*b*n*p + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_208():
    f = x**2*csc(a + b*log(c*x**n))
    F = 2*x**3*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - 3*I/(2*b*n)), (sympy.S(3)/2 - 3*I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(-b*n + 3*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_209():
    f = x*csc(a + b*log(c*x**n))
    F = 2*x**2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - I/(b*n)), (sympy.S(3)/2 - I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(-b*n + 2*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_210():
    f = csc(a + b*log(c*x**n))
    F = 2*x*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half - I/(2*b*n)), (sympy.S(3)/2 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(-b*n + I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_211():
    f = csc(a + b*log(c*x**n))/x
    F = -atanh(cos(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_212():
    f = csc(a + b*log(c*x**n))/x**2
    F = -2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half + I/(2*b*n)), (sympy.S(3)/2 + I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x*(b*n + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_213():
    f = csc(a + b*log(c*x**n))/x**3
    F = -2*(c*x**n)**(I*b)*exp(I*a)*hyper((1, sympy.S.Half + I/(b*n)), (sympy.S(3)/2 + I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(x**2*(b*n + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_214():
    f = csc(a + b*log(c*x**n))**2
    F = -4*x*(c*x**n)**(2*I*b)*exp(2*I*a)*hyper((2, 1 - I/(2*b*n)), (2 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(2*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_215():
    f = csc(a + b*log(c*x**n))**2/x
    F = -cot(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_216():
    f = csc(a + b*log(c*x**n))**3
    F = -8*x*(c*x**n)**(3*I*b)*exp(3*I*a)*hyper((3, sympy.S(3)/2 - I/(2*b*n)), (sympy.S(5)/2 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(-3*b*n + I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_217():
    f = csc(a + b*log(c*x**n))**3/x
    F = -cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n))/(2*b*n) - atanh(cos(a + b*log(c*x**n)))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_218():
    f = csc(a + b*log(c*x**n))**4
    F = 16*x*(c*x**n)**(4*I*b)*exp(4*I*a)*hyper((4, 2 - I/(2*b*n)), (3 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(4*I*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_219():
    f = csc(a + b*log(c*x**n))**4/x
    F = -cot(a + b*log(c*x**n))**3/(3*b*n) - cot(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_220():
    f = 2*b**2*n**2*csc(a + b*log(c*x**n))**3 - (b**2*n**2 + 1)*csc(a + b*log(c*x**n))
    F = -b*n*x*cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n)) - x*csc(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_221():
    f = x**m*csc(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))**3
    F = x**(m + 1)*csc(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))/(2*m + 2) - x**(m + 1)*cot(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))*csc(a + 2*log(c*x**(sqrt(-(m + 1)**2)/2)))/(2*sqrt(-(m + 1)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_222():
    f = x*csc(a + 2*log(c*x**I))**3
    F = -I*x**2*(c*x**I)**(2*I)*exp(I*a)/(-(c*x**I)**(4*I)*exp(2*I*a) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_223():
    f = csc(a + 2*log(c/x**(I/2)))**3
    F = 2*I*x*(c/x**(I/2))**(6*I)*exp(3*I*a)/(-(c/x**(I/2))**(4*I)*exp(2*I*a) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_224():
    f = csc(a + I*log(c*x**n)/(n*(p - 2)))**p
    F = -x*(2 - p)*(-(c*x**n)**(2/(n*(2 - p)))*exp(2*I*a) + 1)*exp(-2*I*a)*csc(a - I*log(c*x**n)/(n*(2 - p)))**p/((c*x**n)**(2/(n*(2 - p)))*(2 - 2*p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_225():
    f = csc(a - I*log(c*x**n)/(n*(p - 2)))**p
    F = x*(1 - exp(2*I*a)/(c*x**n)**(2/(n*(2 - p))))*(2 - p)*csc(a + I*log(c*x**n)/(n*(2 - p)))**p/(2 - 2*p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_226():
    f = sqrt(csc(a + b*log(c*x**n)))
    F = 2*x*sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*sqrt(csc(a + b*log(c*x**n)))*hyper((sympy.S.Half, sympy.S(1)/4 - I/(2*b*n)), (sympy.S(5)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_227():
    f = sqrt(csc(a + b*log(c*x**n)))/x
    F = 2*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_228():
    f = csc(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*csc(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(3)/2, sympy.S(3)/4 - I/(2*b*n)), (sympy.S(7)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(3*I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_229():
    f = csc(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -2*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n) - 2*cos(a + b*log(c*x**n))*sqrt(csc(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_230():
    f = csc(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*csc(a + b*log(c*x**n))**(sympy.S(5)/2)*hyper((sympy.S(5)/2, sympy.S(5)/4 - I/(2*b*n)), (sympy.S(9)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(5*I*b*n + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_231():
    f = csc(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(3*b*n) - 2*cos(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_232():
    f = 1/sqrt(csc(a + b*log(c*x**n)))
    F = 2*x*hyper((sympy.S(-1)/2, -(b*n + 2*I)/(4*b*n)), (sympy.S(3)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2)*sqrt(csc(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_233():
    f = 1/(x*sqrt(csc(a + b*log(c*x**n))))
    F = 2*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_234():
    f = csc(a + b*log(c*x**n))**(sympy.S(-3)/2)
    F = 2*x*hyper((sympy.S(-3)/2, sympy.S(-3)/4 - I/(2*b*n)), (sympy.S(1)/4 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2)*csc(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_235():
    f = 1/(x*csc(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_f(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(3*b*n) - 2*cos(a + b*log(c*x**n))/(3*b*n*sqrt(csc(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_236():
    f = csc(a + b*log(c*x**n))**(sympy.S(-5)/2)
    F = 2*x*hyper((sympy.S(-5)/2, sympy.S(-5)/4 - I/(2*b*n)), (-(b*n + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*(-5*I*b*n + 2)*csc(a + b*log(c*x**n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_237():
    f = 1/(x*csc(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 6*sqrt(sin(a + b*log(c*x**n)))*sqrt(csc(a + b*log(c*x**n)))*elliptic_e(a/2 + b*log(c*x**n)/2 - pi/4, 2)/(5*b*n) - 2*cos(a + b*log(c*x**n))/(5*b*n*csc(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_238():
    f = (e*x)**m*csc(d*(a + b*log(c*x**n)))**3
    F = -8*(c*x**n)**(3*I*b*d)*(e*x)**(m + 1)*exp(3*I*a*d)*hyper((3, -(-3*b*d*n + I*(m + 1))/(2*b*d*n)), (-(-5*b*d*n + I*(m + 1))/(2*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(-3*b*d*n + I*(m + 1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_239():
    f = (e*x)**m*csc(d*(a + b*log(c*x**n)))**2
    F = -4*(c*x**n)**(2*I*b*d)*(e*x)**(m + 1)*exp(2*I*a*d)*hyper((2, -(-2*b*d*n + I*(m + 1))/(2*b*d*n)), (-(-4*b*d*n + I*(m + 1))/(2*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(2*I*b*d*n + m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_240():
    f = (e*x)**m*csc(d*(a + b*log(c*x**n)))
    F = 2*(c*x**n)**(I*b*d)*(e*x)**(m + 1)*exp(I*a*d)*hyper((1, -(-b*d*n + I*m + I)/(2*b*d*n)), (-(-3*b*d*n + I*(m + 1))/(2*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(-b*d*n + I*(m + 1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_241():
    f = x**m*csc(a + b*log(c*x**n))**(sympy.S(5)/2)
    F = 2*x**(m + 1)*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(5)/2)*csc(a + b*log(c*x**n))**(sympy.S(5)/2)*hyper((sympy.S(5)/2, -(-5*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-9*b*n + 2*I*m + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(5*I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_242():
    f = x**m*csc(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*csc(a + b*log(c*x**n))**(sympy.S(3)/2)*hyper((sympy.S(3)/2, -(-3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-7*b*n + 2*I*m + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(3*I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_243():
    f = x**m*sqrt(csc(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*sqrt(csc(a + b*log(c*x**n)))*hyper((sympy.S.Half, -(-b*n + 2*I*m + 2*I)/(4*b*n)), (-(-5*b*n + 2*I*m + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n + 2*m + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_244():
    f = x**m/sqrt(csc(a + b*log(c*x**n)))
    F = 2*x**(m + 1)*hyper((sympy.S(-1)/2, -(b*n + 2*I*m + 2*I)/(4*b*n)), (-(-3*b*n + 2*I*m + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(sqrt(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)*(-I*b*n + 2*m + 2)*sqrt(csc(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_245():
    f = x**m/csc(a + b*log(c*x**n))**(sympy.S(3)/2)
    F = 2*x**(m + 1)*hyper((sympy.S(-3)/2, -(3*b*n + 2*I*m + 2*I)/(4*b*n)), (-(-b*n + 2*I*m + 2*I)/(4*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/((-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**(sympy.S(3)/2)*(-3*I*b*n + 2*m + 2)*csc(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_246():
    f = (e*x)**m*csc(d*(a + b*log(c*x**n)))**p
    F = (e*x)**(m + 1)*(-(c*x**n)**(2*I*b*d)*exp(2*I*a*d) + 1)**p*csc(d*(a + b*log(c*x**n)))**p*hyper((p, -(-b*d*n*p + I*m + I)/(2*b*d*n)), (p/2 + 1 - I*(m + 1)/(2*b*d*n),), (c*x**n)**(2*I*b*d)*exp(2*I*a*d))/(e*(I*b*d*n*p + m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_247():
    f = x*csc(a + b*log(c*x**n))**p
    F = x**2*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*csc(a + b*log(c*x**n))**p*hyper((p, p/2 - I/(b*n)), (p/2 + 1 - I/(b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n*p + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_5_x_pow_m_trig_a_plus_b_log_c_x_pow_n_pow_p_248():
    f = csc(a + b*log(c*x**n))**p
    F = x*(-(c*x**n)**(2*I*b)*exp(2*I*a) + 1)**p*csc(a + b*log(c*x**n))**p*hyper((p, -(-b*n*p + I)/(2*b*n)), (p/2 + 1 - I/(2*b*n),), (c*x**n)**(2*I*b)*exp(2*I*a))/(I*b*n*p + 1)
    assert integrate(f, x) == F

