"""Generated from MathematicaSyntaxTestSuite.

Source: 3 Logarithms/3.1.2 (d x)^m (a+b log(c x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, d1, d2, m, m1, m2, n, p, q, q1, q2 = symbols('a b c d d1 d2 m m1 m2 n p q q1 q2')

def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_1():
    f = x**3*log(c*x)
    F = x**4*log(c*x)/4 - x**4/16
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_2():
    f = x**2*log(c*x)
    F = x**3*log(c*x)/3 - x**3/9
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_3():
    f = x*log(c*x)
    F = x**2*log(c*x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_4():
    f = log(c*x)
    F = x*log(c*x) - x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_5():
    f = log(c*x)/x
    F = log(c*x)**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_6():
    f = log(c*x)/x**2
    F = -log(c*x)/x - 1/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_7():
    f = log(c*x)/x**3
    F = -log(c*x)/(2*x**2) - 1/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_8():
    f = x**3*log(c*x)**2
    F = x**4*log(c*x)**2/4 - x**4*log(c*x)/8 + x**4/32
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_9():
    f = x**2*log(c*x)**2
    F = x**3*log(c*x)**2/3 - 2*x**3*log(c*x)/9 + 2*x**3/27
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_10():
    f = x*log(c*x)**2
    F = x**2*log(c*x)**2/2 - x**2*log(c*x)/2 + x**2/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_11():
    f = log(c*x)**2
    F = x*log(c*x)**2 - 2*x*log(c*x) + 2*x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_12():
    f = log(c*x)**2/x
    F = log(c*x)**3/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_13():
    f = log(c*x)**2/x**2
    F = -log(c*x)**2/x - 2*log(c*x)/x - 2/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_14():
    f = log(c*x)**2/x**3
    F = -log(c*x)**2/(2*x**2) - log(c*x)/(2*x**2) - 1/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_15():
    f = x**3*log(c*x)**3
    F = x**4*log(c*x)**3/4 - 3*x**4*log(c*x)**2/16 + 3*x**4*log(c*x)/32 - 3*x**4/128
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_16():
    f = x**2*log(c*x)**3
    F = x**3*log(c*x)**3/3 - x**3*log(c*x)**2/3 + 2*x**3*log(c*x)/9 - 2*x**3/27
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_17():
    f = x*log(c*x)**3
    F = x**2*log(c*x)**3/2 - 3*x**2*log(c*x)**2/4 + 3*x**2*log(c*x)/4 - 3*x**2/8
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_18():
    f = log(c*x)**3
    F = x*log(c*x)**3 - 3*x*log(c*x)**2 + 6*x*log(c*x) - 6*x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_19():
    f = log(c*x)**3/x
    F = log(c*x)**4/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_20():
    f = log(c*x)**3/x**2
    F = -log(c*x)**3/x - 3*log(c*x)**2/x - 6*log(c*x)/x - 6/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_21():
    f = log(c*x)**3/x**3
    F = -log(c*x)**3/(2*x**2) - 3*log(c*x)**2/(4*x**2) - 3*log(c*x)/(4*x**2) - 3/(8*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_22():
    f = x**3/log(c*x)
    F = sympy.Function('ExpIntegralEi')((Integer(4) * sympy.log((Symbol('c') * x)))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_23():
    f = x**2/log(c*x)
    F = sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('c') * x)))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_24():
    f = x/log(c*x)
    F = sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('c') * x)))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_25():
    f = 1/log(c*x)
    F = sympy.Function('LogIntegral')((Symbol('c') * x)) * (Symbol('c'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_26():
    f = 1/(x*log(c*x))
    F = log(log(c*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_27():
    f = 1/(x**2*log(c*x))
    F = Symbol('c') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log((Symbol('c') * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_28():
    f = 1/(x**3*log(c*x))
    F = (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.log((Symbol('c') * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_29():
    f = x**3/log(c*x)**2
    F = ((Integer(4) * sympy.Function('ExpIntegralEi')((Integer(4) * sympy.log((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * (sympy.log((Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_30():
    f = x**2/log(c*x)**2
    F = ((Integer(3) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * (sympy.log((Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_31():
    f = x/log(c*x)**2
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * (sympy.log((Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_32():
    f = log(c*x)**(-2)
    F = (Integer(-1) * (x * (sympy.log((Symbol('c') * x)))**(Integer(-1)))) + (sympy.Function('LogIntegral')((Symbol('c') * x)) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_33():
    f = 1/(x*log(c*x)**2)
    F = -1/log(c*x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_34():
    f = 1/(x**2*log(c*x)**2)
    F = (Integer(-1) * (Symbol('c') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log((Symbol('c') * x)))))) + (Integer(-1) * ((x * sympy.log((Symbol('c') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_35():
    f = 1/(x**3*log(c*x)**2)
    F = (Integer(-2) * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.log((Symbol('c') * x))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Symbol('c') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_36():
    f = x**3/log(c*x)**3
    F = ((Integer(8) * sympy.Function('ExpIntegralEi')((Integer(4) * sympy.log((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * ((Integer(2) * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * (sympy.log((Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_37():
    f = x**2/log(c*x)**3
    F = ((Integer(9) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(2) * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * sympy.log((Symbol('c') * x))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_38():
    f = x/log(c*x)**3
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(2)) * (sympy.log((Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_39():
    f = log(c*x)**(-3)
    F = (Integer(-1) * (x * ((Integer(2) * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (x * ((Integer(2) * sympy.log((Symbol('c') * x))))**(Integer(-1)))) + (sympy.Function('LogIntegral')((Symbol('c') * x)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_40():
    f = 1/(x*log(c*x)**3)
    F = -1/(2*log(c*x)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_41():
    f = 1/(x**2*log(c*x)**3)
    F = ((Integer(2))**(Integer(-1)) * Symbol('c') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log((Symbol('c') * x))))) + (Integer(-1) * ((Integer(2) * x * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.log((Symbol('c') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_42():
    f = 1/(x**3*log(c*x)**3)
    F = (Integer(2) * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.log((Symbol('c') * x))))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * (sympy.log((Symbol('c') * x)))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.log((Symbol('c') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_43():
    f = x**3*(a + b*log(c*x**n))
    F = -b*n*x**4/16 + x**4*(a + b*log(c*x**n))/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_44():
    f = x**2*(a + b*log(c*x**n))
    F = -b*n*x**3/9 + x**3*(a + b*log(c*x**n))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_45():
    f = x*(a + b*log(c*x**n))
    F = -b*n*x**2/4 + x**2*(a + b*log(c*x**n))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_46():
    f = a + b*log(c*x**n)
    F = a*x - b*n*x + b*x*log(c*x**n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_47():
    f = (a + b*log(c*x**n))/x
    F = (a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_48():
    f = (a + b*log(c*x**n))/x**2
    F = -b*n/x - (a + b*log(c*x**n))/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_49():
    f = (a + b*log(c*x**n))/x**3
    F = -b*n/(4*x**2) - (a + b*log(c*x**n))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_50():
    f = x**3*(a + b*log(c*x**n))**2
    F = b**2*n**2*x**4/32 - b*n*x**4*(a + b*log(c*x**n))/8 + x**4*(a + b*log(c*x**n))**2/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_51():
    f = x**2*(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x**3/27 - 2*b*n*x**3*(a + b*log(c*x**n))/9 + x**3*(a + b*log(c*x**n))**2/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_52():
    f = x*(a + b*log(c*x**n))**2
    F = b**2*n**2*x**2/4 - b*n*x**2*(a + b*log(c*x**n))/2 + x**2*(a + b*log(c*x**n))**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_53():
    f = (a + b*log(c*x**n))**2
    F = -2*a*b*n*x + 2*b**2*n**2*x - 2*b**2*n*x*log(c*x**n) + x*(a + b*log(c*x**n))**2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_54():
    f = (a + b*log(c*x**n))**2/x
    F = (a + b*log(c*x**n))**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_55():
    f = (a + b*log(c*x**n))**2/x**2
    F = -2*b**2*n**2/x - 2*b*n*(a + b*log(c*x**n))/x - (a + b*log(c*x**n))**2/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_56():
    f = (a + b*log(c*x**n))**2/x**3
    F = -b**2*n**2/(4*x**2) - b*n*(a + b*log(c*x**n))/(2*x**2) - (a + b*log(c*x**n))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_57():
    f = x**3*(a + b*log(c*x**n))**3
    F = -3*b**3*n**3*x**4/128 + 3*b**2*n**2*x**4*(a + b*log(c*x**n))/32 - 3*b*n*x**4*(a + b*log(c*x**n))**2/16 + x**4*(a + b*log(c*x**n))**3/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_58():
    f = x**2*(a + b*log(c*x**n))**3
    F = -2*b**3*n**3*x**3/27 + 2*b**2*n**2*x**3*(a + b*log(c*x**n))/9 - b*n*x**3*(a + b*log(c*x**n))**2/3 + x**3*(a + b*log(c*x**n))**3/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_59():
    f = x*(a + b*log(c*x**n))**3
    F = -3*b**3*n**3*x**2/8 + 3*b**2*n**2*x**2*(a + b*log(c*x**n))/4 - 3*b*n*x**2*(a + b*log(c*x**n))**2/4 + x**2*(a + b*log(c*x**n))**3/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_60():
    f = (a + b*log(c*x**n))**3
    F = 6*a*b**2*n**2*x - 6*b**3*n**3*x + 6*b**3*n**2*x*log(c*x**n) - 3*b*n*x*(a + b*log(c*x**n))**2 + x*(a + b*log(c*x**n))**3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_61():
    f = (a + b*log(c*x**n))**3/x
    F = (a + b*log(c*x**n))**4/(4*b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_62():
    f = (a + b*log(c*x**n))**3/x**2
    F = -6*b**3*n**3/x - 6*b**2*n**2*(a + b*log(c*x**n))/x - 3*b*n*(a + b*log(c*x**n))**2/x - (a + b*log(c*x**n))**3/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_63():
    f = (a + b*log(c*x**n))**3/x**3
    F = -3*b**3*n**3/(8*x**2) - 3*b**2*n**2*(a + b*log(c*x**n))/(4*x**2) - 3*b*n*(a + b*log(c*x**n))**2/(4*x**2) - (a + b*log(c*x**n))**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_64():
    f = (a + b*log(c*x**n))**3/x**4
    F = -2*b**3*n**3/(27*x**3) - 2*b**2*n**2*(a + b*log(c*x**n))/(9*x**3) - b*n*(a + b*log(c*x**n))**2/(3*x**3) - (a + b*log(c*x**n))**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_65():
    f = x**3/(a + b*log(c*x**n))
    F = ((x)**(Integer(4)) * sympy.Function('ExpIntegralEi')(((Integer(4) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * (sympy.E)**(((Integer(4) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_66():
    f = x**2/(a + b*log(c*x**n))
    F = ((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_67():
    f = x/(a + b*log(c*x**n))
    F = ((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_68():
    f = 1/(a + b*log(c*x**n))
    F = (x * sympy.Function('ExpIntegralEi')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_69():
    f = 1/(x*(a + b*log(c*x**n)))
    F = log(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_70():
    f = 1/(x**2*(a + b*log(c*x**n)))
    F = ((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Symbol('b') * Symbol('n') * x))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_71():
    f = 1/(x**3*(a + b*log(c*x**n)))
    F = ((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_72():
    f = 1/(x**4*(a + b*log(c*x**n)))
    F = ((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('n') * (x)**(Integer(3))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_73():
    f = x**3/(a + b*log(c*x**n))**2
    F = ((Integer(4) * (x)**(Integer(4)) * sympy.Function('ExpIntegralEi')(((Integer(4) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(4) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_74():
    f = x**2/(a + b*log(c*x**n))**2
    F = ((Integer(3) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_75():
    f = x/(a + b*log(c*x**n))**2
    F = ((Integer(2) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_76():
    f = (a + b*log(c*x**n))**(-2)
    F = ((x * sympy.Function('ExpIntegralEi')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (x * ((Symbol('b') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_77():
    f = 1/(x*(a + b*log(c*x**n))**2)
    F = -1/(b*n*(a + b*log(c*x**n)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_78():
    f = 1/(x**2*(a + b*log(c*x**n))**2)
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_79():
    f = 1/(x**3*(a + b*log(c*x**n))**2)
    F = ((Integer(-2) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_80():
    f = 1/(x**4*(a + b*log(c*x**n))**2)
    F = ((Integer(-3) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_81():
    f = x**3/(a + b*log(c*x**n))**3
    F = ((Integer(8) * (x)**(Integer(4)) * sympy.Function('ExpIntegralEi')(((Integer(4) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(4) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(3)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_82():
    f = x**2/(a + b*log(c*x**n))**3
    F = ((Integer(9) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(3)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_83():
    f = x/(a + b*log(c*x**n))**3
    F = ((Integer(2) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(3)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(2)) * (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_84():
    f = (a + b*log(c*x**n))**(-3)
    F = ((x * sympy.Function('ExpIntegralEi')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(3)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (x * ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (x * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_85():
    f = 1/(x*(a + b*log(c*x**n))**3)
    F = -1/(2*b*n*(a + b*log(c*x**n))**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_86():
    f = 1/(x**2*(a + b*log(c*x**n))**3)
    F = (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * x * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_87():
    f = 1/(x**3*(a + b*log(c*x**n))**3)
    F = ((Integer(2) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_88():
    f = 1/(x**4*(a + b*log(c*x**n))**3)
    F = ((Integer(9) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1))) + (Integer(3) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_89():
    f = (d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))
    F = -4*b*n*(d*x)**(sympy.S(7)/2)/(49*d) + 2*(d*x)**(sympy.S(7)/2)*(a + b*log(c*x**n))/(7*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_90():
    f = (d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))
    F = -4*b*n*(d*x)**(sympy.S(5)/2)/(25*d) + 2*(d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))/(5*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_91():
    f = sqrt(d*x)*(a + b*log(c*x**n))
    F = -4*b*n*(d*x)**(sympy.S(3)/2)/(9*d) + 2*(d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))/(3*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_92():
    f = (a + b*log(c*x**n))/sqrt(d*x)
    F = -4*b*n*sqrt(d*x)/d + 2*sqrt(d*x)*(a + b*log(c*x**n))/d
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_93():
    f = (a + b*log(c*x**n))/(d*x)**(sympy.S(3)/2)
    F = -4*b*n/(d*sqrt(d*x)) - (2*a + 2*b*log(c*x**n))/(d*sqrt(d*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_94():
    f = (a + b*log(c*x**n))/(d*x)**(sympy.S(5)/2)
    F = -4*b*n/(9*d*(d*x)**(sympy.S(3)/2)) - (2*a + 2*b*log(c*x**n))/(3*d*(d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_95():
    f = (d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))**2
    F = 16*b**2*n**2*(d*x)**(sympy.S(7)/2)/(343*d) - 8*b*n*(d*x)**(sympy.S(7)/2)*(a + b*log(c*x**n))/(49*d) + 2*(d*x)**(sympy.S(7)/2)*(a + b*log(c*x**n))**2/(7*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_96():
    f = (d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))**2
    F = 16*b**2*n**2*(d*x)**(sympy.S(5)/2)/(125*d) - 8*b*n*(d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))/(25*d) + 2*(d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))**2/(5*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_97():
    f = sqrt(d*x)*(a + b*log(c*x**n))**2
    F = 16*b**2*n**2*(d*x)**(sympy.S(3)/2)/(27*d) - 8*b*n*(d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))/(9*d) + 2*(d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))**2/(3*d)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_98():
    f = (a + b*log(c*x**n))**2/sqrt(d*x)
    F = 16*b**2*n**2*sqrt(d*x)/d - 8*b*n*sqrt(d*x)*(a + b*log(c*x**n))/d + 2*sqrt(d*x)*(a + b*log(c*x**n))**2/d
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_99():
    f = (a + b*log(c*x**n))**2/(d*x)**(sympy.S(3)/2)
    F = -16*b**2*n**2/(d*sqrt(d*x)) - 8*b*n*(a + b*log(c*x**n))/(d*sqrt(d*x)) - 2*(a + b*log(c*x**n))**2/(d*sqrt(d*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_100():
    f = (a + b*log(c*x**n))**2/(d*x)**(sympy.S(5)/2)
    F = -16*b**2*n**2/(27*d*(d*x)**(sympy.S(3)/2)) - 8*b*n*(a + b*log(c*x**n))/(9*d*(d*x)**(sympy.S(3)/2)) - 2*(a + b*log(c*x**n))**2/(3*d*(d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_101():
    f = (d*x)**(sympy.S(5)/2)/(a + b*log(c*x**n))
    F = (((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(7) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(7) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(7) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_102():
    f = (d*x)**(sympy.S(3)/2)/(a + b*log(c*x**n))
    F = (((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(5) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(5) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(5) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_103():
    f = sqrt(d*x)/(a + b*log(c*x**n))
    F = (((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_104():
    f = 1/(sqrt(d*x)*(a + b*log(c*x**n)))
    F = (sympy.sqrt((Symbol('d') * x)) * sympy.Function('ExpIntegralEi')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('a') * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * Symbol('n') * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('n')))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_105():
    f = 1/((d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n)))
    F = ((sympy.E)**((Symbol('a') * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('n')))**(Integer(-1))) * sympy.Function('ExpIntegralEi')(((Integer(-1) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_106():
    f = 1/((d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n)))
    F = ((sympy.E)**(((Integer(3) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * Symbol('n') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_107():
    f = (d*x)**(sympy.S(5)/2)/(a + b*log(c*x**n))**2
    F = ((Integer(7) * ((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(7) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(7) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(7) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_108():
    f = (d*x)**(sympy.S(3)/2)/(a + b*log(c*x**n))**2
    F = ((Integer(5) * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(5) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(5) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(5) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_109():
    f = sqrt(d*x)/(a + b*log(c*x**n))**2
    F = ((Integer(3) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * ((Integer(2) * Symbol('n')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_110():
    f = 1/(sqrt(d*x)*(a + b*log(c*x**n))**2)
    F = ((sympy.sqrt((Symbol('d') * x)) * sympy.Function('ExpIntegralEi')(((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**((Symbol('a') * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Symbol('n'))**(Integer(2)) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('n')))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt((Symbol('d') * x)) * ((Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_111():
    f = 1/((d*x)**(sympy.S(3)/2)*(a + b*log(c*x**n))**2)
    F = ((Integer(-1) * ((sympy.E)**((Symbol('a') * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('n')))**(Integer(-1))) * sympy.Function('ExpIntegralEi')(((Integer(-1) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * sympy.sqrt((Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt((Symbol('d') * x)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_112():
    f = 1/((d*x)**(sympy.S(5)/2)*(a + b*log(c*x**n))**2)
    F = ((Integer(-3) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')(((Integer(-3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('n') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_113():
    f = sqrt(a + b*log(c*x**n))
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('n'))))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_114():
    f = x**3*sqrt(log(a*x**n))
    F = (((Integer(-1) * (Integer(16))**(Integer(-1))) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_115():
    f = x**2*sqrt(log(a*x**n))
    F = (((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_116():
    f = x*sqrt(log(a*x**n))
    F = (((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_117():
    f = sqrt(log(a*x**n))
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))) + (x * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_118():
    f = sqrt(log(a*x**n))/x
    F = 2*log(a*x**n)**(sympy.S(3)/2)/(3*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_119():
    f = sqrt(log(a*x**n))/x**2
    F = ((sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_120():
    f = sqrt(log(a*x**n))/x**3
    F = ((sympy.sqrt(Symbol('n')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_121():
    f = x**3*log(a*x**n)**(sympy.S(3)/2)
    F = (((Integer(3) * (Integer(128))**(Integer(-1))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Integer(32))**(Integer(-1))) * Symbol('n') * (x)**(Integer(4)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_122():
    f = x**2*log(a*x**n)**(sympy.S(3)/2)
    F = (((Integer(12))**(Integer(-1)) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * Symbol('n') * (x)**(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_123():
    f = x*log(a*x**n)**(sympy.S(3)/2)
    F = (((Integer(3) * (Integer(16))**(Integer(-1))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('n') * (x)**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_124():
    f = log(a*x**n)**(sympy.S(3)/2)
    F = (((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('n') * x * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))))) + (x * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_125():
    f = log(a*x**n)**(sympy.S(3)/2)/x
    F = 2*log(a*x**n)**(sympy.S(5)/2)/(5*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_126():
    f = log(a*x**n)**(sympy.S(3)/2)/x**2
    F = ((Integer(3) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_127():
    f = log(a*x**n)**(sympy.S(3)/2)/x**3
    F = ((Integer(3) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_128():
    f = x**3/sqrt(log(a*x**n))
    F = (sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('n')) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_129():
    f = x**2/sqrt(log(a*x**n))
    F = (sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_130():
    f = x/sqrt(log(a*x**n))
    F = (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_131():
    f = 1/sqrt(log(a*x**n))
    F = (sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_132():
    f = 1/(x*sqrt(log(a*x**n)))
    F = 2*sqrt(log(a*x**n))/n
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_133():
    f = 1/(x**2*sqrt(log(a*x**n)))
    F = (sympy.sqrt(sympy.pi) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * x))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_134():
    f = 1/(x**3*sqrt(log(a*x**n)))
    F = (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * (x)**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_135():
    f = x**3/log(a*x**n)**(sympy.S(3)/2)
    F = ((Integer(4) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * ((Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_136():
    f = x**2/log(a*x**n)**(sympy.S(3)/2)
    F = ((Integer(2) * sympy.sqrt((Integer(3) * sympy.pi)) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_137():
    f = x/log(a*x**n)**(sympy.S(3)/2)
    F = ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2))) * ((Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_138():
    f = log(a*x**n)**(sympy.S(-3)/2)
    F = ((Integer(2) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x) * ((Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_139():
    f = 1/(x*log(a*x**n)**(sympy.S(3)/2))
    F = -2/(n*sqrt(log(a*x**n)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_140():
    f = 1/(x**2*log(a*x**n)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * ((Symbol('n') * x * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_141():
    f = 1/(x**3*log(a*x**n)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * (((Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * ((Symbol('n') * (x)**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_142():
    f = x**3/log(a*x**n)**(sympy.S(5)/2)
    F = ((Integer(32) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * (Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * ((Integer(3) * Symbol('n') * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(4))) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_143():
    f = x**2/log(a*x**n)**(sympy.S(5)/2)
    F = ((Integer(4) * sympy.sqrt((Integer(3) * sympy.pi)) * (x)**(Integer(3)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * Symbol('n') * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3))) * (((Symbol('n'))**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_144():
    f = x/log(a*x**n)**(sympy.S(5)/2)
    F = ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * (x)**(Integer(2)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2))) * ((Integer(3) * Symbol('n') * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2))) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_145():
    f = log(a*x**n)**(sympy.S(-5)/2)
    F = ((Integer(4) * sympy.sqrt(sympy.pi) * x * sympy.Function('Erfi')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x) * ((Integer(3) * Symbol('n') * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * x) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_146():
    f = 1/(x*log(a*x**n)**(sympy.S(5)/2))
    F = -2/(3*n*log(a*x**n)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_147():
    f = 1/(x**2*log(a*x**n)**(sympy.S(5)/2))
    F = ((Integer(4) * sympy.sqrt(sympy.pi) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n'))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (Integer(2) * ((Integer(3) * Symbol('n') * x * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(4) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * x * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_148():
    f = 1/(x**3*log(a*x**n)**(sympy.S(5)/2))
    F = ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * ((Symbol('a') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * ((Integer(3) * Symbol('n') * (x)**(Integer(2)) * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(8) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * (x)**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_149():
    f = (d*x)**m*(a + a*(m + 1)*log(c*x**n)/n)
    F = a*(d*x)**(m + 1)*log(c*x**n)/(d*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_150():
    f = (d*x)**m*(a + b*log(c*x**n))**3
    F = -6*b**3*n**3*(d*x)**(m + 1)/(d*(m + 1)**4) + 6*b**2*n**2*(d*x)**(m + 1)*(a + b*log(c*x**n))/(d*(m + 1)**3) - 3*b*n*(d*x)**(m + 1)*(a + b*log(c*x**n))**2/(d*(m + 1)**2) + (d*x)**(m + 1)*(a + b*log(c*x**n))**3/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_151():
    f = (d*x)**m*(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*(d*x)**(m + 1)/(d*(m + 1)**3) - 2*b*n*(d*x)**(m + 1)*(a + b*log(c*x**n))/(d*(m + 1)**2) + (d*x)**(m + 1)*(a + b*log(c*x**n))**2/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_152():
    f = (d*x)**m*(a + b*log(c*x**n))
    F = -b*n*(d*x)**(m + 1)/(d*(m + 1)**2) + (d*x)**(m + 1)*(a + b*log(c*x**n))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_153():
    f = (d*x)**m/(a + b*log(c*x**n))
    F = (((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('b') * Symbol('d') * Symbol('n'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_154():
    f = (d*x)**m/(a + b*log(c*x**n))**2
    F = (((Integer(1) + Symbol('m')) * ((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('d') * Symbol('n') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_155():
    f = (d*x)**m/(a + b*log(c*x**n))**3
    F = ((((Integer(1) + Symbol('m')))**(Integer(2)) * ((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('n'))**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + Symbol('m')) * ((Symbol('d') * x))**((Integer(1) + Symbol('m')))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('n'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_156():
    f = (d*x)**(n - 1)*log(c*x**n)**3
    F = (d*x)**n*log(c*x**n)**3/(d*n) - 3*(d*x)**n*log(c*x**n)**2/(d*n) + 6*(d*x)**n*log(c*x**n)/(d*n) - 6*(d*x)**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_157():
    f = (d*x)**(n - 1)*log(c*x**n)**2
    F = (d*x)**n*log(c*x**n)**2/(d*n) - 2*(d*x)**n*log(c*x**n)/(d*n) + 2*(d*x)**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_158():
    f = (d*x)**(n - 1)*log(c*x**n)
    F = (d*x)**n*log(c*x**n)/(d*n) - (d*x)**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_159():
    f = (d*x)**(n - 1)/log(c*x**n)
    F = ((x)**((Integer(1) + (Integer(-1) * Symbol('n')))) * ((Symbol('d') * x))**((Integer(-1) + Symbol('n'))) * sympy.Function('LogIntegral')((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('c') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_160():
    f = (d*x)**(n - 1)/log(c*x**n)**2
    F = (Integer(-1) * (((Symbol('d') * x))**(Symbol('n')) * ((Symbol('d') * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))**(Integer(-1)))) + (((x)**((Integer(1) + (Integer(-1) * Symbol('n')))) * ((Symbol('d') * x))**((Integer(-1) + Symbol('n'))) * sympy.Function('LogIntegral')((Symbol('c') * (x)**(Symbol('n'))))) * ((Symbol('c') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_161():
    f = (d*x)**(n - 1)/log(c*x**n)**3
    F = (Integer(-1) * (((Symbol('d') * x))**(Symbol('n')) * ((Integer(2) * Symbol('d') * Symbol('n') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d') * x))**(Symbol('n')) * ((Integer(2) * Symbol('d') * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))**(Integer(-1)))) + (((x)**((Integer(1) + (Integer(-1) * Symbol('n')))) * ((Symbol('d') * x))**((Integer(-1) + Symbol('n'))) * sympy.Function('LogIntegral')((Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('c') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_162():
    f = x**m*log(a*x**n)**(sympy.S(3)/2)
    F = ((Integer(3) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(4) * ((Integer(1) + Symbol('m')))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('n') * (x)**((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * ((Integer(2) * ((Integer(1) + Symbol('m')))**(Integer(2))))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_163():
    f = x**m*sqrt(log(a*x**n))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('n')) * sympy.sqrt(sympy.pi) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * ((Integer(1) + Symbol('m')))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_164():
    f = x**m/sqrt(log(a*x**n))
    F = (sympy.sqrt(sympy.pi) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(Symbol('n')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_165():
    f = x**m/log(a*x**n)**(sympy.S(3)/2)
    F = ((Integer(2) * sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.pi) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('n'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(1) + Symbol('m')))) * ((Symbol('n') * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_166():
    f = x**m/log(a*x**n)**(sympy.S(5)/2)
    F = ((Integer(4) * ((Integer(1) + Symbol('m')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((sympy.sqrt((Integer(1) + Symbol('m'))) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((((Symbol('a') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(3) * (Symbol('n'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(1) + Symbol('m')))) * ((Integer(3) * Symbol('n') * (sympy.log((Symbol('a') * (x)**(Symbol('n')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(1) + Symbol('m')) * (x)**((Integer(1) + Symbol('m')))) * ((Integer(3) * (Symbol('n'))**(Integer(2)) * sympy.sqrt(sympy.log((Symbol('a') * (x)**(Symbol('n')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_167():
    f = (d*x)**m*(a + b*log(c*x**n))**p
    F = (((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_168():
    f = x**2*(a + b*log(c*x**n))**p
    F = ((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_169():
    f = x*(a + b*log(c*x**n))**p
    F = ((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_170():
    f = (a + b*log(c*x**n))**p
    F = (x * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_171():
    f = (a + b*log(c*x**n))**p/x
    F = (a + b*log(c*x**n))**(p + 1)/(b*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_172():
    f = (a + b*log(c*x**n))**p/x**2
    F = Integer(-1) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))**(Symbol('p')) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_173():
    f = (a + b*log(c*x**n))**p/x**3
    F = Integer(-1) * (((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))**(Symbol('p')) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_174():
    f = (a + b*log(c*x**n))**p/x**4
    F = Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))**(Symbol('p')) * (x)**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_175():
    f = (d*x)**m*(a + b*log(c*x))**p
    F = (((Symbol('c') * x))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * ((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_176():
    f = x**2*(a + b*log(c*x))**p
    F = ((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c'))**(Integer(3))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_177():
    f = x*(a + b*log(c*x))**p
    F = ((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_178():
    f = (a + b*log(c*x))**p
    F = (sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_179():
    f = (a + b*log(c*x))**p/x
    F = (a + b*log(c*x))**(p + 1)/(b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_180():
    f = (a + b*log(c*x))**p/x**2
    F = ((Integer(-1) * Symbol('c')) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_181():
    f = (a + b*log(c*x))**p/x**3
    F = ((Integer(-1) * (Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('p'))))) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_182():
    f = (a + b*log(c*x))**p/x**4
    F = ((Integer(-1) * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p'))))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_183():
    f = (d*x)**m*(a + b*log(c*sqrt(x)))**p
    F = (((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(2) * (Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * (((Integer(2))**(Symbol('p')) * (sympy.E)**(((Integer(2) * Symbol('a') * (Integer(1) + Symbol('m'))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('c') * sympy.sqrt(x)))**((Integer(2) * (Integer(1) + Symbol('m')))) * ((Integer(-1) * (((Integer(1) + Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_184():
    f = x**2*(a + b*log(c*sqrt(x)))**p
    F = ((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(6) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * (((Integer(2))**(Symbol('p')) * (sympy.E)**(((Integer(6) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c'))**(Integer(6))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_185():
    f = x*(a + b*log(c*sqrt(x)))**p
    F = ((Integer(2))**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('p'))))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(4) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * (((sympy.E)**(((Integer(4) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c'))**(Integer(4))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_186():
    f = (a + b*log(c*sqrt(x)))**p
    F = (sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * (((Integer(2))**(Symbol('p')) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_187():
    f = (a + b*log(c*sqrt(x)))**p/x
    F = 2*(a + b*log(c*sqrt(x)))**(p + 1)/(b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_188():
    f = (a + b*log(c*sqrt(x)))**p/x**2
    F = ((Integer(-1) * (Integer(2))**((Integer(-1) * Symbol('p')))) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_189():
    f = (a + b*log(c*sqrt(x)))**p/x**3
    F = ((Integer(-1) * (Integer(2))**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('p')))))) * (Symbol('c'))**(Integer(4)) * (sympy.E)**(((Integer(4) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(4) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_190():
    f = (a + b*log(c*sqrt(x)))**p/x**4
    F = ((Integer(-1) * (Integer(2))**((Integer(-1) * Symbol('p')))) * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * (Symbol('c'))**(Integer(6)) * (sympy.E)**(((Integer(6) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Integer(6) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x)))))) * (Symbol('b'))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))))**(Symbol('p'))) * ((((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * sympy.sqrt(x))))) * (Symbol('b'))**(Integer(-1))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_191():
    f = x**(n - 1)*(a + b*log(c*x**n))**p
    F = (sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * (Symbol('b'))**(Integer(-1)))))**(Symbol('p')) * (Symbol('c') * Symbol('n'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_192():
    f = (d*x**q)**m*(a + b*log(c*x**n))**p
    F = (x * ((Symbol('d') * (x)**(Symbol('q'))))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + (Symbol('m') * Symbol('q'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') + (Symbol('a') * Symbol('m') * Symbol('q'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + (Symbol('m') * Symbol('q'))) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + (Symbol('m') * Symbol('q'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + (Symbol('m') * Symbol('q')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_1_2_d_x_pow_m_a_plus_b_log_c_x_pow_n_pow_p_193():
    f = (d1*x**q1)**m1*(d2*x**q2)**m2*(a + b*log(c*x**n))**p
    F = (x * ((Symbol('d1') * (x)**(Symbol('q1'))))**(Symbol('m1')) * ((Symbol('d2') * (x)**(Symbol('q2'))))**(Symbol('m2')) * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * (((Integer(1) + (Symbol('m1') * Symbol('q1')) + (Symbol('m2') * Symbol('q2'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + (Symbol('m1') * Symbol('q1')) + (Symbol('m2') * Symbol('q2')))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + (Symbol('m1') * Symbol('q1')) + (Symbol('m2') * Symbol('q2'))) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + (Symbol('m1') * Symbol('q1')) + (Symbol('m2') * Symbol('q2'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Integer(1) + (Symbol('m1') * Symbol('q1')) + (Symbol('m2') * Symbol('q2')))))**(Integer(-1))
    assert integrate(f, x) == F

