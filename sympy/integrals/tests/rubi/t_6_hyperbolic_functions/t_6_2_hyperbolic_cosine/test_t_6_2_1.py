"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.2 Hyperbolic cosine/6.2.1 (c+d x)^m (a+b cosh)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_1():
    f = (c + d*x)**4*cosh(a + b*x)
    F = (c + d*x)**4*sinh(a + b*x)/b - 4*d*(c + d*x)**3*cosh(a + b*x)/b**2 + 12*d**2*(c + d*x)**2*sinh(a + b*x)/b**3 - 24*d**3*(c + d*x)*cosh(a + b*x)/b**4 + 24*d**4*sinh(a + b*x)/b**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_2():
    f = (c + d*x)**3*cosh(a + b*x)
    F = (c + d*x)**3*sinh(a + b*x)/b - 3*d*(c + d*x)**2*cosh(a + b*x)/b**2 + 6*d**2*(c + d*x)*sinh(a + b*x)/b**3 - 6*d**3*cosh(a + b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_3():
    f = (c + d*x)**2*cosh(a + b*x)
    F = (c + d*x)**2*sinh(a + b*x)/b - 2*d*(c + d*x)*cosh(a + b*x)/b**2 + 2*d**2*sinh(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_4():
    f = (c + d*x)*cosh(a + b*x)
    F = (c + d*x)*sinh(a + b*x)/b - d*cosh(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_5():
    f = cosh(a + b*x)/(c + d*x)
    F = ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * (Symbol('d'))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_6():
    f = cosh(a + b*x)/(c + d*x)**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_7():
    f = cosh(a + b*x)/(c + d*x)**3
    F = ((Integer(-1) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_8():
    f = (c + d*x)**4*cosh(a + b*x)**2
    F = (c + d*x)**5/(10*d) + (c + d*x)**4*sinh(a + b*x)*cosh(a + b*x)/(2*b) - d*(c + d*x)**3*cosh(a + b*x)**2/b**2 + d*(c + d*x)**3/(2*b**2) + 3*d**2*(c + d*x)**2*sinh(a + b*x)*cosh(a + b*x)/(2*b**3) + 3*d**4*x/(4*b**4) - 3*d**3*(c + d*x)*cosh(a + b*x)**2/(2*b**4) + 3*d**4*sinh(a + b*x)*cosh(a + b*x)/(4*b**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_9():
    f = (c + d*x)**3*cosh(a + b*x)**2
    F = (c + d*x)**4/(8*d) + (c + d*x)**3*sinh(a + b*x)*cosh(a + b*x)/(2*b) + 3*c*d**2*x/(4*b**2) + 3*d**3*x**2/(8*b**2) - 3*d*(c + d*x)**2*cosh(a + b*x)**2/(4*b**2) + 3*d**2*(c + d*x)*sinh(a + b*x)*cosh(a + b*x)/(4*b**3) - 3*d**3*cosh(a + b*x)**2/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_10():
    f = (c + d*x)**2*cosh(a + b*x)**2
    F = (c + d*x)**3/(6*d) + (c + d*x)**2*sinh(a + b*x)*cosh(a + b*x)/(2*b) + d**2*x/(4*b**2) - d*(c + d*x)*cosh(a + b*x)**2/(2*b**2) + d**2*sinh(a + b*x)*cosh(a + b*x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_11():
    f = (c + d*x)*cosh(a + b*x)**2
    F = c*x/2 + d*x**2/4 + (c + d*x)*sinh(a + b*x)*cosh(a + b*x)/(2*b) - d*cosh(a + b*x)**2/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_12():
    f = cosh(a + b*x)**2/(c + d*x)
    F = ((sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (sympy.log((Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_13():
    f = cosh(a + b*x)**2/(c + d*x)**2
    F = (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x))) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_14():
    f = cosh(a + b*x)**2/(c + d*x)**3
    F = ((Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_15():
    f = cosh(a + b*x)**2/(c + d*x)**4
    F = ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x))) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_16():
    f = (c + d*x)**4*cosh(a + b*x)**3
    F = (c + d*x)**4*sinh(a + b*x)*cosh(a + b*x)**2/(3*b) + 2*(c + d*x)**4*sinh(a + b*x)/(3*b) - 4*d*(c + d*x)**3*cosh(a + b*x)**3/(9*b**2) - 8*d*(c + d*x)**3*cosh(a + b*x)/(3*b**2) + 4*d**2*(c + d*x)**2*sinh(a + b*x)*cosh(a + b*x)**2/(9*b**3) + 80*d**2*(c + d*x)**2*sinh(a + b*x)/(9*b**3) - 8*d**3*(c + d*x)*cosh(a + b*x)**3/(27*b**4) - 160*d**3*(c + d*x)*cosh(a + b*x)/(9*b**4) + 8*d**4*sinh(a + b*x)**3/(81*b**5) + 488*d**4*sinh(a + b*x)/(27*b**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_17():
    f = (c + d*x)**3*cosh(a + b*x)**3
    F = (c + d*x)**3*sinh(a + b*x)*cosh(a + b*x)**2/(3*b) + 2*(c + d*x)**3*sinh(a + b*x)/(3*b) - d*(c + d*x)**2*cosh(a + b*x)**3/(3*b**2) - 2*d*(c + d*x)**2*cosh(a + b*x)/b**2 + 2*d**2*(c + d*x)*sinh(a + b*x)*cosh(a + b*x)**2/(9*b**3) + 40*d**2*(c + d*x)*sinh(a + b*x)/(9*b**3) - 2*d**3*cosh(a + b*x)**3/(27*b**4) - 40*d**3*cosh(a + b*x)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_18():
    f = (c + d*x)**2*cosh(a + b*x)**3
    F = (c + d*x)**2*sinh(a + b*x)*cosh(a + b*x)**2/(3*b) + 2*(c + d*x)**2*sinh(a + b*x)/(3*b) - 2*d*(c + d*x)*cosh(a + b*x)**3/(9*b**2) - 4*d*(c + d*x)*cosh(a + b*x)/(3*b**2) + 2*d**2*sinh(a + b*x)**3/(27*b**3) + 14*d**2*sinh(a + b*x)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_19():
    f = (c + d*x)*cosh(a + b*x)**3
    F = (c + d*x)*sinh(a + b*x)*cosh(a + b*x)**2/(3*b) + (2*c + 2*d*x)*sinh(a + b*x)/(3*b) - d*cosh(a + b*x)**3/(9*b**2) - 2*d*cosh(a + b*x)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_20():
    f = cosh(a + b*x)**3/(c + d*x)
    F = ((Integer(3) * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_21():
    f = cosh(a + b*x)**3/(c + d*x)**2
    F = (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x))) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_22():
    f = cosh(a + b*x)**3/(c + d*x)**3
    F = ((Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_23():
    f = x**3*cosh(a + b*x)**4
    F = 3*x**4/32 + x**3*sinh(a + b*x)*cosh(a + b*x)**3/(4*b) + 3*x**3*sinh(a + b*x)*cosh(a + b*x)/(8*b) - 3*x**2*cosh(a + b*x)**4/(16*b**2) - 9*x**2*cosh(a + b*x)**2/(16*b**2) + 45*x**2/(128*b**2) + 3*x*sinh(a + b*x)*cosh(a + b*x)**3/(32*b**3) + 45*x*sinh(a + b*x)*cosh(a + b*x)/(64*b**3) - 3*cosh(a + b*x)**4/(128*b**4) - 45*cosh(a + b*x)**2/(128*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_24():
    f = x**2*cosh(a + b*x)**4
    F = x**3/8 + x**2*sinh(a + b*x)*cosh(a + b*x)**3/(4*b) + 3*x**2*sinh(a + b*x)*cosh(a + b*x)/(8*b) - x*cosh(a + b*x)**4/(8*b**2) - 3*x*cosh(a + b*x)**2/(8*b**2) + 15*x/(64*b**2) + sinh(a + b*x)*cosh(a + b*x)**3/(32*b**3) + 15*sinh(a + b*x)*cosh(a + b*x)/(64*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_25():
    f = x*cosh(a + b*x)**4
    F = 3*x**2/16 + x*sinh(a + b*x)*cosh(a + b*x)**3/(4*b) + 3*x*sinh(a + b*x)*cosh(a + b*x)/(8*b) - cosh(a + b*x)**4/(16*b**2) - 3*cosh(a + b*x)**2/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_26():
    f = (c + d*x)**3*sech(a + b*x)
    F = ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * sympy.I) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Integer(3) * sympy.I) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (((Integer(6) * sympy.I) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Integer(6) * sympy.I) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Integer(6) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((Integer(6) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_27():
    f = (c + d*x)**2*sech(a + b*x)
    F = ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((Integer(2) * sympy.I) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Integer(2) * sympy.I) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (((Integer(2) * sympy.I) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Integer(2) * sympy.I) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_28():
    f = (c + d*x)*sech(a + b*x)
    F = ((Integer(2) * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_29():
    f = sech(a + b*x)/(c + d*x)
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_30():
    f = sech(a + b*x)/(c + d*x)**2
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_31():
    f = (c + d*x)**3*sech(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_32():
    f = (c + d*x)**2*sech(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_33():
    f = (c + d*x)*sech(a + b*x)**2
    F = (c + d*x)*tanh(a + b*x)/b - d*log(cosh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_34():
    f = sech(a + b*x)**2/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_35():
    f = sech(a + b*x)**2/(c + d*x)**2
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_36():
    f = (c + d*x)**3*sech(a + b*x)**3
    F = ((Integer(-6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((Integer(3) * sympy.I) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_37():
    f = (c + d*x)**2*sech(a + b*x)**3
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_38():
    f = (c + d*x)*sech(a + b*x)**3
    F = (((Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(2))**(Integer(-1))) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.I * (Integer(2))**(Integer(-1))) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_39():
    f = sech(a + b*x)**3/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_40():
    f = sech(a + b*x)**3/(c + d*x)**2
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_41():
    f = (c + d*x)**(sympy.S(5)/2)*cosh(a + b*x)
    F = ((Integer(-5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_42():
    f = (c + d*x)**(sympy.S(3)/2)*cosh(a + b*x)
    F = ((Integer(-3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_43():
    f = sqrt(c + d*x)*cosh(a + b*x)
    F = ((sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_44():
    f = cosh(a + b*x)/sqrt(c + d*x)
    F = (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_45():
    f = cosh(a + b*x)/(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_46():
    f = cosh(a + b*x)/(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(-2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_47():
    f = cosh(a + b*x)/(c + d*x)**(sympy.S(7)/2)
    F = ((Integer(-2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_48():
    f = (c + d*x)**(sympy.S(5)/2)*cosh(a + b*x)**2
    F = ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(256) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(256) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(64) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_49():
    f = (c + d*x)**(sympy.S(3)/2)*cosh(a + b*x)**2
    F = ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_50():
    f = sqrt(c + d*x)*cosh(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_51():
    f = cosh(a + b*x)**2/sqrt(c + d*x)
    F = (sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_52():
    f = cosh(a + b*x)**2/(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_53():
    f = cosh(a + b*x)**2/(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(-2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_54():
    f = cosh(a + b*x)**2/(c + d*x)**(sympy.S(7)/2)
    F = ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_55():
    f = cosh(a + b*x)**2/(c + d*x)**(sympy.S(9)/2)
    F = ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(105) * (Symbol('d'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(7) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(105) * (Symbol('d'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(105) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(32) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(105) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(35) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(105) * (Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_56():
    f = (c + d*x)**(sympy.S(5)/2)*cosh(a + b*x)**3
    F = ((Integer(-5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(18) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(45) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(576) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(45) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(576) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(45) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Integer(5) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(144) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_57():
    f = (c + d*x)**(sympy.S(3)/2)*cosh(a + b*x)**3
    F = (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(9) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(96) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(9) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(96) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_58():
    f = sqrt(c + d*x)*cosh(a + b*x)**3
    F = ((Integer(3) * sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(48) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('d')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(48) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_59():
    f = cosh(a + b*x)**3/sqrt(c + d*x)
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_60():
    f = cosh(a + b*x)**3/(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(Symbol('b')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_61():
    f = cosh(a + b*x)**3/(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(-2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_62():
    f = cosh(a + b*x)**3/(c + d*x)**(sympy.S(7)/2)
    F = ((Integer(16) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(5) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_63():
    f = (d*x)**(sympy.S(3)/2)*cosh(f*x)
    F = ((Integer(-3) * Symbol('d') * sympy.sqrt((Symbol('d') * x)) * sympy.cosh((Symbol('f') * x))) * ((Integer(2) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('f'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('f'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sinh((Symbol('f') * x))) * (Symbol('f'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_64():
    f = sqrt(d*x)*cosh(f*x)
    F = ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('f'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('f'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('d') * x)) * sympy.sinh((Symbol('f') * x))) * (Symbol('f'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_65():
    f = cosh(f*x)/sqrt(d*x)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('f'))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_66():
    f = cosh(f*x)/(d*x)**(sympy.S(3)/2)
    F = ((Integer(-2) * sympy.cosh((Symbol('f') * x))) * ((Symbol('d') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_67():
    f = cosh(f*x)/(d*x)**(sympy.S(5)/2)
    F = ((Integer(-2) * sympy.cosh((Symbol('f') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('f'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('f'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('f')) * sympy.sqrt((Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('f') * sympy.sinh((Symbol('f') * x))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * x))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_68():
    f = sqrt(c + d*x)*sech(a + b*x)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_69():
    f = sech(a + b*x)/sqrt(c + d*x)
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_70():
    f = cosh(x)**(sympy.S(3)/2)/x**3
    F = ((Integer(-1) * (sympy.cosh(x))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.cosh(x)) * sympy.sinh(x)) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.cosh(x))))**(Integer(-1)), x)) * (Integer(8))**(Integer(-1)))) + ((Integer(9) * sympy.Function('Unintegrable')(((sympy.cosh(x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)) * (Integer(8))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_71():
    f = x*sqrt(cosh(x)) + x/cosh(x)**(sympy.S(3)/2)
    F = 2*x*sinh(x)/sqrt(cosh(x)) - 4*sqrt(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_72():
    f = -x/(3*sqrt(cosh(x))) + x/cosh(x)**(sympy.S(5)/2)
    F = 2*x*sinh(x)/(3*cosh(x)**(sympy.S(3)/2)) + 4/(3*sqrt(cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_73():
    f = 3*x*sqrt(cosh(x))/5 + x/cosh(x)**(sympy.S(7)/2)
    F = 6*x*sinh(x)/(5*sqrt(cosh(x))) + 2*x*sinh(x)/(5*cosh(x)**(sympy.S(5)/2)) - 12*sqrt(cosh(x))/5 + 4/(15*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_74():
    f = x**2*sqrt(cosh(x)) + x**2/cosh(x)**(sympy.S(3)/2)
    F = 2*x**2*sinh(x)/sqrt(cosh(x)) - 8*x*sqrt(cosh(x)) - 16*I*elliptic_e(I*x/2, 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_75():
    f = (b*cosh(e + f*x))**n*(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_76():
    f = (c + d*x)**m*cosh(a + b*x)**3
    F = (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-3) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * ((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(8) * Symbol('b') * ((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * (((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**(((Integer(-3) * Symbol('a')) + ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(3) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * (((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_77():
    f = (c + d*x)**m*cosh(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * ((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * (((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_78():
    f = (c + d*x)**m*cosh(a + b*x)
    F = (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * ((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_79():
    f = (c + d*x)**m*sech(a + b*x)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_80():
    f = (c + d*x)**m*sech(a + b*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_81():
    f = x**(m + 3)*cosh(a + b*x)
    F = ((Integer(-1) * ((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_82():
    f = x**(m + 2)*cosh(a + b*x)
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_83():
    f = x**(m + 1)*cosh(a + b*x)
    F = ((Integer(-1) * ((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_84():
    f = x**m*cosh(a + b*x)
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_85():
    f = x**(m - 1)*cosh(a + b*x)
    F = ((Integer(-1) * ((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Symbol('b') * x))) * ((Integer(2) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_86():
    f = x**(m - 2)*cosh(a + b*x)
    F = ((Symbol('b') * (sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_87():
    f = x**(m - 3)*cosh(a + b*x)
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Symbol('b') * x))) * ((Integer(2) * (sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_88():
    f = x**(m + 3)*cosh(a + b*x)**2
    F = ((x)**((Integer(4) + Symbol('m'))) * ((Integer(2) * (Integer(4) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-6) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-6) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_89():
    f = x**(m + 2)*cosh(a + b*x)**2
    F = ((x)**((Integer(3) + Symbol('m'))) * ((Integer(2) * (Integer(3) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-5) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-5) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_90():
    f = x**(m + 1)*cosh(a + b*x)**2
    F = ((x)**((Integer(2) + Symbol('m'))) * ((Integer(2) * (Integer(2) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_91():
    f = x**m*cosh(a + b*x)**2
    F = ((x)**((Integer(1) + Symbol('m'))) * ((Integer(2) * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * ((Symbol('b') * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * ((Symbol('b') * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_92():
    f = x**(m - 1)*cosh(a + b*x)**2
    F = ((x)**(Symbol('m')) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Integer(-2) * Symbol('b') * x))) * (((Integer(-1) * (Symbol('b') * x)))**(Symbol('m')))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_93():
    f = x**(m - 2)*cosh(a + b*x)**2
    F = ((Integer(-1) * (x)**((Integer(-1) + Symbol('m')))) * ((Integer(2) * (Integer(1) + (Integer(-1) * Symbol('m')))))**(Integer(-1))) + (((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * Symbol('b') * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((Integer(-1) * (Symbol('b') * x)))**(Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * Symbol('b') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_94():
    f = x**(m - 3)*cosh(a + b*x)**2
    F = ((Integer(-1) * (x)**((Integer(-2) + Symbol('m')))) * ((Integer(2) * (Integer(2) + (Integer(-1) * Symbol('m')))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((Integer(2))**(Symbol('m')) * ((Integer(-1) * (Symbol('b') * x)))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((Integer(2))**(Symbol('m')) * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_95():
    f = -x*sqrt(sech(x))/3 + x/sech(x)**(sympy.S(3)/2)
    F = 2*x*sinh(x)/(3*sqrt(sech(x))) - 4/(9*sech(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_96():
    f = -3*x/(5*sqrt(sech(x))) + x/sech(x)**(sympy.S(5)/2)
    F = 2*x*sinh(x)/(5*sech(x)**(sympy.S(3)/2)) - 4/(25*sech(x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_97():
    f = -5*x*sqrt(sech(x))/21 + x/sech(x)**(sympy.S(7)/2)
    F = 10*x*sinh(x)/(21*sqrt(sech(x))) + 2*x*sinh(x)/(7*sech(x)**(sympy.S(5)/2)) - 20/(63*sech(x)**(sympy.S(3)/2)) - 4/(49*sech(x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_98():
    f = -x**2*sqrt(sech(x))/3 + x**2/sech(x)**(sympy.S(3)/2)
    F = 2*x**2*sinh(x)/(3*sqrt(sech(x))) - 8*x/(9*sech(x)**(sympy.S(3)/2)) + 16*sinh(x)/(27*sqrt(sech(x))) - 16*I*sqrt(cosh(x))*elliptic_f(I*x/2, 2)*sqrt(sech(x))/27
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_99():
    f = (c + d*x)**3*(a*cosh(e + f*x) + a)
    F = -6*a*d**3*cosh(e + f*x)/f**4 + 6*a*d**2*(c + d*x)*sinh(e + f*x)/f**3 - 3*a*d*(c + d*x)**2*cosh(e + f*x)/f**2 + a*(c + d*x)**3*sinh(e + f*x)/f + a*(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_100():
    f = (c + d*x)**2*(a*cosh(e + f*x) + a)
    F = 2*a*d**2*sinh(e + f*x)/f**3 - 2*a*d*(c + d*x)*cosh(e + f*x)/f**2 + a*(c + d*x)**2*sinh(e + f*x)/f + a*(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_101():
    f = (c + d*x)*(a*cosh(e + f*x) + a)
    F = -a*d*cosh(e + f*x)/f**2 + a*(c + d*x)*sinh(e + f*x)/f + a*(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_102():
    f = (a*cosh(e + f*x) + a)/(c + d*x)
    F = ((Symbol('a') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('a') * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_103():
    f = (a*cosh(e + f*x) + a)/(c + d*x)**2
    F = (Integer(-1) * (Symbol('a') * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('a') * Symbol('f') * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('a') * Symbol('f') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_104():
    f = (a*cosh(e + f*x) + a)/(c + d*x)**3
    F = (Integer(-1) * (Symbol('a') * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('f'))**(Integer(2)) * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('f') * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('f'))**(Integer(2)) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_105():
    f = (c + d*x)**3*(a*cosh(e + f*x) + a)**2
    F = 3*a**2*c*d**2*x/(4*f**2) + 3*a**2*d**3*x**2/(8*f**2) - 3*a**2*d**3*cosh(e + f*x)**2/(8*f**4) - 12*a**2*d**3*cosh(e + f*x)/f**4 + 3*a**2*d**2*(c + d*x)*sinh(e + f*x)*cosh(e + f*x)/(4*f**3) + 12*a**2*d**2*(c + d*x)*sinh(e + f*x)/f**3 - 3*a**2*d*(c + d*x)**2*cosh(e + f*x)**2/(4*f**2) - 6*a**2*d*(c + d*x)**2*cosh(e + f*x)/f**2 + a**2*(c + d*x)**3*sinh(e + f*x)*cosh(e + f*x)/(2*f) + 2*a**2*(c + d*x)**3*sinh(e + f*x)/f + 3*a**2*(c + d*x)**4/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_106():
    f = (c + d*x)**2*(a*cosh(e + f*x) + a)**2
    F = a**2*d**2*x/(4*f**2) + a**2*d**2*sinh(e + f*x)*cosh(e + f*x)/(4*f**3) + 4*a**2*d**2*sinh(e + f*x)/f**3 - a**2*d*(c + d*x)*cosh(e + f*x)**2/(2*f**2) - 4*a**2*d*(c + d*x)*cosh(e + f*x)/f**2 + a**2*(c + d*x)**2*sinh(e + f*x)*cosh(e + f*x)/(2*f) + 2*a**2*(c + d*x)**2*sinh(e + f*x)/f + a**2*(c + d*x)**3/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_107():
    f = (c + d*x)*(a*cosh(e + f*x) + a)**2
    F = a**2*c*x/2 + a**2*d*x**2/4 - a**2*d*cosh(e + f*x)**2/(4*f**2) - 2*a**2*d*cosh(e + f*x)/f**2 + a**2*(c + d*x)*sinh(e + f*x)*cosh(e + f*x)/(2*f) + 2*a**2*(c + d*x)*sinh(e + f*x)/f + a**2*(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_108():
    f = (a*cosh(e + f*x) + a)**2/(c + d*x)
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_109():
    f = (a*cosh(e + f*x) + a)**2/(c + d*x)**2
    F = ((Integer(-4) * (Symbol('a'))**(Integer(2)) * (sympy.cosh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(4))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x))) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_110():
    f = (a*cosh(e + f*x) + a)**2/(c + d*x)**3
    F = ((Integer(-2) * (Symbol('a'))**(Integer(2)) * (sympy.cosh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(4))) * ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('f') * (sympy.cosh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(3)) * sympy.sinh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_111():
    f = (c + d*x)**3/(a*cosh(e + f*x) + a)
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Symbol('a') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(12) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_112():
    f = (c + d*x)**2/(a*cosh(e + f*x) + a)
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_113():
    f = (c + d*x)/(a*cosh(e + f*x) + a)
    F = -2*d*log(cosh(e/2 + f*x/2))/(a*f**2) + (c + d*x)*tanh(e/2 + f*x/2)/(a*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_114():
    f = 1/((c + d*x)*(a*cosh(e + f*x) + a))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * (Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_115():
    f = 1/((c + d*x)**2*(a*cosh(e + f*x) + a))
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_116():
    f = (c + d*x)**3/(a*cosh(e + f*x) + a)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(3)) * sympy.log(sympy.cosh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.sech(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.sech(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_117():
    f = (c + d*x)**2/(a*cosh(e + f*x) + a)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * (sympy.sech(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.sech(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2)) * sympy.tanh(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_118():
    f = (c + d*x)/(a*cosh(e + f*x) + a)**2
    F = -2*d*log(cosh(e/2 + f*x/2))/(3*a**2*f**2) + d*sech(e/2 + f*x/2)**2/(6*a**2*f**2) + (c + d*x)*tanh(e/2 + f*x/2)*sech(e/2 + f*x/2)**2/(6*a**2*f) + (c + d*x)*tanh(e/2 + f*x/2)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_119():
    f = 1/((c + d*x)*(a*cosh(e + f*x) + a)**2)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_120():
    f = 1/((c + d*x)**2*(a*cosh(e + f*x) + a)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_121():
    f = x**3*sqrt(a*cosh(c + d*x) + a)
    F = 2*x**3*sqrt(a*cosh(c + d*x) + a)*tanh(c/2 + d*x/2)/d - 12*x**2*sqrt(a*cosh(c + d*x) + a)/d**2 + 48*x*sqrt(a*cosh(c + d*x) + a)*tanh(c/2 + d*x/2)/d**3 - 96*sqrt(a*cosh(c + d*x) + a)/d**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_122():
    f = x**2*sqrt(a*cosh(c + d*x) + a)
    F = 2*x**2*sqrt(a*cosh(c + d*x) + a)*tanh(c/2 + d*x/2)/d - 8*x*sqrt(a*cosh(c + d*x) + a)/d**2 + 16*sqrt(a*cosh(c + d*x) + a)*tanh(c/2 + d*x/2)/d**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_123():
    f = x*sqrt(a*cosh(c + d*x) + a)
    F = 2*x*sqrt(a*cosh(c + d*x) + a)*tanh(c/2 + d*x/2)/d - 4*sqrt(a*cosh(c + d*x) + a)/d**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_124():
    f = sqrt(a*cosh(c + d*x) + a)/x
    F = (sympy.cosh((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CoshIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) + (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sinh((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_125():
    f = sqrt(a*cosh(c + d*x) + a)/x**2
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * (x)**(Integer(-1)))) + ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CoshIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sinh((Symbol('c') * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1))) + ((Symbol('d') * sympy.cosh((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('SinhIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_126():
    f = sqrt(a*cosh(c + d*x) + a)/x**3
    F = ((Integer(-1) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.cosh((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CoshIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * (Integer(8))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.sech(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sinh((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))) * sympy.tanh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(4) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_127():
    f = x**3*sqrt(a*cosh(x) + a)
    F = 2*x**3*sqrt(a*cosh(x) + a)*tanh(x/2) - 12*x**2*sqrt(a*cosh(x) + a) + 48*x*sqrt(a*cosh(x) + a)*tanh(x/2) - 96*sqrt(a*cosh(x) + a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_128():
    f = x**2*sqrt(a*cosh(x) + a)
    F = 2*x**2*sqrt(a*cosh(x) + a)*tanh(x/2) - 8*x*sqrt(a*cosh(x) + a) + 16*sqrt(a*cosh(x) + a)*tanh(x/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_129():
    f = x*sqrt(a*cosh(x) + a)
    F = 2*x*sqrt(a*cosh(x) + a)*tanh(x/2) - 4*sqrt(a*cosh(x) + a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_130():
    f = sqrt(a*cosh(x) + a)/x
    F = sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_131():
    f = sqrt(a*cosh(x) + a)/x**2
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.sech((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinhIntegral')((x * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_132():
    f = sqrt(a*cosh(x) + a)/x**3
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(8))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.tanh((x * (Integer(2))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_133():
    f = x**3*(a*cosh(x) + a)**(sympy.S(3)/2)
    F = 4*a*x**3*sqrt(a*cosh(x) + a)*sinh(x/2)*cosh(x/2)/3 + 8*a*x**3*sqrt(a*cosh(x) + a)*tanh(x/2)/3 - 8*a*x**2*sqrt(a*cosh(x) + a)*cosh(x/2)**2/3 - 16*a*x**2*sqrt(a*cosh(x) + a) + 32*a*x*sqrt(a*cosh(x) + a)*sinh(x/2)*cosh(x/2)/9 + 640*a*x*sqrt(a*cosh(x) + a)*tanh(x/2)/9 - 64*a*sqrt(a*cosh(x) + a)*cosh(x/2)**2/27 - 1280*a*sqrt(a*cosh(x) + a)/9
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_134():
    f = x**2*(a*cosh(x) + a)**(sympy.S(3)/2)
    F = 4*a*x**2*sqrt(a*cosh(x) + a)*sinh(x/2)*cosh(x/2)/3 + 8*a*x**2*sqrt(a*cosh(x) + a)*tanh(x/2)/3 - 16*a*x*sqrt(a*cosh(x) + a)*cosh(x/2)**2/9 - 32*a*x*sqrt(a*cosh(x) + a)/3 + 32*a*sqrt(a*cosh(x) + a)*sinh(x/2)**2*tanh(x/2)/27 + 224*a*sqrt(a*cosh(x) + a)*tanh(x/2)/9
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_135():
    f = x*(a*cosh(x) + a)**(sympy.S(3)/2)
    F = 4*a*x*sqrt(a*cosh(x) + a)*sinh(x/2)*cosh(x/2)/3 + 8*a*x*sqrt(a*cosh(x) + a)*tanh(x/2)/3 - 8*a*sqrt(a*cosh(x) + a)*cosh(x/2)**2/9 - 16*a*sqrt(a*cosh(x) + a)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_136():
    f = (a*cosh(x) + a)**(sympy.S(3)/2)/x
    F = ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_137():
    f = (a*cosh(x) + a)**(sympy.S(3)/2)/x**2
    F = ((Integer(-2) * Symbol('a') * (sympy.cosh((x * (Integer(2))**(Integer(-1)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))) * (x)**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.sech((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinhIntegral')((x * (Integer(2))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.sech((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_138():
    f = (a*cosh(x) + a)**(sympy.S(3)/2)/x**3
    F = (Integer(-1) * ((Symbol('a') * (sympy.cosh((x * (Integer(2))**(Integer(-1)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))) * ((x)**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + ((Integer(9) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.Function('CoshIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1)))) * sympy.sech((x * (Integer(2))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x)))) * sympy.sinh((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_139():
    f = x**3/sqrt(a*cosh(c + d*x) + a)
    F = ((Integer(4) * (x)**(Integer(3)) * sympy.atan((sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Integer(12) * sympy.I) * (x)**(Integer(2)) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (((Integer(12) * sympy.I) * (x)**(Integer(2)) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Integer(48) * sympy.I) * x * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Integer(48) * sympy.I) * x * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(96) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (((Integer(96) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_140():
    f = x**2/sqrt(a*cosh(c + d*x) + a)
    F = ((Integer(4) * (x)**(Integer(2)) * sympy.atan((sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Integer(8) * sympy.I) * x * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (((Integer(8) * sympy.I) * x * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Integer(16) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Integer(16) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_141():
    f = x/sqrt(a*cosh(c + d*x) + a)
    F = ((Integer(4) * x * sympy.atan((sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Integer(4) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (((Integer(4) * sympy.I) * sympy.cosh(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_142():
    f = 1/(x*sqrt(a*cosh(c + d*x) + a))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_143():
    f = 1/(x**2*sqrt(a*cosh(c + d*x) + a))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_144():
    f = x**3/(a*cosh(x) + a)**(sympy.S(3)/2)
    F = ((Integer(3) * (x)**(Integer(2))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * x * sympy.atan((sympy.E)**((x * (Integer(2))**(Integer(-1))))) * sympy.cosh((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.atan((sympy.E)**((x * (Integer(2))**(Integer(-1))))) * sympy.cosh((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (((Integer(24) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * sympy.I) * (x)**(Integer(2)) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(24) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (((Integer(3) * sympy.I) * (x)**(Integer(2)) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (((Integer(12) * sympy.I) * x * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * (((Integer(12) * sympy.I) * x * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(24) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (((Integer(24) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.tanh((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_145():
    f = x**2/(a*cosh(x) + a)**(sympy.S(3)/2)
    F = ((Integer(2) * x) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.atan((sympy.E)**((x * (Integer(2))**(Integer(-1))))) * sympy.cosh((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.atan(sympy.sinh((x * (Integer(2))**(Integer(-1))))) * sympy.cosh((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2) * sympy.I) * x * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (((Integer(2) * sympy.I) * x * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (((Integer(4) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * (((Integer(4) * sympy.I) * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.tanh((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_146():
    f = x/(a*cosh(x) + a)**(sympy.S(3)/2)
    F = ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)) + ((x * sympy.atan((sympy.E)**((x * (Integer(2))**(Integer(-1))))) * sympy.cosh((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))) + ((sympy.I * sympy.cosh((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((x * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1))) + ((x * sympy.tanh((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cosh(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_147():
    f = 1/(x*(a*cosh(x) + a)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('a') * sympy.cosh(x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_148():
    f = 1/(x**2*(a*cosh(x) + a)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * ((Symbol('a') + (Symbol('a') * sympy.cosh(x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_149():
    f = (a*cosh(c + d*x) + a)**(sympy.S(1)/3)/x
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_150():
    f = (c + d*x)**m*(a*cosh(e + f*x) + a)**n
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_151():
    f = (c + d*x)**m*(a*cosh(e + f*x) + a)**3
    F = ((Integer(5) * (Symbol('a'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('e')) + (Integer(-1) * ((Integer(3) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-3) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(3) * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(15) * (Symbol('a'))**(Integer(3)) * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(8) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-2) * Symbol('e')) + ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-3) * Symbol('e')) + ((Integer(3) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(3) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_152():
    f = (c + d*x)**m*(a*cosh(e + f*x) + a)**2
    F = ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('a'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('e')) + ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_153():
    f = (c + d*x)**m*(a*cosh(e + f*x) + a)
    F = ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_154():
    f = (c + d*x)**m/(a*cosh(e + f*x) + a)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_155():
    f = (c + d*x)**m/(a*cosh(e + f*x) + a)**2
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (((Symbol('a') + (Symbol('a') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_156():
    f = (a + b*cosh(e + f*x))*(c + d*x)**3
    F = a*(c + d*x)**4/(4*d) - 6*b*d**3*cosh(e + f*x)/f**4 + 6*b*d**2*(c + d*x)*sinh(e + f*x)/f**3 - 3*b*d*(c + d*x)**2*cosh(e + f*x)/f**2 + b*(c + d*x)**3*sinh(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_157():
    f = (a + b*cosh(e + f*x))*(c + d*x)**2
    F = a*(c + d*x)**3/(3*d) + 2*b*d**2*sinh(e + f*x)/f**3 - 2*b*d*(c + d*x)*cosh(e + f*x)/f**2 + b*(c + d*x)**2*sinh(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_158():
    f = (a + b*cosh(e + f*x))*(c + d*x)
    F = a*(c + d*x)**2/(2*d) - b*d*cosh(e + f*x)/f**2 + b*(c + d*x)*sinh(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_159():
    f = (a + b*cosh(e + f*x))/(c + d*x)
    F = ((Symbol('b') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_160():
    f = (a + b*cosh(e + f*x))/(c + d*x)**2
    F = (Integer(-1) * (Symbol('a') * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * Symbol('f') * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('f') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_161():
    f = (a + b*cosh(e + f*x))/(c + d*x)**3
    F = ((Integer(-1) * Symbol('a')) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('f') * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_162():
    f = (a + b*cosh(e + f*x))**2*(c + d*x)**3
    F = a**2*(c + d*x)**4/(4*d) - 12*a*b*d**3*cosh(e + f*x)/f**4 + 12*a*b*d**2*(c + d*x)*sinh(e + f*x)/f**3 - 6*a*b*d*(c + d*x)**2*cosh(e + f*x)/f**2 + 2*a*b*(c + d*x)**3*sinh(e + f*x)/f + 3*b**2*c*d**2*x/(4*f**2) + 3*b**2*d**3*x**2/(8*f**2) - 3*b**2*d**3*cosh(e + f*x)**2/(8*f**4) + 3*b**2*d**2*(c + d*x)*sinh(e + f*x)*cosh(e + f*x)/(4*f**3) - 3*b**2*d*(c + d*x)**2*cosh(e + f*x)**2/(4*f**2) + b**2*(c + d*x)**3*sinh(e + f*x)*cosh(e + f*x)/(2*f) + b**2*(c + d*x)**4/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_163():
    f = (a + b*cosh(e + f*x))**2*(c + d*x)**2
    F = a**2*(c + d*x)**3/(3*d) + 4*a*b*d**2*sinh(e + f*x)/f**3 - 4*a*b*d*(c + d*x)*cosh(e + f*x)/f**2 + 2*a*b*(c + d*x)**2*sinh(e + f*x)/f + b**2*d**2*x/(4*f**2) + b**2*d**2*sinh(e + f*x)*cosh(e + f*x)/(4*f**3) - b**2*d*(c + d*x)*cosh(e + f*x)**2/(2*f**2) + b**2*(c + d*x)**2*sinh(e + f*x)*cosh(e + f*x)/(2*f) + b**2*(c + d*x)**3/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_164():
    f = (a + b*cosh(e + f*x))**2*(c + d*x)
    F = a**2*(c + d*x)**2/(2*d) - 2*a*b*d*cosh(e + f*x)/f**2 + 2*a*b*(c + d*x)*sinh(e + f*x)/f + b**2*c*x/2 + b**2*d*x**2/4 - b**2*d*cosh(e + f*x)**2/(4*f**2) + b**2*(c + d*x)*sinh(e + f*x)*cosh(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_165():
    f = (a + b*cosh(e + f*x))**2/(c + d*x)
    F = ((Integer(2) * Symbol('a') * Symbol('b') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_166():
    f = (a + b*cosh(e + f*x))**2/(c + d*x)**2
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x))) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_167():
    f = (a + b*cosh(e + f*x))**2/(c + d*x)**3
    F = ((Integer(-1) * (Symbol('a'))**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.cosh((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * Symbol('b') * (Symbol('f'))**(Integer(2)) * sympy.cosh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * Symbol('f') * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.cosh((Symbol('e') + (Symbol('f') * x))) * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('a') * Symbol('b') * (Symbol('f'))**(Integer(2)) * sympy.sinh((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_168():
    f = (c + d*x)**3/(a + b*cosh(e + f*x))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_169():
    f = (c + d*x)**2/(a + b*cosh(e + f*x))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_170():
    f = (c + d*x)/(a + b*cosh(e + f*x))
    F = (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('f')))**(Integer(-1)))) + ((Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_171():
    f = 1/((a + b*cosh(e + f*x))*(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_172():
    f = 1/((a + b*cosh(e + f*x))*(c + d*x)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_173():
    f = (c + d*x)**3/(a + b*cosh(e + f*x))**2
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * Symbol('a') * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_174():
    f = (c + d*x)**2/(a + b*cosh(e + f*x))**2
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_175():
    f = (c + d*x)/(a + b*cosh(e + f*x))**2
    F = ((Symbol('a') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)))) + ((Symbol('d') * sympy.log((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x)) * sympy.sinh((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_176():
    f = 1/((a + b*cosh(e + f*x))**2*(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_177():
    f = 1/((a + b*cosh(e + f*x))**2*(c + d*x)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_178():
    f = (a + b*cosh(e + f*x))**n*(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_179():
    f = (a + b*cosh(e + f*x))**3*(c + d*x)**m
    F = (((Symbol('a'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('e')) + (Integer(-1) * ((Integer(3) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-3) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(3) * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(8) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('e')) + ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-3) * Symbol('e')) + ((Integer(3) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(3) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(8) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_180():
    f = (a + b*cosh(e + f*x))**2*(c + d*x)**m
    F = (((Symbol('a'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + ((Symbol('a') * Symbol('b') * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('e')) + ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_181():
    f = (a + b*cosh(e + f*x))*(c + d*x)**m
    F = ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * ((Integer(-1) * ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('e')) + ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * (((Symbol('f') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_182():
    f = (c + d*x)**m/(a + b*cosh(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_1_c_plus_d_x_pow_m_a_plus_b_cosh_pow_n_183():
    f = (c + d*x)**m/(a + b*cosh(e + f*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F

