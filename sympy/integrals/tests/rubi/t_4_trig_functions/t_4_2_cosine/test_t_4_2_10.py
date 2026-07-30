"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.10 (c+d x)^m (a+b cos)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_1():
    f = (c + d*x)**4*cos(a + b*x)
    F = (c + d*x)**4*sin(a + b*x)/b + 4*d*(c + d*x)**3*cos(a + b*x)/b**2 - 12*d**2*(c + d*x)**2*sin(a + b*x)/b**3 - 24*d**3*(c + d*x)*cos(a + b*x)/b**4 + 24*d**4*sin(a + b*x)/b**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_2():
    f = (c + d*x)**3*cos(a + b*x)
    F = (c + d*x)**3*sin(a + b*x)/b + 3*d*(c + d*x)**2*cos(a + b*x)/b**2 - 6*d**2*(c + d*x)*sin(a + b*x)/b**3 - 6*d**3*cos(a + b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_3():
    f = (c + d*x)**2*cos(a + b*x)
    F = (c + d*x)**2*sin(a + b*x)/b + 2*d*(c + d*x)*cos(a + b*x)/b**2 - 2*d**2*sin(a + b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_4():
    f = (c + d*x)*cos(a + b*x)
    F = (c + d*x)*sin(a + b*x)/b + d*cos(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_5():
    f = cos(a + b*x)/(c + d*x)
    F = ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_6():
    f = cos(a + b*x)/(c + d*x)**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_7():
    f = cos(a + b*x)/(c + d*x)**3
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_8():
    f = cos(a + b*x)/(c + d*x)**4
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_9():
    f = (c + d*x)**4*cos(a + b*x)**2
    F = (c + d*x)**5/(10*d) + (c + d*x)**4*sin(a + b*x)*cos(a + b*x)/(2*b) + d*(c + d*x)**3*cos(a + b*x)**2/b**2 - d*(c + d*x)**3/(2*b**2) - 3*d**2*(c + d*x)**2*sin(a + b*x)*cos(a + b*x)/(2*b**3) + 3*d**4*x/(4*b**4) - 3*d**3*(c + d*x)*cos(a + b*x)**2/(2*b**4) + 3*d**4*sin(a + b*x)*cos(a + b*x)/(4*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_10():
    f = (c + d*x)**3*cos(a + b*x)**2
    F = (c + d*x)**4/(8*d) + (c + d*x)**3*sin(a + b*x)*cos(a + b*x)/(2*b) - 3*c*d**2*x/(4*b**2) - 3*d**3*x**2/(8*b**2) + 3*d*(c + d*x)**2*cos(a + b*x)**2/(4*b**2) - 3*d**2*(c + d*x)*sin(a + b*x)*cos(a + b*x)/(4*b**3) - 3*d**3*cos(a + b*x)**2/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_11():
    f = (c + d*x)**2*cos(a + b*x)**2
    F = (c + d*x)**3/(6*d) + (c + d*x)**2*sin(a + b*x)*cos(a + b*x)/(2*b) - d**2*x/(4*b**2) + d*(c + d*x)*cos(a + b*x)**2/(2*b**2) - d**2*sin(a + b*x)*cos(a + b*x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_12():
    f = (c + d*x)*cos(a + b*x)**2
    F = c*x/2 + d*x**2/4 + (c + d*x)*sin(a + b*x)*cos(a + b*x)/(2*b) + d*cos(a + b*x)**2/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_13():
    f = cos(a + b*x)**2/(c + d*x)
    F = ((sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (sympy.log((Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_14():
    f = cos(a + b*x)**2/(c + d*x)**2
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CosIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_15():
    f = cos(a + b*x)**2/(c + d*x)**3
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_16():
    f = (c + d*x)**4*cos(a + b*x)**3
    F = (c + d*x)**4*sin(a + b*x)*cos(a + b*x)**2/(3*b) + 2*(c + d*x)**4*sin(a + b*x)/(3*b) + 4*d*(c + d*x)**3*cos(a + b*x)**3/(9*b**2) + 8*d*(c + d*x)**3*cos(a + b*x)/(3*b**2) - 4*d**2*(c + d*x)**2*sin(a + b*x)*cos(a + b*x)**2/(9*b**3) - 80*d**2*(c + d*x)**2*sin(a + b*x)/(9*b**3) - 8*d**3*(c + d*x)*cos(a + b*x)**3/(27*b**4) - 160*d**3*(c + d*x)*cos(a + b*x)/(9*b**4) - 8*d**4*sin(a + b*x)**3/(81*b**5) + 488*d**4*sin(a + b*x)/(27*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_17():
    f = (c + d*x)**3*cos(a + b*x)**3
    F = (c + d*x)**3*sin(a + b*x)*cos(a + b*x)**2/(3*b) + 2*(c + d*x)**3*sin(a + b*x)/(3*b) + d*(c + d*x)**2*cos(a + b*x)**3/(3*b**2) + 2*d*(c + d*x)**2*cos(a + b*x)/b**2 - 2*d**2*(c + d*x)*sin(a + b*x)*cos(a + b*x)**2/(9*b**3) - 40*d**2*(c + d*x)*sin(a + b*x)/(9*b**3) - 2*d**3*cos(a + b*x)**3/(27*b**4) - 40*d**3*cos(a + b*x)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_18():
    f = (c + d*x)**2*cos(a + b*x)**3
    F = (c + d*x)**2*sin(a + b*x)*cos(a + b*x)**2/(3*b) + 2*(c + d*x)**2*sin(a + b*x)/(3*b) + 2*d*(c + d*x)*cos(a + b*x)**3/(9*b**2) + 4*d*(c + d*x)*cos(a + b*x)/(3*b**2) + 2*d**2*sin(a + b*x)**3/(27*b**3) - 14*d**2*sin(a + b*x)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_19():
    f = (c + d*x)*cos(a + b*x)**3
    F = (c + d*x)*sin(a + b*x)*cos(a + b*x)**2/(3*b) + (2*c + 2*d*x)*sin(a + b*x)/(3*b) + d*cos(a + b*x)**3/(9*b**2) + 2*d*cos(a + b*x)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_20():
    f = cos(a + b*x)**3/(c + d*x)
    F = ((Integer(3) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_21():
    f = cos(a + b*x)**3/(c + d*x)**2
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('CosIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_22():
    f = cos(a + b*x)**3/(c + d*x)**3
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_23():
    f = x**3*cos(a + b*x)**4
    F = 3*x**4/32 + x**3*sin(a + b*x)*cos(a + b*x)**3/(4*b) + 3*x**3*sin(a + b*x)*cos(a + b*x)/(8*b) + 3*x**2*cos(a + b*x)**4/(16*b**2) + 9*x**2*cos(a + b*x)**2/(16*b**2) - 45*x**2/(128*b**2) - 3*x*sin(a + b*x)*cos(a + b*x)**3/(32*b**3) - 45*x*sin(a + b*x)*cos(a + b*x)/(64*b**3) - 3*cos(a + b*x)**4/(128*b**4) - 45*cos(a + b*x)**2/(128*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_24():
    f = x**2*cos(a + b*x)**4
    F = x**3/8 + x**2*sin(a + b*x)*cos(a + b*x)**3/(4*b) + 3*x**2*sin(a + b*x)*cos(a + b*x)/(8*b) + x*cos(a + b*x)**4/(8*b**2) + 3*x*cos(a + b*x)**2/(8*b**2) - 15*x/(64*b**2) - sin(a + b*x)*cos(a + b*x)**3/(32*b**3) - 15*sin(a + b*x)*cos(a + b*x)/(64*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_25():
    f = x*cos(a + b*x)**4
    F = 3*x**2/16 + x*sin(a + b*x)*cos(a + b*x)**3/(4*b) + 3*x*sin(a + b*x)*cos(a + b*x)/(8*b) + cos(a + b*x)**4/(16*b**2) + 3*cos(a + b*x)**2/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_26():
    f = cos(a + b*x)**4/x
    F = ((Integer(2))**(Integer(-1)) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(8))**(Integer(-1)) * sympy.cos((Integer(4) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(4) * Symbol('b') * x))) + ((Integer(3) * sympy.log(x)) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sin((Integer(4) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(4) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_27():
    f = cos(a + b*x)**4/x**2
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(4)) * (x)**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sin((Integer(2) * Symbol('a'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CosIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sin((Integer(4) * Symbol('a'))))) + (Integer(-1) * (Symbol('b') * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cos((Integer(4) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(4) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_28():
    f = cos(a + b*x)**4/x**3
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(4)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.cos((Integer(4) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(4) * Symbol('b') * x)))) + ((Integer(2) * Symbol('b') * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) + ((Symbol('b'))**(Integer(2)) * sympy.sin((Integer(4) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_29():
    f = (c + d*x)**3*sec(a + b*x)
    F = (Integer(-1) * ((Integer(2) * sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_30():
    f = (c + d*x)**2*sec(a + b*x)
    F = (Integer(-1) * ((Integer(2) * sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_31():
    f = (c + d*x)*sec(a + b*x)
    F = (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_32():
    f = sec(a + b*x)/(c + d*x)
    F = sympy.Function('Unintegrable')((sympy.sec((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_33():
    f = (c + d*x)**3*sec(a + b*x)**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tan((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_34():
    f = (c + d*x)**2*sec(a + b*x)**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tan((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_35():
    f = (c + d*x)*sec(a + b*x)**2
    F = (c + d*x)*tan(a + b*x)/b + d*log(cos(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_36():
    f = sec(a + b*x)**2/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_37():
    f = (c + d*x)**3*sec(a + b*x)**3
    F = (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sec((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.sec((Symbol('a') + (Symbol('b') * x))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_38():
    f = (c + d*x)**2*sec(a + b*x)**3
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.atanh(sympy.sin((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.sec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sec((Symbol('a') + (Symbol('b') * x))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_39():
    f = (c + d*x)*sec(a + b*x)**3
    F = (Integer(-1) * ((sympy.I * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.sec((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sec((Symbol('a') + (Symbol('b') * x))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_40():
    f = sec(a + b*x)**2/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_41():
    f = (c + d*x)**(sympy.S(5)/2)*cos(a + b*x)
    F = ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_42():
    f = (c + d*x)**(sympy.S(3)/2)*cos(a + b*x)
    F = ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_43():
    f = sqrt(c + d*x)*cos(a + b*x)
    F = (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_44():
    f = cos(a + b*x)/sqrt(c + d*x)
    F = ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_45():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_46():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_47():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_48():
    f = (c + d*x)**(sympy.S(5)/2)*cos(a + b*x)**2
    F = (Integer(-1) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * Symbol('d')))**(Integer(-1))) + ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(128) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(15) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(128) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(64) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_49():
    f = (c + d*x)**(sympy.S(3)/2)*cos(a + b*x)**2
    F = (Integer(-1) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_50():
    f = sqrt(c + d*x)*cos(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_51():
    f = cos(a + b*x)**2/sqrt(c + d*x)
    F = (sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_52():
    f = cos(a + b*x)**2/(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_53():
    f = cos(a + b*x)**2/(c + d*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(3) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_54():
    f = cos(a + b*x)**2/(c + d*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(15) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(32) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(15) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(15) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_55():
    f = cos(a + b*x)**2/(c + d*x)**(sympy.S(9)/2)
    F = (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(105) * (Symbol('d'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(7) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(105) * (Symbol('d'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(128) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(105) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(105) * (Symbol('d'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(35) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(105) * (Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_56():
    f = (c + d*x)**(sympy.S(5)/2)*cos(a + b*x)**3
    F = ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(18) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(45) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(144) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(144) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(45) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(45) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(144) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_57():
    f = (c + d*x)**(sympy.S(3)/2)*cos(a + b*x)**3
    F = ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(24) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(24) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(9) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_58():
    f = sqrt(c + d*x)*cos(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(12) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(12) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_59():
    f = cos(a + b*x)**3/sqrt(c + d*x)
    F = ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_60():
    f = cos(a + b*x)**3/(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_61():
    f = cos(a + b*x)**3/(c + d*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_62():
    f = cos(a + b*x)**3/(c + d*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(24) * (Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(5) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.cos(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(5) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_63():
    f = x**(sympy.S(3)/2)*cos(x)
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt(x) * sympy.cos(x)) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(x))))) + ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_64():
    f = sqrt(x)*cos(x)
    F = ((Integer(-1) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1))))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(x)))) + (sympy.sqrt(x) * sympy.sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_65():
    f = cos(x)/sqrt(x)
    F = sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_66():
    f = cos(x)/x**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos(x)) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(x)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_67():
    f = (c + d*x)**(sympy.S(4)/3)*cos(a + b*x)
    F = ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('d'))**(Integer(2)) * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(9) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('d'))**(Integer(2)) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(9) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_68():
    f = (c + d*x)**(sympy.S(2)/3)*cos(a + b*x)
    F = ((Symbol('d') * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('d') * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_69():
    f = (c + d*x)**(sympy.S(1)/3)*cos(a + b*x)
    F = ((Symbol('d') * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('d') * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(6) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_70():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(1)/3)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(2) * Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_71():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(2)/3)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(2) * Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_72():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(4)/3)
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_73():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(5)/3)
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_74():
    f = cos(a + b*x)/(c + d*x)**(sympy.S(7)/3)
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(9) * sympy.I * Symbol('b') * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.I * Symbol('b') * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (Integer(8) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(9) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_75():
    f = x*sqrt(cos(a + b*x))
    F = sympy.Function('Unintegrable')((x * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_76():
    f = sqrt(cos(a + b*x))
    F = 2*elliptic_e(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_77():
    f = sqrt(cos(a + b*x))/x
    F = sympy.Function('Unintegrable')((sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_78():
    f = x*cos(a + b*x)**(sympy.S(3)/2)
    F = ((Integer(4) * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * sympy.Function('Unintegrable')((x * (sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)) + ((Integer(2) * x * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x)))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_79():
    f = cos(a + b*x)**(sympy.S(3)/2)
    F = 2*sin(a + b*x)*sqrt(cos(a + b*x))/(3*b) + 2*elliptic_f(a/2 + b*x/2, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_80():
    f = cos(a + b*x)**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_81():
    f = x*cos(a + b*x)**(sympy.S(3)/2) - x/(3*sqrt(cos(a + b*x)))
    F = 2*x*sin(a + b*x)*sqrt(cos(a + b*x))/(3*b) + 4*cos(a + b*x)**(sympy.S(3)/2)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_82():
    f = cos(x)**(sympy.S(3)/2)/x**3
    F = (Integer(-1) * ((sympy.cos(x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.cos(x))))**(Integer(-1)), x)) + (Integer(-1) * ((Integer(9) * (Integer(8))**(Integer(-1))) * sympy.Function('Unintegrable')(((sympy.cos(x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x))) + ((Integer(3) * sympy.sqrt(sympy.cos(x)) * sympy.sin(x)) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_83():
    f = x/sqrt(cos(a + b*x))
    F = sympy.Function('Unintegrable')((x * (sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_84():
    f = 1/sqrt(cos(a + b*x))
    F = 2*elliptic_f(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_85():
    f = 1/(x*sqrt(cos(a + b*x)))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_86():
    f = x/cos(a + b*x)**(sympy.S(3)/2)
    F = ((Integer(4) * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))) + (Integer(-1) * sympy.Function('Unintegrable')((x * sympy.sqrt(sympy.cos((Symbol('a') + (Symbol('b') * x))))), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_87():
    f = cos(a + b*x)**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)/(b*sqrt(cos(a + b*x))) - 2*elliptic_e(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_88():
    f = 1/(x*cos(a + b*x)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.cos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_89():
    f = x*sqrt(cos(a + b*x)) + x/cos(a + b*x)**(sympy.S(3)/2)
    F = 2*x*sin(a + b*x)/(b*sqrt(cos(a + b*x))) + 4*sqrt(cos(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_90():
    f = x*sqrt(cos(x)) + x/cos(x)**(sympy.S(3)/2)
    F = 2*x*sin(x)/sqrt(cos(x)) + 4*sqrt(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_91():
    f = -x/(3*sqrt(cos(x))) + x/cos(x)**(sympy.S(5)/2)
    F = 2*x*sin(x)/(3*cos(x)**(sympy.S(3)/2)) - 4/(3*sqrt(cos(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_92():
    f = 3*x*sqrt(cos(x))/5 + x/cos(x)**(sympy.S(7)/2)
    F = 6*x*sin(x)/(5*sqrt(cos(x))) + 2*x*sin(x)/(5*cos(x)**(sympy.S(5)/2)) + 12*sqrt(cos(x))/5 - 4/(15*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_93():
    f = x**2*sqrt(cos(x)) + x**2/cos(x)**(sympy.S(3)/2)
    F = 2*x**2*sin(x)/sqrt(cos(x)) + 8*x*sqrt(cos(x)) - 16*elliptic_e(x/2, 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_94():
    f = -x*sqrt(sec(x))/3 + x/sec(x)**(sympy.S(3)/2)
    F = 2*x*sin(x)/(3*sqrt(sec(x))) + 4/(9*sec(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_95():
    f = -3*x/(5*sqrt(sec(x))) + x/sec(x)**(sympy.S(5)/2)
    F = 2*x*sin(x)/(5*sec(x)**(sympy.S(3)/2)) + 4/(25*sec(x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_96():
    f = -5*x*sqrt(sec(x))/21 + x/sec(x)**(sympy.S(7)/2)
    F = 10*x*sin(x)/(21*sqrt(sec(x))) + 2*x*sin(x)/(7*sec(x)**(sympy.S(5)/2)) + 20/(63*sec(x)**(sympy.S(3)/2)) + 4/(49*sec(x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_97():
    f = -x**2*sqrt(sec(x))/3 + x**2/sec(x)**(sympy.S(3)/2)
    F = 2*x**2*sin(x)/(3*sqrt(sec(x))) + 8*x/(9*sec(x)**(sympy.S(3)/2)) - 16*sin(x)/(27*sqrt(sec(x))) - 16*sqrt(cos(x))*elliptic_f(x/2, 2)*sqrt(sec(x))/27
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_98():
    f = (b*cos(e + f*x))**n*(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_99():
    f = (c + d*x)**m*cos(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + ((Integer(3) * sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(3) * sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + ((sympy.I * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(3) * sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((Integer(3) * sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_100():
    f = (c + d*x)**m*cos(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * ((Integer(2) * Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Integer(2) * sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(2) * sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_101():
    f = (c + d*x)**m*cos(a + b*x)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))) + ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * (((sympy.I * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_102():
    f = (c + d*x)**m*sec(a + b*x)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.sec((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_103():
    f = (c + d*x)**m*sec(a + b*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.sec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_104():
    f = x**(m + 3)*cos(a + b*x)
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(4)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_105():
    f = x**(m + 2)*cos(a + b*x)
    F = ((sympy.I * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(3)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_106():
    f = x**(m + 1)*cos(a + b*x)
    F = (((sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_107():
    f = x**m*cos(a + b*x)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))) + ((sympy.I * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_108():
    f = x**(m - 1)*cos(a + b*x)
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * ((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_109():
    f = x**(m - 2)*cos(a + b*x)
    F = (((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * ((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_110():
    f = x**(m - 3)*cos(a + b*x)
    F = (((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * ((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')))**(Integer(-1))) + (((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (sympy.I * Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_111():
    f = x**(m + 3)*cos(a + b*x)**2
    F = ((x)**((Integer(4) + Symbol('m'))) * ((Integer(2) * (Integer(4) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-6) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-6) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_112():
    f = x**(m + 2)*cos(a + b*x)**2
    F = ((x)**((Integer(3) + Symbol('m'))) * ((Integer(2) * (Integer(3) + Symbol('m'))))**(Integer(-1))) + ((sympy.I * (Integer(2))**((Integer(-5) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(2))**((Integer(-5) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_113():
    f = x**(m + 1)*cos(a + b*x)**2
    F = ((x)**((Integer(2) + Symbol('m'))) * ((Integer(2) * (Integer(2) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_114():
    f = x**m*cos(a + b*x)**2
    F = ((x)**((Integer(1) + Symbol('m'))) * ((Integer(2) * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_115():
    f = x**(m - 1)*cos(a + b*x)**2
    F = ((x)**(Symbol('m')) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Integer(-2) * sympy.I * Symbol('b') * x))) * ((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_116():
    f = x**(m - 2)*cos(a + b*x)**2
    F = (Integer(-1) * ((x)**((Integer(-1) + Symbol('m'))) * ((Integer(2) * (Integer(1) + (Integer(-1) * Symbol('m')))))**(Integer(-1)))) + ((sympy.I * (Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * ((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * Symbol('b') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_117():
    f = x**(m - 3)*cos(a + b*x)**2
    F = (Integer(-1) * ((x)**((Integer(-2) + Symbol('m'))) * ((Integer(2) * (Integer(2) + (Integer(-1) * Symbol('m')))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Integer(-2) * sympy.I * Symbol('b') * x))) * (((Integer(2))**(Symbol('m')) * (((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-2) + Symbol('m')), (Integer(2) * sympy.I * Symbol('b') * x))) * (((Integer(2))**(Symbol('m')) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * x))**(Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_118():
    f = (c + d*x)**3*(a*cos(e + f*x) + a)
    F = -6*a*d**3*cos(e + f*x)/f**4 - 6*a*d**2*(c + d*x)*sin(e + f*x)/f**3 + 3*a*d*(c + d*x)**2*cos(e + f*x)/f**2 + a*(c + d*x)**3*sin(e + f*x)/f + a*(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_119():
    f = (c + d*x)**2*(a*cos(e + f*x) + a)
    F = -2*a*d**2*sin(e + f*x)/f**3 + 2*a*d*(c + d*x)*cos(e + f*x)/f**2 + a*(c + d*x)**2*sin(e + f*x)/f + a*(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_120():
    f = (c + d*x)*(a*cos(e + f*x) + a)
    F = a*d*cos(e + f*x)/f**2 + a*(c + d*x)*sin(e + f*x)/f + a*(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_121():
    f = (a*cos(e + f*x) + a)/(c + d*x)
    F = ((Symbol('a') * sympy.cos((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sin((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_122():
    f = (a*cos(e + f*x) + a)/(c + d*x)**2
    F = (Integer(-1) * (Symbol('a') * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('f') * sympy.Function('CosIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sin((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('f') * sympy.cos((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_123():
    f = (c + d*x)**3*(a*cos(e + f*x) + a)**2
    F = -3*a**2*c*d**2*x/(4*f**2) - 3*a**2*d**3*x**2/(8*f**2) - 3*a**2*d**3*cos(e + f*x)**2/(8*f**4) - 12*a**2*d**3*cos(e + f*x)/f**4 - 3*a**2*d**2*(c + d*x)*sin(e + f*x)*cos(e + f*x)/(4*f**3) - 12*a**2*d**2*(c + d*x)*sin(e + f*x)/f**3 + 3*a**2*d*(c + d*x)**2*cos(e + f*x)**2/(4*f**2) + 6*a**2*d*(c + d*x)**2*cos(e + f*x)/f**2 + a**2*(c + d*x)**3*sin(e + f*x)*cos(e + f*x)/(2*f) + 2*a**2*(c + d*x)**3*sin(e + f*x)/f + 3*a**2*(c + d*x)**4/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_124():
    f = (c + d*x)**2*(a*cos(e + f*x) + a)**2
    F = -a**2*d**2*x/(4*f**2) - a**2*d**2*sin(e + f*x)*cos(e + f*x)/(4*f**3) - 4*a**2*d**2*sin(e + f*x)/f**3 + a**2*d*(c + d*x)*cos(e + f*x)**2/(2*f**2) + 4*a**2*d*(c + d*x)*cos(e + f*x)/f**2 + a**2*(c + d*x)**2*sin(e + f*x)*cos(e + f*x)/(2*f) + 2*a**2*(c + d*x)**2*sin(e + f*x)/f + a**2*(c + d*x)**3/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_125():
    f = (c + d*x)*(a*cos(e + f*x) + a)**2
    F = a**2*c*x/2 + a**2*d*x**2/4 + a**2*d*cos(e + f*x)**2/(4*f**2) + 2*a**2*d*cos(e + f*x)/f**2 + a**2*(c + d*x)*sin(e + f*x)*cos(e + f*x)/(2*f) + 2*a**2*(c + d*x)*sin(e + f*x)/f + a**2*(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_126():
    f = (a*cos(e + f*x) + a)**2/(c + d*x)
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.cos((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.cos(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_127():
    f = (a*cos(e + f*x) + a)**2/(c + d*x)**2
    F = (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (sympy.cos(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(4))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.Function('CosIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x))) * sympy.sin(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.Function('CosIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x))) * sympy.sin((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.cos((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * sympy.cos(((Integer(2) * Symbol('e')) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(2) * Symbol('f') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_128():
    f = (c + d*x)**3/(a*cos(e + f*x) + a)
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(12) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_129():
    f = (c + d*x)**2/(a*cos(e + f*x) + a)
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + ((Integer(4) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_130():
    f = (c + d*x)/(a*cos(e + f*x) + a)
    F = 2*d*log(cos(e/2 + f*x/2))/(a*f**2) + (c + d*x)*tan(e/2 + f*x/2)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_131():
    f = 1/((c + d*x)*(a*cos(e + f*x) + a))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * (Symbol('a') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_132():
    f = 1/((c + d*x)**2*(a*cos(e + f*x) + a))
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_133():
    f = (c + d*x)**3/(a*cos(e + f*x) + a)**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * (Symbol('d'))**(Integer(3)) * sympy.log(sympy.cos(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.sec(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.sec(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_134():
    f = (c + d*x)**2/(a*cos(e + f*x) + a)**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((Integer(4) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * (sympy.sec(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.sec(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1))))))**(Integer(2)) * sympy.tan(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_135():
    f = (c + d*x)/(a*cos(e + f*x) + a)**2
    F = 2*d*log(cos(e/2 + f*x/2))/(3*a**2*f**2) - d*sec(e/2 + f*x/2)**2/(6*a**2*f**2) + (c + d*x)*tan(e/2 + f*x/2)*sec(e/2 + f*x/2)**2/(6*a**2*f) + (c + d*x)*tan(e/2 + f*x/2)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_136():
    f = 1/((c + d*x)*(a*cos(e + f*x) + a)**2)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_137():
    f = 1/((c + d*x)**2*(a*cos(e + f*x) + a)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_138():
    f = (c + d*x)**3/(-a*cos(e + f*x) + a)
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.cot(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x)))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1)))) + ((Integer(12) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_139():
    f = (c + d*x)**2/(-a*cos(e + f*x) + a)
    F = (Integer(-1) * ((sympy.I * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.cot(((Symbol('e') * (Integer(2))**(Integer(-1))) + ((Symbol('f') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + ((Integer(4) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x)))))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_140():
    f = (c + d*x)/(-a*cos(e + f*x) + a)
    F = 2*d*log(sin(e/2 + f*x/2))/(a*f**2) - (c + d*x)*cot(e/2 + f*x/2)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_141():
    f = 1/((c + d*x)*(-a*cos(e + f*x) + a))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)) * (Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_142():
    f = 1/((c + d*x)**2*(-a*cos(e + f*x) + a))
    F = sympy.Function('Unintegrable')(((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_143():
    f = x**3*sqrt(a*cos(c + d*x) + a)
    F = 2*x**3*sqrt(a*cos(c + d*x) + a)*tan(c/2 + d*x/2)/d + 12*x**2*sqrt(a*cos(c + d*x) + a)/d**2 - 48*x*sqrt(a*cos(c + d*x) + a)*tan(c/2 + d*x/2)/d**3 - 96*sqrt(a*cos(c + d*x) + a)/d**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_144():
    f = x**2*sqrt(a*cos(c + d*x) + a)
    F = 2*x**2*sqrt(a*cos(c + d*x) + a)*tan(c/2 + d*x/2)/d + 8*x*sqrt(a*cos(c + d*x) + a)/d**2 - 16*sqrt(a*cos(c + d*x) + a)*tan(c/2 + d*x/2)/d**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_145():
    f = x*sqrt(a*cos(c + d*x) + a)
    F = 2*x*sqrt(a*cos(c + d*x) + a)*tan(c/2 + d*x/2)/d + 4*sqrt(a*cos(c + d*x) + a)/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_146():
    f = sqrt(a*cos(c + d*x) + a)
    F = 2*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_147():
    f = sqrt(a*cos(c + d*x) + a)/x
    F = (sympy.cos((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CosIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_148():
    f = sqrt(a*cos(c + d*x) + a)/x**2
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CosIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((Symbol('c') * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('d') * sympy.cos((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('SinIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_149():
    f = sqrt(a*cos(c + d*x) + a)/x**3
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.Function('CosIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1)))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))))) + ((Integer(8))**(Integer(-1)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((Symbol('c') * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) + ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_150():
    f = x**3*sqrt(a*cos(x) + a)
    F = 2*x**3*sqrt(a*cos(x) + a)*tan(x/2) + 12*x**2*sqrt(a*cos(x) + a) - 48*x*sqrt(a*cos(x) + a)*tan(x/2) - 96*sqrt(a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_151():
    f = x**2*sqrt(a*cos(x) + a)
    F = 2*x**2*sqrt(a*cos(x) + a)*tan(x/2) + 8*x*sqrt(a*cos(x) + a) - 16*sqrt(a*cos(x) + a)*tan(x/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_152():
    f = x*sqrt(a*cos(x) + a)
    F = 2*x*sqrt(a*cos(x) + a)*tan(x/2) + 4*sqrt(a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_153():
    f = sqrt(a*cos(x) + a)
    F = 2*a*sin(x)/sqrt(a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_154():
    f = sqrt(a*cos(x) + a)/x
    F = sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_155():
    f = sqrt(a*cos(x) + a)/x**2
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')((x * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_156():
    f = sqrt(a*cos(x) + a)/x**3
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))))) + ((sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.tan((x * (Integer(2))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_157():
    f = x**3*sqrt(-a*cos(x) + a)
    F = -2*x**3*sqrt(-a*cos(x) + a)*cot(x/2) + 12*x**2*sqrt(-a*cos(x) + a) + 48*x*sqrt(-a*cos(x) + a)*cot(x/2) - 96*sqrt(-a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_158():
    f = x**2*sqrt(-a*cos(x) + a)
    F = -2*x**2*sqrt(-a*cos(x) + a)*cot(x/2) + 8*x*sqrt(-a*cos(x) + a) + 16*sqrt(-a*cos(x) + a)*cot(x/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_159():
    f = x*sqrt(-a*cos(x) + a)
    F = -2*x*sqrt(-a*cos(x) + a)*cot(x/2) + 4*sqrt(-a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_160():
    f = sqrt(-a*cos(x) + a)
    F = -2*a*sin(x)/sqrt(-a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_161():
    f = sqrt(-a*cos(x) + a)/x
    F = sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * sympy.csc((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')((x * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_162():
    f = sqrt(-a*cos(x) + a)/x**2
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * sympy.Function('CosIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.csc((x * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_163():
    f = sqrt(-a*cos(x) + a)/x**3
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * sympy.cot((x * (Integer(2))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))) * sympy.csc((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')((x * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_164():
    f = x**3*(a*cos(x) + a)**(sympy.S(3)/2)
    F = 4*a*x**3*sqrt(a*cos(x) + a)*sin(x/2)*cos(x/2)/3 + 8*a*x**3*sqrt(a*cos(x) + a)*tan(x/2)/3 + 8*a*x**2*sqrt(a*cos(x) + a)*cos(x/2)**2/3 + 16*a*x**2*sqrt(a*cos(x) + a) - 32*a*x*sqrt(a*cos(x) + a)*sin(x/2)*cos(x/2)/9 - 640*a*x*sqrt(a*cos(x) + a)*tan(x/2)/9 - 64*a*sqrt(a*cos(x) + a)*cos(x/2)**2/27 - 1280*a*sqrt(a*cos(x) + a)/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_165():
    f = x**2*(a*cos(x) + a)**(sympy.S(3)/2)
    F = 4*a*x**2*sqrt(a*cos(x) + a)*sin(x/2)*cos(x/2)/3 + 8*a*x**2*sqrt(a*cos(x) + a)*tan(x/2)/3 + 16*a*x*sqrt(a*cos(x) + a)*cos(x/2)**2/9 + 32*a*x*sqrt(a*cos(x) + a)/3 + 32*a*sqrt(a*cos(x) + a)*sin(x/2)**2*tan(x/2)/27 - 224*a*sqrt(a*cos(x) + a)*tan(x/2)/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_166():
    f = x*(a*cos(x) + a)**(sympy.S(3)/2)
    F = 4*a*x*sqrt(a*cos(x) + a)*sin(x/2)*cos(x/2)/3 + 8*a*x*sqrt(a*cos(x) + a)*tan(x/2)/3 + 8*a*sqrt(a*cos(x) + a)*cos(x/2)**2/9 + 16*a*sqrt(a*cos(x) + a)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_167():
    f = (a*cos(x) + a)**(sympy.S(3)/2)/x
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_168():
    f = (a*cos(x) + a)**(sympy.S(3)/2)/x**2
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * (sympy.cos((x * (Integer(2))**(Integer(-1)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')((x * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_169():
    f = (a*cos(x) + a)**(sympy.S(3)/2)/x**3
    F = (Integer(-1) * ((Symbol('a') * (sympy.cos((x * (Integer(2))**(Integer(-1)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))) * ((x)**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')((x * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * ((Integer(9) * (Integer(16))**(Integer(-1))) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.Function('CosIntegral')(((Integer(3) * x) * (Integer(2))**(Integer(-1)))) * sympy.sec((x * (Integer(2))**(Integer(-1)))))) + ((Integer(3) * Symbol('a') * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x)))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_170():
    f = x**3/sqrt(a*cos(c + d*x) + a)
    F = (Integer(-1) * ((Integer(4) * sympy.I * (x)**(Integer(3)) * sympy.atan((sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(12) * sympy.I * (x)**(Integer(2)) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * sympy.I * (x)**(Integer(2)) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(48) * x * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(48) * x * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(96) * sympy.I * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(96) * sympy.I * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_171():
    f = x**2/sqrt(a*cos(c + d*x) + a)
    F = (Integer(-1) * ((Integer(4) * sympy.I * (x)**(Integer(2)) * sympy.atan((sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(8) * sympy.I * x * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * sympy.I * x * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(16) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_172():
    f = x/sqrt(a*cos(c + d*x) + a)
    F = (Integer(-1) * ((Integer(4) * sympy.I * x * sympy.atan((sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1)))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(4) * sympy.I * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * sympy.cos(((Symbol('c') * (Integer(2))**(Integer(-1))) + ((Symbol('d') * x) * (Integer(2))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_173():
    f = 1/sqrt(a*cos(c + d*x) + a)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_174():
    f = 1/(x*sqrt(a*cos(c + d*x) + a))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_175():
    f = x**3/sqrt(-a*cos(x) + a)
    F = (Integer(-1) * ((Integer(4) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(12) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(48) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(48) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(96) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(96) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_176():
    f = x**2/sqrt(-a*cos(x) + a)
    F = (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(8) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(16) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_177():
    f = x/sqrt(-a*cos(x) + a)
    F = (Integer(-1) * ((Integer(4) * x * sympy.atanh((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(4) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1)))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.sin((x * (Integer(2))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_178():
    f = 1/sqrt(-a*cos(x) + a)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(x)/(2*sqrt(-a*cos(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_179():
    f = 1/(x*sqrt(-a*cos(x) + a))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.cos(x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_180():
    f = x**3/(a*cos(x) + a)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24) * sympy.I * x * sympy.atan((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.cos((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (x)**(Integer(3)) * sympy.atan((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.cos((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(24) * sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * x * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(12) * x * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(24) * sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.tan((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_181():
    f = x**2/(a*cos(x) + a)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * x) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.atan((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.cos((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(4) * sympy.atanh(sympy.sin((x * (Integer(2))**(Integer(-1))))) * sympy.cos((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + ((Integer(2) * sympy.I * x * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * x * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((Integer(4) * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.tan((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_182():
    f = x/(a*cos(x) + a)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.atan((sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))) * sympy.cos((x * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.cos((x * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((sympy.I * x) * (Integer(2))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))) + ((x * sympy.tan((x * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('a') * sympy.cos(x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_183():
    f = 1/(x*(a*cos(x) + a)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('a') * sympy.cos(x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_184():
    f = (a*cos(c + d*x) + a)**(sympy.S(1)/3)/x
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_185():
    f = x**3/(a + b*cos(x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * x))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_186():
    f = x**2/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_187():
    f = x/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((sympy.I * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_188():
    f = 1/(x*(a + b*cos(x)))
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + (Symbol('b') * sympy.cos(x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_10_c_plus_d_x_pow_m_a_plus_b_cos_pow_n_189():
    f = (e + f*x)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((sympy.I * Symbol('a') * (Symbol('e') + (Symbol('f') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('a') * (Symbol('e') + (Symbol('f') * x)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * sympy.log((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('f') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * Symbol('f') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e') + (Symbol('f') * x)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F

