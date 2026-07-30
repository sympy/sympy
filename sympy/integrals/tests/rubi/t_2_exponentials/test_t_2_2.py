"""Generated from MathematicaSyntaxTestSuite.

Source: 2 Exponentials/2.2 (c+d x)^m (F^(g (e+f x)))^n (a+b (F^(g (e+f x)))^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

F, a, b, c, d, e, f, g, m, n, p = symbols('F a b c d e f g m n p')

def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_1():
    f = x**3/(a + b*exp(c + d*x))
    F = ((x)**(Integer(4)) * ((Integer(4) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_2():
    f = x**2/(a + b*exp(c + d*x))
    F = ((x)**(Integer(3)) * ((Integer(3) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_3():
    f = x/(a + b*exp(c + d*x))
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1))))) * ((Symbol('a') * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_4():
    f = 1/(a + b*exp(c + d*x))
    F = x/a - log(a + b*exp(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_5():
    f = 1/(x*(a + b*exp(c + d*x)))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))) * x))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_6():
    f = 1/(x**2*(a + b*exp(c + d*x)))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))) * (x)**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_7():
    f = 1/(a + b*exp(c - d*x))
    F = x/a + log(a + b*exp(c - d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_8():
    f = 1/(a + b*exp(-c - d*x))
    F = x/a + log(a + b*exp(-c - d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_9():
    f = x**3/(a + b*exp(c + d*x))**2
    F = (Integer(-1) * ((x)**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Symbol('a') * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_10():
    f = x**2/(a + b*exp(c + d*x))**2
    F = (Integer(-1) * ((x)**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Symbol('a') * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_11():
    f = x/(a + b*exp(c + d*x))**2
    F = (Integer(-1) * (x * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (x * ((Symbol('a') * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_12():
    f = (a + b*exp(c + d*x))**(-2)
    F = 1/(a*d*(a + b*exp(c + d*x))) + x/a**2 - log(a + b*exp(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_13():
    f = 1/(x*(a + b*exp(c + d*x))**2)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * x))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_14():
    f = 1/(x**2*(a + b*exp(c + d*x))**2)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_15():
    f = (a + b*exp(c - d*x))**(-2)
    F = -1/(a*d*(a + b*exp(c - d*x))) + x/a**2 + log(a + b*exp(c - d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_16():
    f = (a + b*exp(-c - d*x))**(-2)
    F = -1/(a*d*(a + b*exp(-c - d*x))) + x/a**2 + log(a + b*exp(-c - d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_17():
    f = x**3/(a + b*exp(c + d*x))**3
    F = ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(9) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_18():
    f = x**2/(a + b*exp(c + d*x))**3
    F = (x * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (x * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_19():
    f = x/(a + b*exp(c + d*x))**3
    F = (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (x * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (x * (((Symbol('a'))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_20():
    f = (a + b*exp(c + d*x))**(-3)
    F = 1/(2*a*d*(a + b*exp(c + d*x))**2) + 1/(a**2*d*(a + b*exp(c + d*x))) + x/a**3 - log(a + b*exp(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_21():
    f = 1/(x*(a + b*exp(c + d*x))**3)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * x))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_22():
    f = 1/(x**2*(a + b*exp(c + d*x))**3)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * (x)**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_23():
    f = (a + b*exp(c - d*x))**(-3)
    F = -1/(2*a*d*(a + b*exp(c - d*x))**2) - 1/(a**2*d*(a + b*exp(c - d*x))) + x/a**3 + log(a + b*exp(c - d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_24():
    f = (a + b*exp(-c - d*x))**(-3)
    F = -1/(2*a*d*(a + b*exp(-c - d*x))**2) - 1/(a**2*d*(a + b*exp(-c - d*x))) + x/a**3 + log(a + b*exp(-c - d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_25():
    f = (a + b*(F**(g*(e + f*x)))**n)*(c + d*x)**3
    F = a*(c + d*x)**4/(4*d) - 6*b*d**3*(F**(e*g + f*g*x))**n/(f**4*g**4*n**4*log(F)**4) + 6*b*d**2*(c + d*x)*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 3*b*d*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + b*(c + d*x)**3*(F**(e*g + f*g*x))**n/(f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_26():
    f = (a + b*(F**(g*(e + f*x)))**n)*(c + d*x)**2
    F = a*(c + d*x)**3/(3*d) + 2*b*d**2*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 2*b*d*(c + d*x)*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + b*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_27():
    f = (a + b*(F**(g*(e + f*x)))**n)*(c + d*x)
    F = a*(c + d*x)**2/(2*d) - b*d*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + b*(c + d*x)*(F**(e*g + f*g*x))**n/(f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_28():
    f = a + b*(F**(g*(e + f*x)))**n
    F = a*x + b*(F**(g*(e + f*x)))**n/(f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_29():
    f = (a + b*(F**(g*(e + f*x)))**n)/(c + d*x)
    F = ((Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_30():
    f = (a + b*(F**(g*(e + f*x)))**n)/(c + d*x)**2
    F = (Integer(-1) * (Symbol('a') * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * Symbol('f') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_31():
    f = (a + b*(F**(g*(e + f*x)))**n)/(c + d*x)**3
    F = (Integer(-1) * (Symbol('a') * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_32():
    f = (a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x)**3
    F = a**2*(c + d*x)**4/(4*d) - 12*a*b*d**3*(F**(e*g + f*g*x))**n/(f**4*g**4*n**4*log(F)**4) + 12*a*b*d**2*(c + d*x)*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 6*a*b*d*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 2*a*b*(c + d*x)**3*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) - 3*b**2*d**3*(F**(e*g + f*g*x))**(2*n)/(8*f**4*g**4*n**4*log(F)**4) + 3*b**2*d**2*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(4*f**3*g**3*n**3*log(F)**3) - 3*b**2*d*(c + d*x)**2*(F**(e*g + f*g*x))**(2*n)/(4*f**2*g**2*n**2*log(F)**2) + b**2*(c + d*x)**3*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_33():
    f = (a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x)**2
    F = a**2*(c + d*x)**3/(3*d) + 4*a*b*d**2*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 4*a*b*d*(c + d*x)*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 2*a*b*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) + b**2*d**2*(F**(e*g + f*g*x))**(2*n)/(4*f**3*g**3*n**3*log(F)**3) - b**2*d*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(2*f**2*g**2*n**2*log(F)**2) + b**2*(c + d*x)**2*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_34():
    f = (a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x)
    F = a**2*(c + d*x)**2/(2*d) - 2*a*b*d*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 2*a*b*(c + d*x)*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) - b**2*d*(F**(e*g + f*g*x))**(2*n)/(4*f**2*g**2*n**2*log(F)**2) + b**2*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_35():
    f = (a + b*(F**(g*(e + f*x)))**n)**2
    F = a**2*x + 2*a*b*(F**(g*(e + f*x)))**n/(f*g*n*log(F)) + b**2*(F**(g*(e + f*x)))**(2*n)/(2*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_36():
    f = (a + b*(F**(g*(e + f*x)))**n)**2/(c + d*x)
    F = ((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_37():
    f = (a + b*(F**(g*(e + f*x)))**n)**2/(c + d*x)**2
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_38():
    f = (a + b*(F**(g*(e + f*x)))**n)**2/(c + d*x)**3
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n')))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('a') * Symbol('b') * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_39():
    f = (a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x)**3
    F = a**3*(c + d*x)**4/(4*d) - 18*a**2*b*d**3*(F**(e*g + f*g*x))**n/(f**4*g**4*n**4*log(F)**4) + 18*a**2*b*d**2*(c + d*x)*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 9*a**2*b*d*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 3*a**2*b*(c + d*x)**3*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) - 9*a*b**2*d**3*(F**(e*g + f*g*x))**(2*n)/(8*f**4*g**4*n**4*log(F)**4) + 9*a*b**2*d**2*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(4*f**3*g**3*n**3*log(F)**3) - 9*a*b**2*d*(c + d*x)**2*(F**(e*g + f*g*x))**(2*n)/(4*f**2*g**2*n**2*log(F)**2) + 3*a*b**2*(c + d*x)**3*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F)) - 2*b**3*d**3*(F**(e*g + f*g*x))**(3*n)/(27*f**4*g**4*n**4*log(F)**4) + 2*b**3*d**2*(c + d*x)*(F**(e*g + f*g*x))**(3*n)/(9*f**3*g**3*n**3*log(F)**3) - b**3*d*(c + d*x)**2*(F**(e*g + f*g*x))**(3*n)/(3*f**2*g**2*n**2*log(F)**2) + b**3*(c + d*x)**3*(F**(e*g + f*g*x))**(3*n)/(3*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_40():
    f = (a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x)**2
    F = a**3*(c + d*x)**3/(3*d) + 6*a**2*b*d**2*(F**(e*g + f*g*x))**n/(f**3*g**3*n**3*log(F)**3) - 6*a**2*b*d*(c + d*x)*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 3*a**2*b*(c + d*x)**2*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) + 3*a*b**2*d**2*(F**(e*g + f*g*x))**(2*n)/(4*f**3*g**3*n**3*log(F)**3) - 3*a*b**2*d*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(2*f**2*g**2*n**2*log(F)**2) + 3*a*b**2*(c + d*x)**2*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F)) + 2*b**3*d**2*(F**(e*g + f*g*x))**(3*n)/(27*f**3*g**3*n**3*log(F)**3) - 2*b**3*d*(c + d*x)*(F**(e*g + f*g*x))**(3*n)/(9*f**2*g**2*n**2*log(F)**2) + b**3*(c + d*x)**2*(F**(e*g + f*g*x))**(3*n)/(3*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_41():
    f = (a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x)
    F = a**3*(c + d*x)**2/(2*d) - 3*a**2*b*d*(F**(e*g + f*g*x))**n/(f**2*g**2*n**2*log(F)**2) + 3*a**2*b*(c + d*x)*(F**(e*g + f*g*x))**n/(f*g*n*log(F)) - 3*a*b**2*d*(F**(e*g + f*g*x))**(2*n)/(4*f**2*g**2*n**2*log(F)**2) + 3*a*b**2*(c + d*x)*(F**(e*g + f*g*x))**(2*n)/(2*f*g*n*log(F)) - b**3*d*(F**(e*g + f*g*x))**(3*n)/(9*f**2*g**2*n**2*log(F)**2) + b**3*(c + d*x)*(F**(e*g + f*g*x))**(3*n)/(3*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_42():
    f = (a + b*(F**(g*(e + f*x)))**n)**3
    F = a**3*x + 3*a**2*b*(F**(g*(e + f*x)))**n/(f*g*n*log(F)) + 3*a*b**2*(F**(g*(e + f*x)))**(2*n)/(2*f*g*n*log(F)) + b**3*(F**(g*(e + f*x)))**(3*n)/(3*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_43():
    f = (a + b*(F**(g*(e + f*x)))**n)**3/(c + d*x)
    F = ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(((Integer(3) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(3) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n'))) * sympy.Function('ExpIntegralEi')(((Integer(3) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_44():
    f = (a + b*(F**(g*(e + f*x)))**n)**3/(c + d*x)**2
    F = (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('f') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('f') * (Symbol('F'))**(((Integer(3) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(3) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.Function('ExpIntegralEi')(((Integer(3) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_45():
    f = (a + b*(F**(g*(e + f*x)))**n)**3/(c + d*x)**3
    F = (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n'))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n')))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n')))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * (((Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('f') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n'))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(9) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**(((Integer(3) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(3) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n'))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(3) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_46():
    f = (c + d*x)**3/(a + b*(F**(g*(e + f*x)))**n)
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_47():
    f = (c + d*x)**2/(a + b*(F**(g*(e + f*x)))**n)
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_48():
    f = (c + d*x)/(a + b*(F**(g*(e + f*x)))**n)
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_49():
    f = 1/(a + b*(F**(g*(e + f*x)))**n)
    F = x/a - log(a + b*(F**(g*(e + f*x)))**n)/(a*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_50():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)*(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_51():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)*(c + d*x)**2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_52():
    f = (c + d*x)**3/(a + b*(F**(g*(e + f*x)))**n)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Symbol('a') * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_53():
    f = (c + d*x)**2/(a + b*(F**(g*(e + f*x)))**n)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_54():
    f = (c + d*x)/(a + b*(F**(g*(e + f*x)))**n)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Symbol('d') * sympy.log((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_55():
    f = (a + b*(F**(g*(e + f*x)))**n)**(-2)
    F = 1/(a*f*g*n*(a + b*(F**(g*(e + f*x)))**n)*log(F)) + x/a**2 - log(a + b*(F**(g*(e + f*x)))**n)/(a**2*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_56():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x))
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_57():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_58():
    f = (c + d*x)**3/(a + b*(F**(g*(e + f*x)))**n)**3
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Integer(2) * Symbol('a') * Symbol('f') * ((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Integer(2)) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(9) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(4)) * (Symbol('g'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_59():
    f = (c + d*x)**2/(a + b*(F**(g*(e + f*x)))**n)**3
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * x) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))) * (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('f') * ((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Integer(2)) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('g'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_60():
    f = (c + d*x)/(a + b*(F**(g*(e + f*x)))**n)**3
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('d') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * x) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)) * ((Integer(2) * Symbol('a') * Symbol('f') * ((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Integer(2)) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Symbol('c') + (Symbol('d') * x)) * (((Symbol('a'))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))) * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * sympy.log((Symbol('a') + (Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('g') * (Symbol('e') + (Symbol('f') * x)))))**(Symbol('n'))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('g'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_61():
    f = (a + b*(F**(g*(e + f*x)))**n)**(-3)
    F = 1/(2*a*f*g*n*(a + b*(F**(g*(e + f*x)))**n)**2*log(F)) + 1/(a**2*f*g*n*(a + b*(F**(g*(e + f*x)))**n)*log(F)) + x/a**3 - log(a + b*(F**(g*(e + f*x)))**n)/(a**3*f*g*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_62():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x))
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_63():
    f = 1/((a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_64():
    f = (a + b*exp(x))*sqrt(c + d*x)
    F = (Symbol('b') * (sympy.E)**(x) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(2) * Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**((Symbol('c') * (Symbol('d'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_65():
    f = (a + b*exp(x))**2*sqrt(c + d*x)
    F = (Integer(2) * Symbol('a') * Symbol('b') * (sympy.E)**(x) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * x)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**((Symbol('c') * (Symbol('d'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**(((Integer(2) * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_66():
    f = (a + b*exp(x))**3*sqrt(c + d*x)
    F = (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (sympy.E)**(x) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * x)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(3) * x)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**((Symbol('c') * (Symbol('d'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**(((Integer(2) * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('d')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((sympy.E)**(((Integer(3) * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_67():
    f = sqrt(c + d*x)/(a + b*exp(x))
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_68():
    f = sqrt(c + d*x)/(a + b*exp(x))**2
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_69():
    f = sqrt(c + d*x)/(a + b*exp(x))**3
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * (((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_70():
    f = (a + b*(F**(g*(e + f*x)))**n)**3*(c + d*x)**m
    F = (((Symbol('a'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(3)) * (Symbol('F'))**(((Integer(3) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(3) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(3) * Symbol('n'))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Integer(3) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1))) + ((Integer(3) * (Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_71():
    f = (a + b*(F**(g*(e + f*x)))**n)**2*(c + d*x)**m
    F = (((Symbol('a'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Integer(2))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**(((Integer(2) * (Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Integer(2) * Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**((Integer(2) * Symbol('n'))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Integer(2) * Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_72():
    f = (a + b*(F**(g*(e + f*x)))**n)*(c + d*x)**m
    F = ((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((((Symbol('e') + (Integer(-1) * ((Symbol('c') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * Symbol('g') * Symbol('n')) + (Integer(-1) * (Symbol('g') * Symbol('n') * (Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('f') * Symbol('g') * Symbol('n') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('f') * Symbol('g') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_73():
    f = (c + d*x)**m/(a + b*(F**(g*(e + f*x)))**n)
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_74():
    f = (c + d*x)**m/(a + b*(F**(g*(e + f*x)))**n)**2
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_75():
    f = (a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('F'))**(((Symbol('e') * Symbol('g')) + (Symbol('f') * Symbol('g') * x))))**(Symbol('n')))))**(Symbol('p')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_76():
    f = F**(c + d*x)*x**3/(F**(c + d*x)*b + a)
    F = (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_77():
    f = F**(c + d*x)*x**2/(F**(c + d*x)*b + a)
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_78():
    f = F**(c + d*x)*x/(F**(c + d*x)*b + a)
    F = ((x * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_79():
    f = F**(c + d*x)/(F**(c + d*x)*b + a)
    F = log(F**(c + d*x)*b + a)/(b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_80():
    f = F**(c + d*x)/(x*(F**(c + d*x)*b + a))
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('c') + (Symbol('d') * x))) * (((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_81():
    f = F**(c + d*x)/(x**2*(F**(c + d*x)*b + a))
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('c') + (Symbol('d') * x))) * (((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (x)**(Integer(2))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_82():
    f = F**(c + d*x)*x**3/(F**(c + d*x)*b + a)**2
    F = ((x)**(Integer(3)) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_83():
    f = F**(c + d*x)*x**2/(F**(c + d*x)*b + a)**2
    F = ((x)**(Integer(2)) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * sympy.log(Symbol('F'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_84():
    f = F**(c + d*x)*x/(F**(c + d*x)*b + a)**2
    F = -x/(b*d*(F**(c + d*x)*b + a)*log(F)) + x/(a*b*d*log(F)) - log(F**(c + d*x)*b + a)/(a*b*d**2*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_85():
    f = F**(c + d*x)/(F**(c + d*x)*b + a)**2
    F = -1/(b*d*(F**(c + d*x)*b + a)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_86():
    f = F**(c + d*x)/(x*(F**(c + d*x)*b + a)**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * x * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (x)**(Integer(2))))**(Integer(-1)), x) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_87():
    f = F**(c + d*x)/(x**2*(F**(c + d*x)*b + a)**2)
    F = (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (x)**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (x)**(Integer(3))))**(Integer(-1)), x)) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_88():
    f = F**(c + d*x)*x**3/(F**(c + d*x)*b + a)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_89():
    f = F**(c + d*x)*x**2/(F**(c + d*x)*b + a)**3
    F = (Integer(-1) * (x * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (x * ((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a'))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_90():
    f = F**(c + d*x)*x/(F**(c + d*x)*b + a)**3
    F = -x/(2*b*d*(F**(c + d*x)*b + a)**2*log(F)) + 1/(2*a*b*d**2*(F**(c + d*x)*b + a)*log(F)**2) + x/(2*a**2*b*d*log(F)) - log(F**(c + d*x)*b + a)/(2*a**2*b*d**2*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_91():
    f = F**(c + d*x)/(F**(c + d*x)*b + a)**3
    F = -1/(2*b*d*(F**(c + d*x)*b + a)**2*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_92():
    f = F**(c + d*x)/(x*(F**(c + d*x)*b + a)**3)
    F = (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * x * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)), x) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_2_c_plus_d_x_pow_m_F_pow_g_e_plus_f_x_pow_n_a_plus_b_F_pow_g_e_plus_f_x_pow_n_pow_p_93():
    f = F**(c + d*x)/(x**2*(F**(c + d*x)*b + a)**3)
    F = (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * (x)**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (Symbol('F'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1)), x) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F

