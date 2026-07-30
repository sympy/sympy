"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.2 Hyperbolic cosine/6.2.4 (d+e x)^m cosh(a+b x+c x^2)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e = symbols('a b c d e')

def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_1():
    f = x**2*cosh(a + b*x + c*x**2)
    F = (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_2():
    f = x*cosh(a + b*x + c*x**2)
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_3():
    f = cosh(a + b*x + c*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_4():
    f = cosh(a + b*x + c*x**2)/x
    F = sympy.Function('Unintegrable')((sympy.cosh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_5():
    f = -b*sinh(a + b*x + c*x**2)/x + cosh(a + b*x + c*x**2)/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('c')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('c')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_6():
    f = x**2*cosh(a + b*x - c*x**2)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_7():
    f = x*cosh(a + b*x - c*x**2)
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_8():
    f = cosh(a + b*x - c*x**2)
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_9():
    f = cosh(a + b*x - c*x**2)/x
    F = sympy.Function('Unintegrable')((sympy.cosh((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_10():
    f = -b*sinh(a + b*x - c*x**2)/x + cosh(a + b*x - c*x**2)/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('c')) * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('c')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_11():
    f = x**2*cosh(x**2 + x + sympy.S(1)/4)
    F = ((Integer(-1) * (Integer(3) * (Integer(16))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * (Integer(2) * x)))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(2) * x)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.sinh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * x * sympy.sinh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_12():
    f = x*cosh(x**2 + x + sympy.S(1)/4)
    F = ((Integer(8))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * (Integer(2) * x)))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(2) * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.sinh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_13():
    f = cosh(x**2 + x + sympy.S(1)/4)
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * (Integer(2) * x)))))) + ((Integer(4))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_14():
    f = cosh(x**2 + x + sympy.S(1)/4)/x
    F = sympy.Function('Unintegrable')((sympy.cosh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_15():
    f = cosh(x**2 + x + sympy.S(1)/4)/x**2
    F = (Integer(-1) * (sympy.cosh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * (Integer(2) * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(2) * x))))) + sympy.Function('Unintegrable')((sympy.sinh(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_16():
    f = x**2*cosh(a + b*x + c*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_17():
    f = x*cosh(a + b*x + c*x**2)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_18():
    f = cosh(a + b*x + c*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_19():
    f = cosh(a + b*x + c*x**2)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_20():
    f = x**2*cosh(a + b*x - c*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(-2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_21():
    f = x*cosh(a + b*x - c*x**2)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_22():
    f = cosh(a + b*x - c*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_23():
    f = cosh(a + b*x - c*x**2)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_24():
    f = x**2*cosh(x**2 + x + sympy.S(1)/4)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + ((Integer(16))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(2)))**(Integer(-1))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sinh(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))) + ((Integer(8))**(Integer(-1)) * x * sympy.sinh(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_25():
    f = x*cosh(x**2 + x + sympy.S(1)/4)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(2)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(2)))**(Integer(-1)))))) + ((Integer(8))**(Integer(-1)) * sympy.sinh(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_26():
    f = cosh(x**2 + x + sympy.S(1)/4)**2
    F = (x * (Integer(2))**(Integer(-1))) + ((Integer(8))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(2)))**(Integer(-1))))) + ((Integer(8))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(2)))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_27():
    f = cosh(x**2 + x + sympy.S(1)/4)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cosh(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_28():
    f = (d + e*x)**2*cosh(a + b*x + c*x**2)
    F = (((Symbol('e'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * (Symbol('d') + (Symbol('e') * x)) * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_29():
    f = (d + e*x)*cosh(a + b*x + c*x**2)
    F = ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * sympy.sinh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_30():
    f = cosh(a + b*x + c*x**2)/(d + e*x)
    F = sympy.Function('Unintegrable')((sympy.cosh((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_31():
    f = (d + e*x)**2*cosh(a + b*x + c*x**2)**2
    F = (((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * ((Integer(6) * Symbol('e')))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * (Symbol('d') + (Symbol('e') * x)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_32():
    f = (d + e*x)*cosh(a + b*x + c*x**2)**2
    F = (((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * ((Integer(4) * Symbol('e')))**(Integer(-1))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.E)**(((Integer(-2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('e') * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_4_d_plus_e_x_pow_m_cosh_a_plus_b_x_plus_c_x_pow_2_pow_n_33():
    f = cosh(a + b*x + c*x**2)**2/(d + e*x)
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)) + (sympy.log((Symbol('d') + (Symbol('e') * x))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F

