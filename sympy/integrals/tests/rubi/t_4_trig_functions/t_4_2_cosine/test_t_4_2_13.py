"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.13 (d+e x)^m cos(a+b x+c x^2)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e = symbols('a b c d e')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_1():
    f = x**2*cos(a + b*x + c*x**2)
    F = (((Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_2():
    f = x*cos(a + b*x + c*x**2)
    F = (Integer(-1) * ((Symbol('b') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_3():
    f = cos(a + b*x + c*x**2)
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_4():
    f = cos(a + b*x + c*x**2)/x
    F = sympy.Function('Unintegrable')((sympy.cos((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_5():
    f = b*sin(a + b*x + c*x**2)/x + cos(a + b*x + c*x**2)/x**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))))) + (Integer(-1) * (sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_6():
    f = x**2*cos(a + b*x - c*x**2)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_7():
    f = x*cos(a + b*x - c*x**2)
    F = (Integer(-1) * ((Symbol('b') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_8():
    f = cos(a + b*x - c*x**2)
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_9():
    f = cos(a + b*x - c*x**2)/x
    F = sympy.Function('Unintegrable')((sympy.cos((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_10():
    f = b*sin(a + b*x - c*x**2)/x + cos(a + b*x - c*x**2)/x**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) + (Integer(-1) * (sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_11():
    f = x**2*cos(x**2 + x + sympy.S(1)/4)
    F = ((Integer(4))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt((Integer(2) * sympy.pi)))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt((Integer(2) * sympy.pi)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.sin(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * x * sympy.sin(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_12():
    f = x*cos(x**2 + x + sympy.S(1)/4)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt((Integer(2) * sympy.pi)))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.sin(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_13():
    f = cos(x**2 + x + sympy.S(1)/4)
    F = sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt((Integer(2) * sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_14():
    f = cos(x**2 + x + sympy.S(1)/4)/x
    F = sympy.Function('Unintegrable')((sympy.cos(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_15():
    f = cos(x**2 + x + sympy.S(1)/4)/x**2
    F = (Integer(-1) * (sympy.cos(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt((Integer(2) * sympy.pi)))**(Integer(-1)))))) + (Integer(-1) * sympy.Function('Unintegrable')((sympy.sin(((Integer(4))**(Integer(-1)) + x + (x)**(Integer(2)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_16():
    f = x**2*cos(a + b*x + c*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_17():
    f = x*cos(a + b*x + c*x**2)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_18():
    f = cos(a + b*x + c*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_19():
    f = cos(a + b*x + c*x**2)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_20():
    f = x**2*cos(a + b*x - c*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_21():
    f = x*cos(a + b*x - c*x**2)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_22():
    f = cos(a + b*x - c*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_23():
    f = cos(a + b*x - c*x**2)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(-1) * (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_24():
    f = x**2*cos(x**2 + x + sympy.S(1)/4)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + ((Integer(16))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.sin(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))) + ((Integer(8))**(Integer(-1)) * x * sympy.sin(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_25():
    f = x*cos(x**2 + x + sympy.S(1)/4)**2
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + ((Integer(8))**(Integer(-1)) * sympy.sin(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_26():
    f = cos(x**2 + x + sympy.S(1)/4)**2
    F = (x * (Integer(2))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_27():
    f = cos(x**2 + x + sympy.S(1)/4)**2/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cos(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_28():
    f = cos(x**2 + x + sympy.S(1)/4)**2/x**2
    F = (Integer(-1) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.cos(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + (Integer(-1) * sympy.Function('Unintegrable')((sympy.sin(((Integer(2))**(Integer(-1)) + (Integer(2) * x) + (Integer(2) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_29():
    f = (d + e*x)**2*cos(a + b*x + c*x**2)
    F = (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * (Symbol('d') + (Symbol('e') * x)) * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_30():
    f = (d + e*x)*cos(a + b*x + c*x**2)
    F = ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * sympy.sin((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_31():
    f = cos(a + b*x + c*x**2)/(d + e*x)
    F = sympy.Function('Unintegrable')((sympy.cos((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_32():
    f = (d + e*x)**2*cos(a + b*x + c*x**2)**2
    F = (((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * ((Integer(6) * Symbol('e')))**(Integer(-1))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(16) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(16) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Symbol('e') * (Symbol('d') + (Symbol('e') * x)) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_33():
    f = (d + e*x)*cos(a + b*x + c*x**2)**2
    F = (((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * ((Integer(4) * Symbol('e')))**(Integer(-1))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('FresnelC')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')))**(Integer(-1))))))) * ((Integer(8) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_13_d_plus_e_x_pow_m_cos_a_plus_b_x_plus_c_x_pow_2_pow_n_34():
    f = cos(a + b*x + c*x**2)**2/(d + e*x)
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x) + (Integer(2) * Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)) + (sympy.log((Symbol('d') + (Symbol('e') * x))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F

