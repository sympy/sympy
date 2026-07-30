"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.12 (e x)^m (a+b cos(c+d x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, m, n, p = symbols('a b c d e m n p')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_1():
    f = x**3*cos(a + b*x**2)
    F = x**2*sin(a + b*x**2)/(2*b) + cos(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_2():
    f = x**2*cos(a + b*x**2)
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a'))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((x * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_3():
    f = x*cos(a + b*x**2)
    F = sin(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_4():
    f = cos(a + b*x**2)
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_5():
    f = cos(a + b*x**2)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(2))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_6():
    f = cos(a + b*x**2)/x**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a'))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_7():
    f = cos(a + b*x**2)/x**3
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(2)))) * sympy.sin(Symbol('a')))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_8():
    f = x**3*cos(a + b*x**2)**2
    F = x**4/8 + x**2*sin(a + b*x**2)*cos(a + b*x**2)/(4*b) + cos(a + b*x**2)**2/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_9():
    f = x**2*cos(a + b*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((x * sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_10():
    f = x*cos(a + b*x**2)**2
    F = x**2/4 + sin(a + b*x**2)*cos(a + b*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_11():
    f = cos(a + b*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_12():
    f = cos(a + b*x**2)**2/x
    F = ((Integer(4))**(Integer(-1)) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2))))) + (sympy.log(x) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_13():
    f = cos(a + b*x**2)**2/x**2
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_14():
    f = cos(a + b*x**2)**2/x**3
    F = (Integer(-1) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.cos((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2)))) * sympy.sin((Integer(2) * Symbol('a'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_15():
    f = x**3*cos(a + b*x**2)**3
    F = x**2*sin(a + b*x**2)*cos(a + b*x**2)**2/(6*b) + x**2*sin(a + b*x**2)/(3*b) + cos(a + b*x**2)**3/(18*b**2) + cos(a + b*x**2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_16():
    f = x**2*cos(a + b*x**2)**3
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x))) * ((Integer(24) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a'))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin((Integer(3) * Symbol('a')))) * ((Integer(24) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * x * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((x * sympy.sin(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * (x)**(Integer(2)))))) * ((Integer(24) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_17():
    f = x*cos(a + b*x**2)**3
    F = -sin(a + b*x**2)**3/(6*b) + sin(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_18():
    f = cos(a + b*x**2)**3
    F = ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a'))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin((Integer(3) * Symbol('a')))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_19():
    f = cos(a + b*x**2)**3/x
    F = ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(2))))) + ((Integer(8))**(Integer(-1)) * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2))))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sin((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_20():
    f = cos(a + b*x**2)**3/x**2
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin(Symbol('a')))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * x)) * sympy.sin((Integer(3) * Symbol('a')))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_21():
    f = cos(a + b*x**2)**3/x**3
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.cos((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(2)))) * sympy.sin(Symbol('a')))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2)))) * sympy.sin((Integer(3) * Symbol('a'))))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_22():
    f = x*cos(a + b*x**2)**7
    F = -sin(a + b*x**2)**7/(14*b) + 3*sin(a + b*x**2)**5/(10*b) - sin(a + b*x**2)**3/(2*b) + sin(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_23():
    f = x**(sympy.S(5)/2)*cos(a + b*x**2)
    F = (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((Integer(16) * Symbol('b') * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * (Integer(16) * Symbol('b') * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))))**(Integer(-1))) + (((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_24():
    f = x**(sympy.S(3)/2)*cos(a + b*x**2)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((Integer(16) * Symbol('b') * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * (Integer(16) * Symbol('b') * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt(x) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_25():
    f = sqrt(x)*cos(a + b*x**2)
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((Integer(4) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * (Integer(4) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_26():
    f = cos(a + b*x**2)/sqrt(x)
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((Integer(4) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * (Integer(4) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_27():
    f = cos(a + b*x**2)/x**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_28():
    f = cos(a + b*x**2)/x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))) * ((Integer(3) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * (Integer(3) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_29():
    f = x**(sympy.S(5)/2)*cos(a + b*x**2)**2
    F = ((x)**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Integer(7))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * ((Integer(64) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('b') * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (Integer(64) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('b') * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))))**(Integer(-1))) + (((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_30():
    f = x**(sympy.S(3)/2)*cos(a + b*x**2)**2
    F = ((x)**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Integer(5))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * ((Integer(64) * (Integer(2))**((Integer(4))**(Integer(-1))) * Symbol('b') * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (Integer(64) * (Integer(2))**((Integer(4))**(Integer(-1))) * Symbol('b') * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt(x) * sympy.sin((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_31():
    f = sqrt(x)*cos(a + b*x**2)**2
    F = ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * ((Integer(8) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (Integer(8) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_32():
    f = cos(a + b*x**2)**2/sqrt(x)
    F = sympy.sqrt(x) + (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * ((Integer(8) * (Integer(2))**((Integer(4))**(Integer(-1))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (Integer(8) * (Integer(2))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_33():
    f = cos(a + b*x**2)**2/x**(sympy.S(3)/2)
    F = (Integer(-1) * (sympy.sqrt(x))**(Integer(-1))) + (Integer(-1) * (sympy.cos((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(3) * (Integer(4))**(Integer(-1))), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_34():
    f = cos(a + b*x**2)**2/x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(2))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * ((Integer(3) * (Integer(2))**((Integer(4))**(Integer(-1))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(4))**(Integer(-1)), (Integer(2) * sympy.I * Symbol('b') * (x)**(Integer(2))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (Integer(3) * (Integer(2))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('b') * (x)**(Integer(2))))**((Integer(4))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_35():
    f = cos(a + b/x)
    F = (x * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Integer(-1)))))) + (Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(-1)))) * sympy.sin(Symbol('a'))) + (Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_36():
    f = cos(a + b/x)/x
    F = ((Integer(-1) * sympy.cos(Symbol('a'))) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Integer(-1))))) + (sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_37():
    f = cos(a + b/x)/x**2
    F = -sin(a + b/x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_38():
    f = cos(a + b/x)/x**3
    F = -sin(a + b/x)/(b*x) - cos(a + b/x)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_39():
    f = cos(a + b/x)/x**4
    F = -sin(a + b/x)/(b*x**2) - 2*cos(a + b/x)/(b**2*x) + 2*sin(a + b/x)/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_40():
    f = cos(a + b/x**2)
    F = (x * sympy.cos((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1))))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1)))) * sympy.sin(Symbol('a')))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_41():
    f = cos(a + b/x**2)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_42():
    f = cos(a + b/x**2)/x**2
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1)))) * sympy.sin(Symbol('a'))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_43():
    f = cos(a + b/x**2)/x**3
    F = -sin(a + b/x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_44():
    f = cos(a + b/x**2)/x**4
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1))))) * (x)**(Integer(-1)))) * sympy.sin(Symbol('a'))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_45():
    f = cos(sqrt(x))**2/sqrt(x)
    F = sqrt(x) + sin(sqrt(x))*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_46():
    f = cos(sqrt(x))/sqrt(x)
    F = 2*sin(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_47():
    f = cos(sqrt(x))
    F = 2*sqrt(x)*sin(sqrt(x)) + 2*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_48():
    f = cos(sqrt(x))**2
    F = sqrt(x)*sin(sqrt(x))*cos(sqrt(x)) + x/2 + cos(sqrt(x))**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_49():
    f = x**(sympy.S(3)/2)*cos(a + b*x**(sympy.S(1)/3))
    F = ((Integer(135135) * sympy.sqrt(x) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(32) * (Symbol('b'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(3861) * (x)**((Integer(7) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(39) * (x)**((Integer(11) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(405405) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1)))))) * ((Integer(64) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(405405) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a'))) * ((Integer(64) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(405405) * (x)**((Integer(6))**(Integer(-1))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(64) * (Symbol('b'))**(Integer(7))))**(Integer(-1)))) + ((Integer(27027) * (x)**((Integer(5) * (Integer(6))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(16) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(429) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(13) * (Integer(6))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_50():
    f = sqrt(x)*cos(a + b*x**(sympy.S(1)/3))
    F = (Integer(-1) * ((Integer(315) * (x)**((Integer(6))**(Integer(-1))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(21) * (x)**((Integer(5) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(315) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(315) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a'))) * ((Integer(8) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(105) * sympy.sqrt(x) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(7) * (Integer(6))**(Integer(-1)))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_51():
    f = cos(a + b*x**(sympy.S(1)/3))/sqrt(x)
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1)))))) * ((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a'))) * ((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(6))**(Integer(-1))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_52():
    f = cos(a + b*x**(sympy.S(1)/3))/x**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))))) + (Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a'))) + ((Integer(4) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((x)**((Integer(6))**(Integer(-1))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_53():
    f = cos(a + b*x**(sympy.S(1)/3))/x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(105) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(315) * (x)**((Integer(6))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * (Integer(315))**(Integer(-1))) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Symbol('a')) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))))) + (Integer(-1) * ((Integer(32) * (Integer(315))**(Integer(-1))) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a')))) + ((Integer(4) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(21) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(315) * sympy.sqrt(x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_54():
    f = cos(a + b*x**(sympy.S(1)/3))/x**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(5) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(715) * (x)**((Integer(11) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(45045) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(128) * (Symbol('b'))**(Integer(6)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(675675) * sympy.sqrt(x)))**(Integer(-1))) + ((Integer(256) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Symbol('a')) * sympy.Function('FresnelC')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1)))))) * (Integer(675675))**(Integer(-1))) + (Integer(-1) * ((Integer(256) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (x)**((Integer(6))**(Integer(-1))))) * sympy.sin(Symbol('a'))) * (Integer(675675))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(65) * (x)**((Integer(13) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(6435) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(64) * (Symbol('b'))**(Integer(5)) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(225225) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(256) * (Symbol('b'))**(Integer(7)) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(675675) * (x)**((Integer(6))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_55():
    f = x**(sympy.S(3)/2)*cos(a + b*x**(sympy.S(1)/3))**2
    F = (Integer(-1) * ((Integer(135135) * sympy.sqrt(x)) * ((Integer(4096) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + ((Integer(3861) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))) * ((Integer(256) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(39) * (x)**((Integer(11) * (Integer(6))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Integer(5))**(Integer(-1))) + ((Integer(135135) * sympy.sqrt(x) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(2048) * (Symbol('b'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(3861) * (x)**((Integer(7) * (Integer(6))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(128) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(39) * (x)**((Integer(11) * (Integer(6))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(405405) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(32768) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(405405) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(32768) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(27027) * (x)**((Integer(5) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(512) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(429) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(32) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(13) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(405405) * (x)**((Integer(6))**(Integer(-1))) * sympy.sin((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))) * ((Integer(16384) * (Symbol('b'))**(Integer(7))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_56():
    f = sqrt(x)*cos(a + b*x**(sympy.S(1)/3))**2
    F = ((Integer(315) * (x)**((Integer(6))**(Integer(-1)))) * ((Integer(256) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(21) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(315) * (x)**((Integer(6))**(Integer(-1))) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(128) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(21) * (x)**((Integer(5) * (Integer(6))**(Integer(-1)))) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(315) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(512) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(315) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(512) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(105) * sympy.sqrt(x) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(32) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(7) * (Integer(6))**(Integer(-1)))) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_57():
    f = cos(a + b*x**(sympy.S(1)/3))**2/sqrt(x)
    F = sympy.sqrt(x) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (x)**((Integer(6))**(Integer(-1))) * sympy.sin((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_58():
    f = cos(a + b*x**(sympy.S(1)/3))**2/x**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + (Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((x)**((Integer(6))**(Integer(-1))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_59():
    f = cos(a + b*x**(sympy.S(1)/3))**2/x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(105) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(256) * (Symbol('b'))**(Integer(4))) * ((Integer(315) * (x)**((Integer(6))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(105) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(512) * (Symbol('b'))**(Integer(4)) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(315) * (x)**((Integer(6))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(512) * (Integer(315))**(Integer(-1))) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))))) + (Integer(-1) * ((Integer(512) * (Integer(315))**(Integer(-1))) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a'))))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(21) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(315) * sympy.sqrt(x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_60():
    f = cos(a + b*x**(sympy.S(1)/3))**2/x**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2))) * ((Integer(715) * (x)**((Integer(11) * (Integer(6))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(256) * (Symbol('b'))**(Integer(4))) * ((Integer(45045) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4096) * (Symbol('b'))**(Integer(6))) * ((Integer(675675) * sympy.sqrt(x)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(5) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(715) * (x)**((Integer(11) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(512) * (Symbol('b'))**(Integer(4)) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(45045) * (x)**((Integer(7) * (Integer(6))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8192) * (Symbol('b'))**(Integer(6)) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))))**(Integer(2))) * ((Integer(675675) * sympy.sqrt(x)))**(Integer(-1))) + ((Integer(32768) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * (Integer(675675))**(Integer(-1))) + (Integer(-1) * ((Integer(32768) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(Symbol('b')) * (x)**((Integer(6))**(Integer(-1)))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Integer(2) * Symbol('a')))) * (Integer(675675))**(Integer(-1)))) + ((Integer(8) * Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(65) * (x)**((Integer(13) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(6435) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2048) * (Symbol('b'))**(Integer(5)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(225225) * (x)**((Integer(5) * (Integer(6))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(32768) * (Symbol('b'))**(Integer(7)) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Symbol('b') * (x)**((Integer(3))**(Integer(-1))))))) * ((Integer(675675) * (x)**((Integer(6))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_61():
    f = cos(x**(sympy.S(1)/3))**3
    F = x**(sympy.S(2)/3)*sin(x**(sympy.S(1)/3))*cos(x**(sympy.S(1)/3))**2 + 2*x**(sympy.S(2)/3)*sin(x**(sympy.S(1)/3)) + 2*x**(sympy.S(1)/3)*cos(x**(sympy.S(1)/3))**3/3 + 4*x**(sympy.S(1)/3)*cos(x**(sympy.S(1)/3)) + 2*sin(x**(sympy.S(1)/3))**3/9 - 14*sin(x**(sympy.S(1)/3))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_62():
    f = cos(x**(sympy.S(1)/6))/x**(sympy.S(5)/6)
    F = 6*sin(x**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_63():
    f = (b*cos(c + d*x**n))**p*(e*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') * x))**(Symbol('m')) * ((Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * (x)**(Symbol('n')))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_64():
    f = (e*x)**m*(a + b*cos(c + d*x**n))**p
    F = sympy.Function('Unintegrable')((((Symbol('e') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * (x)**(Symbol('n'))))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_65():
    f = (b*cos(c + d*x**n))**p*(e*x)**(n - 1)
    F = -(b*cos(c + d*x**n))**(p + 1)*(e*x)**n*sin(c + d*x**n)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), cos(c + d*x**n)**2)/(b*d*e*n*x**n*(p + 1)*sqrt(sin(c + d*x**n)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_66():
    f = (b*cos(c + d*x**n))**p*(e*x)**(2*n - 1)
    F = (((Symbol('e') * x))**((Integer(2) * Symbol('n'))) * sympy.Function('Unintegrable')(((x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * (x)**(Symbol('n')))))))**(Symbol('p'))), x)) * (((x)**((Integer(2) * Symbol('n'))) * Symbol('e')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_67():
    f = (e*x)**(n - 1)*(a + b*cos(c + d*x**n))**p
    F = sqrt(2)*(e*x)**n*(a + b*cos(c + d*x**n))**p*sin(c + d*x**n)*appellf1(sympy.S.Half, sympy.S.Half, -p, sympy.S(3)/2, sympy.S.Half - cos(c + d*x**n)/2, b*(1 - cos(c + d*x**n))/(a + b))/(d*e*n*x**n*((a + b*cos(c + d*x**n))/(a + b))**p*sqrt(cos(c + d*x**n) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_68():
    f = (e*x)**(2*n - 1)*(a + b*cos(c + d*x**n))**p
    F = (((Symbol('e') * x))**((Integer(2) * Symbol('n'))) * sympy.Function('Unintegrable')(((x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * (x)**(Symbol('n'))))))))**(Symbol('p'))), x)) * (((x)**((Integer(2) * Symbol('n'))) * Symbol('e')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_69():
    f = cos(a + b*x**n)/x
    F = ((sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_70():
    f = cos(a + b*x**n)**2/x
    F = ((sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (sympy.log(x) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_71():
    f = cos(a + b*x**n)**3/x
    F = ((Integer(3) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_72():
    f = cos(a + b*x**n)**4/x
    F = ((sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((sympy.cos((Integer(4) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(4) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1))) + ((Integer(3) * sympy.log(x)) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Integer(4) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(4) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_73():
    f = cos(a + b*x**n)
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_74():
    f = cos(a + b*x**n)**2
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * (Symbol('n'))**(Integer(-1))))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * (Symbol('n'))**(Integer(-1))))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(2) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_75():
    f = cos(a + b*x**n)**3
    F = (Integer(-1) * ((Integer(3) * (sympy.E)**((sympy.I * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Integer(3) * sympy.I * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(-3) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**((Symbol('n'))**(Integer(-1))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(3) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**((Symbol('n'))**(Integer(-1))) * (sympy.E)**((Integer(3) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_76():
    f = x**m*cos(a + b*x**n)
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_77():
    f = x**m*cos(a + b*x**n)**2
    F = ((x)**((Integer(1) + Symbol('m'))) * ((Integer(2) * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-2) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(2))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('n'))) * (Symbol('n'))**(Integer(-1)))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(2) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(2))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('n'))) * (Symbol('n'))**(Integer(-1)))) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_78():
    f = x**m*cos(a + b*x**n)**3
    F = (Integer(-1) * ((Integer(3) * (sympy.E)**((sympy.I * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Integer(3) * sympy.I * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-3) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (((Integer(-1) * sympy.I) * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(3) * sympy.I * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (sympy.E)**((Integer(3) * sympy.I * Symbol('a'))) * ((sympy.I * Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_79():
    f = x**(-n - 1)*cos(a + b*x**n)
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n')))) * sympy.sin(Symbol('a'))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_80():
    f = x**(-n - 1)*cos(a + b*x**n)**2
    F = (Integer(-1) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * (sympy.cos((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n')))) * sympy.sin((Integer(2) * Symbol('a')))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_81():
    f = x**(-n - 1)*cos(a + b*x**n)**3
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (sympy.cos((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n')))) * sympy.sin(Symbol('a'))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n')))) * sympy.sin((Integer(3) * Symbol('a')))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_82():
    f = x**(-2*n - 1)*cos(a + b*x**n)
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_83():
    f = x**(-2*n - 1)*cos(a + b*x**n)**2
    F = (Integer(-1) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(4) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * (sympy.cos((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))) + ((Symbol('b') * sympy.sin((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_84():
    f = x**(-2*n - 1)*cos(a + b*x**n)**3
    F = (Integer(-1) * ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (sympy.cos((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(8) * Symbol('n'))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.sin((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * (Integer(8) * Symbol('n'))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1))) + ((Integer(9) * (Symbol('b'))**(Integer(2)) * sympy.sin((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(8) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_85():
    f = x**2*cos((a + b*x)**2)
    F = (((Symbol('a'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sin(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sin(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_86():
    f = x*cos((a + b*x)**2)
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.sin(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_87():
    f = cos((a + b*x)**2)
    F = (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_88():
    f = cos((a + b*x)**2)/x
    F = sympy.Function('Unintegrable')((sympy.cos(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_89():
    f = cos((a + b*x)**2)/x**2
    F = sympy.Function('Unintegrable')((sympy.cos(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_90():
    f = x**2*cos(a + b*sqrt(c + d*x))
    F = 2*c**2*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b*d**3) - 4*c*(c + d*x)**(sympy.S(3)/2)*sin(a + b*sqrt(c + d*x))/(b*d**3) + 2*(c + d*x)**(sympy.S(5)/2)*sin(a + b*sqrt(c + d*x))/(b*d**3) + 2*c**2*cos(a + b*sqrt(c + d*x))/(b**2*d**3) - 12*c*(c + d*x)*cos(a + b*sqrt(c + d*x))/(b**2*d**3) + 10*(c + d*x)**2*cos(a + b*sqrt(c + d*x))/(b**2*d**3) + 24*c*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b**3*d**3) - 40*(c + d*x)**(sympy.S(3)/2)*sin(a + b*sqrt(c + d*x))/(b**3*d**3) + 24*c*cos(a + b*sqrt(c + d*x))/(b**4*d**3) - (120*c + 120*d*x)*cos(a + b*sqrt(c + d*x))/(b**4*d**3) + 240*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b**5*d**3) + 240*cos(a + b*sqrt(c + d*x))/(b**6*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_91():
    f = x*cos(a + b*sqrt(c + d*x))
    F = -2*c*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b*d**2) + 2*(c + d*x)**(sympy.S(3)/2)*sin(a + b*sqrt(c + d*x))/(b*d**2) - 2*c*cos(a + b*sqrt(c + d*x))/(b**2*d**2) + (6*c + 6*d*x)*cos(a + b*sqrt(c + d*x))/(b**2*d**2) - 12*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b**3*d**2) - 12*cos(a + b*sqrt(c + d*x))/(b**4*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_92():
    f = cos(a + b*sqrt(c + d*x))
    F = 2*sqrt(c + d*x)*sin(a + b*sqrt(c + d*x))/(b*d) + 2*cos(a + b*sqrt(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_93():
    f = cos(a + b*sqrt(c + d*x))/x
    F = (sympy.cos((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('CosIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) + (sympy.cos((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('CosIntegral')(((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('SinIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))) + (sympy.sin((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('SinIntegral')(((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_94():
    f = cos(a + b*sqrt(c + d*x))/x**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))) * (x)**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('CosIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c'))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('CosIntegral')(((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) * sympy.sin((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c')))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.cos((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('SinIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * sympy.cos((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('SinIntegral')(((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_95():
    f = x**2*cos(a + b*(c + d*x)**(sympy.S(1)/3))
    F = 3*c**2*(c + d*x)**(sympy.S(2)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) - 6*c*(c + d*x)**(sympy.S(5)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) + 3*(c + d*x)**(sympy.S(8)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) + 6*c**2*(c + d*x)**(sympy.S(1)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) - 30*c*(c + d*x)**(sympy.S(4)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) + 24*(c + d*x)**(sympy.S(7)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) - 6*c**2*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) + 120*c*(c + d*x)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) - 168*(c + d*x)**2*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) + 360*c*(c + d*x)**(sympy.S(2)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**3) - 1008*(c + d*x)**(sympy.S(5)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**3) - 720*c*(c + d*x)**(sympy.S(1)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**3) + 5040*(c + d*x)**(sympy.S(4)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**3) - 720*c*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**3) + (20160*c + 20160*d*x)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**3) - 60480*(c + d*x)**(sympy.S(2)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**7*d**3) - 120960*(c + d*x)**(sympy.S(1)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**8*d**3) + 120960*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**9*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_96():
    f = x*cos(a + b*(c + d*x)**(sympy.S(1)/3))
    F = -3*c*(c + d*x)**(sympy.S(2)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**2) + 3*(c + d*x)**(sympy.S(5)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**2) - 6*c*(c + d*x)**(sympy.S(1)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**2) + 15*(c + d*x)**(sympy.S(4)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**2) + 6*c*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**2) - (60*c + 60*d*x)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**2) - 180*(c + d*x)**(sympy.S(2)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**2) + 360*(c + d*x)**(sympy.S(1)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**2) + 360*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_97():
    f = cos(a + b*(c + d*x)**(sympy.S(1)/3))
    F = 3*(c + d*x)**(sympy.S(2)/3)*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d) + 6*(c + d*x)**(sympy.S(1)/3)*cos(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d) - 6*sin(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_98():
    f = cos(a + b*(c + d*x)**(sympy.S(1)/3))/x
    F = (sympy.cos((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('CosIntegral')(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (sympy.cos((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (sympy.cos((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('CosIntegral')((((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) + (sympy.sin((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinIntegral')(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (sympy.sin((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('SinIntegral')((((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_12_e_x_pow_m_a_plus_b_cos_c_plus_d_x_pow_n_pow_p_99():
    f = cos(a + b*(c + d*x)**(sympy.S(1)/3))/x**2
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('CosIntegral')(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * sympy.sin((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('d') * sympy.Function('CosIntegral')((((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * sympy.Function('CosIntegral')((((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * sympy.sin((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.cos((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinIntegral')(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * sympy.cos((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('d') * sympy.cos((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('SinIntegral')((((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F

