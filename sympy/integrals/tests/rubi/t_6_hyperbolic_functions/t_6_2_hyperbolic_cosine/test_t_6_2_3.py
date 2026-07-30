"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.2 Hyperbolic cosine/6.2.3 (e x)^m (a+b cosh(c+d x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, m, n, p = symbols('a b c d e m n p')

def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_1():
    f = x**3*cosh(a + b*x**2)
    F = x**2*sinh(a + b*x**2)/(2*b) - cosh(a + b*x**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_2():
    f = x**2*cosh(a + b*x**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_3():
    f = x*cosh(a + b*x**2)
    F = sinh(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_4():
    f = cosh(a + b*x**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * sympy.sqrt(Symbol('b')))))**(Integer(-1))) + (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_5():
    f = cosh(a + b*x**2)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(2))))) + ((Integer(2))**(Integer(-1)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_6():
    f = cosh(a + b*x**2)/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_7():
    f = cosh(a + b*x**2)/x**3
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(2)))) * sympy.sinh(Symbol('a'))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_8():
    f = x**3*cosh(a + b*x**2)**2
    F = x**4/8 + x**2*sinh(a + b*x**2)*cosh(a + b*x**2)/(4*b) - cosh(a + b*x**2)**2/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_9():
    f = x**2*cosh(a + b*x**2)**2
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * (Integer(32) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x))) * ((Integer(32) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((x * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_10():
    f = x*cosh(a + b*x**2)**2
    F = x**2/4 + sinh(a + b*x**2)*cosh(a + b*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_11():
    f = cosh(a + b*x**2)**2
    F = (x * (Integer(2))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * (Integer(8) * sympy.sqrt(Symbol('b')))))**(Integer(-1))) + (((sympy.E)**((Integer(2) * Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x))) * ((Integer(8) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_12():
    f = cosh(a + b*x**2)**2/x
    F = ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2))))) + (sympy.log(x) * (Integer(2))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_13():
    f = cosh(a + b*x**2)**2/x**2
    F = (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x))) * ((sympy.E)**((Integer(2) * Symbol('a'))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (sympy.E)**((Integer(2) * Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('b')) * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_14():
    f = cosh(a + b*x**2)**2/x**3
    F = (Integer(-1) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2)))) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_15():
    f = x**3*cosh(a + b*x**2)**3
    F = x**2*sinh(a + b*x**2)*cosh(a + b*x**2)**2/(6*b) + x**2*sinh(a + b*x**2)/(3*b) - cosh(a + b*x**2)**3/(18*b**2) - cosh(a + b*x**2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_16():
    f = x**2*cosh(a + b*x**2)**3
    F = ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(32) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * (Integer(96) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x))) * ((Integer(32) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Integer(3) * Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x))) * ((Integer(96) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * x * sympy.sinh((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((x * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * (x)**(Integer(2)))))) * ((Integer(24) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_17():
    f = x*cosh(a + b*x**2)**3
    F = sinh(a + b*x**2)**3/(6*b) + sinh(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_18():
    f = cosh(a + b*x**2)**3
    F = ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(16) * sympy.sqrt(Symbol('b')))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * (Integer(16) * sympy.sqrt(Symbol('b')))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x))) * ((Integer(16) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + (((sympy.E)**((Integer(3) * Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x))) * ((Integer(16) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_19():
    f = cosh(a + b*x**2)**3/x
    F = ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(2))))) + ((Integer(8))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2))))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(2))))) + ((Integer(8))**(Integer(-1)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_20():
    f = cosh(a + b*x**2)**3/x**2
    F = (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * (((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (((Integer(8))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x))) * ((sympy.E)**((Integer(3) * Symbol('a'))))**(Integer(-1)))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * (sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x))) + ((Integer(8))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (sympy.E)**((Integer(3) * Symbol('a'))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('b')) * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_21():
    f = cosh(a + b*x**2)**3/x**3
    F = (Integer(-1) * ((Integer(3) * sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(2)))) * sympy.sinh(Symbol('a'))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2)))) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(2))))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_22():
    f = x*cosh(a + b*x**2)**7
    F = sinh(a + b*x**2)**7/(14*b) + 3*sinh(a + b*x**2)**5/(10*b) + sinh(a + b*x**2)**3/(2*b) + sinh(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_23():
    f = x**2*cosh(x**3)
    F = sinh(x**3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_24():
    f = cosh(x**(-5))/x**6
    F = -sinh(x**(-5))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_25():
    f = cosh(a + b/x)
    F = (x * sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Integer(-1)))))) + (Integer(-1) * (Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(-1)))) * sympy.sinh(Symbol('a')))) + (Integer(-1) * (Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_26():
    f = cosh(a + b/x)/x
    F = ((Integer(-1) * sympy.cosh(Symbol('a'))) * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Integer(-1))))) + (Integer(-1) * (sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_27():
    f = cosh(a + b/x)/x**2
    F = -sinh(a + b/x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_28():
    f = cosh(a + b/x)/x**3
    F = -sinh(a + b/x)/(b*x) + cosh(a + b/x)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_29():
    f = cosh(a + b/x)/x**4
    F = -sinh(a + b/x)/(b*x**2) + 2*cosh(a + b/x)/(b**2*x) - 2*sinh(a + b/x)/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_30():
    f = cosh(a + b/x**2)
    F = (x * sympy.cosh((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) + (((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(Symbol('b')) * (sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_31():
    f = cosh(a + b/x**2)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_32():
    f = cosh(a + b/x**2)/x**2
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * sympy.sqrt(Symbol('b')))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('b'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_33():
    f = cosh(a + b/x**2)/x**3
    F = -sinh(a + b/x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_34():
    f = cosh(a + b/x**2)/x**4
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))) * (((sympy.E)**(Symbol('a')) * (Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_35():
    f = cosh(a + b*x**n)
    F = (Integer(-1) * (((sympy.E)**(Symbol('a')) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_36():
    f = cosh(a + b*x**n)/x
    F = ((sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))) + ((sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_37():
    f = cosh(a + b*x**n)**2
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * (Symbol('n'))**(Integer(-1))))) * (sympy.E)**((Integer(2) * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(-2) * Symbol('b') * (x)**(Symbol('n'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * (Symbol('n'))**(Integer(-1))))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_38():
    f = cosh(a + b*x**n)**2/x
    F = ((sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (sympy.log(x) * (Integer(2))**(Integer(-1))) + ((sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_39():
    f = cosh(a + b*x**n)**3
    F = (Integer(-1) * (((sympy.E)**((Integer(3) * Symbol('a'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(-3) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**((Symbol('n'))**(Integer(-1))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('a')) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**((Symbol('n'))**(Integer(-1))) * (sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(8) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_40():
    f = cosh(a + b*x**n)**3/x
    F = ((Integer(3) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((Integer(3) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_41():
    f = (b*cosh(c + d*x**n))**p*(e*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') * x))**(Symbol('m')) * ((Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * (x)**(Symbol('n')))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_42():
    f = (e*x)**m*(a + b*cosh(c + d*x**n))**p
    F = sympy.Function('Unintegrable')((((Symbol('e') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * (x)**(Symbol('n'))))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_43():
    f = (b*cosh(c + d*x**n))**p*(e*x)**(n - 1)
    F = -(b*cosh(c + d*x**n))**(p + 1)*(e*x)**n*sinh(c + d*x**n)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), cosh(c + d*x**n)**2)/(b*d*e*n*x**n*sqrt(-sinh(c + d*x**n)**2)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_44():
    f = (b*cosh(c + d*x**n))**p*(e*x)**(2*n - 1)
    F = (((Symbol('e') * x))**((Integer(2) * Symbol('n'))) * sympy.Function('Unintegrable')(((x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * (x)**(Symbol('n')))))))**(Symbol('p'))), x)) * (((x)**((Integer(2) * Symbol('n'))) * Symbol('e')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_45():
    f = (e*x)**(n - 1)*(a + b*cosh(c + d*x**n))**p
    F = sqrt(2)*(e*x)**n*(a + b*cosh(c + d*x**n))**p*sinh(c + d*x**n)*appellf1(sympy.S.Half, sympy.S.Half, -p, sympy.S(3)/2, sympy.S.Half - cosh(c + d*x**n)/2, b*(1 - cosh(c + d*x**n))/(a + b))/(d*e*n*x**n*((a + b*cosh(c + d*x**n))/(a + b))**p*sqrt(cosh(c + d*x**n) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_46():
    f = (e*x)**(2*n - 1)*(a + b*cosh(c + d*x**n))**p
    F = (((Symbol('e') * x))**((Integer(2) * Symbol('n'))) * sympy.Function('Unintegrable')(((x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * (x)**(Symbol('n'))))))))**(Symbol('p'))), x)) * (((x)**((Integer(2) * Symbol('n'))) * Symbol('e')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_47():
    f = x**m*cosh(a + b*x**n)
    F = (Integer(-1) * (((sympy.E)**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_48():
    f = x**m*cosh(a + b*x**n)**2
    F = ((x)**((Integer(1) + Symbol('m'))) * ((Integer(2) * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-2) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(2))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('n'))) * (Symbol('n'))**(Integer(-1)))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(2))**(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('n'))) * (Symbol('n'))**(Integer(-1)))) * (sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_49():
    f = x**m*cosh(a + b*x**n)**3
    F = (Integer(-1) * (((sympy.E)**((Integer(3) * Symbol('a'))) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-3) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Symbol('b') * (x)**(Symbol('n'))))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * (((Integer(3))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(8) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_50():
    f = x**(-n - 1)*cosh(a + b*x**n)
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Symbol('n')))) * sympy.sinh(Symbol('a'))) * (Symbol('n'))**(Integer(-1))) + ((Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_51():
    f = x**(-n - 1)*cosh(a + b*x**n)**2
    F = (Integer(-1) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh((Integer(2) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n')))) * sympy.sinh((Integer(2) * Symbol('a')))) * (Symbol('n'))**(Integer(-1))) + ((Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_52():
    f = x**(-n - 1)*cosh(a + b*x**n)**3
    F = (Integer(-1) * ((Integer(3) * sympy.cosh((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh((Integer(3) * (Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * (x)**(Symbol('n')))) * sympy.sinh(Symbol('a'))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n')))) * sympy.sinh((Integer(3) * Symbol('a')))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * (x)**(Symbol('n'))))) * ((Integer(4) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_53():
    f = x**(n/2 - 1)*cosh(a + b*x**n)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1))))))) * (((sympy.E)**(Symbol('a')) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('n'))))**(Integer(-1))) + (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_54():
    f = x**2*cosh((a + b*x)**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sinh(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sinh(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_55():
    f = x*cosh((a + b*x)**2)
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.sinh(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_56():
    f = cosh((a + b*x)**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_57():
    f = cosh((a + b*x)**2)/x
    F = Symbol('b') * sympy.Function('CannotIntegrate')((sympy.cosh(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b') * x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_58():
    f = x**2*cosh(a + b*sqrt(c + d*x))
    F = 2*c**2*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b*d**3) - 4*c*(c + d*x)**(sympy.S(3)/2)*sinh(a + b*sqrt(c + d*x))/(b*d**3) + 2*(c + d*x)**(sympy.S(5)/2)*sinh(a + b*sqrt(c + d*x))/(b*d**3) - 2*c**2*cosh(a + b*sqrt(c + d*x))/(b**2*d**3) + 12*c*(c + d*x)*cosh(a + b*sqrt(c + d*x))/(b**2*d**3) - 10*(c + d*x)**2*cosh(a + b*sqrt(c + d*x))/(b**2*d**3) - 24*c*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b**3*d**3) + 40*(c + d*x)**(sympy.S(3)/2)*sinh(a + b*sqrt(c + d*x))/(b**3*d**3) + 24*c*cosh(a + b*sqrt(c + d*x))/(b**4*d**3) - (120*c + 120*d*x)*cosh(a + b*sqrt(c + d*x))/(b**4*d**3) + 240*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b**5*d**3) - 240*cosh(a + b*sqrt(c + d*x))/(b**6*d**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_59():
    f = x*cosh(a + b*sqrt(c + d*x))
    F = -2*c*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b*d**2) + 2*(c + d*x)**(sympy.S(3)/2)*sinh(a + b*sqrt(c + d*x))/(b*d**2) + 2*c*cosh(a + b*sqrt(c + d*x))/(b**2*d**2) - (6*c + 6*d*x)*cosh(a + b*sqrt(c + d*x))/(b**2*d**2) + 12*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b**3*d**2) - 12*cosh(a + b*sqrt(c + d*x))/(b**4*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_60():
    f = cosh(a + b*sqrt(c + d*x))
    F = 2*sqrt(c + d*x)*sinh(a + b*sqrt(c + d*x))/(b*d) - 2*cosh(a + b*sqrt(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_61():
    f = cosh(a + b*sqrt(c + d*x))/x
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('CoshIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))) + (sympy.cosh((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('CoshIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('SinhIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))))) + (sympy.sinh((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('SinhIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_62():
    f = cosh(a + b*sqrt(c + d*x))/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('CoshIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))) * sympy.sinh((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c'))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('CoshIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) * sympy.sinh((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c')))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.cosh((Symbol('a') + (Symbol('b') * sympy.sqrt(Symbol('c'))))) * sympy.Function('SinhIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.cosh((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('c')))))) * sympy.Function('SinhIntegral')((Symbol('b') * (sympy.sqrt(Symbol('c')) + sympy.sqrt((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_63():
    f = x**2*cosh(a + b*(c + d*x)**(sympy.S(1)/3))
    F = 3*c**2*(c + d*x)**(sympy.S(2)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) - 6*c*(c + d*x)**(sympy.S(5)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) + 3*(c + d*x)**(sympy.S(8)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**3) - 6*c**2*(c + d*x)**(sympy.S(1)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) + 30*c*(c + d*x)**(sympy.S(4)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) - 24*(c + d*x)**(sympy.S(7)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**3) + 6*c**2*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) - 120*c*(c + d*x)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) + 168*(c + d*x)**2*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**3) + 360*c*(c + d*x)**(sympy.S(2)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**3) - 1008*(c + d*x)**(sympy.S(5)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**3) - 720*c*(c + d*x)**(sympy.S(1)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**3) + 5040*(c + d*x)**(sympy.S(4)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**3) + 720*c*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**3) - (20160*c + 20160*d*x)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**3) + 60480*(c + d*x)**(sympy.S(2)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**7*d**3) - 120960*(c + d*x)**(sympy.S(1)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**8*d**3) + 120960*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**9*d**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_64():
    f = x*cosh(a + b*(c + d*x)**(sympy.S(1)/3))
    F = -3*c*(c + d*x)**(sympy.S(2)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**2) + 3*(c + d*x)**(sympy.S(5)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d**2) + 6*c*(c + d*x)**(sympy.S(1)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**2) - 15*(c + d*x)**(sympy.S(4)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d**2) - 6*c*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**2) + (60*c + 60*d*x)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d**2) - 180*(c + d*x)**(sympy.S(2)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**4*d**2) + 360*(c + d*x)**(sympy.S(1)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**5*d**2) - 360*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**6*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_65():
    f = cosh(a + b*(c + d*x)**(sympy.S(1)/3))
    F = 3*(c + d*x)**(sympy.S(2)/3)*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b*d) - 6*(c + d*x)**(sympy.S(1)/3)*cosh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**2*d) + 6*sinh(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_66():
    f = cosh(a + b*(c + d*x)**(sympy.S(1)/3))/x
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((Symbol('b') * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (sympy.cosh((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('CoshIntegral')(((Integer(-1) * Symbol('b')) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) + (sympy.cosh((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('CoshIntegral')((Symbol('b') * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((Symbol('b') * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))))) + (Integer(-1) * (sympy.sinh((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((Symbol('b') * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))))) + (sympy.sinh((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('SinhIntegral')((Symbol('b') * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_3_e_x_pow_m_a_plus_b_cosh_c_plus_d_x_pow_n_pow_p_67():
    f = cosh(a + b*(c + d*x)**(sympy.S(1)/3))/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))) * (x)**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('CoshIntegral')((Symbol('b') * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * sympy.sinh((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('d') * sympy.Function('CoshIntegral')((Symbol('b') * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * sympy.Function('CoshIntegral')(((Integer(-1) * Symbol('b')) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * sympy.sinh((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.cosh((Symbol('a') + (Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((Symbol('b') * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * Symbol('d') * sympy.cosh((Symbol('a') + ((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((Symbol('b') * (((Integer(-1))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1)))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('d') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Integer(-1))**((Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))))))) * sympy.Function('SinhIntegral')((Symbol('b') * (((Integer(-1))**((Integer(3))**(Integer(-1))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('c') + (Symbol('d') * x)))**((Integer(3))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F

