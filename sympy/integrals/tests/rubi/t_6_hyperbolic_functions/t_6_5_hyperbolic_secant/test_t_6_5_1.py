"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.5 Hyperbolic secant/6.5.1 (c+d x)^m (a+b sech)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_1():
    f = (c + d*x)**3*sech(a + b*x)
    F = ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_2():
    f = (c + d*x)**2*sech(a + b*x)
    F = ((Integer(2) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_3():
    f = (c + d*x)*sech(a + b*x)
    F = ((Integer(2) * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_4():
    f = sech(a + b*x)/(c + d*x)
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_5():
    f = (c + d*x)**3*sech(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_6():
    f = (c + d*x)**2*sech(a + b*x)**2
    F = (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_7():
    f = (c + d*x)*sech(a + b*x)**2
    F = (c + d*x)*tanh(a + b*x)/b - d*log(cosh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_8():
    f = sech(a + b*x)**2/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_9():
    f = (c + d*x)**3*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_10():
    f = (c + d*x)**2*sech(a + b*x)**3
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_11():
    f = (c + d*x)*sech(a + b*x)**3
    F = (((Symbol('c') + (Symbol('d') * x)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('d') * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_12():
    f = sech(a + b*x)**3/(c + d*x)
    F = sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_13():
    f = -x*sqrt(sech(x))/3 + x/sech(x)**(sympy.S(3)/2)
    F = 2*x*sinh(x)/(3*sqrt(sech(x))) - 4/(9*sech(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_14():
    f = -3*x/(5*sqrt(sech(x))) + x/sech(x)**(sympy.S(5)/2)
    F = 2*x*sinh(x)/(5*sech(x)**(sympy.S(3)/2)) - 4/(25*sech(x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_15():
    f = -5*x*sqrt(sech(x))/21 + x/sech(x)**(sympy.S(7)/2)
    F = 10*x*sinh(x)/(21*sqrt(sech(x))) + 2*x*sinh(x)/(7*sech(x)**(sympy.S(5)/2)) - 20/(63*sech(x)**(sympy.S(3)/2)) - 4/(49*sech(x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_1_c_plus_d_x_pow_m_a_plus_b_sech_pow_n_16():
    f = -x**2*sqrt(sech(x))/3 + x**2/sech(x)**(sympy.S(3)/2)
    F = 2*x**2*sinh(x)/(3*sqrt(sech(x))) - 8*x/(9*sech(x)**(sympy.S(3)/2)) + 16*sinh(x)/(27*sqrt(sech(x))) - 16*I*sqrt(cosh(x))*elliptic_f(I*x/2, 2)*sqrt(sech(x))/27
    assert integrate(f, x) == F

