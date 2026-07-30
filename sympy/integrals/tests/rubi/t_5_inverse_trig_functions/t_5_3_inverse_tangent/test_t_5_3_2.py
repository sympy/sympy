"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.3 Inverse tangent/5.3.2 (d x)^m (a+b arctan(c x^n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n, p = symbols('a b c d m n p')

def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_1():
    f = x**5*(a + b*atan(c*x))
    F = -b*x**5/(30*c) + b*x**3/(18*c**3) - b*x/(6*c**5) + b*atan(c*x)/(6*c**6) + x**6*(a + b*atan(c*x))/6
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_2():
    f = x**4*(a + b*atan(c*x))
    F = -b*x**4/(20*c) + b*x**2/(10*c**3) - b*log(c**2*x**2 + 1)/(10*c**5) + x**5*(a + b*atan(c*x))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_3():
    f = x**3*(a + b*atan(c*x))
    F = -b*x**3/(12*c) + b*x/(4*c**3) - b*atan(c*x)/(4*c**4) + x**4*(a + b*atan(c*x))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_4():
    f = x**2*(a + b*atan(c*x))
    F = -b*x**2/(6*c) + b*log(c**2*x**2 + 1)/(6*c**3) + x**3*(a + b*atan(c*x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_5():
    f = x*(a + b*atan(c*x))
    F = -b*x/(2*c) + b*atan(c*x)/(2*c**2) + x**2*(a + b*atan(c*x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_6():
    f = a + b*atan(c*x)
    F = a*x + b*x*atan(c*x) - b*log(c**2*x**2 + 1)/(2*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_7():
    f = (a + b*atan(c*x))/x
    F = (Symbol('a') * sympy.log(x)) + ((sympy.I * (Integer(2))**(Integer(-1))) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * x))) + (Integer(-1) * ((sympy.I * (Integer(2))**(Integer(-1))) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * x))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_8():
    f = (a + b*atan(c*x))/x**2
    F = b*c*log(x) - b*c*log(c**2*x**2 + 1)/2 - (a + b*atan(c*x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_9():
    f = (a + b*atan(c*x))/x**3
    F = -b*c**2*atan(c*x)/2 - b*c/(2*x) - (a + b*atan(c*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_10():
    f = (a + b*atan(c*x))/x**4
    F = -b*c**3*log(x)/3 + b*c**3*log(c**2*x**2 + 1)/6 - b*c/(6*x**2) - (a + b*atan(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_11():
    f = (a + b*atan(c*x))/x**5
    F = b*c**4*atan(c*x)/4 + b*c**3/(4*x) - b*c/(12*x**3) - (a + b*atan(c*x))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_12():
    f = (a + b*atan(c*x))/x**6
    F = b*c**5*log(x)/5 - b*c**5*log(c**2*x**2 + 1)/10 + b*c**3/(10*x**2) - b*c/(20*x**4) - (a + b*atan(c*x))/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_13():
    f = x**5*(a + b*atan(c*x))**2
    F = -a*b*x/(3*c**5) + b**2*x**4/(60*c**2) - 4*b**2*x**2/(45*c**4) - b**2*x*atan(c*x)/(3*c**5) + 23*b**2*log(c**2*x**2 + 1)/(90*c**6) - b*x**5*(a + b*atan(c*x))/(15*c) + b*x**3*(a + b*atan(c*x))/(9*c**3) + x**6*(a + b*atan(c*x))**2/6 + (a + b*atan(c*x))**2/(6*c**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_14():
    f = x**4*(a + b*atan(c*x))**2
    F = ((Integer(-3) * (Symbol('b'))**(Integer(2)) * x) * ((Integer(10) * (Symbol('c'))**(Integer(4))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(30) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.atan((Symbol('c') * x))) * ((Integer(10) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(5) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(10) * Symbol('c')))**(Integer(-1)))) + (((sympy.I * (Integer(5))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('c'))**(Integer(5)))**(Integer(-1))) + (((x)**(Integer(5)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * (Integer(5))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Integer(5) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + (((sympy.I * (Integer(5))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_15():
    f = x**3*(a + b*atan(c*x))**2
    F = a*b*x/(2*c**3) + b**2*x**2/(12*c**2) + b**2*x*atan(c*x)/(2*c**3) - b**2*log(c**2*x**2 + 1)/(3*c**4) - b*x**3*(a + b*atan(c*x))/(6*c) + x**4*(a + b*atan(c*x))**2/4 - (a + b*atan(c*x))**2/(4*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_16():
    f = x**2*(a + b*atan(c*x))**2
    F = (((Symbol('b'))**(Integer(2)) * x) * ((Integer(3) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.atan((Symbol('c') * x))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(3) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((sympy.I * (Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_17():
    f = x*(a + b*atan(c*x))**2
    F = -a*b*x/c - b**2*x*atan(c*x)/c + b**2*log(c**2*x**2 + 1)/(2*c**2) + x**2*(a + b*atan(c*x))**2/2 + (a + b*atan(c*x))**2/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_18():
    f = (a + b*atan(c*x))**2
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * (Symbol('c'))**(Integer(-1))) + (x * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) + ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_19():
    f = (a + b*atan(c*x))**2/x
    F = (Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))))) + (sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * (Integer(2))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_20():
    f = (a + b*atan(c*x))**2/x**2
    F = ((Integer(-1) * sympy.I) * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(2) * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_21():
    f = (a + b*atan(c*x))**2/x**3
    F = b**2*c**2*log(x) - b**2*c**2*log(c**2*x**2 + 1)/2 - b*c*(a + b*atan(c*x))/x - c**2*(a + b*atan(c*x))**2/2 - (a + b*atan(c*x))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_22():
    f = (a + b*atan(c*x))**2/x**4
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)))) * ((Integer(3) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.atan((Symbol('c') * x))) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) * (Integer(3))**(Integer(-1)))) + ((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_23():
    f = (a + b*atan(c*x))**2/x**5
    F = -2*b**2*c**4*log(x)/3 + b**2*c**4*log(c**2*x**2 + 1)/3 - b**2*c**2/(12*x**2) + b*c**3*(a + b*atan(c*x))/(2*x) - b*c*(a + b*atan(c*x))/(6*x**3) + c**4*(a + b*atan(c*x))**2/4 - (a + b*atan(c*x))**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_24():
    f = x**5*(a + b*atan(c*x))**3
    F = ((Integer(19) * (Symbol('b'))**(Integer(3)) * x) * ((Integer(60) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (x)**(Integer(3))) * ((Integer(60) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(19) * (Symbol('b'))**(Integer(3)) * sympy.atan((Symbol('c') * x))) * ((Integer(60) * (Symbol('c'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(15) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(20) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Integer(23) * sympy.I) * (Integer(30))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('c'))**(Integer(6)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * x * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(5))))**(Integer(-1)))) + ((Symbol('b') * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(5)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(10) * Symbol('c')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(6) * (Symbol('c'))**(Integer(6))))**(Integer(-1))) + (((x)**(Integer(6)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(6))**(Integer(-1))) + (Integer(-1) * ((Integer(23) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Integer(15) * (Symbol('c'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(23) * sympy.I) * (Integer(30))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(6)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_25():
    f = x**4*(a + b*atan(c*x))**3
    F = ((Integer(-9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * x) * ((Integer(10) * (Symbol('c'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (x)**(Integer(2))) * ((Integer(20) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('b'))**(Integer(3)) * x * sympy.atan((Symbol('c') * x))) * ((Integer(10) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(10) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(9) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(20) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(10) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (x)**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(20) * Symbol('c')))**(Integer(-1)))) + (((sympy.I * (Integer(5))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * ((Symbol('c'))**(Integer(5)))**(Integer(-1))) + (((x)**(Integer(5)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(5))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Integer(5) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('c'))**(Integer(5))))**(Integer(-1))) + ((((Integer(3) * sympy.I) * (Integer(5))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(5)))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Integer(10) * (Symbol('c'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_26():
    f = x**3*(a + b*atan(c*x))**3
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(3)) * x)) * ((Integer(4) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.atan((Symbol('c') * x))) * ((Integer(4) * (Symbol('c'))**(Integer(4))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(4) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(4))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_27():
    f = x**2*(a + b*atan(c*x))**3
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * x) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * x * sympy.atan((Symbol('c') * x))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((sympy.I * (Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_28():
    f = x*(a + b*atan(c*x))**3
    F = ((((Integer(-3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_29():
    f = (a + b*atan(c*x))**3
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Symbol('c'))**(Integer(-1))) + (x * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))) + (((Integer(3) * sympy.I) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_30():
    f = (a + b*atan(c*x))**3/x
    F = (Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))))) + (((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) * (Integer(2))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))) * (Integer(2))**(Integer(-1))) + (((Integer(3) * sympy.I) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1))))))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_31():
    f = (a + b*atan(c*x))**3/x**2
    F = ((Integer(-1) * sympy.I) * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * sympy.I) * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_32():
    f = (a + b*atan(c*x))**3/x**3
    F = (((Integer(-3) * sympy.I) * (Integer(2))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) + (Integer(-1) * (((Integer(3) * sympy.I) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_33():
    f = (a + b*atan(c*x))**3/x**4
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.log(x)) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))))) + (sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_34():
    f = (a + b*atan(c*x))**3/x**5
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * sympy.atan((Symbol('c') * x))) * (Integer(4))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (sympy.I * Symbol('b') * (Symbol('c'))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * x))**(Integer(-1))) + (((Symbol('c'))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))))) + (sympy.I * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_35():
    f = x/atan(a*x)
    F = sympy.Function('Unintegrable')((x * (sympy.atan((Symbol('a') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_36():
    f = 1/atan(a*x)
    F = sympy.Function('Unintegrable')((sympy.atan((Symbol('a') * x)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_37():
    f = 1/(x*atan(a*x))
    F = sympy.Function('Unintegrable')(((x * sympy.atan((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_38():
    f = x/atan(a*x)**2
    F = sympy.Function('Unintegrable')((x * ((sympy.atan((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_39():
    f = atan(a*x)**(-2)
    F = sympy.Function('Unintegrable')(((sympy.atan((Symbol('a') * x)))**(Integer(2)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_40():
    f = 1/(x*atan(a*x)**2)
    F = sympy.Function('Unintegrable')(((x * (sympy.atan((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_41():
    f = x*sqrt(atan(a*x))
    F = sympy.Function('Unintegrable')((x * sympy.sqrt(sympy.atan((Symbol('a') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_42():
    f = sqrt(atan(a*x))
    F = sympy.Function('Unintegrable')(sympy.sqrt(sympy.atan((Symbol('a') * x))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_43():
    f = sqrt(atan(a*x))/x
    F = sympy.Function('Unintegrable')((sympy.sqrt(sympy.atan((Symbol('a') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_44():
    f = x*atan(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((x * (sympy.atan((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_45():
    f = atan(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((sympy.atan((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_46():
    f = atan(a*x)**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.atan((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_47():
    f = x/sqrt(atan(a*x))
    F = sympy.Function('Unintegrable')((x * (sympy.sqrt(sympy.atan((Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_48():
    f = 1/sqrt(atan(a*x))
    F = sympy.Function('Unintegrable')((sympy.sqrt(sympy.atan((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_49():
    f = 1/(x*sqrt(atan(a*x)))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.atan((Symbol('a') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_50():
    f = x/atan(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((x * ((sympy.atan((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_51():
    f = atan(a*x)**(sympy.S(-3)/2)
    F = sympy.Function('Unintegrable')((sympy.atan((Symbol('a') * x)))**((Integer(-3) * (Integer(2))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_52():
    f = 1/(x*atan(a*x)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.atan((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_53():
    f = sqrt(x)*atan(x)
    F = 2*x**(sympy.S(3)/2)*atan(x)/3 - 4*sqrt(x)/3 - sqrt(2)*log(-sqrt(2)*sqrt(x) + x + 1)/6 + sqrt(2)*log(sqrt(2)*sqrt(x) + x + 1)/6 + sqrt(2)*atan(sqrt(2)*sqrt(x) - 1)/3 + sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_54():
    f = (d*x)**m*(a + b*atan(c*x))**3
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_55():
    f = (d*x)**m*(a + b*atan(c*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_56():
    f = (d*x)**m*(a + b*atan(c*x))
    F = -b*c*(d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -c**2*x**2)/(d**2*(m + 1)*(m + 2)) + (d*x)**(m + 1)*(a + b*atan(c*x))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_57():
    f = (d*x)**m/(a + b*atan(c*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_58():
    f = (a + b*atan(c*x))**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_59():
    f = (d*x)**m*(a + b*atan(c*x))**p
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * x)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_60():
    f = x**7*(a + b*atan(c*x**2))
    F = -b*x**6/(24*c) + b*x**2/(8*c**3) - b*atan(c*x**2)/(8*c**4) + x**8*(a + b*atan(c*x**2))/8
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_61():
    f = x**5*(a + b*atan(c*x**2))
    F = -b*x**4/(12*c) + b*log(c**2*x**4 + 1)/(12*c**3) + x**6*(a + b*atan(c*x**2))/6
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_62():
    f = x**3*(a + b*atan(c*x**2))
    F = -b*x**2/(4*c) + b*atan(c*x**2)/(4*c**2) + x**4*(a + b*atan(c*x**2))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_63():
    f = x*(a + b*atan(c*x**2))
    F = -b*log(c**2*x**4 + 1)/(4*c) + x**2*(a + b*atan(c*x**2))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_64():
    f = (a + b*atan(c*x**2))/x
    F = (Symbol('a') * sympy.log(x)) + ((Integer(4))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (x)**(Integer(2))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_65():
    f = (a + b*atan(c*x**2))/x**3
    F = b*c*log(x) - b*c*log(c**2*x**4 + 1)/4 - (a + b*atan(c*x**2))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_66():
    f = (a + b*atan(c*x**2))/x**5
    F = -b*c**2*atan(c*x**2)/4 - b*c/(4*x**2) - (a + b*atan(c*x**2))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_67():
    f = (a + b*atan(c*x**2))/x**7
    F = -b*c**3*log(x)/3 + b*c**3*log(c**2*x**4 + 1)/12 - b*c/(12*x**4) - (a + b*atan(c*x**2))/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_68():
    f = x**4*(a + b*atan(c*x**2))
    F = -2*b*x**3/(15*c) + sqrt(2)*b*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(20*c**(sympy.S(5)/2)) - sqrt(2)*b*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(20*c**(sympy.S(5)/2)) + sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x - 1)/(10*c**(sympy.S(5)/2)) + sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x + 1)/(10*c**(sympy.S(5)/2)) + x**5*(a + b*atan(c*x**2))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_69():
    f = x**2*(a + b*atan(c*x**2))
    F = -2*b*x/(3*c) - sqrt(2)*b*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(12*c**(sympy.S(3)/2)) + sqrt(2)*b*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(12*c**(sympy.S(3)/2)) + sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x - 1)/(6*c**(sympy.S(3)/2)) + sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x + 1)/(6*c**(sympy.S(3)/2)) + x**3*(a + b*atan(c*x**2))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_70():
    f = a + b*atan(c*x**2)
    F = a*x + b*x*atan(c*x**2) - sqrt(2)*b*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(4*sqrt(c)) + sqrt(2)*b*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/(4*sqrt(c)) - sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x - 1)/(2*sqrt(c)) - sqrt(2)*b*atan(sqrt(2)*sqrt(c)*x + 1)/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_71():
    f = (a + b*atan(c*x**2))/x**2
    F = -sqrt(2)*b*sqrt(c)*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/4 + sqrt(2)*b*sqrt(c)*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/4 + sqrt(2)*b*sqrt(c)*atan(sqrt(2)*sqrt(c)*x - 1)/2 + sqrt(2)*b*sqrt(c)*atan(sqrt(2)*sqrt(c)*x + 1)/2 - (a + b*atan(c*x**2))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_72():
    f = (a + b*atan(c*x**2))/x**4
    F = -sqrt(2)*b*c**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/12 + sqrt(2)*b*c**(sympy.S(3)/2)*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/12 - sqrt(2)*b*c**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(c)*x - 1)/6 - sqrt(2)*b*c**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(c)*x + 1)/6 - 2*b*c/(3*x) - (a + b*atan(c*x**2))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_73():
    f = (a + b*atan(c*x**2))/x**6
    F = sqrt(2)*b*c**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(c)*x + c*x**2 + 1)/20 - sqrt(2)*b*c**(sympy.S(5)/2)*log(sqrt(2)*sqrt(c)*x + c*x**2 + 1)/20 - sqrt(2)*b*c**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(c)*x - 1)/10 - sqrt(2)*b*c**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(c)*x + 1)/10 - 2*b*c/(15*x**3) - (a + b*atan(c*x**2))/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_74():
    f = x**7*(a + b*atan(c*x**2))**2
    F = a*b*x**2/(4*c**3) + b**2*x**4/(24*c**2) + b**2*x**2*atan(c*x**2)/(4*c**3) - b**2*log(c**2*x**4 + 1)/(6*c**4) - b*x**6*(a + b*atan(c*x**2))/(12*c) + x**8*(a + b*atan(c*x**2))**2/8 - (a + b*atan(c*x**2))**2/(8*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_75():
    f = x**5*(a + b*atan(c*x**2))**2
    F = (((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(6) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.atan((Symbol('c') * (x)**(Integer(2))))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2))))))) * ((Integer(6) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) + (Integer(-1) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_76():
    f = x**3*(a + b*atan(c*x**2))**2
    F = -a*b*x**2/(2*c) - b**2*x**2*atan(c*x**2)/(2*c) + b**2*log(c**2*x**4 + 1)/(4*c**2) + x**4*(a + b*atan(c*x**2))**2/4 + (a + b*atan(c*x**2))**2/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_77():
    f = x*(a + b*atan(c*x**2))**2
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) + ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_78():
    f = (a + b*atan(c*x**2))**2/x
    F = (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_79():
    f = (a + b*atan(c*x**2))**2/x**3
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_80():
    f = (a + b*atan(c*x**2))**2/x**5
    F = b**2*c**2*log(x) - b**2*c**2*log(c**2*x**4 + 1)/4 - b*c*(a + b*atan(c*x**2))/(2*x**2) - c**2*(a + b*atan(c*x**2))**2/4 - (a + b*atan(c*x**2))**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_81():
    f = x**2*(a + b*atan(c*x**2))**2
    F = (Integer(-1) * ((Integer(4) * Symbol('a') * Symbol('b') * x) * ((Integer(3) * Symbol('c')))**(Integer(-1)))) + ((Integer(2) * (Integer(9))**(Integer(-1))) * sympy.I * Symbol('a') * Symbol('b') * (x)**(Integer(3))) + ((Integer(4) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('a') * Symbol('b') * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log(((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(9))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('b') * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(12))**(Integer(-1)) * (x)**(Integer(3)) * (((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))**(Integer(2))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * Symbol('a') * Symbol('b') * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(12))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * (sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(2)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(6) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) * ((Integer(6) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(6) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(6) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_82():
    f = (a + b*atan(c*x**2))**2
    F = ((Symbol('a'))**(Integer(2)) * x) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('a') * Symbol('b') * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + ((Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('a') * Symbol('b') * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + ((Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log(((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + ((Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (sympy.I * Symbol('a') * Symbol('b') * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * x * (sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)))) + (Integer(-1) * (sympy.I * Symbol('a') * Symbol('b') * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * x * (sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(2)))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_83():
    f = (a + b*atan(c*x**2))**2/x**2
    F = ((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * (sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('a') * Symbol('b') * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * (sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2)))) + (Integer(-1) * (Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1)))))) + (Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log(((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + (Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) + (Integer(-1) * (Integer(2) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('b') * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))) + (Integer(-1) * ((((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))**(Integer(2)) * ((Integer(4) * x))**(Integer(-1)))) + ((sympy.I * Symbol('a') * Symbol('b') * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(2))) * ((Integer(4) * x))**(Integer(-1))) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_84():
    f = (a + b*atan(c*x**2))**2/x**4
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) + (Integer(-1) * ((Integer(4) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log(((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('b') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))) + (Integer(-1) * ((((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))**(Integer(2)) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * Symbol('b') * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(6) * (x)**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(2))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))))) + ((Integer(3))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_85():
    f = (a + b*atan(c*x**2))**2/x**6
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('a') * Symbol('b') * (Symbol('c'))**(Integer(2))) * ((Integer(5) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) * ((Integer(15) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(15))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2)))) + ((Integer(2) * (Integer(5))**(Integer(-1))) * (Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) + ((Integer(4) * (Integer(15))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) + ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(2))) + ((Integer(2) * (Integer(5))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * (Integer(5))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log(((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * (Integer(5))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1)))))) + ((Integer(2) * (Integer(5))**(Integer(-1))) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) * ((Integer(5) * x))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))))) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))))) * ((Integer(5) * x))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * ((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))))) + (Integer(-1) * ((((Integer(2) * Symbol('a')) + (sympy.I * Symbol('b') * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))))))**(Integer(2)) * ((Integer(20) * (x)**(Integer(5))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * Symbol('b') * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))) * ((Integer(10) * (x)**(Integer(5))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.log((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(2))) * ((Integer(20) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + ((Integer(10))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(4))**(Integer(-1))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))) + ((Integer(10))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((sympy.sqrt(Integer(2)) * ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))) + ((Integer(10))**(Integer(-1)) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * (Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) + ((Integer(10))**(Integer(-1)) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(1) + ((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * x))) * ((Integer(1) + ((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_86():
    f = x**3*(a + b*atan(c*x**2))**3
    F = (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3)) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_87():
    f = x*(a + b*atan(c*x**2))**3
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_88():
    f = (a + b*atan(c*x**2))**3/x
    F = (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))) + ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_89():
    f = (a + b*atan(c*x**2))**3/x**3
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_90():
    f = (a + b*atan(c*x**2))**3/x**5
    F = ((Integer(-1) * (Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(2))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_91():
    f = (d*x)**m*(a + b*atan(c*x**2))**3
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_92():
    f = (d*x)**m*(a + b*atan(c*x**2))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_93():
    f = (d*x)**m*(a + b*atan(c*x**2))
    F = -2*b*c*(d*x)**(m + 3)*hyper((1, m/4 + sympy.S(3)/4), (m/4 + sympy.S(7)/4,), -c**2*x**4)/(d**3*(m + 1)*(m + 3)) + (d*x)**(m + 1)*(a + b*atan(c*x**2))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_94():
    f = (d*x)**m/(a + b*atan(c*x**2))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_95():
    f = (d*x)**m/(a + b*atan(c*x**2))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(2)))))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_96():
    f = x**11*(a + b*atan(c*x**3))
    F = -b*x**9/(36*c) + b*x**3/(12*c**3) - b*atan(c*x**3)/(12*c**4) + x**12*(a + b*atan(c*x**3))/12
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_97():
    f = x**8*(a + b*atan(c*x**3))
    F = -b*x**6/(18*c) + b*log(c**2*x**6 + 1)/(18*c**3) + x**9*(a + b*atan(c*x**3))/9
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_98():
    f = x**5*(a + b*atan(c*x**3))
    F = -b*x**3/(6*c) + b*atan(c*x**3)/(6*c**2) + x**6*(a + b*atan(c*x**3))/6
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_99():
    f = x**2*(a + b*atan(c*x**3))
    F = -b*log(c**2*x**6 + 1)/(6*c) + x**3*(a + b*atan(c*x**3))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_100():
    f = (a + b*atan(c*x**3))/x
    F = (Symbol('a') * sympy.log(x)) + ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (x)**(Integer(3))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (x)**(Integer(3))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_101():
    f = (a + b*atan(c*x**3))/x**4
    F = b*c*log(x) - b*c*log(c**2*x**6 + 1)/6 - (a + b*atan(c*x**3))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_102():
    f = (a + b*atan(c*x**3))/x**7
    F = -b*c**2*atan(c*x**3)/6 - b*c/(6*x**3) - (a + b*atan(c*x**3))/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_103():
    f = (a + b*atan(c*x**3))/x**10
    F = -b*c**3*log(x)/3 + b*c**3*log(c**2*x**6 + 1)/18 - b*c/(18*x**6) - (a + b*atan(c*x**3))/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_104():
    f = x**3*(a + b*atan(c*x**3))
    F = -3*b*x/(4*c) - sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 - sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(16*c**(sympy.S(4)/3)) + sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 + sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(16*c**(sympy.S(4)/3)) + b*atan(c**(sympy.S(1)/3)*x)/(4*c**(sympy.S(4)/3)) + b*atan(2*c**(sympy.S(1)/3)*x - sqrt(3))/(8*c**(sympy.S(4)/3)) + b*atan(2*c**(sympy.S(1)/3)*x + sqrt(3))/(8*c**(sympy.S(4)/3)) + x**4*(a + b*atan(c*x**3))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_105():
    f = a + b*atan(c*x**3)
    F = a*x + b*x*atan(c*x**3) + b*log(c**(sympy.S(2)/3)*x**2 + 1)/(2*c**(sympy.S(1)/3)) - b*log(c**(sympy.S(4)/3)*x**4 - c**(sympy.S(2)/3)*x**2 + 1)/(4*c**(sympy.S(1)/3)) + sqrt(3)*b*atan(sqrt(3)*(-2*c**(sympy.S(2)/3)*x**2 + 1)/3)/(2*c**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_106():
    f = (a + b*atan(c*x**3))/x**3
    F = -sqrt(3)*b*c**(sympy.S(2)/3)*log(c**(sympy.S(2)/3)*x**2 - sqrt(3)*c**(sympy.S(1)/3)*x + 1)/8 + sqrt(3)*b*c**(sympy.S(2)/3)*log(c**(sympy.S(2)/3)*x**2 + sqrt(3)*c**(sympy.S(1)/3)*x + 1)/8 + b*c**(sympy.S(2)/3)*atan(c**(sympy.S(1)/3)*x)/2 + b*c**(sympy.S(2)/3)*atan(2*c**(sympy.S(1)/3)*x - sqrt(3))/4 + b*c**(sympy.S(2)/3)*atan(2*c**(sympy.S(1)/3)*x + sqrt(3))/4 - (a + b*atan(c*x**3))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_107():
    f = (a + b*atan(c*x**3))/x**6
    F = b*c**(sympy.S(5)/3)*log(c**(sympy.S(2)/3)*x**2 + 1)/10 - b*c**(sympy.S(5)/3)*log(c**(sympy.S(4)/3)*x**4 - c**(sympy.S(2)/3)*x**2 + 1)/20 + sqrt(3)*b*c**(sympy.S(5)/3)*atan(sqrt(3)*(-2*c**(sympy.S(2)/3)*x**2 + 1)/3)/10 - 3*b*c/(10*x**2) - (a + b*atan(c*x**3))/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_108():
    f = x**7*(a + b*atan(c*x**3))
    F = -3*b*x**5/(40*c) + sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 - sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(32*c**(sympy.S(8)/3)) - sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 + sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(32*c**(sympy.S(8)/3)) + b*atan(c**(sympy.S(1)/3)*x)/(8*c**(sympy.S(8)/3)) + b*atan(2*c**(sympy.S(1)/3)*x - sqrt(3))/(16*c**(sympy.S(8)/3)) + b*atan(2*c**(sympy.S(1)/3)*x + sqrt(3))/(16*c**(sympy.S(8)/3)) + x**8*(a + b*atan(c*x**3))/8
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_109():
    f = x**4*(a + b*atan(c*x**3))
    F = -3*b*x**2/(10*c) + b*log(c**(sympy.S(2)/3)*x**2 + 1)/(10*c**(sympy.S(5)/3)) - b*log(c**(sympy.S(4)/3)*x**4 - c**(sympy.S(2)/3)*x**2 + 1)/(20*c**(sympy.S(5)/3)) - sqrt(3)*b*atan(sqrt(3)*(-2*c**(sympy.S(2)/3)*x**2 + 1)/3)/(10*c**(sympy.S(5)/3)) + x**5*(a + b*atan(c*x**3))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_110():
    f = x*(a + b*atan(c*x**3))
    F = -sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 - sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(8*c**(sympy.S(2)/3)) + sqrt(3)*b*log(c**(sympy.S(2)/3)*x**2 + sqrt(3)*c**(sympy.S(1)/3)*x + 1)/(8*c**(sympy.S(2)/3)) - b*atan(c**(sympy.S(1)/3)*x)/(2*c**(sympy.S(2)/3)) - b*atan(2*c**(sympy.S(1)/3)*x - sqrt(3))/(4*c**(sympy.S(2)/3)) - b*atan(2*c**(sympy.S(1)/3)*x + sqrt(3))/(4*c**(sympy.S(2)/3)) + x**2*(a + b*atan(c*x**3))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_111():
    f = (a + b*atan(c*x**3))/x**2
    F = b*c**(sympy.S(1)/3)*log(c**(sympy.S(2)/3)*x**2 + 1)/2 - b*c**(sympy.S(1)/3)*log(c**(sympy.S(4)/3)*x**4 - c**(sympy.S(2)/3)*x**2 + 1)/4 - sqrt(3)*b*c**(sympy.S(1)/3)*atan(sqrt(3)*(-2*c**(sympy.S(2)/3)*x**2 + 1)/3)/2 - (a + b*atan(c*x**3))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_112():
    f = (a + b*atan(c*x**3))/x**5
    F = -sqrt(3)*b*c**(sympy.S(4)/3)*log(c**(sympy.S(2)/3)*x**2 - sqrt(3)*c**(sympy.S(1)/3)*x + 1)/16 + sqrt(3)*b*c**(sympy.S(4)/3)*log(c**(sympy.S(2)/3)*x**2 + sqrt(3)*c**(sympy.S(1)/3)*x + 1)/16 - b*c**(sympy.S(4)/3)*atan(c**(sympy.S(1)/3)*x)/4 - b*c**(sympy.S(4)/3)*atan(2*c**(sympy.S(1)/3)*x - sqrt(3))/8 - b*c**(sympy.S(4)/3)*atan(2*c**(sympy.S(1)/3)*x + sqrt(3))/8 - 3*b*c/(4*x) - (a + b*atan(c*x**3))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_113():
    f = x**11*(a + b*atan(c*x**3))**2
    F = a*b*x**3/(6*c**3) + b**2*x**6/(36*c**2) + b**2*x**3*atan(c*x**3)/(6*c**3) - b**2*log(c**2*x**6 + 1)/(9*c**4) - b*x**9*(a + b*atan(c*x**3))/(18*c) + x**12*(a + b*atan(c*x**3))**2/12 - (a + b*atan(c*x**3))**2/(12*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_114():
    f = x**8*(a + b*atan(c*x**3))**2
    F = (((Symbol('b'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(9) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.atan((Symbol('c') * (x)**(Integer(3))))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(6)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3))))))) * ((Integer(9) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9))**(Integer(-1)) * (x)**(Integer(9)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_115():
    f = x**5*(a + b*atan(c*x**3))**2
    F = -a*b*x**3/(3*c) - b**2*x**3*atan(c*x**3)/(3*c) + b**2*log(c**2*x**6 + 1)/(6*c**2) + x**6*(a + b*atan(c*x**3))**2/6 + (a + b*atan(c*x**3))**2/(6*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_116():
    f = x**2*(a + b*atan(c*x**3))**2
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(3) * Symbol('c')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) + ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))) * ((Integer(3) * Symbol('c')))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(3) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_117():
    f = (a + b*atan(c*x**3))**2/x
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))))) + ((Integer(3))**(Integer(-1)) * sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))))) + ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_118():
    f = (a + b*atan(c*x**3))**2/x**4
    F = ((Integer(-1) * (Integer(3))**(Integer(-1))) * sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_119():
    f = (a + b*atan(c*x**3))**2/x**7
    F = b**2*c**2*log(x) - b**2*c**2*log(c**2*x**6 + 1)/6 - b*c*(a + b*atan(c*x**3))/(3*x**3) - c**2*(a + b*atan(c*x**3))**2/6 - (a + b*atan(c*x**3))**2/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_120():
    f = (a + b*atan(c*x**3))**2/x**10
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3))))))) * ((Integer(9) * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(9))**(Integer(-1)) * sympy.I * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * ((Integer(9) * (x)**(Integer(9))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Integer(9))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1)))))))) + ((Integer(9))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_121():
    f = x**8*(a + b*atan(c*x**3))**3
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') * (x)**(Integer(3))))) * ((Integer(3) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(6)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(6) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9))**(Integer(-1)) * (x)**(Integer(9)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(6)))))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_122():
    f = x**5*(a + b*atan(c*x**3))**3
    F = (Integer(-1) * ((sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3)) * ((Integer(6) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_123():
    f = x**2*(a + b*atan(c*x**3))**3
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) * ((Integer(3) * Symbol('c')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) + ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_124():
    f = (a + b*atan(c*x**3))**3/x
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + (sympy.I * Symbol('c') * (x)**(Integer(3)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_125():
    f = (a + b*atan(c*x**3))**3/x**4
    F = ((Integer(-1) * (Integer(3))**(Integer(-1))) * sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_126():
    f = (a + b*atan(c*x**3))**3/x**7
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (x)**(Integer(3))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_127():
    f = (d*x)**m*(a + b*atan(c*x**3))**3
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_128():
    f = (d*x)**m*(a + b*atan(c*x**3))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_129():
    f = (d*x)**m*(a + b*atan(c*x**3))
    F = -3*b*c*(d*x)**(m + 4)*hyper((1, m/6 + sympy.S(2)/3), (m/6 + sympy.S(5)/3,), -c**2*x**6)/(d**4*(m + 1)*(m + 4)) + (d*x)**(m + 1)*(a + b*atan(c*x**3))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_130():
    f = (d*x)**m/(a + b*atan(c*x**3))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_131():
    f = (d*x)**m/(a + b*atan(c*x**3))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.atan((Symbol('c') * (x)**(Integer(3)))))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_132():
    f = x**3*(a + b*atan(c/x))
    F = b*c**4*atan(x/c)/4 - b*c**3*x/4 + b*c*x**3/12 + x**4*(a + b*atan(c/x))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_133():
    f = x**2*(a + b*atan(c/x))
    F = -b*c**3*log(c**2 + x**2)/6 + b*c*x**2/6 + x**3*(a + b*atan(c/x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_134():
    f = x*(a + b*atan(c/x))
    F = -b*c**2*atan(x/c)/2 + b*c*x/2 + x**2*(a + b*atan(c/x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_135():
    f = a + b*atan(c/x)
    F = a*x + b*c*log(c**2 + x**2)/2 + b*x*atan(c/x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_136():
    f = (a + b*atan(c/x))/x
    F = (Symbol('a') * sympy.log(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_137():
    f = (a + b*atan(c/x))/x**2
    F = b*log(c**2/x**2 + 1)/(2*c) - (a + b*atan(c/x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_138():
    f = (a + b*atan(c/x))/x**3
    F = b/(2*c*x) + b*atan(x/c)/(2*c**2) - (a + b*atan(c/x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_139():
    f = (a + b*atan(c/x))/x**4
    F = b/(6*c*x**2) + b*log(x)/(3*c**3) - b*log(c**2 + x**2)/(6*c**3) - (a + b*atan(c/x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_140():
    f = x**3*(a + b*atan(c/x))**2
    F = -2*b**2*c**4*log(x)/3 - b**2*c**4*log(c**2/x**2 + 1)/3 + b**2*c**2*x**2/12 - b*c**3*x*(a + b*acot(x/c))/2 + b*c*x**3*(a + b*acot(x/c))/6 - c**4*(a + b*acot(x/c))**2/4 + x**4*(a + b*acot(x/c))**2/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_141():
    f = x**2*(a + b*atan(c/x))**2
    F = ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * x) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.acot((x * (Symbol('c'))**(Integer(-1))))) + ((Integer(3))**(Integer(-1)) * Symbol('b') * Symbol('c') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_142():
    f = x*(a + b*atan(c/x))**2
    F = b**2*c**2*log(x) + b**2*c**2*log(c**2/x**2 + 1)/2 + b*c*x*(a + b*acot(x/c)) + c**2*(a + b*acot(x/c))**2/2 + x**2*(a + b*acot(x/c))**2/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_143():
    f = (a + b*atan(c/x))**2
    F = (sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log(((Integer(2) * Symbol('c')) * ((Symbol('c') + (sympy.I * x)))**(Integer(-1)))))) + (sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c')) * ((Symbol('c') + (sympy.I * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_144():
    f = (a + b*atan(c/x))**2/x
    F = (Integer(-2) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_145():
    f = (a + b*atan(c/x))**2/x**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log((Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_146():
    f = (a + b*atan(c/x))**2/x**3
    F = a*b/(c*x) + b**2*acot(x/c)/(c*x) - b**2*log(c**2/x**2 + 1)/(2*c**2) - (a + b*acot(x/c))**2/(2*x**2) - (a + b*acot(x/c))**2/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_147():
    f = x**3*(a + b*atan(c/x))**3
    F = ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * x) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * sympy.acot((x * (Symbol('c'))**(Integer(-1))))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * Symbol('b') * (Symbol('c'))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(3)) * x * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * Symbol('c') * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('c'))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_148():
    f = x**2*(a + b*atan(c/x))**3
    F = ((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + (Symbol('b') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1)))))) + ((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.log(x)) + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_149():
    f = x*(a + b*atan(c/x))**3
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('b') * Symbol('c') * x * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) + ((Integer(2))**(Integer(-1)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('c')) * (x)**(Integer(-1))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_150():
    f = (a + b*atan(c/x))**3
    F = (sympy.I * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) + (Integer(-1) * (Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('c')) * ((Symbol('c') + (sympy.I * x)))**(Integer(-1)))))) + (Integer(3) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c')) * ((Symbol('c') + (sympy.I * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c')) * ((Symbol('c') + (sympy.I * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_151():
    f = (a + b*atan(c/x))**3/x
    F = (Integer(-2) * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1)))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_152():
    f = (a + b*atan(c/x))**3/x**2
    F = (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_153():
    f = (a + b*atan(c/x))**3/x**3
    F = ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(2))) * ((Integer(2) * Symbol('c') * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((x * (Symbol('c'))**(Integer(-1)))))) * sympy.log((Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * Symbol('c')) * (x)**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_154():
    f = x**2*atan(sqrt(x))
    F = -x**(sympy.S(5)/2)/15 + x**(sympy.S(3)/2)/9 - sqrt(x)/3 + x**3*atan(sqrt(x))/3 + atan(sqrt(x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_155():
    f = x*atan(sqrt(x))
    F = -x**(sympy.S(3)/2)/6 + sqrt(x)/2 + x**2*atan(sqrt(x))/2 - atan(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_156():
    f = atan(sqrt(x))
    F = -sqrt(x) + x*atan(sqrt(x)) + atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_157():
    f = atan(sqrt(x))/x
    F = (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * sympy.sqrt(x)))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * sympy.sqrt(x)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_158():
    f = atan(sqrt(x))/x**2
    F = -atan(sqrt(x)) - atan(sqrt(x))/x - 1/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_159():
    f = atan(sqrt(x))/x**3
    F = atan(sqrt(x))/2 - atan(sqrt(x))/(2*x**2) + 1/(2*sqrt(x)) - 1/(6*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_160():
    f = x**(sympy.S(3)/2)*atan(sqrt(x))
    F = 2*x**(sympy.S(5)/2)*atan(sqrt(x))/5 - x**2/10 + x/5 - log(x + 1)/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_161():
    f = sqrt(x)*atan(sqrt(x))
    F = 2*x**(sympy.S(3)/2)*atan(sqrt(x))/3 - x/3 + log(x + 1)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_162():
    f = atan(sqrt(x))/sqrt(x)
    F = 2*sqrt(x)*atan(sqrt(x)) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_163():
    f = atan(sqrt(x))/x**(sympy.S(3)/2)
    F = log(x) - log(x + 1) - 2*atan(sqrt(x))/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_164():
    f = atan(sqrt(x))/x**(sympy.S(5)/2)
    F = -log(x)/3 + log(x + 1)/3 - 1/(3*x) - 2*atan(sqrt(x))/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_165():
    f = atan(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('a') * (x)**(Integer(5))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('a') * (x)**(Integer(5))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_2_d_x_pow_m_a_plus_b_arctan_c_x_pow_n_pow_p_166():
    f = atan(a*x**n)/x
    F = ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('a') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('a') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F

