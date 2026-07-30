"""Generated from MathematicaSyntaxTestSuite.

Source: 2 Exponentials/2.1 u (F^(c (a+b x)))^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

F, a, b, c, d, e, f, g, h, m, n = symbols('F a b c d e f g h m n')

def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_1():
    f = F**(c*(a + b*x))*(d + e*x)**m
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_2():
    f = F**(c*(a + b*x))*(d + e*x)**4
    F = F**(c*(a + b*x))*(d + e*x)**4/(b*c*log(F)) - 4*F**(c*(a + b*x))*e*(d + e*x)**3/(b**2*c**2*log(F)**2) + 12*F**(c*(a + b*x))*e**2*(d + e*x)**2/(b**3*c**3*log(F)**3) - 24*F**(c*(a + b*x))*e**3*(d + e*x)/(b**4*c**4*log(F)**4) + 24*F**(c*(a + b*x))*e**4/(b**5*c**5*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_3():
    f = F**(c*(a + b*x))*(d + e*x)**3
    F = F**(c*(a + b*x))*(d + e*x)**3/(b*c*log(F)) - 3*F**(c*(a + b*x))*e*(d + e*x)**2/(b**2*c**2*log(F)**2) + 6*F**(c*(a + b*x))*e**2*(d + e*x)/(b**3*c**3*log(F)**3) - 6*F**(c*(a + b*x))*e**3/(b**4*c**4*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_4():
    f = F**(c*(a + b*x))*(d + e*x)**2
    F = F**(c*(a + b*x))*(d + e*x)**2/(b*c*log(F)) - 2*F**(c*(a + b*x))*e*(d + e*x)/(b**2*c**2*log(F)**2) + 2*F**(c*(a + b*x))*e**2/(b**3*c**3*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_5():
    f = F**(c*(a + b*x))*(d + e*x)
    F = F**(c*(a + b*x))*(d + e*x)/(b*c*log(F)) - F**(c*(a + b*x))*e/(b**2*c**2*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_6():
    f = F**(c*(a + b*x))
    F = F**(c*(a + b*x))/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_7():
    f = F**(c*(a + b*x))/(d + e*x)
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1))))) * (Symbol('e'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_8():
    f = F**(c*(a + b*x))/(d + e*x)**2
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_9():
    f = F**(c*(a + b*x))/(d + e*x)**3
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_10():
    f = F**(c*(a + b*x))/(d + e*x)**4
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_11():
    f = F**(c*(a + b*x))/(d + e*x)**5
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(12) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(24) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(24) * (Symbol('e'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(4)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(24) * (Symbol('e'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_12():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(4)))
    F = F**(c*(a + b*x))*(d + e*x)**4/(b*c*log(F)) - 4*F**(c*(a + b*x))*e*(d + e*x)**3/(b**2*c**2*log(F)**2) + 12*F**(c*(a + b*x))*e**2*(d + e*x)**2/(b**3*c**3*log(F)**3) - 24*F**(c*(a + b*x))*e**3*(d + e*x)/(b**4*c**4*log(F)**4) + 24*F**(c*(a + b*x))*e**4/(b**5*c**5*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_13():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(3)))
    F = F**(c*(a + b*x))*(d + e*x)**3/(b*c*log(F)) - 3*F**(c*(a + b*x))*e*(d + e*x)**2/(b**2*c**2*log(F)**2) + 6*F**(c*(a + b*x))*e**2*(d + e*x)/(b**3*c**3*log(F)**3) - 6*F**(c*(a + b*x))*e**3/(b**4*c**4*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_14():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(2)))
    F = F**(c*(a + b*x))*(d + e*x)**2/(b*c*log(F)) - 2*F**(c*(a + b*x))*e*(d + e*x)/(b**2*c**2*log(F)**2) + 2*F**(c*(a + b*x))*e**2/(b**3*c**3*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_15():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_16():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_17():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(6) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_18():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(5))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(12) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(24) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(24) * (Symbol('e'))**(Integer(4)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(4)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(24) * (Symbol('e'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_19():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**(Symbol('m'))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * (((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Symbol('m') * Symbol('n'))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Symbol('m') * Symbol('n'))) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_20():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(4))))**(Symbol('m'))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * (((Symbol('d') + (Symbol('e') * x)))**(Integer(4)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Integer(4) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(4) * Symbol('m'))) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_21():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Symbol('m'))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * (((Symbol('d') + (Symbol('e') * x)))**(Integer(3)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Integer(3) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(3) * Symbol('m'))) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_22():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Symbol('m'))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * (((Symbol('d') + (Symbol('e') * x)))**(Integer(2)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Integer(2) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(2) * Symbol('m'))) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_23():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(1))))**(Symbol('m'))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * ((((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_24():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(1))))**(Symbol('m')))**(Integer(-1))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**(Symbol('m'))) * ((((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_25():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Symbol('m')))**(Integer(-1))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(2) * Symbol('m')))) * (((((Symbol('d') + (Symbol('e') * x)))**(Integer(2)))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_26():
    f = (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((sympy.Function('Expand')(((Symbol('d') + (Symbol('e') * x)))**(Integer(3))))**(Symbol('m')))**(Integer(-1))
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * (Integer(3) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(3) * Symbol('m')))) * (((((Symbol('d') + (Symbol('e') * x)))**(Integer(3)))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_27():
    f = F**(5*x + 2)
    F = F**(5*x + 2)/(5*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_28():
    f = F**(a + b*x)
    F = F**(a + b*x)/(b*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_29():
    f = 10**(5*x + 2)
    F = 2**(5*x + 2)*5**(5*x + 1)/log(10)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_30():
    f = F**(a + b*x)*x**(sympy.S(7)/2)
    F = ((Integer(105) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(16) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(105) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.sqrt(x)) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(35) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_31():
    f = F**(a + b*x)*x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(15) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.sqrt(x)) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_32():
    f = F**(a + b*x)*x**(sympy.S(3)/2)
    F = ((Integer(3) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.sqrt(x)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_33():
    f = F**(a + b*x)*sqrt(x)
    F = (Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.sqrt(x)) * ((Symbol('b') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_34():
    f = F**(a + b*x)/sqrt(x)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_35():
    f = F**(a + b*x)/x**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x)))) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F'))))) * sympy.sqrt(sympy.log(Symbol('F'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_36():
    f = F**(a + b*x)/x**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.log(Symbol('F'))) * ((Integer(3) * sympy.sqrt(x)))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_37():
    f = F**(a + b*x)/x**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.log(Symbol('F'))) * ((Integer(15) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(15) * sympy.sqrt(x)))**(Integer(-1)))) + ((Integer(8) * (Integer(15))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_38():
    f = F**(a + b*x)/x**(sympy.S(9)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(7) * (x)**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * sympy.log(Symbol('F'))) * ((Integer(35) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(105) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * x))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(105) * sympy.sqrt(x)))**(Integer(-1)))) + ((Integer(16) * (Integer(105))**(Integer(-1))) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_39():
    f = F**(c*(a + b*x))*(d + e*x)**(sympy.S(7)/2)
    F = ((Integer(105) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(105) * (Symbol('e'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(35) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('e') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_40():
    f = F**(c*(a + b*x))*(d + e*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(15) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('e') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_41():
    f = F**(c*(a + b*x))*(d + e*x)**(sympy.S(3)/2)
    F = ((Integer(3) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('e') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_42():
    f = F**(c*(a + b*x))*sqrt(d + e*x)
    F = (Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_43():
    f = F**(c*(a + b*x))/sqrt(d + e*x)
    F = ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('e')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_44():
    f = F**(c*(a + b*x))/(d + e*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_45():
    f = F**(c*(a + b*x))/(d + e*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_46():
    f = F**(c*(a + b*x))/(d + e*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(5) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(15) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(15) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_47():
    f = F**(c*(a + b*x))/(d + e*x)**(sympy.S(9)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(7) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('c') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * sympy.log(Symbol('F'))) * ((Integer(35) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(105) * (Symbol('e'))**(Integer(3)) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(105) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt(sympy.log(Symbol('F')))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(105) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_48():
    f = x**(sympy.S(13)/2)*exp(-b*x)
    F = (Integer(-1) * ((Integer(135135) * sympy.sqrt(x)) * (((sympy.E)**((Symbol('b') * x)) * (Integer(64) * (Symbol('b'))**(Integer(7)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(45045) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((sympy.E)**((Symbol('b') * x)) * (Integer(32) * (Symbol('b'))**(Integer(6)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9009) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((sympy.E)**((Symbol('b') * x)) * (Integer(16) * (Symbol('b'))**(Integer(5)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(1287) * (x)**((Integer(7) * (Integer(2))**(Integer(-1))))) * (((sympy.E)**((Symbol('b') * x)) * (Integer(8) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(143) * (x)**((Integer(9) * (Integer(2))**(Integer(-1))))) * (((sympy.E)**((Symbol('b') * x)) * (Integer(4) * (Symbol('b'))**(Integer(3)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(13) * (x)**((Integer(11) * (Integer(2))**(Integer(-1))))) * (((sympy.E)**((Symbol('b') * x)) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((x)**((Integer(13) * (Integer(2))**(Integer(-1)))) * (((sympy.E)**((Symbol('b') * x)) * Symbol('b')))**(Integer(-1)))) + ((Integer(135135) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * ((Integer(128) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_49():
    f = F**(c*(a + b*x))*(d + e*x)**(sympy.S(4)/3)
    F = Integer(-1) * ((Symbol('e') * (Symbol('F'))**((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(7) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2)) * ((Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_50():
    f = (d + e*x)**(sympy.S(4)/3)*(F**(c*(a + b*x)))**n
    F = Integer(-1) * ((Symbol('e') * (Symbol('F'))**(((Symbol('c') * (Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))) * Symbol('n')) + (Integer(-1) * (Symbol('c') * Symbol('n') * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('F'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(7) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * Symbol('c') * Symbol('n') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2)) * ((Integer(-1) * ((Symbol('b') * Symbol('c') * Symbol('n') * (Symbol('d') + (Symbol('e') * x)) * sympy.log(Symbol('F'))) * (Symbol('e'))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_51():
    f = F**(c*(a + b*x))*(d + e*x)
    F = F**(c*(a + b*x))*(d + e*x)/(b*c*log(F)) - F**(c*(a + b*x))*e/(b**2*c**2*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_52():
    f = F**(c*(a + b*x))*(d + e*x + f*x**2)
    F = F**(c*(a + b*x))*d/(b*c*log(F)) + F**(c*(a + b*x))*e*x/(b*c*log(F)) + F**(c*(a + b*x))*f*x**2/(b*c*log(F)) - F**(c*(a + b*x))*e/(b**2*c**2*log(F)**2) - 2*F**(c*(a + b*x))*f*x/(b**2*c**2*log(F)**2) + 2*F**(c*(a + b*x))*f/(b**3*c**3*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_53():
    f = F**(c*(a + b*x))*(d + e*x + f*x**2 + g*x**3)
    F = F**(c*(a + b*x))*d/(b*c*log(F)) + F**(c*(a + b*x))*e*x/(b*c*log(F)) + F**(c*(a + b*x))*f*x**2/(b*c*log(F)) + F**(c*(a + b*x))*g*x**3/(b*c*log(F)) - F**(c*(a + b*x))*e/(b**2*c**2*log(F)**2) - 2*F**(c*(a + b*x))*f*x/(b**2*c**2*log(F)**2) - 3*F**(c*(a + b*x))*g*x**2/(b**2*c**2*log(F)**2) + 2*F**(c*(a + b*x))*f/(b**3*c**3*log(F)**3) + 6*F**(c*(a + b*x))*g*x/(b**3*c**3*log(F)**3) - 6*F**(c*(a + b*x))*g/(b**4*c**4*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_54():
    f = F**(c*(a + b*x))*(d + e*x + f*x**2 + g*x**3 + h*x**4)
    F = F**(c*(a + b*x))*d/(b*c*log(F)) + F**(c*(a + b*x))*e*x/(b*c*log(F)) + F**(c*(a + b*x))*f*x**2/(b*c*log(F)) + F**(c*(a + b*x))*g*x**3/(b*c*log(F)) + F**(c*(a + b*x))*h*x**4/(b*c*log(F)) - F**(c*(a + b*x))*e/(b**2*c**2*log(F)**2) - 2*F**(c*(a + b*x))*f*x/(b**2*c**2*log(F)**2) - 3*F**(c*(a + b*x))*g*x**2/(b**2*c**2*log(F)**2) - 4*F**(c*(a + b*x))*h*x**3/(b**2*c**2*log(F)**2) + 2*F**(c*(a + b*x))*f/(b**3*c**3*log(F)**3) + 6*F**(c*(a + b*x))*g*x/(b**3*c**3*log(F)**3) + 12*F**(c*(a + b*x))*h*x**2/(b**3*c**3*log(F)**3) - 6*F**(c*(a + b*x))*g/(b**4*c**4*log(F)**4) - 24*F**(c*(a + b*x))*h*x/(b**4*c**4*log(F)**4) + 24*F**(c*(a + b*x))*h/(b**5*c**5*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_55():
    f = x**m*(a + b*x)**3*exp(-a - b*x)
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_56():
    f = x**3*(a + b*x)**3*exp(-a - b*x)
    F = -a**3*x**3*exp(-a - b*x)/b - 3*a**3*x**2*exp(-a - b*x)/b**2 - 6*a**3*x*exp(-a - b*x)/b**3 - 6*a**3*exp(-a - b*x)/b**4 - 3*a**2*x**4*exp(-a - b*x) - 12*a**2*x**3*exp(-a - b*x)/b - 36*a**2*x**2*exp(-a - b*x)/b**2 - 72*a**2*x*exp(-a - b*x)/b**3 - 72*a**2*exp(-a - b*x)/b**4 - 3*a*b*x**5*exp(-a - b*x) - 15*a*x**4*exp(-a - b*x) - 60*a*x**3*exp(-a - b*x)/b - 180*a*x**2*exp(-a - b*x)/b**2 - 360*a*x*exp(-a - b*x)/b**3 - 360*a*exp(-a - b*x)/b**4 - b**2*x**6*exp(-a - b*x) - 6*b*x**5*exp(-a - b*x) - 30*x**4*exp(-a - b*x) - 120*x**3*exp(-a - b*x)/b - 360*x**2*exp(-a - b*x)/b**2 - 720*x*exp(-a - b*x)/b**3 - 720*exp(-a - b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_57():
    f = x**2*(a + b*x)**3*exp(-a - b*x)
    F = -a**3*x**2*exp(-a - b*x)/b - 2*a**3*x*exp(-a - b*x)/b**2 - 2*a**3*exp(-a - b*x)/b**3 - 3*a**2*x**3*exp(-a - b*x) - 9*a**2*x**2*exp(-a - b*x)/b - 18*a**2*x*exp(-a - b*x)/b**2 - 18*a**2*exp(-a - b*x)/b**3 - 3*a*b*x**4*exp(-a - b*x) - 12*a*x**3*exp(-a - b*x) - 36*a*x**2*exp(-a - b*x)/b - 72*a*x*exp(-a - b*x)/b**2 - 72*a*exp(-a - b*x)/b**3 - b**2*x**5*exp(-a - b*x) - 5*b*x**4*exp(-a - b*x) - 20*x**3*exp(-a - b*x) - 60*x**2*exp(-a - b*x)/b - 120*x*exp(-a - b*x)/b**2 - 120*exp(-a - b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_58():
    f = x*(a + b*x)**3*exp(-a - b*x)
    F = a*(a + b*x)**3*exp(-a - b*x)/b**2 + 3*a*(a + b*x)**2*exp(-a - b*x)/b**2 + 6*a*(a + b*x)*exp(-a - b*x)/b**2 + 6*a*exp(-a - b*x)/b**2 - (a + b*x)**4*exp(-a - b*x)/b**2 - 4*(a + b*x)**3*exp(-a - b*x)/b**2 - 12*(a + b*x)**2*exp(-a - b*x)/b**2 - 24*(a + b*x)*exp(-a - b*x)/b**2 - 24*exp(-a - b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_59():
    f = (a + b*x)**3*exp(-a - b*x)
    F = -(a + b*x)**3*exp(-a - b*x)/b - 3*(a + b*x)**2*exp(-a - b*x)/b - 6*(a + b*x)*exp(-a - b*x)/b - 6*exp(-a - b*x)/b
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_60():
    f = (a + b*x)**3*exp(-a - b*x)/x
    F = (Integer(-2) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) + (Integer(-1) * (Integer(3) * Symbol('a') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) + (Integer(-1) * (Integer(2) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * x)) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * x)) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (x)**(Integer(2)))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_61():
    f = (a + b*x)**3*exp(-a - b*x)/x**2
    F = ((Integer(-1) * Symbol('b')) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * x)) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * Symbol('b') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_62():
    f = (a + b*x)**3*exp(-a - b*x)/x**3
    F = ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (x)**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * x))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))) + (((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_63():
    f = (a + b*x)**3*exp(-a - b*x)/x**4
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (x)**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1)))) + (((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((Integer(6))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((sympy.E)**(Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_64():
    f = F**(a + b*(c + d*x))*x**m*(e + f*x)**2
    F = (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))**(Symbol('m')) * ((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), ((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))**(Symbol('m')) * ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2)))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * Symbol('d') * x * sympy.log(Symbol('F'))))**(Symbol('m')) * (Symbol('b') * Symbol('d') * sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_65():
    f = F**(a + b*(c + d*x))*x**3*(e + f*x)**2
    F = F**(a + b*c + b*d*x)*e**2*x**3/(b*d*log(F)) + 2*F**(a + b*c + b*d*x)*e*f*x**4/(b*d*log(F)) + F**(a + b*c + b*d*x)*f**2*x**5/(b*d*log(F)) - 3*F**(a + b*c + b*d*x)*e**2*x**2/(b**2*d**2*log(F)**2) - 8*F**(a + b*c + b*d*x)*e*f*x**3/(b**2*d**2*log(F)**2) - 5*F**(a + b*c + b*d*x)*f**2*x**4/(b**2*d**2*log(F)**2) + 6*F**(a + b*c + b*d*x)*e**2*x/(b**3*d**3*log(F)**3) + 24*F**(a + b*c + b*d*x)*e*f*x**2/(b**3*d**3*log(F)**3) + 20*F**(a + b*c + b*d*x)*f**2*x**3/(b**3*d**3*log(F)**3) - 6*F**(a + b*c + b*d*x)*e**2/(b**4*d**4*log(F)**4) - 48*F**(a + b*c + b*d*x)*e*f*x/(b**4*d**4*log(F)**4) - 60*F**(a + b*c + b*d*x)*f**2*x**2/(b**4*d**4*log(F)**4) + 48*F**(a + b*c + b*d*x)*e*f/(b**5*d**5*log(F)**5) + 120*F**(a + b*c + b*d*x)*f**2*x/(b**5*d**5*log(F)**5) - 120*F**(a + b*c + b*d*x)*f**2/(b**6*d**6*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_66():
    f = F**(a + b*(c + d*x))*x**2*(e + f*x)**2
    F = F**(a + b*c + b*d*x)*e**2*x**2/(b*d*log(F)) + 2*F**(a + b*c + b*d*x)*e*f*x**3/(b*d*log(F)) + F**(a + b*c + b*d*x)*f**2*x**4/(b*d*log(F)) - 2*F**(a + b*c + b*d*x)*e**2*x/(b**2*d**2*log(F)**2) - 6*F**(a + b*c + b*d*x)*e*f*x**2/(b**2*d**2*log(F)**2) - 4*F**(a + b*c + b*d*x)*f**2*x**3/(b**2*d**2*log(F)**2) + 2*F**(a + b*c + b*d*x)*e**2/(b**3*d**3*log(F)**3) + 12*F**(a + b*c + b*d*x)*e*f*x/(b**3*d**3*log(F)**3) + 12*F**(a + b*c + b*d*x)*f**2*x**2/(b**3*d**3*log(F)**3) - 12*F**(a + b*c + b*d*x)*e*f/(b**4*d**4*log(F)**4) - 24*F**(a + b*c + b*d*x)*f**2*x/(b**4*d**4*log(F)**4) + 24*F**(a + b*c + b*d*x)*f**2/(b**5*d**5*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_67():
    f = F**(a + b*(c + d*x))*x*(e + f*x)**2
    F = F**(a + b*c + b*d*x)*e**2*x/(b*d*log(F)) + 2*F**(a + b*c + b*d*x)*e*f*x**2/(b*d*log(F)) + F**(a + b*c + b*d*x)*f**2*x**3/(b*d*log(F)) - F**(a + b*c + b*d*x)*e**2/(b**2*d**2*log(F)**2) - 4*F**(a + b*c + b*d*x)*e*f*x/(b**2*d**2*log(F)**2) - 3*F**(a + b*c + b*d*x)*f**2*x**2/(b**2*d**2*log(F)**2) + 4*F**(a + b*c + b*d*x)*e*f/(b**3*d**3*log(F)**3) + 6*F**(a + b*c + b*d*x)*f**2*x/(b**3*d**3*log(F)**3) - 6*F**(a + b*c + b*d*x)*f**2/(b**4*d**4*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_68():
    f = F**(a + b*(c + d*x))*(e + f*x)**2
    F = F**(a + b*c + b*d*x)*(e + f*x)**2/(b*d*log(F)) - 2*F**(a + b*c + b*d*x)*f*(e + f*x)/(b**2*d**2*log(F)**2) + 2*F**(a + b*c + b*d*x)*f**2/(b**3*d**3*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_69():
    f = F**(a + b*(c + d*x))*(e + f*x)**2/x
    F = ((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F'))))) + (Integer(-1) * (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * x) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_70():
    f = F**(a + b*(c + d*x))*(e + f*x)**2/x**2
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * (x)**(Integer(-1)))) + (Integer(2) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F'))))) + (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1))) + (Symbol('b') * Symbol('d') * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * sympy.log(Symbol('F')))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_71():
    f = F**(a + b*(c + d*x))*(e + f*x)**2/x**3
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * (x)**(Integer(-1)))) + ((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F'))))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(2) * Symbol('b') * Symbol('d') * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * sympy.log(Symbol('F'))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_72():
    f = F**(a + b*(c + d*x))*(e + f*x)**2/x**4
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((x)**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * (x)**(Integer(-1)))) + (Symbol('b') * Symbol('d') * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * sympy.log(Symbol('F'))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * x))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2))) + ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_73():
    f = F**(a + b*(c + d*x))*(e + f*x)**2/x**5
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * sympy.log(Symbol('F'))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(24) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')) + (Symbol('b') * Symbol('d') * x))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(24) * x))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3)) * Symbol('e') * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(3))) + ((Integer(24))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * (Symbol('d'))**(Integer(4)) * (Symbol('e'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * Symbol('c')))) * sympy.Function('ExpIntegralEi')((Symbol('b') * Symbol('d') * x * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_74():
    f = (a + b*x)**4*(c + d*x)**3*exp(-a - b*x)
    F = -d**3*(a + b*x)**7*exp(-a - b*x)/b**4 - 7*d**3*(a + b*x)**6*exp(-a - b*x)/b**4 - 42*d**3*(a + b*x)**5*exp(-a - b*x)/b**4 - 210*d**3*(a + b*x)**4*exp(-a - b*x)/b**4 - 840*d**3*(a + b*x)**3*exp(-a - b*x)/b**4 - 2520*d**3*(a + b*x)**2*exp(-a - b*x)/b**4 - 5040*d**3*(a + b*x)*exp(-a - b*x)/b**4 - 5040*d**3*exp(-a - b*x)/b**4 - 3*d**2*(a + b*x)**6*(-a*d + b*c)*exp(-a - b*x)/b**4 - 18*d**2*(a + b*x)**5*(-a*d + b*c)*exp(-a - b*x)/b**4 - 90*d**2*(a + b*x)**4*(-a*d + b*c)*exp(-a - b*x)/b**4 - 360*d**2*(a + b*x)**3*(-a*d + b*c)*exp(-a - b*x)/b**4 - 1080*d**2*(a + b*x)**2*(-a*d + b*c)*exp(-a - b*x)/b**4 - 2160*d**2*(a + b*x)*(-a*d + b*c)*exp(-a - b*x)/b**4 - 2160*d**2*(-a*d + b*c)*exp(-a - b*x)/b**4 - 3*d*(a + b*x)**5*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - 15*d*(a + b*x)**4*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - 60*d*(a + b*x)**3*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - 180*d*(a + b*x)**2*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - 360*d*(a + b*x)*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - 360*d*(-a*d + b*c)**2*exp(-a - b*x)/b**4 - (a + b*x)**4*(-a*d + b*c)**3*exp(-a - b*x)/b**4 - 4*(a + b*x)**3*(-a*d + b*c)**3*exp(-a - b*x)/b**4 - 12*(a + b*x)**2*(-a*d + b*c)**3*exp(-a - b*x)/b**4 - 24*(a + b*x)*(-a*d + b*c)**3*exp(-a - b*x)/b**4 - 24*(-a*d + b*c)**3*exp(-a - b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_75():
    f = (a + b*x)**4*(c + d*x)**2*exp(-a - b*x)
    F = -d**2*(a + b*x)**6*exp(-a - b*x)/b**3 - 6*d**2*(a + b*x)**5*exp(-a - b*x)/b**3 - 30*d**2*(a + b*x)**4*exp(-a - b*x)/b**3 - 120*d**2*(a + b*x)**3*exp(-a - b*x)/b**3 - 360*d**2*(a + b*x)**2*exp(-a - b*x)/b**3 - 720*d**2*(a + b*x)*exp(-a - b*x)/b**3 - 720*d**2*exp(-a - b*x)/b**3 - 2*d*(a + b*x)**5*(-a*d + b*c)*exp(-a - b*x)/b**3 - 10*d*(a + b*x)**4*(-a*d + b*c)*exp(-a - b*x)/b**3 - 40*d*(a + b*x)**3*(-a*d + b*c)*exp(-a - b*x)/b**3 - 120*d*(a + b*x)**2*(-a*d + b*c)*exp(-a - b*x)/b**3 - 240*d*(a + b*x)*(-a*d + b*c)*exp(-a - b*x)/b**3 - 240*d*(-a*d + b*c)*exp(-a - b*x)/b**3 - (a + b*x)**4*(-a*d + b*c)**2*exp(-a - b*x)/b**3 - 4*(a + b*x)**3*(-a*d + b*c)**2*exp(-a - b*x)/b**3 - 12*(a + b*x)**2*(-a*d + b*c)**2*exp(-a - b*x)/b**3 - 24*(a + b*x)*(-a*d + b*c)**2*exp(-a - b*x)/b**3 - 24*(-a*d + b*c)**2*exp(-a - b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_76():
    f = (a + b*x)**4*(c + d*x)*exp(-a - b*x)
    F = -d*(a + b*x)**5*exp(-a - b*x)/b**2 - 5*d*(a + b*x)**4*exp(-a - b*x)/b**2 - 20*d*(a + b*x)**3*exp(-a - b*x)/b**2 - 60*d*(a + b*x)**2*exp(-a - b*x)/b**2 - 120*d*(a + b*x)*exp(-a - b*x)/b**2 - 120*d*exp(-a - b*x)/b**2 - (a + b*x)**4*(-a*d + b*c)*exp(-a - b*x)/b**2 - (a + b*x)**3*(-4*a*d + 4*b*c)*exp(-a - b*x)/b**2 - (a + b*x)**2*(-12*a*d + 12*b*c)*exp(-a - b*x)/b**2 - (a + b*x)*(-24*a*d + 24*b*c)*exp(-a - b*x)/b**2 - (-24*a*d + 24*b*c)*exp(-a - b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_77():
    f = (a + b*x)**4*exp(-a - b*x)
    F = -(a + b*x)**4*exp(-a - b*x)/b - 4*(a + b*x)**3*exp(-a - b*x)/b - 12*(a + b*x)**2*exp(-a - b*x)/b - 24*(a + b*x)*exp(-a - b*x)/b - 24*exp(-a - b*x)/b
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_78():
    f = (a + b*x)**4*exp(-a - b*x)/(c + d*x)
    F = (Integer(-1) * ((Integer(6) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('a') + (Symbol('b') * x))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('d'))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * (Symbol('d'))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_79():
    f = (a + b*x)**4*exp(-a - b*x)/(c + d*x)**2
    F = (Integer(-1) * ((Integer(2) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_80():
    f = (a + b*x)**4*exp(-a - b*x)/(c + d*x)**3
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(6)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(7))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_81():
    f = (a + b*x)**4*exp(-a - b*x)/(c + d*x)**4
    F = (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(6)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(6)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(7)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(7)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(8))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_82():
    f = (a + b*x)**4*exp(-a - b*x)/(c + d*x)**5
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(4) * (Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(12) * (Symbol('d'))**(Integer(6)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('d'))**(Integer(6)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(24) * (Symbol('d'))**(Integer(7)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(5)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('d'))**(Integer(6)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('d'))**(Integer(7)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(24) * (Symbol('d'))**(Integer(8)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(4)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(6)))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(4)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(7)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(4)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(8))))**(Integer(-1))) + (((Symbol('b'))**(Integer(4)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(24) * (Symbol('d'))**(Integer(9))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_83():
    f = F**(c*(a + b*x))*x**m*(e*n + e*(b*c*x*log(F) + m + 1)*log(d*x) + e)*log(d*x)**n
    F = F**(c*(a + b*x))*e*x**(m + 1)*log(d*x)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_84():
    f = F**(c*(a + b*x))*x**2*(e*n + e*(b*c*x*log(F) + 3)*log(d*x) + e)*log(d*x)**n
    F = F**(c*(a + b*x))*e*x**3*log(d*x)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_85():
    f = F**(c*(a + b*x))*x*(e*n + e*(b*c*x*log(F) + 2)*log(d*x) + e)*log(d*x)**n
    F = F**(c*(a + b*x))*e*x**2*log(d*x)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_86():
    f = F**(c*(a + b*x))*(e*n + e*(b*c*x*log(F) + 1)*log(d*x) + e)*log(d*x)**n
    F = F**(c*(a + b*x))*e*x*log(d*x)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_87():
    f = F**(c*(a + b*x))*(b*c*e*x*log(F)*log(d*x) + e*n + e)*log(d*x)**n/x
    F = F**(c*(a + b*x))*e*log(d*x)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_88():
    f = F**(c*(a + b*x))*(e*n + e*(b*c*x*log(F) - 1)*log(d*x) + e)*log(d*x)**n/x**2
    F = F**(c*(a + b*x))*e*log(d*x)**(n + 1)/x
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_89():
    f = F**(c*(a + b*x))*(e*n + e*(b*c*x*log(F) - 2)*log(d*x) + e)*log(d*x)**n/x**3
    F = F**(c*(a + b*x))*e*log(d*x)**(n + 1)/x**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_90():
    f = x**4*sqrt(exp(a + b*x))
    F = 2*x**4*sqrt(exp(a + b*x))/b - 16*x**3*sqrt(exp(a + b*x))/b**2 + 96*x**2*sqrt(exp(a + b*x))/b**3 - 384*x*sqrt(exp(a + b*x))/b**4 + 768*sqrt(exp(a + b*x))/b**5
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_91():
    f = x**3*sqrt(exp(a + b*x))
    F = 2*x**3*sqrt(exp(a + b*x))/b - 12*x**2*sqrt(exp(a + b*x))/b**2 + 48*x*sqrt(exp(a + b*x))/b**3 - 96*sqrt(exp(a + b*x))/b**4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_92():
    f = x**2*sqrt(exp(a + b*x))
    F = 2*x**2*sqrt(exp(a + b*x))/b - 8*x*sqrt(exp(a + b*x))/b**2 + 16*sqrt(exp(a + b*x))/b**3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_93():
    f = x*sqrt(exp(a + b*x))
    F = 2*x*sqrt(exp(a + b*x))/b - 4*sqrt(exp(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_94():
    f = sqrt(exp(a + b*x))
    F = 2*sqrt(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_95():
    f = sqrt(exp(a + b*x))/x
    F = (sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * x) * (Integer(2))**(Integer(-1))))) * ((sympy.E)**(((Symbol('b') * x) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_96():
    f = sqrt(exp(a + b*x))/x**2
    F = (Integer(-1) * (sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1)))) + (((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * x) * (Integer(2))**(Integer(-1))))) * ((sympy.E)**(((Symbol('b') * x) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_97():
    f = sqrt(exp(a + b*x))/x**3
    F = (Integer(-1) * (sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * x))**(Integer(-1)))) + (((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * x) * (Integer(2))**(Integer(-1))))) * ((sympy.E)**(((Symbol('b') * x) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_1_u_F_pow_c_a_plus_b_x_pow_n_98():
    f = sqrt(exp(a + b*x))/x**4
    F = (Integer(-1) * (sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(24) * x))**(Integer(-1)))) + (((Integer(48))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * x) * (Integer(2))**(Integer(-1))))) * ((sympy.E)**(((Symbol('b') * x) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F

