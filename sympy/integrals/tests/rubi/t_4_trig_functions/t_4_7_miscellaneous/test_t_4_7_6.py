"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.7 Miscellaneous/4.7.6 f^(a+b x+c x^2) trig(d+e x+f x^2)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

F, a, b, c, d, e, f, m, n = symbols('F a b c d e f m n')

def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_1():
    f = F**(c*(a + b*x))*sin(d + e*x)**n
    F = -F**(c*(a + b*x))*sin(d + e*x)**n*hyper((-n, -(I*b*c*log(F) + e*n)/(2*e)), (-I*b*c*log(F)/(2*e) - n/2 + 1,), exp(2*I*(d + e*x)))/((1 - exp(2*I*(d + e*x)))**n*(-b*c*log(F) + I*e*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_2():
    f = F**(c*(a + b*x))*sin(d + e*x)**3
    F = 6*F**(c*(a + b*x))*b*c*e**2*log(F)*sin(d + e*x)/(b**4*c**4*log(F)**4 + 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) + F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**3/(b**2*c**2*log(F)**2 + 9*e**2) - 6*F**(c*(a + b*x))*e**3*cos(d + e*x)/(b**4*c**4*log(F)**4 + 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) - 3*F**(c*(a + b*x))*e*sin(d + e*x)**2*cos(d + e*x)/(b**2*c**2*log(F)**2 + 9*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_3():
    f = F**(c*(a + b*x))*sin(d + e*x)**2
    F = F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**2/(b**2*c**2*log(F)**2 + 4*e**2) - 2*F**(c*(a + b*x))*e*sin(d + e*x)*cos(d + e*x)/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e**2/(b*c*(b**2*c**2*log(F)**2 + 4*e**2)*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_4():
    f = F**(c*(a + b*x))*sin(d + e*x)
    F = F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2) - F**(c*(a + b*x))*e*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_5():
    f = F**(c*(a + b*x))*csc(d + e*x)
    F = -2*F**(c*(a + b*x))*exp(I*(d + e*x))*hyper((1, (-I*b*c*log(F) + e)/(2*e)), (-I*b*c*log(F)/(2*e) + sympy.S(3)/2,), exp(2*I*(d + e*x)))/(-I*b*c*log(F) + e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_6():
    f = F**(c*(a + b*x))*csc(d + e*x)**2
    F = -4*F**(c*(a + b*x))*exp(2*I*(d + e*x))*hyper((2, -I*b*c*log(F)/(2*e) + 1), (-I*b*c*log(F)/(2*e) + 2,), exp(2*I*(d + e*x)))/(b*c*log(F) + 2*I*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_7():
    f = F**(c*(a + b*x))*csc(d + e*x)**3
    F = -F**(c*(a + b*x))*b*c*log(F)*csc(d + e*x)/(2*e**2) - F**(c*(a + b*x))*cot(d + e*x)*csc(d + e*x)/(2*e) - F**(c*(a + b*x))*(I*b*c*log(F) + e)*exp(I*(d + e*x))*hyper((1, (-I*b*c*log(F) + e)/(2*e)), (-I*b*c*log(F)/(2*e) + sympy.S(3)/2,), exp(2*I*(d + e*x)))/e**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_8():
    f = F**(c*(a + b*x))*csc(d + e*x)**4
    F = -F**(c*(a + b*x))*b*c*log(F)*csc(d + e*x)**2/(6*e**2) - F**(c*(a + b*x))*cot(d + e*x)*csc(d + e*x)**2/(3*e) + 2*F**(c*(a + b*x))*(-b*c*log(F) + 2*I*e)*exp(2*I*(d + e*x))*hyper((2, -I*b*c*log(F)/(2*e) + 1), (-I*b*c*log(F)/(2*e) + 2,), exp(2*I*(d + e*x)))/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_9():
    f = exp(x)*sin(x)**4
    F = exp(x)*sin(x)**4/17 - 4*exp(x)*sin(x)**3*cos(x)/17 + 12*exp(x)*sin(x)**2/85 - 24*exp(x)*sin(x)*cos(x)/85 + 24*exp(x)/85
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_10():
    f = F**(c*(a + b*x))*cos(d + e*x)**n
    F = -F**(c*(a + b*x))*cos(d + e*x)**n*hyper((-n, -(I*b*c*log(F) + e*n)/(2*e)), (-I*b*c*log(F)/(2*e) - n/2 + 1,), -exp(2*I*(d + e*x)))/((-b*c*log(F) + I*e*n)*(exp(2*I*(d + e*x)) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_11():
    f = F**(c*(a + b*x))*cos(d + e*x)**3
    F = 6*F**(c*(a + b*x))*b*c*e**2*log(F)*cos(d + e*x)/(b**4*c**4*log(F)**4 + 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) + F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**3/(b**2*c**2*log(F)**2 + 9*e**2) + 6*F**(c*(a + b*x))*e**3*sin(d + e*x)/(b**4*c**4*log(F)**4 + 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) + 3*F**(c*(a + b*x))*e*sin(d + e*x)*cos(d + e*x)**2/(b**2*c**2*log(F)**2 + 9*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_12():
    f = F**(c*(a + b*x))*cos(d + e*x)**2
    F = F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**2/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e*sin(d + e*x)*cos(d + e*x)/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e**2/(b*c*(b**2*c**2*log(F)**2 + 4*e**2)*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_13():
    f = F**(c*(a + b*x))*cos(d + e*x)
    F = F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + F**(c*(a + b*x))*e*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_14():
    f = F**(c*(a + b*x))*sec(d + e*x)
    F = 2*F**(c*(a + b*x))*exp(I*(d + e*x))*hyper((1, (-I*b*c*log(F) + e)/(2*e)), (-I*b*c*log(F)/(2*e) + sympy.S(3)/2,), -exp(2*I*(d + e*x)))/(b*c*log(F) + I*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_15():
    f = F**(c*(a + b*x))*sec(d + e*x)**2
    F = 4*F**(c*(a + b*x))*exp(2*I*(d + e*x))*hyper((2, -I*b*c*log(F)/(2*e) + 1), (-I*b*c*log(F)/(2*e) + 2,), -exp(2*I*(d + e*x)))/(b*c*log(F) + 2*I*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_16():
    f = F**(c*(a + b*x))*sec(d + e*x)**3
    F = -F**(c*(a + b*x))*b*c*log(F)*sec(d + e*x)/(2*e**2) + F**(c*(a + b*x))*tan(d + e*x)*sec(d + e*x)/(2*e) - F**(c*(a + b*x))*(-b*c*log(F) + I*e)*exp(I*(d + e*x))*hyper((1, (-I*b*c*log(F) + e)/(2*e)), (-I*b*c*log(F)/(2*e) + sympy.S(3)/2,), -exp(2*I*(d + e*x)))/e**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_17():
    f = F**(c*(a + b*x))*sec(d + e*x)**4
    F = -F**(c*(a + b*x))*b*c*log(F)*sec(d + e*x)**2/(6*e**2) + F**(c*(a + b*x))*tan(d + e*x)*sec(d + e*x)**2/(3*e) - 2*F**(c*(a + b*x))*(-b*c*log(F) + 2*I*e)*exp(2*I*(d + e*x))*hyper((2, -I*b*c*log(F)/(2*e) + 1), (-I*b*c*log(F)/(2*e) + 2,), -exp(2*I*(d + e*x)))/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_18():
    f = exp(x)*cos(x)**4
    F = 4*exp(x)*sin(x)*cos(x)**3/17 + 24*exp(x)*sin(x)*cos(x)/85 + exp(x)*cos(x)**4/17 + 12*exp(x)*cos(x)**2/85 + 24*exp(x)/85
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_19():
    f = exp(c*(a + b*x))*tan(d + e*x)**3
    F = -6*I*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) + 12*I*exp(c*(a + b*x))*hyper((2, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) - 8*I*exp(c*(a + b*x))*hyper((3, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) + I*exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_20():
    f = exp(c*(a + b*x))*tan(d + e*x)**2
    F = 4*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) - 4*exp(c*(a + b*x))*hyper((2, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) - exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_21():
    f = exp(c*(a + b*x))*tan(d + e*x)
    F = 2*I*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), -exp(2*I*(d + e*x)))/(b*c) - I*exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_22():
    f = exp(c*(a + b*x))*cot(d + e*x)
    F = -2*I*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) + I*exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_23():
    f = exp(c*(a + b*x))*cot(d + e*x)**2
    F = 4*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) - 4*exp(c*(a + b*x))*hyper((2, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) - exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_24():
    f = exp(c*(a + b*x))*cot(d + e*x)**3
    F = 6*I*exp(c*(a + b*x))*hyper((1, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) - 12*I*exp(c*(a + b*x))*hyper((2, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) + 8*I*exp(c*(a + b*x))*hyper((3, -I*b*c/(2*e)), (-I*b*c/(2*e) + 1,), exp(2*I*(d + e*x)))/(b*c) - I*exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_25():
    f = F**(a + b*x)*cot(c/2 + d*x/2 + pi/4)
    F = -2*I*F**(a + b*x)*hyper((1, -I*b*log(F)/d), (-I*b*log(F)/d + 1,), I*exp(I*(c + d*x)))/(b*log(F)) + I*F**(a + b*x)/(b*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_26():
    f = F**(c*(a + b*x))*sec(d + e*x)**n
    F = F**(a*c + b*c*x)*(exp(2*I*(d + e*x)) + 1)**n*hyper((n, (-I*b*c*log(F) + e*n)/(2*e)), (-I*b*c*log(F)/(2*e) + n/2 + 1,), -exp(2*I*(d + e*x)))*sec(d + e*x)**n/(b*c*log(F) + I*e*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_27():
    f = F**(c*(a + b*x))*csc(d + e*x)**n
    F = -F**(a*c + b*c*x)*(1 - exp(-2*I*(d + e*x)))**n*csc(d + e*x)**n*hyper((n, (I*b*c*log(F) + e*n)/(2*e)), (I*b*c*log(F)/(2*e) + n/2 + 1,), exp(-2*I*(d + e*x)))/(-b*c*log(F) + I*e*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_28():
    f = F**(c*(a + b*x))*(f*x)**m*sin(d + e*x)
    F = (Integer(-1) * (((Symbol('F'))**((Symbol('a') * Symbol('c'))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (x * ((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))))) * (((sympy.E)**((sympy.I * Symbol('d'))) * ((x * ((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))))**(Symbol('m')) * (Integer(2) * (Symbol('e') + (sympy.I * Symbol('b') * Symbol('c') * sympy.log(Symbol('F')))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('d'))) * (Symbol('F'))**((Symbol('a') * Symbol('c'))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * x) * ((sympy.I * Symbol('e')) + (Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))))) * (((((Integer(-1) * x) * ((sympy.I * Symbol('e')) + (Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))))**(Symbol('m')) * (Integer(2) * (Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('c') * sympy.log(Symbol('F'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_29():
    f = F**(c*(a + b*x))*(f*x)**m/sin(d + e*x)
    F = sympy.Function('CannotIntegrate')(((Symbol('F'))**(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))) * ((Symbol('f') * x))**(Symbol('m')) * sympy.csc((Symbol('d') + (Symbol('e') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_30():
    f = F**(c*(a + b*x))*(f*x)**m/sin(d + e*x)**2
    F = sympy.Function('CannotIntegrate')(((Symbol('F'))**(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))) * ((Symbol('f') * x))**(Symbol('m')) * (sympy.csc((Symbol('d') + (Symbol('e') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_31():
    f = F**(c*(a + b*x))*f*(f*x)**(m - 2)*(e*x*cos(d + e*x) + (b*c*x*log(F) + m - 1)*sin(d + e*x))
    F = F**(a*c + b*c*x)*(f*x)**(m - 1)*sin(d + e*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_32():
    f = F**(c*(a + b*x))*f*(f*x)**m*(e*x*cos(d + e*x) + (b*c*x*log(F) + m + 1)*sin(d + e*x))
    F = F**(c*(a + b*x))*f*x*(f*x)**m*sin(d + e*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_33():
    f = F**(c*(a + b*x))*(f*x)**m*(e*x*cos(d + e*x) + (b*c*x*log(F) + m)*sin(d + e*x))/x
    F = F**(a*c + b*c*x)*(f*x)**m*sin(d + e*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_34():
    f = F**(c*(a + b*x))*(b*c*log(F)*sin(d + e*x) + e*cos(d + e*x))
    F = F**(c*(a + b*x))*sin(d + e*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_35():
    f = F**(c*(a + b*x))*(e*x*cos(d + e*x) + (b*c*x*log(F) - 1)*sin(d + e*x))/x**2
    F = F**(a*c + b*c*x)*sin(d + e*x)/x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_36():
    f = F**(c*(a + b*x))*(e*x*cos(d + e*x) + (b*c*x*log(F) - 2)*sin(d + e*x))/x**3
    F = F**(a*c + b*c*x)*sin(d + e*x)/x**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_37():
    f = exp(a + b*x)*sin(c + d*x)*cos(c + d*x)
    F = b*exp(a + b*x)*sin(2*c + 2*d*x)/(2*b**2 + 8*d**2) - d*exp(a + b*x)*cos(2*c + 2*d*x)/(b**2 + 4*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_38():
    f = exp(a + b*x)*sin(c + d*x)**2*cos(c + d*x)
    F = -b*exp(a + b*x)*cos(3*c + 3*d*x)/(4*b**2 + 36*d**2) + b*exp(a + b*x)*cos(c + d*x)/(4*b**2 + 4*d**2) - 3*d*exp(a + b*x)*sin(3*c + 3*d*x)/(4*b**2 + 36*d**2) + d*exp(a + b*x)*sin(c + d*x)/(4*b**2 + 4*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_39():
    f = exp(a + b*x)*sin(c + d*x)**3*cos(c + d*x)
    F = -b*exp(a + b*x)*sin(4*c + 4*d*x)/(8*b**2 + 128*d**2) + b*exp(a + b*x)*sin(2*c + 2*d*x)/(4*b**2 + 16*d**2) + d*exp(a + b*x)*cos(4*c + 4*d*x)/(2*b**2 + 32*d**2) - d*exp(a + b*x)*cos(2*c + 2*d*x)/(2*b**2 + 8*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_40():
    f = exp(a + b*x)*sin(c + d*x)*cos(c + d*x)**2
    F = b*exp(a + b*x)*sin(3*c + 3*d*x)/(4*b**2 + 36*d**2) + b*exp(a + b*x)*sin(c + d*x)/(4*b**2 + 4*d**2) - 3*d*exp(a + b*x)*cos(3*c + 3*d*x)/(4*b**2 + 36*d**2) - d*exp(a + b*x)*cos(c + d*x)/(4*b**2 + 4*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_41():
    f = exp(a + b*x)*sin(c + d*x)**2*cos(c + d*x)**2
    F = -b*exp(a + b*x)*cos(4*c + 4*d*x)/(8*b**2 + 128*d**2) - d*exp(a + b*x)*sin(4*c + 4*d*x)/(2*b**2 + 32*d**2) + exp(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_42():
    f = exp(a + b*x)*sin(c + d*x)**3*cos(c + d*x)**2
    F = -b*exp(a + b*x)*sin(5*c + 5*d*x)/(16*b**2 + 400*d**2) + b*exp(a + b*x)*sin(3*c + 3*d*x)/(16*b**2 + 144*d**2) + b*exp(a + b*x)*sin(c + d*x)/(8*b**2 + 8*d**2) + 5*d*exp(a + b*x)*cos(5*c + 5*d*x)/(16*b**2 + 400*d**2) - 3*d*exp(a + b*x)*cos(3*c + 3*d*x)/(16*b**2 + 144*d**2) - d*exp(a + b*x)*cos(c + d*x)/(8*b**2 + 8*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_43():
    f = exp(a + b*x)*sin(c + d*x)*cos(c + d*x)**3
    F = b*exp(a + b*x)*sin(4*c + 4*d*x)/(8*b**2 + 128*d**2) + b*exp(a + b*x)*sin(2*c + 2*d*x)/(4*b**2 + 16*d**2) - d*exp(a + b*x)*cos(4*c + 4*d*x)/(2*b**2 + 32*d**2) - d*exp(a + b*x)*cos(2*c + 2*d*x)/(2*b**2 + 8*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_44():
    f = exp(a + b*x)*sin(c + d*x)**2*cos(c + d*x)**3
    F = -b*exp(a + b*x)*cos(5*c + 5*d*x)/(16*b**2 + 400*d**2) - b*exp(a + b*x)*cos(3*c + 3*d*x)/(16*b**2 + 144*d**2) + b*exp(a + b*x)*cos(c + d*x)/(8*b**2 + 8*d**2) - 5*d*exp(a + b*x)*sin(5*c + 5*d*x)/(16*b**2 + 400*d**2) - 3*d*exp(a + b*x)*sin(3*c + 3*d*x)/(16*b**2 + 144*d**2) + d*exp(a + b*x)*sin(c + d*x)/(8*b**2 + 8*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_45():
    f = exp(a + b*x)*sin(c + d*x)**3*cos(c + d*x)**3
    F = -b*exp(a + b*x)*sin(6*c + 6*d*x)/(32*b**2 + 1152*d**2) + 3*b*exp(a + b*x)*sin(2*c + 2*d*x)/(32*b**2 + 128*d**2) + 3*d*exp(a + b*x)*cos(6*c + 6*d*x)/(16*b**2 + 576*d**2) - 3*d*exp(a + b*x)*cos(2*c + 2*d*x)/(16*b**2 + 64*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_46():
    f = x*exp(x)*sin(x)
    F = x*exp(x)*sin(x)/2 - x*exp(x)*cos(x)/2 + exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_47():
    f = x**2*exp(x)*sin(x)
    F = x**2*exp(x)*sin(x)/2 - x**2*exp(x)*cos(x)/2 + x*exp(x)*cos(x) - exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_48():
    f = x*exp(x)*cos(x)
    F = x*exp(x)*sin(x)/2 + x*exp(x)*cos(x)/2 - exp(x)*sin(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_49():
    f = x**2*exp(x)*cos(x)
    F = x**2*exp(x)*sin(x)/2 + x**2*exp(x)*cos(x)/2 - x*exp(x)*sin(x) + exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_50():
    f = (2*sin(4*x) - 5*cos(4*x))*exp(3*x)
    F = -14*exp(3*x)*sin(4*x)/25 - 23*exp(3*x)*cos(4*x)/25
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_51():
    f = exp(x)*sin(x) + exp(-x)*sin(x)
    F = exp(x)*sin(x)/2 - exp(x)*cos(x)/2 - exp(-x)*sin(x)/2 - exp(-x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_52():
    f = F**(a + b*x)*cos(c + d*x)/(e*sin(c + d*x) + e)
    F = -2*I*F**(a + b*x)*hyper((1, -I*b*log(F)/d), (-I*b*log(F)/d + 1,), I*exp(I*(c + d*x)))/(b*e*log(F)) + I*F**(a + b*x)/(b*e*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_53():
    f = F**(a + b*x)*cos(c + d*x)/(-e*sin(c + d*x) + e)
    F = 2*I*F**(a + b*x)*hyper((1, -I*b*log(F)/d), (-I*b*log(F)/d + 1,), -I*exp(I*(c + d*x)))/(b*e*log(F)) - I*F**(a + b*x)/(b*e*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_54():
    f = F**(a + b*x)*sin(c + d*x)/(e*cos(c + d*x) + e)
    F = 2*I*F**(a + b*x)*hyper((1, -I*b*log(F)/d), (-I*b*log(F)/d + 1,), -exp(I*(c + d*x)))/(b*e*log(F)) - I*F**(a + b*x)/(b*e*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_55():
    f = F**(a + b*x)*sin(c + d*x)/(-e*cos(c + d*x) + e)
    F = -2*I*F**(a + b*x)*hyper((1, -I*b*log(F)/d), (-I*b*log(F)/d + 1,), exp(I*(c + d*x)))/(b*e*log(F)) + I*F**(a + b*x)/(b*e*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_56():
    f = exp(x**2)*sin(b*x)
    F = ((Integer(4))**(Integer(-1)) * sympy.I * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((sympy.I * Symbol('b')) + (Integer(2) * x))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_57():
    f = exp(x**2)*cos(b*x)
    F = ((Integer(4))**(Integer(-1)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x))))) + ((Integer(4))**(Integer(-1)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((sympy.I * Symbol('b')) + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_58():
    f = exp(x**2)*sin(a + b*x)
    F = ((Integer(4))**(Integer(-1)) * sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((sympy.I * Symbol('b')) + (Integer(2) * x))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_59():
    f = exp(x**2)*cos(a + b*x)
    F = ((Integer(4))**(Integer(-1)) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x))))) + ((Integer(4))**(Integer(-1)) * (sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((sympy.I * Symbol('b')) + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_60():
    f = x*exp(2*x**2)*cos(2*x**2)
    F = exp(2*x**2)*sin(2*x**2)/8 + exp(2*x**2)*cos(2*x**2)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_61():
    f = exp(x)*sin(exp(x))
    F = -cos(exp(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_62():
    f = exp(x)*csc(exp(x))*sec(exp(x))
    F = log(tan(exp(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_63():
    f = exp(x)*cos(exp(x))
    F = sin(exp(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_64():
    f = exp(2*x)*cos(exp(2*x))
    F = sin(exp(2*x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_65():
    f = exp(-2*x)*cos(exp(-2*x))
    F = -sin(exp(-2*x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_66():
    f = exp(2*x)*cos(exp(x))
    F = exp(x)*sin(exp(x)) + cos(exp(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_67():
    f = exp(3*x - 1)*sin(exp(3*x - 1) + 1)*cos(exp(3*x - 1))
    F = exp(3*x - 1)*sin(1)/6 - cos(2*exp(3*x - 1) + 1)/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_68():
    f = exp(x)*tan(exp(x))
    F = -log(cos(exp(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_69():
    f = exp(x)*sec(exp(x))
    F = atanh(sin(exp(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_70():
    f = exp(x)*tan(exp(x))*sec(exp(x))
    F = sec(exp(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_71():
    f = exp(x)*csc(exp(x))**2
    F = -cot(exp(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_72():
    f = exp(x)*sin(a + b*x)
    F = -b*exp(x)*cos(a + b*x)/(b**2 + 1) + exp(x)*sin(a + b*x)/(b**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_73():
    f = exp(x)*sin(a + c*x**2)
    F = (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('a')) + (Symbol('c'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + ((Integer(2) * sympy.I) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * Symbol('c') * x)))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('a')) + (Symbol('c'))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_74():
    f = exp(x)*sin(a + b*x + c*x**2)
    F = (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((Integer(4))**(Integer(-1)) * sympy.I * ((Integer(4) * Symbol('a')) + (((Integer(1) + (sympy.I * Symbol('b'))))**(Integer(2)) * (Symbol('c'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.I * Symbol('b')) + (Integer(2) * sympy.I * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((sympy.I * ((sympy.I + Symbol('b')))**(Integer(2))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('b'))) + (Integer(-1) * (Integer(2) * sympy.I * Symbol('c') * x)))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_75():
    f = exp(x**2)*sin(a + b*x)
    F = ((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x)) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(2) * x)) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_76():
    f = exp(x**2)*sin(a + c*x**2)
    F = (((sympy.I * (Integer(4))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))) * x))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))) * (sympy.E)**((sympy.I * Symbol('a')))))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**((sympy.I * Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (sympy.I * Symbol('c')))) * x))) * (sympy.sqrt((Integer(1) + (sympy.I * Symbol('c')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_77():
    f = exp(x**2)*sin(a + b*x + c*x**2)
    F = (Integer(-1) * ((sympy.I * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(-1) * (Integer(2) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * x))) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))))))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (((Integer(4) * sympy.I) + (Integer(4) * Symbol('c'))))**(Integer(-1))))))) * (Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + (sympy.I * Symbol('c')))))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(2) * (Integer(1) + (sympy.I * Symbol('c'))) * x)) * ((Integer(2) * sympy.sqrt((Integer(1) + (sympy.I * Symbol('c'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + (sympy.I * Symbol('c'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_78():
    f = f**(a + b*x)*sin(d + f*x**2)
    F = (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_79():
    f = f**(a + b*x)*sin(d + f*x**2)**2
    F = (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(4) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) + ((((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(4) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) * ((sympy.E)**(((sympy.I * (Integer(8))**(Integer(-1))) * ((Integer(16) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_80():
    f = f**(a + b*x)*sin(d + f*x**2)**3
    F = ((Integer(3) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((Integer(16))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erf')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(6) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((Integer(16) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(16))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erfi')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(6) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((sympy.E)**(((sympy.I * (Integer(12))**(Integer(-1))) * ((Integer(36) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_81():
    f = f**(a + b*x)*sin(d + e*x + f*x**2)
    F = (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + ((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(4))**(Integer(-1))) * ((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_82():
    f = f**(a + b*x)*sin(d + e*x + f*x**2)**2
    F = (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * ((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(2) * sympy.I) * Symbol('e')) + ((Integer(4) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) + (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * (((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(2) * sympy.I) * Symbol('e')) + ((Integer(4) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_83():
    f = f**(a + b*x)*sin(d + e*x + f*x**2)**3
    F = ((Integer(3) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + ((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((Integer(16))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * ((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erf')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(3) * sympy.I) * Symbol('e')) + ((Integer(6) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(4))**(Integer(-1))) * ((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * (((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erfi')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(3) * sympy.I) * Symbol('e')) + ((Integer(6) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_84():
    f = f**(a + c*x**2)*sin(d + e*x)
    F = ((((Integer(-1) * sympy.I) * (Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_85():
    f = f**(a + c*x**2)*sin(d + e*x)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * x * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('c') * x * sympy.log(Symbol('f')))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_86():
    f = f**(a + c*x**2)*sin(d + e*x)**3
    F = ((((Integer(-3) * sympy.I) * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.I * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * sympy.I) * (Integer(16))**(Integer(-1))) * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.I * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_87():
    f = f**(a + c*x**2)*sin(d + f*x**2)
    F = (((sympy.I * (Integer(4))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * (((sympy.E)**((sympy.I * Symbol('d'))) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**((sympy.I * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * (sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_88():
    f = f**(a + c*x**2)*sin(d + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(8) * (sympy.E)**(((Integer(2) * sympy.I) * Symbol('d'))) * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * sympy.I) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(8) * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_89():
    f = f**(a + c*x**2)*sin(d + f*x**2)**3
    F = ((((Integer(3) * sympy.I) * (Integer(16))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * (((sympy.E)**((sympy.I * Symbol('d'))) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(16))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * (((sympy.E)**(((Integer(3) * sympy.I) * Symbol('d'))) * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * sympy.I) * (Integer(16))**(Integer(-1))) * (sympy.E)**((sympy.I * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * (sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))) + (((sympy.I * (Integer(16))**(Integer(-1))) * (sympy.E)**(((Integer(3) * sympy.I) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * (sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_90():
    f = f**(a + c*x**2)*sin(d + e*x + f*x**2)
    F = ((sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_91():
    f = f**(a + c*x**2)*sin(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * (sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * (sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_92():
    f = f**(a + c*x**2)*sin(d + e*x + f*x**2)**3
    F = ((Integer(3) * sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + ((sympy.I * (sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_93():
    f = f**(a + b*x + c*x**2)*sin(d + e*x)
    F = ((((Integer(-1) * sympy.I) * (Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.I * (Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_94():
    f = f**(a + b*x + c*x**2)*sin(d + e*x)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(2) * sympy.I) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (Integer(-1) * (((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_95():
    f = f**(a + b*x + c*x**2)*sin(d + e*x)**3
    F = ((((Integer(-3) * sympy.I) * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.I * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * sympy.I) * (Integer(16))**(Integer(-1))) * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.I * (Integer(16))**(Integer(-1))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (Integer(-1) * (((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_96():
    f = f**(a + b*x + c*x**2)*sin(d + f*x**2)
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_97():
    f = f**(a + b*x + c*x**2)*sin(d + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_98():
    f = f**(a + b*x + c*x**2)*sin(d + f*x**2)**3
    F = (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + ((sympy.I * (sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(12) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + ((sympy.I * (sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_99():
    f = f**(a + b*x + c*x**2)*sin(d + e*x + f*x**2)
    F = ((sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_100():
    f = f**(a + b*x + c*x**2)*sin(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_101():
    f = f**(a + b*x + c*x**2)*sin(d + e*x + f*x**2)**3
    F = ((Integer(3) * sympy.I * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + ((sympy.I * (sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(3) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_102():
    f = f**(a + b*x + c*x**2)*sin(a + b*x + e*x**2)
    F = ((sympy.I * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * (sympy.I + (Integer(-1) * sympy.log(Symbol('f'))))) + (Integer(2) * x * ((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * (((sympy.E)**(((sympy.I + (Integer(-1) * sympy.log(Symbol('f')))) * (Symbol('a') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.I + (Integer(-1) * sympy.log(Symbol('f'))))) * (((Integer(4) * sympy.I * Symbol('e')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))))) * (Integer(4) * sympy.sqrt(((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.E)**(((sympy.I + sympy.log(Symbol('f'))) * (Symbol('a') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.I + sympy.log(Symbol('f')))) * (((Integer(4) * sympy.I * Symbol('e')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * (sympy.I + sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_103():
    f = exp(x)*cos(a + b*x)
    F = b*exp(x)*sin(a + b*x)/(b**2 + 1) + exp(x)*cos(a + b*x)/(b**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_104():
    f = exp(x)*cos(a + c*x**2)
    F = ((Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('a')) + (Symbol('c'))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + ((Integer(2) * sympy.I) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * Symbol('c') * x)))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('a')) + (Symbol('c'))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_105():
    f = exp(x)*cos(a + b*x + c*x**2)
    F = (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((Integer(4))**(Integer(-1)) * sympy.I * ((Integer(4) * Symbol('a')) + (((Integer(1) + (sympy.I * Symbol('b'))))**(Integer(2)) * (Symbol('c'))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.I * Symbol('b')) + (Integer(2) * sympy.I * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((sympy.I * ((sympy.I + Symbol('b')))**(Integer(2))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('b'))) + (Integer(-1) * (Integer(2) * sympy.I * Symbol('c') * x)))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_106():
    f = exp(x**2)*cos(a + b*x)
    F = (((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(-1) * sympy.I) * Symbol('b')) + (Integer(2) * x)) * (Integer(2))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (((sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(2) * x)) * (Integer(2))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_107():
    f = exp(x**2)*cos(a + c*x**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))) * x))) * ((Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))) * (sympy.E)**((sympy.I * Symbol('a')))))**(Integer(-1))) + (((sympy.E)**((sympy.I * Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (sympy.I * Symbol('c')))) * x))) * ((Integer(4) * sympy.sqrt((Integer(1) + (sympy.I * Symbol('c'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_108():
    f = exp(x**2)*cos(a + b*x + c*x**2)
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(-1) * (Integer(2) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * x))) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))))))**(Integer(-1))))) * (((sympy.E)**((sympy.I * (Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (((Integer(4) * sympy.I) + (Integer(4) * Symbol('c'))))**(Integer(-1))))))) * (Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c'))))))))**(Integer(-1)))) + (((sympy.E)**(((sympy.I * Symbol('a')) + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + (sympy.I * Symbol('c')))))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('b')) + (Integer(2) * (Integer(1) + (sympy.I * Symbol('c'))) * x)) * ((Integer(2) * sympy.sqrt((Integer(1) + (sympy.I * Symbol('c'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + (sympy.I * Symbol('c'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_109():
    f = f**(a + b*x)*cos(d + f*x**2)
    F = ((Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_110():
    f = f**(a + b*x)*cos(d + f*x**2)**2
    F = (((Integer(-1) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(4) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) + (Integer(-1) * ((((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(4) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) * ((sympy.E)**(((sympy.I * (Integer(8))**(Integer(-1))) * ((Integer(16) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_111():
    f = f**(a + b*x)*cos(d + f*x**2)**3
    F = ((Integer(-3) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erf')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(6) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * (((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((Integer(16) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erfi')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(6) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * ((sympy.E)**(((sympy.I * (Integer(12))**(Integer(-1))) * ((Integer(36) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_112():
    f = f**(a + b*x)*cos(d + e*x + f*x**2)
    F = ((Integer(-1) * ((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + ((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(4))**(Integer(-1))) * ((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_113():
    f = f**(a + b*x)*cos(d + e*x + f*x**2)**2
    F = (((Integer(-1) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (sympy.I * (Integer(16))**(Integer(-1))))) * (sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * ((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(2) * sympy.I) * Symbol('e')) + ((Integer(4) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(8))**(Integer(-1))) * (((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(4))**(Integer(-1)) + (sympy.I * (Integer(4))**(Integer(-1)))) * (((Integer(2) * sympy.I) * Symbol('e')) + ((Integer(4) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_114():
    f = f**(a + b*x)*cos(d + e*x + f*x**2)**3
    F = ((Integer(-3) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**(((sympy.I * (Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('d')) + ((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (Symbol('f'))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * ((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erf')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(3) * sympy.I) * Symbol('e')) + ((Integer(6) * sympy.I) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(-1))**((Integer(4))**(Integer(-1))) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(4))**(Integer(-1))) * ((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(-1))**((Integer(4))**(Integer(-1))) * ((sympy.I * Symbol('e')) + ((Integer(2) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (Integer(-1) * (((Integer(16))**(Integer(-1)) + (sympy.I * (Integer(16))**(Integer(-1)))) * (sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + (((sympy.I * (Integer(12))**(Integer(-1))) * (((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('Erfi')(((((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (((Integer(3) * sympy.I) * Symbol('e')) + ((Integer(6) * sympy.I) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))))) * ((sympy.sqrt(Integer(6)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_115():
    f = f**(a + c*x**2)*cos(d + e*x)
    F = ((Integer(-1) * ((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_116():
    f = f**(a + c*x**2)*cos(d + e*x)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * x * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('c') * x * sympy.log(Symbol('f')))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_117():
    f = f**(a + c*x**2)*cos(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_118():
    f = f**(a + c*x**2)*cos(d + f*x**2)
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(4) * (sympy.E)**((sympy.I * Symbol('d'))) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((sympy.I * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_119():
    f = f**(a + c*x**2)*cos(d + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(8) * (sympy.E)**(((Integer(2) * sympy.I) * Symbol('d'))) * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * sympy.I) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(8) * sympy.sqrt((((Integer(2) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_120():
    f = f**(a + c*x**2)*cos(d + f*x**2)**3
    F = ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**((sympy.I * Symbol('d'))) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**(((Integer(3) * sympy.I) * Symbol('d'))) * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((sympy.I * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * sympy.I) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt((((Integer(3) * sympy.I) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_121():
    f = f**(a + c*x**2)*cos(d + e*x + f*x**2)
    F = (((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_122():
    f = f**(a + c*x**2)*cos(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * (sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * (sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_123():
    f = f**(a + c*x**2)*cos(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(((sympy.I * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_124():
    f = f**(a + b*x + c*x**2)*cos(d + e*x)
    F = ((Integer(-1) * ((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_125():
    f = f**(a + b*x + c*x**2)*cos(d + e*x)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(-2) * sympy.I) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(2) * sympy.I) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**((((Integer(2) * sympy.I) * Symbol('d')) + (Integer(-1) * (((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(2) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_126():
    f = f**(a + b*x + c*x**2)*cos(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((((Integer(-3) * sympy.I) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((((Integer(3) * sympy.I) * Symbol('d')) + (Integer(-1) * (((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((((Integer(3) * sympy.I) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_127():
    f = f**(a + b*x + c*x**2)*cos(d + f*x**2)
    F = (Integer(-1) * (((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (((sympy.E)**(((sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_128():
    f = f**(a + b*x + c*x**2)*cos(d + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_129():
    f = f**(a + b*x + c*x**2)*cos(d + f*x**2)**3
    F = (Integer(-1) * ((Integer(3) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(12) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**(((sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_130():
    f = f**(a + b*x + c*x**2)*cos(d + e*x + f*x**2)
    F = (((sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_131():
    f = f**(a + b*x + c*x**2)*cos(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * sympy.I * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(8) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_132():
    f = f**(a + b*x + c*x**2)*cos(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**((((Integer(-1) * sympy.I) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (sympy.I * Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(((sympy.I * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (sympy.I * Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(4) * sympy.I * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * sympy.I * Symbol('d')) + (Integer(-1) * ((((Integer(3) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * sympy.I * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * sympy.I * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_133():
    f = f**(a + b*x + c*x**2)*cos(a + b*x + e*x**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * (sympy.I + (Integer(-1) * sympy.log(Symbol('f'))))) + (Integer(2) * x * ((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * (((sympy.E)**(((sympy.I + (Integer(-1) * sympy.log(Symbol('f')))) * (Symbol('a') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.I + (Integer(-1) * sympy.log(Symbol('f'))))) * (((Integer(4) * sympy.I * Symbol('e')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))))) * (Integer(4) * sympy.sqrt(((sympy.I * Symbol('e')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))))**(Integer(-1))) + (((sympy.E)**(((sympy.I + sympy.log(Symbol('f'))) * (Symbol('a') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.I + sympy.log(Symbol('f')))) * (((Integer(4) * sympy.I * Symbol('e')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * (sympy.I + sympy.log(Symbol('f')))) + (Integer(2) * x * ((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((sympy.I * Symbol('e')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_134():
    f = F**(c*(a + b*x))*(f*sin(d + e*x) + f)**2
    F = F**(a*c + b*c*x)*b*c*f**2*log(F)*sin(d + e*x)**2/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(a*c + b*c*x)*b*c*f**2*log(F)*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2) - 2*F**(a*c + b*c*x)*e*f**2*sin(d + e*x)*cos(d + e*x)/(b**2*c**2*log(F)**2 + 4*e**2) - 2*F**(a*c + b*c*x)*e*f**2*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e**2*f**2/(b*c*(b**2*c**2*log(F)**2 + 4*e**2)*log(F)) + F**(a*c + b*c*x)*f**2/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_135():
    f = F**(c*(a + b*x))*(f*sin(d + e*x) + f)
    F = F**(a*c + b*c*x)*b*c*f*log(F)*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2) - F**(a*c + b*c*x)*e*f*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*f/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_136():
    f = F**(c*(a + b*x))/(f*sin(d + e*x) + f)
    F = -2*F**(c*(a + b*x))*exp(I*(d + e*x))*hyper((2, -I*b*c*log(F)/e + 1), (-I*b*c*log(F)/e + 2,), I*exp(I*(d + e*x)))/(f*(-I*b*c*log(F) + e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_137():
    f = F**(c*(a + b*x))/(f*sin(d + e*x) + f)**2
    F = -F**(c*(a + b*x))*b*c*log(F)*csc(d/2 + e*x/2 + pi/4)**2/(6*e**2*f**2) - F**(c*(a + b*x))*cot(d/2 + e*x/2 + pi/4)*csc(d/2 + e*x/2 + pi/4)**2/(6*e*f**2) - 2*F**(c*(a + b*x))*(I*b*c*log(F) + e)*exp(I*(d + e*x))*hyper((2, -I*b*c*log(F)/e + 1), (-I*b*c*log(F)/e + 2,), I*exp(I*(d + e*x)))/(3*e**2*f**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_138():
    f = F**(c*(a + b*x))*(f*cos(d + e*x) + f)**2
    F = F**(a*c + b*c*x)*b*c*f**2*log(F)*cos(d + e*x)**2/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(a*c + b*c*x)*b*c*f**2*log(F)*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e*f**2*sin(d + e*x)*cos(d + e*x)/(b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(a*c + b*c*x)*e*f**2*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e**2*f**2/(b*c*(b**2*c**2*log(F)**2 + 4*e**2)*log(F)) + F**(a*c + b*c*x)*f**2/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_139():
    f = F**(c*(a + b*x))*(f*cos(d + e*x) + f)
    F = F**(a*c + b*c*x)*b*c*f*log(F)*cos(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*e*f*sin(d + e*x)/(b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*f/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_140():
    f = F**(c*(a + b*x))/(f*cos(d + e*x) + f)
    F = 2*F**(c*(a + b*x))*exp(I*(d + e*x))*hyper((2, -I*b*c*log(F)/e + 1), (-I*b*c*log(F)/e + 2,), -exp(I*(d + e*x)))/(f*(b*c*log(F) + I*e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_6_f_pow_a_plus_b_x_plus_c_x_pow_2_trig_d_plus_e_x_plus_f_x_pow_2_pow_n_141():
    f = F**(c*(a + b*x))/(f*cos(d + e*x) + f)**2
    F = -F**(c*(a + b*x))*b*c*log(F)*sec(d/2 + e*x/2)**2/(6*e**2*f**2) + F**(c*(a + b*x))*tan(d/2 + e*x/2)*sec(d/2 + e*x/2)**2/(6*e*f**2) - 2*F**(c*(a + b*x))*(-b*c*log(F) + I*e)*exp(I*(d + e*x))*hyper((2, -I*b*c*log(F)/e + 1), (-I*b*c*log(F)/e + 2,), -exp(I*(d + e*x)))/(3*e**2*f**2)
    assert integrate(f, x) == F

