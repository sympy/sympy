"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.2.1 (a+b cot)^m (c+d cot)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_1():
    f = (I*a*cot(c + d*x) + a)**n
    F = I*(I*a*cot(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*cot(c + d*x)/2 + sympy.S.Half)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_2():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a)
    F = -sqrt(2)*a*e**(sympy.S(5)/2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d + 2*a*e**2*sqrt(e*cot(c + d*x))/d - 2*a*e*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - 2*a*(e*cot(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_3():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)
    F = -sqrt(2)*a*e**(sympy.S(3)/2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d - 2*a*e*sqrt(e*cot(c + d*x))/d - 2*a*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_4():
    f = sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)
    F = sqrt(2)*a*sqrt(e)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d - 2*a*sqrt(e*cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_5():
    f = (a*cot(c + d*x) + a)/sqrt(e*cot(c + d*x))
    F = sqrt(2)*a*atan(sqrt(2)*sqrt(e)*(1 - cot(c + d*x))/(2*sqrt(e*cot(c + d*x))))/(d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_6():
    f = (a*cot(c + d*x) + a)/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = 2*a/(d*e*sqrt(e*cot(c + d*x))) - sqrt(2)*a*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_7():
    f = (a*cot(c + d*x) + a)/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 2*a/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) + 2*a/(d*e**2*sqrt(e*cot(c + d*x))) - sqrt(2)*a*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_8():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a)**2
    F = sqrt(2)*a**2*e**(sympy.S(5)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) - sqrt(2)*a**2*e**(sympy.S(5)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) + sqrt(2)*a**2*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d - sqrt(2)*a**2*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d + 4*a**2*e**2*sqrt(e*cot(c + d*x))/d - 4*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)/(5*d) - 2*a**2*(e*cot(c + d*x))**(sympy.S(7)/2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_9():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)**2
    F = sqrt(2)*a**2*e**(sympy.S(3)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) - sqrt(2)*a**2*e**(sympy.S(3)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) - sqrt(2)*a**2*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d + sqrt(2)*a**2*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d - 4*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - 2*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_10():
    f = sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)**2
    F = -sqrt(2)*a**2*sqrt(e)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) + sqrt(2)*a**2*sqrt(e)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d) - sqrt(2)*a**2*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d + sqrt(2)*a**2*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/d - 4*a**2*sqrt(e*cot(c + d*x))/d - 2*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_11():
    f = (a*cot(c + d*x) + a)**2/sqrt(e*cot(c + d*x))
    F = -2*a**2*sqrt(e*cot(c + d*x))/(d*e) - sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*sqrt(e)) + sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*sqrt(e)) + sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*sqrt(e)) - sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_12():
    f = (a*cot(c + d*x) + a)**2/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = 2*a**2/(d*e*sqrt(e*cot(c + d*x))) + sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(3)/2)) - sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2)) - sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_13():
    f = (a*cot(c + d*x) + a)**2/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 2*a**2/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) + 4*a**2/(d*e**2*sqrt(e*cot(c + d*x))) + sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(5)/2)) - sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(5)/2)) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2)) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_14():
    f = (a*cot(c + d*x) + a)**2/(e*cot(c + d*x))**(sympy.S(7)/2)
    F = 2*a**2/(5*d*e*(e*cot(c + d*x))**(sympy.S(5)/2)) + 4*a**2/(3*d*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(7)/2)) + sqrt(2)*a**2*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(2*d*e**(sympy.S(7)/2)) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(7)/2)) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(d*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_15():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a)**3
    F = 2*sqrt(2)*a**3*e**(sympy.S(5)/2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d + 4*a**3*e**2*sqrt(e*cot(c + d*x))/d + 4*a**3*e*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - 4*a**3*(e*cot(c + d*x))**(sympy.S(5)/2)/(5*d) - 40*a**3*(e*cot(c + d*x))**(sympy.S(7)/2)/(63*d*e) - 2*(e*cot(c + d*x))**(sympy.S(7)/2)*(a**3*cot(c + d*x) + a**3)/(9*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_16():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)**3
    F = -2*sqrt(2)*a**3*e**(sympy.S(3)/2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d + 4*a**3*e*sqrt(e*cot(c + d*x))/d - 4*a**3*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - 32*a**3*(e*cot(c + d*x))**(sympy.S(5)/2)/(35*d*e) - 2*(e*cot(c + d*x))**(sympy.S(5)/2)*(a**3*cot(c + d*x) + a**3)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_17():
    f = sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)**3
    F = -2*sqrt(2)*a**3*sqrt(e)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/d - 4*a**3*sqrt(e*cot(c + d*x))/d - 8*a**3*(e*cot(c + d*x))**(sympy.S(3)/2)/(5*d*e) - 2*(e*cot(c + d*x))**(sympy.S(3)/2)*(a**3*cot(c + d*x) + a**3)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_18():
    f = (a*cot(c + d*x) + a)**3/sqrt(e*cot(c + d*x))
    F = -16*a**3*sqrt(e*cot(c + d*x))/(3*d*e) + 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*sqrt(e)) - 2*sqrt(e*cot(c + d*x))*(a**3*cot(c + d*x) + a**3)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_19():
    f = (a*cot(c + d*x) + a)**3/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = -4*a**3*sqrt(e*cot(c + d*x))/(d*e**2) + 2*sqrt(2)*a**3*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(3)/2)) + (2*a**3*cot(c + d*x) + 2*a**3)/(d*e*sqrt(e*cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_20():
    f = (a*cot(c + d*x) + a)**3/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 16*a**3/(3*d*e**2*sqrt(e*cot(c + d*x))) - 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(5)/2)) + (2*a**3*cot(c + d*x) + 2*a**3)/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_21():
    f = (a*cot(c + d*x) + a)**3/(e*cot(c + d*x))**(sympy.S(7)/2)
    F = 8*a**3/(5*d*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)) + 4*a**3/(d*e**3*sqrt(e*cot(c + d*x))) - 2*sqrt(2)*a**3*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(7)/2)) + (2*a**3*cot(c + d*x) + 2*a**3)/(5*d*e*(e*cot(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_22():
    f = (a*cot(c + d*x) + a)**3/(e*cot(c + d*x))**(sympy.S(9)/2)
    F = 32*a**3/(35*d*e**2*(e*cot(c + d*x))**(sympy.S(5)/2)) + 4*a**3/(3*d*e**3*(e*cot(c + d*x))**(sympy.S(3)/2)) - 4*a**3/(d*e**4*sqrt(e*cot(c + d*x))) + 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(d*e**(sympy.S(9)/2)) + (2*a**3*cot(c + d*x) + 2*a**3)/(7*d*e*(e*cot(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_23():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a*cot(c + d*x) + a)
    F = e**(sympy.S(5)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d) - sqrt(2)*e**(sympy.S(5)/2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(2*a*d) - 2*e**2*sqrt(e*cot(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_24():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a*cot(c + d*x) + a)
    F = -e**(sympy.S(3)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d) + sqrt(2)*e**(sympy.S(3)/2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_25():
    f = sqrt(e*cot(c + d*x))/(a*cot(c + d*x) + a)
    F = sqrt(e)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d) + sqrt(2)*sqrt(e)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_26():
    f = 1/(sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a))
    F = -atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d*sqrt(e)) - sqrt(2)*atanh(sqrt(2)*sqrt(e)*(cot(c + d*x) + 1)/(2*sqrt(e*cot(c + d*x))))/(2*a*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_27():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a))
    F = 2/(a*d*e*sqrt(e*cot(c + d*x))) + atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d*e**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(2*a*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_28():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a))
    F = 2/(3*a*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) - 2/(a*d*e**2*sqrt(e*cot(c + d*x))) - atan(sqrt(e*cot(c + d*x))/sqrt(e))/(a*d*e**(sympy.S(5)/2)) + sqrt(2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(2*a*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_29():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a*cot(c + d*x) + a)**2
    F = e**2*sqrt(e*cot(c + d*x))/(2*d*(a**2*cot(c + d*x) + a**2)) - sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) + sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) - 3*e**(sympy.S(5)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d) - sqrt(2)*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d) + sqrt(2)*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_30():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a*cot(c + d*x) + a)**2
    F = -e*sqrt(e*cot(c + d*x))/(2*d*(a**2*cot(c + d*x) + a**2)) - sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) + sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) + e**(sympy.S(3)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d) + sqrt(2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d) - sqrt(2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_31():
    f = sqrt(e*cot(c + d*x))/(a*cot(c + d*x) + a)**2
    F = sqrt(e*cot(c + d*x))/(2*d*(a**2*cot(c + d*x) + a**2)) + sqrt(2)*sqrt(e)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) - sqrt(2)*sqrt(e)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d) + sqrt(e)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d) + sqrt(2)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d) - sqrt(2)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_32():
    f = 1/(sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)**2)
    F = -sqrt(e*cot(c + d*x))/(2*d*e*(a**2*cot(c + d*x) + a**2)) + sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*sqrt(e)) - sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*sqrt(e)) - 3*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d*sqrt(e)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*sqrt(e)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_33():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)**2)
    F = -1/(2*d*e*sqrt(e*cot(c + d*x))*(a**2*cot(c + d*x) + a**2)) + 5/(2*a**2*d*e*sqrt(e*cot(c + d*x))) - sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*e**(sympy.S(3)/2)) + sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*e**(sympy.S(3)/2)) + 5*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d*e**(sympy.S(3)/2)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*e**(sympy.S(3)/2)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_34():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a)**2)
    F = -1/(2*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)*(a**2*cot(c + d*x) + a**2)) + 7/(6*a**2*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) - 9/(2*a**2*d*e**2*sqrt(e*cot(c + d*x))) - sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*e**(sympy.S(5)/2)) + sqrt(2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(8*a**2*d*e**(sympy.S(5)/2)) - 7*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(2*a**2*d*e**(sympy.S(5)/2)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*e**(sympy.S(5)/2)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(4*a**2*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_35():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a*cot(c + d*x) + a)**3
    F = e**2*sqrt(e*cot(c + d*x))/(4*a*d*(a*cot(c + d*x) + a)**2) - e**(sympy.S(5)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d) + sqrt(2)*e**(sympy.S(5)/2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d) - 5*e**2*sqrt(e*cot(c + d*x))/(8*a**3*d*(cot(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_36():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a*cot(c + d*x) + a)**3
    F = e*sqrt(e*cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + a**3)) - e*sqrt(e*cot(c + d*x))/(4*a*d*(a*cot(c + d*x) + a)**2) + 5*e**(sympy.S(3)/2)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d) + sqrt(2)*e**(sympy.S(3)/2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_37():
    f = sqrt(e*cot(c + d*x))/(a*cot(c + d*x) + a)**3
    F = 3*sqrt(e*cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + a**3)) + sqrt(e*cot(c + d*x))/(4*a*d*(a*cot(c + d*x) + a)**2) - sqrt(e)*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d) - sqrt(2)*sqrt(e)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_38():
    f = 1/(sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)**3)
    F = -sqrt(e*cot(c + d*x))/(4*a*d*e*(a*cot(c + d*x) + a)**2) - 7*sqrt(e*cot(c + d*x))/(8*a**3*d*e*(cot(c + d*x) + 1)) - 11*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d*sqrt(e)) - sqrt(2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_39():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)**3)
    F = -1/(4*a*d*e*sqrt(e*cot(c + d*x))*(a*cot(c + d*x) + a)**2) + 27/(8*a**3*d*e*sqrt(e*cot(c + d*x))) - 9/(8*a**3*d*e*sqrt(e*cot(c + d*x))*(cot(c + d*x) + 1)) + 31*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d*e**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*(sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_40():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a*cot(c + d*x) + a)**3)
    F = -1/(4*a*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)*(a*cot(c + d*x) + a)**2) + 55/(24*a**3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) - 11/(8*a**3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)*(cot(c + d*x) + 1)) - 63/(8*a**3*d*e**2*sqrt(e*cot(c + d*x))) - 59*atan(sqrt(e*cot(c + d*x))/sqrt(e))/(8*a**3*d*e**(sympy.S(5)/2)) + sqrt(2)*atan(sqrt(2)*(-sqrt(e)*cot(c + d*x) + sqrt(e))/(2*sqrt(e*cot(c + d*x))))/(4*a**3*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_41():
    f = sqrt(cot(x) + 1)*cot(x)**2
    F = -2*(cot(x) + 1)**(sympy.S(3)/2)/3 + log(-sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(2*sqrt(2 + 2*sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(2*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2))) + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_42():
    f = sqrt(cot(x) + 1)*cot(x)
    F = -2*sqrt(cot(x) + 1) + sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*cot(x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(cot(x) + 1))) + sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*cot(x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(cot(x) + 1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_43():
    f = (cot(x) + 1)**(sympy.S(3)/2)*cot(x)**2
    F = -2*(cot(x) + 1)**(sympy.S(5)/2)/5 + 2*sqrt(cot(x) + 1) - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*cot(x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(cot(x) + 1))) - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*cot(x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(cot(x) + 1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_44():
    f = (cot(x) + 1)**(sympy.S(3)/2)*cot(x)
    F = -2*(cot(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(cot(x) + 1) - log(-sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(2*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(2*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2))) + sqrt(1 + sqrt(2))*atan((2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_45():
    f = cot(x)**2/sqrt(cot(x) + 1)
    F = -2*sqrt(cot(x) + 1) - log(-sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(4*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(4*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2 + sqrt(1 + sqrt(2))*atan((2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_46():
    f = cot(x)/sqrt(cot(x) + 1)
    F = sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*cot(x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(cot(x) + 1)))/2 + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*cot(x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(cot(x) + 1)))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_47():
    f = cot(x)**2/(cot(x) + 1)**(sympy.S(3)/2)
    F = sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*cot(x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(cot(x) + 1)))/2 + sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*cot(x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(cot(x) + 1)))/2 + 1/sqrt(cot(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_48():
    f = cot(x)/(cot(x) + 1)**(sympy.S(3)/2)
    F = -log(-sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) + sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2 - sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2 - 1/sqrt(cot(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_49():
    f = cot(x)**2/(cot(x) + 1)**(sympy.S(5)/2)
    F = sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*cot(x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(cot(x) + 1)))/4 + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*cot(x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(cot(x) + 1)))/4 - 1/sqrt(cot(x) + 1) + 1/(3*(cot(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_50():
    f = cot(x)/(cot(x) + 1)**(sympy.S(5)/2)
    F = log(-sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(8*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(cot(x) + 1) + cot(x) + 1 + sqrt(2))/(8*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/4 - sqrt(1 + sqrt(2))*atan((2*sqrt(cot(x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/4 - 1/(3*(cot(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_51():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))
    F = -2*a*e*sqrt(e*cot(c + d*x))/d - 2*b*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - sqrt(2)*e**(sympy.S(3)/2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*e**(sympy.S(3)/2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) - sqrt(2)*e**(sympy.S(3)/2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) + sqrt(2)*e**(sympy.S(3)/2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_52():
    f = sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))
    F = -2*b*sqrt(e*cot(c + d*x))/d + sqrt(2)*sqrt(e)*(a - b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*sqrt(e)*(a - b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*sqrt(e)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*sqrt(e)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_53():
    f = (a + b*cot(c + d*x))/sqrt(e*cot(c + d*x))
    F = sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)) + sqrt(2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)) - sqrt(2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_54():
    f = (a + b*cot(c + d*x))/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = 2*a/(d*e*sqrt(e*cot(c + d*x))) - sqrt(2)*(a - b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*(a - b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)) - sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_55():
    f = (a + b*cot(c + d*x))/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 2*a/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) + 2*b/(d*e**2*sqrt(e*cot(c + d*x))) - sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)) + sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)) - sqrt(2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) + sqrt(2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_56():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))**2
    F = -4*a*b*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d) - 2*b**2*(e*cot(c + d*x))**(sympy.S(5)/2)/(5*d*e) - sqrt(2)*e**(sympy.S(3)/2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*e**(sympy.S(3)/2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) - sqrt(2)*e**(sympy.S(3)/2)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) + sqrt(2)*e**(sympy.S(3)/2)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - e*sqrt(e*cot(c + d*x))*(2*a**2 - 2*b**2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_57():
    f = sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))**2
    F = -4*a*b*sqrt(e*cot(c + d*x))/d - 2*b**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(3*d*e) + sqrt(2)*sqrt(e)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*sqrt(e)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*sqrt(e)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*sqrt(e)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_58():
    f = (a + b*cot(c + d*x))**2/sqrt(e*cot(c + d*x))
    F = -2*b**2*sqrt(e*cot(c + d*x))/(d*e) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_59():
    f = (a + b*cot(c + d*x))**2/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = 2*a**2/(d*e*sqrt(e*cot(c + d*x))) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_60():
    f = (a + b*cot(c + d*x))**2/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 2*a**2/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) + 4*a*b/(d*e**2*sqrt(e*cot(c + d*x))) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_61():
    f = (a + b*cot(c + d*x))**2/(e*cot(c + d*x))**(sympy.S(7)/2)
    F = 2*a**2/(5*d*e*(e*cot(c + d*x))**(sympy.S(5)/2)) + 4*a*b/(3*d*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)) - (2*a**2 - 2*b**2)/(d*e**3*sqrt(e*cot(c + d*x))) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(7)/2)) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_62():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))**3
    F = -32*a*b**2*(e*cot(c + d*x))**(sympy.S(5)/2)/(35*d*e) - 2*a*e*sqrt(e*cot(c + d*x))*(a**2 - 3*b**2)/d - 2*b**2*(e*cot(c + d*x))**(sympy.S(5)/2)*(a + b*cot(c + d*x))/(7*d*e) - 2*b*(e*cot(c + d*x))**(sympy.S(3)/2)*(3*a**2 - b**2)/(3*d) - sqrt(2)*e**(sympy.S(3)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) + sqrt(2)*e**(sympy.S(3)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*e**(sympy.S(3)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*e**(sympy.S(3)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_63():
    f = sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))**3
    F = -8*a*b**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(5*d*e) - 2*b**2*(e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))/(5*d*e) - 2*b*sqrt(e*cot(c + d*x))*(3*a**2 - b**2)/d - sqrt(2)*sqrt(e)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*sqrt(e)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d) + sqrt(2)*sqrt(e)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*sqrt(e)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_64():
    f = (a + b*cot(c + d*x))**3/sqrt(e*cot(c + d*x))
    F = -16*a*b**2*sqrt(e*cot(c + d*x))/(3*d*e) - 2*b**2*sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))/(3*d*e) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_65():
    f = (a + b*cot(c + d*x))**3/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = 2*a**2*(a + b*cot(c + d*x))/(d*e*sqrt(e*cot(c + d*x))) - 2*b*sqrt(e*cot(c + d*x))*(a**2 + b**2)/(d*e**2) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_66():
    f = (a + b*cot(c + d*x))**3/(e*cot(c + d*x))**(sympy.S(5)/2)
    F = 16*a**2*b/(3*d*e**2*sqrt(e*cot(c + d*x))) + 2*a**2*(a + b*cot(c + d*x))/(3*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_67():
    f = (a + b*cot(c + d*x))**3/(e*cot(c + d*x))**(sympy.S(7)/2)
    F = 8*a**2*b/(5*d*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)) + 2*a**2*(a + b*cot(c + d*x))/(5*d*e*(e*cot(c + d*x))**(sympy.S(5)/2)) - 2*a*(a**2 - 3*b**2)/(d*e**3*sqrt(e*cot(c + d*x))) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(7)/2)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(7)/2)) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_68():
    f = (a + b*cot(c + d*x))**3/(e*cot(c + d*x))**(sympy.S(9)/2)
    F = 32*a**2*b/(35*d*e**2*(e*cot(c + d*x))**(sympy.S(5)/2)) + 2*a**2*(a + b*cot(c + d*x))/(7*d*e*(e*cot(c + d*x))**(sympy.S(7)/2)) - 2*a*(a**2 - 3*b**2)/(3*d*e**3*(e*cot(c + d*x))**(sympy.S(3)/2)) - 2*b*(3*a**2 - b**2)/(d*e**4*sqrt(e*cot(c + d*x))) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(9)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(9)/2)) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(9)/2)) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_69():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a + b*cot(c + d*x))
    F = 2*a**(sympy.S(5)/2)*e**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(5)/2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)) - sqrt(2)*e**(sympy.S(5)/2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)) - sqrt(2)*e**(sympy.S(5)/2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(5)/2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)) - 2*e**2*sqrt(e*cot(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_70():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a + b*cot(c + d*x))
    F = -2*a**(sympy.S(3)/2)*e**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(sqrt(b)*d*(a**2 + b**2)) - sqrt(2)*e**(sympy.S(3)/2)*(a - b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(3)/2)*(a - b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)) - sqrt(2)*e**(sympy.S(3)/2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(3)/2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_71():
    f = sqrt(e*cot(c + d*x))/(a + b*cot(c + d*x))
    F = 2*sqrt(a)*sqrt(b)*sqrt(e)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(d*(a**2 + b**2)) - sqrt(2)*sqrt(e)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)) + sqrt(2)*sqrt(e)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)) + sqrt(2)*sqrt(e)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)) - sqrt(2)*sqrt(e)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_72():
    f = 1/(sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x)))
    F = sqrt(2)*(a - b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)) - sqrt(2)*(a - b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)) + sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)) - 2*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(sqrt(a)*d*sqrt(e)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_73():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x)))
    F = sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)) - sqrt(2)*(a - b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)) + 2/(a*d*e*sqrt(e*cot(c + d*x))) + 2*b**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(a**(sympy.S(3)/2)*d*e**(sympy.S(3)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_74():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a + b*cot(c + d*x)))
    F = -sqrt(2)*(a - b)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)*(a**2 + b**2)) + sqrt(2)*(a + b)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(5)/2)*(a**2 + b**2)) + 2/(3*a*d*e*(e*cot(c + d*x))**(sympy.S(3)/2)) - 2*b/(a**2*d*e**2*sqrt(e*cot(c + d*x))) - 2*b**(sympy.S(7)/2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(a**(sympy.S(5)/2)*d*e**(sympy.S(5)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_75():
    f = (e*cot(c + d*x))**(sympy.S(7)/2)/(a + b*cot(c + d*x))**2
    F = a**(sympy.S(5)/2)*e**(sympy.S(7)/2)*(3*a**2 + 7*b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)**2) + a**2*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(b*d*(a + b*cot(c + d*x))*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(7)/2)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(7)/2)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(7)/2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(7)/2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) - e**3*sqrt(e*cot(c + d*x))*(3*a**2 + 2*b**2)/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_76():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a + b*cot(c + d*x))**2
    F = -a**(sympy.S(3)/2)*e**(sympy.S(5)/2)*(a**2 + 5*b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)**2) + a**2*e**2*sqrt(e*cot(c + d*x))/(b*d*(a + b*cot(c + d*x))*(a**2 + b**2)) + sqrt(2)*e**(sympy.S(5)/2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(5)/2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(5)/2)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(5)/2)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_77():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a + b*cot(c + d*x))**2
    F = -sqrt(a)*e**(sympy.S(3)/2)*(a**2 - 3*b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(sqrt(b)*d*(a**2 + b**2)**2) - a*e*sqrt(e*cot(c + d*x))/(d*(a + b*cot(c + d*x))*(a**2 + b**2)) - sqrt(2)*e**(sympy.S(3)/2)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(3)/2)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(3)/2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(3)/2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_78():
    f = sqrt(e*cot(c + d*x))/(a + b*cot(c + d*x))**2
    F = b*sqrt(e*cot(c + d*x))/(d*(a + b*cot(c + d*x))*(a**2 + b**2)) - sqrt(2)*sqrt(e)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) + sqrt(2)*sqrt(e)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**2) + sqrt(2)*sqrt(e)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) - sqrt(2)*sqrt(e)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**2) + sqrt(b)*sqrt(e)*(3*a**2 - b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(sqrt(a)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_79():
    f = 1/(sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))**2)
    F = sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)**2) - b**2*sqrt(e*cot(c + d*x))/(a*d*e*(a + b*cot(c + d*x))*(a**2 + b**2)) - b**(sympy.S(3)/2)*(5*a**2 + b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(a**(sympy.S(3)/2)*d*sqrt(e)*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_80():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))**2)
    F = sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)**2) - b**2/(a*d*e*sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))*(a**2 + b**2)) + (2*a**2 + 3*b**2)/(a**2*d*e*sqrt(e*cot(c + d*x))*(a**2 + b**2)) + b**(sympy.S(5)/2)*(7*a**2 + 3*b**2)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(a**(sympy.S(5)/2)*d*e**(sympy.S(3)/2)*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_81():
    f = (e*cot(c + d*x))**(sympy.S(9)/2)/(a + b*cot(c + d*x))**3
    F = a**(sympy.S(5)/2)*e**(sympy.S(9)/2)*(15*a**4 + 46*a**2*b**2 + 63*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*b**(sympy.S(7)/2)*d*(a**2 + b**2)**3) + a**2*e**2*(e*cot(c + d*x))**(sympy.S(5)/2)/(2*b*d*(a + b*cot(c + d*x))**2*(a**2 + b**2)) + a**2*e**3*(e*cot(c + d*x))**(sympy.S(3)/2)*(5*a**2 + 13*b**2)/(4*b**2*d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(9)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(9)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(9)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(9)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) - e**4*sqrt(e*cot(c + d*x))*(15*a**4 + 31*a**2*b**2 + 8*b**4)/(4*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_82():
    f = (e*cot(c + d*x))**(sympy.S(7)/2)/(a + b*cot(c + d*x))**3
    F = -a**(sympy.S(3)/2)*e**(sympy.S(7)/2)*(3*a**4 + 6*a**2*b**2 + 35*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*b**(sympy.S(5)/2)*d*(a**2 + b**2)**3) + a**2*e**2*(e*cot(c + d*x))**(sympy.S(3)/2)/(2*b*d*(a + b*cot(c + d*x))**2*(a**2 + b**2)) + a**2*e**3*sqrt(e*cot(c + d*x))*(3*a**2 + 11*b**2)/(4*b**2*d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) + sqrt(2)*e**(sympy.S(7)/2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(7)/2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(7)/2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(7)/2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_83():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)/(a + b*cot(c + d*x))**3
    F = -sqrt(a)*e**(sympy.S(5)/2)*(a**4 + 18*a**2*b**2 - 15*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*b**(sympy.S(3)/2)*d*(a**2 + b**2)**3) + a**2*e**2*sqrt(e*cot(c + d*x))/(2*b*d*(a + b*cot(c + d*x))**2*(a**2 + b**2)) - a*e**2*sqrt(e*cot(c + d*x))*(a**2 + 9*b**2)/(4*b*d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) - sqrt(2)*e**(sympy.S(5)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(5)/2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(5)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(5)/2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_84():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a + b*cot(c + d*x))**3
    F = -a*e*sqrt(e*cot(c + d*x))/(d*(a + b*cot(c + d*x))**2*(2*a**2 + 2*b**2)) - sqrt(2)*e**(sympy.S(3)/2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(3)/2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) - sqrt(2)*e**(sympy.S(3)/2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) + sqrt(2)*e**(sympy.S(3)/2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - e*sqrt(e*cot(c + d*x))*(3*a**2 - 5*b**2)/(4*d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) - e**(sympy.S(3)/2)*(3*a**4 - 26*a**2*b**2 + 3*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*sqrt(a)*sqrt(b)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_85():
    f = sqrt(e*cot(c + d*x))/(a + b*cot(c + d*x))**3
    F = b*sqrt(e*cot(c + d*x))/(d*(a + b*cot(c + d*x))**2*(2*a**2 + 2*b**2)) + sqrt(2)*sqrt(e)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - sqrt(2)*sqrt(e)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*(a**2 + b**2)**3) - sqrt(2)*sqrt(e)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) + sqrt(2)*sqrt(e)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*(a**2 + b**2)**3) + b*sqrt(e*cot(c + d*x))*(7*a**2 - b**2)/(4*a*d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) + sqrt(b)*sqrt(e)*(15*a**4 - 18*a**2*b**2 - b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*a**(sympy.S(3)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_86():
    f = 1/(sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))**3)
    F = sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*sqrt(e)*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*sqrt(e)*(a**2 + b**2)**3) - b**2*sqrt(e*cot(c + d*x))/(2*a*d*e*(a + b*cot(c + d*x))**2*(a**2 + b**2)) - b**2*sqrt(e*cot(c + d*x))*(11*a**2 + 3*b**2)/(4*a**2*d*e*(a + b*cot(c + d*x))*(a**2 + b**2)**2) - b**(sympy.S(3)/2)*(35*a**4 + 6*a**2*b**2 + 3*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*a**(sympy.S(5)/2)*d*sqrt(e)*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_87():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a + b*cot(c + d*x))**3)
    F = -sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 - sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(1 + sqrt(2)*sqrt(e*cot(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(e)*cot(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*cot(c + d*x)))/(4*d*e**(sympy.S(3)/2)*(a**2 + b**2)**3) - b**2/(2*a*d*e*sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))**2*(a**2 + b**2)) - b**2*(13*a**2 + 5*b**2)/(4*a**2*d*e*sqrt(e*cot(c + d*x))*(a + b*cot(c + d*x))*(a**2 + b**2)**2) + (8*a**4 + 31*a**2*b**2 + 15*b**4)/(4*a**3*d*e*sqrt(e*cot(c + d*x))*(a**2 + b**2)**2) + b**(sympy.S(5)/2)*(63*a**4 + 46*a**2*b**2 + 15*b**4)*atan(sqrt(b)*sqrt(e*cot(c + d*x))/(sqrt(a)*sqrt(e)))/(4*a**(sympy.S(7)/2)*d*e**(sympy.S(3)/2)*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_88():
    f = (a + b*cot(c + d*x))**n
    F = b*(a + b*cot(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*cot(c + d*x))/(a + sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))*(n + 1)) - b*(a + b*cot(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*cot(c + d*x))/(a - sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_89():
    f = (d*tan(e + f*x))**n*(a + b*cot(e + f*x))**m
    F = -(d*tan(e + f*x))**n*(a + b*cot(e + f*x))**m*cot(e + f*x)*appellf1(1 - n, 1, -m, 2 - n, -I*cot(e + f*x), -b*cot(e + f*x)/a)/(2*f*(1 - n)*(1 + b*cot(e + f*x)/a)**m) - (d*tan(e + f*x))**n*(a + b*cot(e + f*x))**m*cot(e + f*x)*appellf1(1 - n, 1, -m, 2 - n, I*cot(e + f*x), -b*cot(e + f*x)/a)/(2*f*(1 - n)*(1 + b*cot(e + f*x)/a)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_90():
    f = (I*cot(c + d*x) + 1)/sqrt(a + b*cot(c + d*x))
    F = 2*I*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_91():
    f = (-I*cot(c + d*x) + 1)/sqrt(a + b*cot(c + d*x))
    F = -2*I*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_92():
    f = (A + B*cot(c + d*x))/(a + b*cot(c + d*x))
    F = x*(A*a + B*b)/(a**2 + b**2) - (A*b - B*a)*log(a*sin(c + d*x) + b*cos(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_93():
    f = (A + B*cot(c + d*x))/(a + b*cot(c + d*x))**2
    F = x*(A*a**2 - A*b**2 + 2*B*a*b)/(a**2 + b**2)**2 - (2*A*a*b - B*a**2 + B*b**2)*log(a*sin(c + d*x) + b*cos(c + d*x))/(d*(a**2 + b**2)**2) + (A*b - B*a)/(d*(a + b*cot(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_94():
    f = (A + B*cot(c + d*x))/(a + b*cot(c + d*x))**3
    F = x*(A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(a**2 + b**2)**3 - (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*log(a*sin(c + d*x) + b*cos(c + d*x))/(d*(a**2 + b**2)**3) + (2*A*a*b - B*a**2 + B*b**2)/(d*(a + b*cot(c + d*x))*(a**2 + b**2)**2) + (A*b - B*a)/(d*(a + b*cot(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_95():
    f = (A + B*cot(c + d*x))*(a + b*cot(c + d*x))**(sympy.S(5)/2)
    F = -2*B*(a + b*cot(c + d*x))**(sympy.S(5)/2)/(5*d) + (a - I*b)**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*(I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/d - (a + b*cot(c + d*x))**(sympy.S(3)/2)*(2*A*b + 2*B*a)/(3*d) - sqrt(a + b*cot(c + d*x))*(4*A*a*b + 2*B*a**2 - 2*B*b**2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_96():
    f = (A + B*cot(c + d*x))*(a + b*cot(c + d*x))**(sympy.S(3)/2)
    F = -2*B*(a + b*cot(c + d*x))**(sympy.S(3)/2)/(3*d) + (a - I*b)**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(3)/2)*(I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*cot(c + d*x))*(2*A*b + 2*B*a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_97():
    f = (A + B*cot(c + d*x))*sqrt(a + b*cot(c + d*x))
    F = -2*B*sqrt(a + b*cot(c + d*x))/d + sqrt(a - I*b)*(I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/d - sqrt(a + I*b)*(I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_98():
    f = (-a + b*cot(c + d*x))*(a + b*cot(c + d*x))**(sympy.S(5)/2)
    F = -2*b*(a + b*cot(c + d*x))**(sympy.S(5)/2)/(5*d) + 2*b*sqrt(a + b*cot(c + d*x))*(a**2 + b**2)/d - (a - I*b)**(sympy.S(5)/2)*(I*a - b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(5)/2)*(I*a + b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_99():
    f = (-a + b*cot(c + d*x))*(a + b*cot(c + d*x))**(sympy.S(3)/2)
    F = -2*b*(a + b*cot(c + d*x))**(sympy.S(3)/2)/(3*d) + sqrt(2)*b*(a**2 + b**2)*log(a + b*cot(c + d*x) - sqrt(2)*sqrt(a + b*cot(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*(a**2 + b**2)*log(a + b*cot(c + d*x) + sqrt(2)*sqrt(a + b*cot(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*(a**2 + b**2)*atanh((-sqrt(2)*sqrt(a + b*cot(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(2)*b*(a**2 + b**2)*atanh((sqrt(2)*sqrt(a + b*cot(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_100():
    f = (-a + b*cot(c + d*x))*sqrt(a + b*cot(c + d*x))
    F = -2*b*sqrt(a + b*cot(c + d*x))/d - sqrt(2)*b*sqrt(a**2 + b**2)*log(a + b*cot(c + d*x) - sqrt(2)*sqrt(a + b*cot(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*sqrt(a**2 + b**2)*log(a + b*cot(c + d*x) + sqrt(2)*sqrt(a + b*cot(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*sqrt(a**2 + b**2)*atanh((-sqrt(2)*sqrt(a + b*cot(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(2)*b*sqrt(a**2 + b**2)*atanh((sqrt(2)*sqrt(a + b*cot(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_101():
    f = (A + B*cot(c + d*x))/sqrt(a + b*cot(c + d*x))
    F = -(I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + (I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_102():
    f = (A + B*cot(c + d*x))/(a + b*cot(c + d*x))**(sympy.S(3)/2)
    F = (2*A*b - 2*B*a)/(d*sqrt(a + b*cot(c + d*x))*(a**2 + b**2)) - (I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + (I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_103():
    f = (A + B*cot(c + d*x))/(a + b*cot(c + d*x))**(sympy.S(5)/2)
    F = (4*A*a*b - 2*B*a**2 + 2*B*b**2)/(d*sqrt(a + b*cot(c + d*x))*(a**2 + b**2)**2) + (2*A*b - 2*B*a)/(d*(a + b*cot(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) - (I*A - B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (I*A + B)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_104():
    f = (-a + b*cot(c + d*x))/sqrt(a + b*cot(c + d*x))
    F = (I*a + b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) - (I*a - b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_105():
    f = (-a + b*cot(c + d*x))/(a + b*cot(c + d*x))**(sympy.S(3)/2)
    F = -4*a*b/(d*sqrt(a + b*cot(c + d*x))*(a**2 + b**2)) + (I*a + b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - (I*a - b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_2_1_a_plus_b_cot_pow_m_c_plus_d_cot_pow_n_106():
    f = (-a + b*cot(c + d*x))/(a + b*cot(c + d*x))**(sympy.S(5)/2)
    F = -4*a*b/(d*(a + b*cot(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) - 2*b*(3*a**2 - b**2)/(d*sqrt(a + b*cot(c + d*x))*(a**2 + b**2)**2) + (I*a + b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - (I*a - b)*atanh(sqrt(a + b*cot(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F

