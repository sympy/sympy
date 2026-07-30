"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.3.1 (a+b csc)^m (d csc)^n (A+B csc).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, a, c, d = symbols('A a c d')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_1():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)**3
    F = -2*A*a*cot(c + d*x)**3/(3*d) - A*a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 7*A*a*cot(c + d*x)*csc(c + d*x)/(8*d) - 2*A*a*cot(c + d*x)/d - 7*A*a*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_2():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)**2
    F = -A*a*cot(c + d*x)*csc(c + d*x)**2/(3*d) - A*a*cot(c + d*x)*csc(c + d*x)/d - 5*A*a*cot(c + d*x)/(3*d) - A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_3():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)
    F = -A*a*cot(c + d*x)*csc(c + d*x)/(2*d) - 2*A*a*cot(c + d*x)/d - 3*A*a*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_4():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*sin(c + d*x)
    F = 2*A*a*x - A*a*cos(c + d*x)/d - A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_5():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*sin(c + d*x)**2
    F = 3*A*a*x/2 - A*a*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*A*a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_6():
    f = (A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*sin(c + d*x)**3
    F = A*a*x - A*a*sin(c + d*x)*cos(c + d*x)/d + A*a*cos(c + d*x)**3/(3*d) - 2*A*a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_7():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)**3
    F = A*a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - A*a*cot(c + d*x)*csc(c + d*x)/(8*d) - A*a*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_8():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)**2
    F = A*a*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_9():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)
    F = A*a*cot(c + d*x)*csc(c + d*x)/(2*d) - A*a*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_10():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*sin(c + d*x)
    F = -A*a*cos(c + d*x)/d + A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_11():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*sin(c + d*x)**2
    F = -A*a*x/2 - A*a*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_12():
    f = (A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*sin(c + d*x)**3
    F = A*a*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_13():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)**3
    F = A*a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - A*a*cot(c + d*x)*csc(c + d*x)/(8*d) - A*a*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_14():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)**2
    F = A*a*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_15():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)*csc(c + d*x)
    F = A*a*cot(c + d*x)*csc(c + d*x)/(2*d) - A*a*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_16():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)/csc(c + d*x)
    F = -A*a*cos(c + d*x)/d + A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_17():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)/csc(c + d*x)**2
    F = -A*a*x/2 - A*a*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_18():
    f = (-A*csc(c + d*x) + A)*(a*csc(c + d*x) + a)/csc(c + d*x)**3
    F = A*a*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_19():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)**3
    F = 2*A*a*cot(c + d*x)**3/(3*d) - A*a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 7*A*a*cot(c + d*x)*csc(c + d*x)/(8*d) + 2*A*a*cot(c + d*x)/d - 7*A*a*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_20():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)**2
    F = -A*a*cot(c + d*x)*csc(c + d*x)**2/(3*d) + A*a*cot(c + d*x)*csc(c + d*x)/d - 5*A*a*cot(c + d*x)/(3*d) + A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_21():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)*csc(c + d*x)
    F = -A*a*cot(c + d*x)*csc(c + d*x)/(2*d) + 2*A*a*cot(c + d*x)/d - 3*A*a*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_22():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)/csc(c + d*x)
    F = -2*A*a*x - A*a*cos(c + d*x)/d - A*a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_23():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)/csc(c + d*x)**2
    F = 3*A*a*x/2 - A*a*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*A*a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_3_1_a_plus_b_csc_pow_m_d_csc_pow_n_A_plus_B_csc_24():
    f = (-A*csc(c + d*x) + A)*(-a*csc(c + d*x) + a)/csc(c + d*x)**3
    F = -A*a*x + A*a*sin(c + d*x)*cos(c + d*x)/d + A*a*cos(c + d*x)**3/(3*d) - 2*A*a*cos(c + d*x)/d
    assert integrate(f, x) == F

