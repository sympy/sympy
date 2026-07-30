"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.3.1 (a+b tan)^m (c+d tan)^n (A+B tan).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_1():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*tan(c + d*x)**2
    F = I*B*a*tan(c + d*x)**3/(3*d) - a*x*(A - I*B) + a*(A - I*B)*tan(c + d*x)/d + a*(I*A + B)*log(cos(c + d*x))/d + a*(I*A + B)*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_2():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*tan(c + d*x)
    F = I*B*a*tan(c + d*x)**2/(2*d) - a*x*(I*A + B) - a*(A - I*B)*log(cos(c + d*x))/d + a*(I*A + B)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_3():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)
    F = I*B*a*tan(c + d*x)/d + a*x*(A - I*B) - a*(I*A + B)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_4():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)
    F = A*a*log(sin(c + d*x))/d - I*B*a*log(cos(c + d*x))/d + a*x*(I*A + B)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_5():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**2
    F = -A*a*cot(c + d*x)/d - a*x*(A - I*B) + a*(I*A + B)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_6():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**3
    F = -A*a*cot(c + d*x)**2/(2*d) - a*x*(I*A + B) - a*(A - I*B)*log(sin(c + d*x))/d - a*(I*A + B)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_7():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**4
    F = -A*a*cot(c + d*x)**3/(3*d) + a*x*(A - I*B) + a*(A - I*B)*cot(c + d*x)/d - a*(I*A + B)*log(sin(c + d*x))/d - a*(I*A + B)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_8():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**5
    F = -A*a*cot(c + d*x)**4/(4*d) + a*x*(I*A + B) + a*(A - I*B)*log(sin(c + d*x))/d + a*(A - I*B)*cot(c + d*x)**2/(2*d) - a*(I*A + B)*cot(c + d*x)**3/(3*d) + a*(I*A + B)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_9():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**2
    F = I*B*(I*a**2*tan(c + d*x) + a**2)*tan(c + d*x)**3/(4*d) - 2*a**2*x*(A - I*B) + 2*a**2*(A - I*B)*tan(c + d*x)/d - a**2*(4*A - 5*I*B)*tan(c + d*x)**3/(12*d) + 2*a**2*(I*A + B)*log(cos(c + d*x))/d + a**2*(I*A + B)*tan(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_10():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)
    F = A*(I*a*tan(c + d*x) + a)**2/(2*d) - I*B*(I*a*tan(c + d*x) + a)**3/(3*a*d) - 2*a**2*x*(I*A + B) - 2*a**2*(A - I*B)*log(cos(c + d*x))/d + a**2*(I*A + B)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_11():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2
    F = B*(I*a*tan(c + d*x) + a)**2/(2*d) + 2*a**2*x*(A - I*B) - a**2*(A - I*B)*tan(c + d*x)/d - 2*a**2*(I*A + B)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_12():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)
    F = A*a**2*log(sin(c + d*x))/d + I*B*(I*a**2*tan(c + d*x) + a**2)/d + 2*a**2*x*(I*A + B) + a**2*(A - 2*I*B)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_13():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**2
    F = -A*(I*a**2*tan(c + d*x) + a**2)*cot(c + d*x)/d + B*a**2*log(cos(c + d*x))/d - 2*a**2*x*(A - I*B) + a**2*(2*I*A + B)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_14():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**3
    F = -A*(I*a**2*tan(c + d*x) + a**2)*cot(c + d*x)**2/(2*d) - 2*a**2*x*(I*A + B) - 2*a**2*(A - I*B)*log(sin(c + d*x))/d - a**2*(3*I*A + 2*B)*cot(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_15():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**4
    F = -A*(I*a**2*tan(c + d*x) + a**2)*cot(c + d*x)**3/(3*d) + 2*a**2*x*(A - I*B) + 2*a**2*(A - I*B)*cot(c + d*x)/d - 2*a**2*(I*A + B)*log(sin(c + d*x))/d - a**2*(4*I*A + 3*B)*cot(c + d*x)**2/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_16():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**5
    F = -A*(I*a**2*tan(c + d*x) + a**2)*cot(c + d*x)**4/(4*d) + 2*a**2*x*(I*A + B) + 2*a**2*(A - I*B)*log(sin(c + d*x))/d + a**2*(A - I*B)*cot(c + d*x)**2/d + 2*a**2*(I*A + B)*cot(c + d*x)/d - a**2*(5*I*A + 4*B)*cot(c + d*x)**3/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_17():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**2
    F = I*B*a*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**3/(5*d) - 4*a**3*x*(A - I*B) + 4*a**3*(A - I*B)*tan(c + d*x)/d - a**3*(45*A - 47*I*B)*tan(c + d*x)**3/(60*d) + 4*a**3*(I*A + B)*log(cos(c + d*x))/d + 2*a**3*(I*A + B)*tan(c + d*x)**2/d - (5*A - 7*I*B)*(I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_18():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)
    F = A*(I*a*tan(c + d*x) + a)**3/(3*d) - I*B*(I*a*tan(c + d*x) + a)**4/(4*a*d) - 4*a**3*x*(I*A + B) - 4*a**3*(A - I*B)*log(cos(c + d*x))/d + 2*a**3*(I*A + B)*tan(c + d*x)/d + a*(A - I*B)*(I*a*tan(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_19():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3
    F = B*(I*a*tan(c + d*x) + a)**3/(3*d) + 4*a**3*x*(A - I*B) - 2*a**3*(A - I*B)*tan(c + d*x)/d - 4*a**3*(I*A + B)*log(cos(c + d*x))/d + a*(I*A + B)*(I*a*tan(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_20():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)
    F = A*a**3*log(sin(c + d*x))/d + I*B*a*(I*a*tan(c + d*x) + a)**2/(2*d) + 4*a**3*x*(I*A + B) + a**3*(3*A - 4*I*B)*log(cos(c + d*x))/d - (A - 2*I*B)*(I*a**3*tan(c + d*x) + a**3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_21():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**2
    F = -A*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)/d - 4*a**3*x*(A - I*B) + a**3*(I*A + 3*B)*log(cos(c + d*x))/d + a**3*(3*I*A + B)*log(sin(c + d*x))/d + (I*A - B)*(I*a**3*tan(c + d*x) + a**3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_22():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3
    F = -A*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**2/(2*d) + I*B*a**3*log(cos(c + d*x))/d - 4*a**3*x*(I*A + B) - a**3*(4*A - 3*I*B)*log(sin(c + d*x))/d - (2*I*A + B)*(I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_23():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**4
    F = -A*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**3/(3*d) + 4*a**3*x*(A - I*B) + a**3*(17*A - 15*I*B)*cot(c + d*x)/(6*d) - 4*a**3*(I*A + B)*log(sin(c + d*x))/d - (5*I*A + 3*B)*(I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)**2/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_24():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**5
    F = -A*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**4/(4*d) + 4*a**3*x*(I*A + B) + 4*a**3*(A - I*B)*log(sin(c + d*x))/d + a**3*(15*A - 14*I*B)*cot(c + d*x)**2/(12*d) + 4*a**3*(I*A + B)*cot(c + d*x)/d - (3*I*A + 2*B)*(I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)**3/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_25():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**6
    F = -A*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**5/(5*d) - 4*a**3*x*(A - I*B) - 4*a**3*(A - I*B)*cot(c + d*x)/d + a**3*(47*A - 45*I*B)*cot(c + d*x)**3/(60*d) + 4*a**3*(I*A + B)*log(sin(c + d*x))/d + 2*a**3*(I*A + B)*cot(c + d*x)**2/d - (7*I*A + 5*B)*(I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)**4/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_26():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*tan(c + d*x)**2
    F = I*B*a*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**3/(6*d) - 8*a**4*x*(A - I*B) + 8*a**4*(A - I*B)*tan(c + d*x)/d - a**4*(92*A - 93*I*B)*tan(c + d*x)**3/(60*d) + 8*a**4*(I*A + B)*log(cos(c + d*x))/d + 4*a**4*(I*A + B)*tan(c + d*x)**2/d - (2*A - 3*I*B)*(I*a**2*tan(c + d*x) + a**2)**2*tan(c + d*x)**3/(10*d) - (12*A - 13*I*B)*(I*a**4*tan(c + d*x) + a**4)*tan(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_27():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*tan(c + d*x)
    F = A*(I*a*tan(c + d*x) + a)**4/(4*d) - I*B*(I*a*tan(c + d*x) + a)**5/(5*a*d) - 8*a**4*x*(I*A + B) - 8*a**4*(A - I*B)*log(cos(c + d*x))/d + 4*a**4*(I*A + B)*tan(c + d*x)/d + a*(A - I*B)*(I*a*tan(c + d*x) + a)**3/(3*d) + (A - I*B)*(I*a**2*tan(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_28():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4
    F = B*(I*a*tan(c + d*x) + a)**4/(4*d) + 8*a**4*x*(A - I*B) - 4*a**4*(A - I*B)*tan(c + d*x)/d - 8*a**4*(I*A + B)*log(cos(c + d*x))/d + a*(I*A + B)*(I*a*tan(c + d*x) + a)**3/(3*d) + (I*A + B)*(I*a**2*tan(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_29():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)
    F = A*a**4*log(sin(c + d*x))/d + I*B*a*(I*a*tan(c + d*x) + a)**3/(3*d) + 8*a**4*x*(I*A + B) + a**4*(7*A - 8*I*B)*log(cos(c + d*x))/d - (A - 2*I*B)*(I*a**2*tan(c + d*x) + a**2)**2/(2*d) - (3*A - 4*I*B)*(I*a**4*tan(c + d*x) + a**4)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_30():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**2
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)/d - 3*B*(I*a**4*tan(c + d*x) + a**4)/d - 8*a**4*x*(A - I*B) + a**4*(4*I*A + B)*log(sin(c + d*x))/d + a**4*(4*I*A + 7*B)*log(cos(c + d*x))/d + (2*I*A - B)*(I*a**2*tan(c + d*x) + a**2)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_31():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**3
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**2/(2*d) - 3*A*(I*a**4*tan(c + d*x) + a**4)/d - 8*a**4*x*(I*A + B) - a**4*(A - 4*I*B)*log(cos(c + d*x))/d - a**4*(7*A - 4*I*B)*log(sin(c + d*x))/d - (5*I*A + 2*B)*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_32():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**4
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3/(3*d) - B*a**4*log(cos(c + d*x))/d + 8*a**4*x*(A - I*B) - a**4*(8*I*A + 7*B)*log(sin(c + d*x))/d + (4*A - 3*I*B)*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)/d - (2*I*A + B)*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_33():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**5
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**4/(4*d) + 8*a**4*x*(I*A + B) + 8*a**4*(A - I*B)*log(sin(c + d*x))/d + a**4*(67*I*A + 64*B)*cot(c + d*x)/(12*d) + (19*A - 16*I*B)*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)**2/(12*d) - (7*I*A + 4*B)*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**3/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_34():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**6
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**5/(5*d) - 8*a**4*x*(A - I*B) - 8*a**4*(A - I*B)*cot(c + d*x)/d + 8*a**4*(I*A + B)*log(sin(c + d*x))/d + a**4*(148*I*A + 145*B)*cot(c + d*x)**2/(60*d) + (28*A - 25*I*B)*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)**3/(30*d) - (8*I*A + 5*B)*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**4/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_35():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*cot(c + d*x)**7
    F = -A*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**6/(6*d) - 8*a**4*x*(I*A + B) - 8*a**4*(A - I*B)*log(sin(c + d*x))/d - 4*a**4*(A - I*B)*cot(c + d*x)**2/d - 8*a**4*(I*A + B)*cot(c + d*x)/d + a**4*(93*I*A + 92*B)*cot(c + d*x)**3/(60*d) + (13*A - 12*I*B)*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)**4/(20*d) - (3*I*A + 2*B)*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**5/(10*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_36():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**3/(2*d*(I*a*tan(c + d*x) + a)) + x*(3*I*A - 3*B)/(2*a) - (A + 2*I*B)*log(cos(c + d*x))/(a*d) - (A + 2*I*B)*tan(c + d*x)**2/(2*a*d) - (3*I*A - 3*B)*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_37():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**2/(2*d*(I*a*tan(c + d*x) + a)) + x*(A + 3*I*B)/(2*a) - (A + 3*I*B)*tan(c + d*x)/(2*a*d) + (I*A - B)*log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_38():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)
    F = I*B*log(cos(c + d*x))/(a*d) - x*(I*A - B)/(2*a) - (A + I*B)/(2*a*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_39():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)
    F = (I*A - B)/(2*d*(I*a*tan(c + d*x) + a)) + x*(A - I*B)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_40():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)
    F = A*log(sin(c + d*x))/(a*d) + (A + I*B)/(2*d*(I*a*tan(c + d*x) + a)) - x*(I*A - B)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_41():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)/(2*d*(I*a*tan(c + d*x) + a)) - x*(3*A + I*B)/(2*a) - (3*A + I*B)*cot(c + d*x)/(2*a*d) - (I*A - B)*log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_42():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**2/(2*d*(I*a*tan(c + d*x) + a)) + x*(3*I*A - 3*B)/(2*a) - (2*A + I*B)*log(sin(c + d*x))/(a*d) - (2*A + I*B)*cot(c + d*x)**2/(2*a*d) + (3*I*A - 3*B)*cot(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_43():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**4/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**3/(2*d*(I*a*tan(c + d*x) + a)) + x*(5*A + 3*I*B)/(2*a) - (5*A + 3*I*B)*cot(c + d*x)**3/(6*a*d) + (5*A + 3*I*B)*cot(c + d*x)/(2*a*d) + (I*A - B)*cot(c + d*x)**2/(a*d) + (2*I*A - 2*B)*log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_44():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = (I*A - B)*tan(c + d*x)**3/(4*d*(I*a*tan(c + d*x) + a)**2) - x*(3*I*A - 9*B)/(4*a**2) + (A + 2*I*B)*log(cos(c + d*x))/(a**2*d) + (A + 2*I*B)*tan(c + d*x)**2/(2*a**2*d*(I*tan(c + d*x) + 1)) + (3*I*A - 9*B)*tan(c + d*x)/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_45():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = B*log(cos(c + d*x))/(a**2*d) + (I*A - B)*tan(c + d*x)**2/(4*d*(I*a*tan(c + d*x) + a)**2) - x*(A + 3*I*B)/(4*a**2) + (I*A - 3*B)/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_46():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = -(A + I*B)/(4*d*(I*a*tan(c + d*x) + a)**2) - x*(I*A + B)/(4*a**2) + (A + 3*I*B)/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_47():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**2
    F = (I*A - B)/(4*d*(I*a*tan(c + d*x) + a)**2) + (I*A + B)/(4*d*(I*a**2*tan(c + d*x) + a**2)) + x*(A - I*B)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_48():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = A*log(sin(c + d*x))/(a**2*d) + (A + I*B)/(4*d*(I*a*tan(c + d*x) + a)**2) - x*(3*I*A - B)/(4*a**2) + (3*A + I*B)/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_49():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = (A + I*B)*cot(c + d*x)/(4*d*(I*a*tan(c + d*x) + a)**2) - x*(9*A + 3*I*B)/(4*a**2) + (2*A + I*B)*cot(c + d*x)/(2*a**2*d*(I*tan(c + d*x) + 1)) - (9*A + 3*I*B)*cot(c + d*x)/(4*a**2*d) - (2*I*A - B)*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_50():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = (A + I*B)*cot(c + d*x)**2/(4*d*(I*a*tan(c + d*x) + a)**2) + x*(15*I*A - 9*B)/(4*a**2) - (2*A + I*B)*cot(c + d*x)**2/(a**2*d) - (4*A + 2*I*B)*log(sin(c + d*x))/(a**2*d) + (5*A + 3*I*B)*cot(c + d*x)**2/(4*a**2*d*(I*tan(c + d*x) + 1)) + (15*I*A - 9*B)*cot(c + d*x)/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_51():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**3
    F = -(I*A - 3*B)*tan(c + d*x)**2/(2*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - B)*tan(c + d*x)**4/(6*d*(I*a*tan(c + d*x) + a)**3) + (5*A + 11*I*B)*tan(c + d*x)**3/(24*a*d*(I*a*tan(c + d*x) + a)**2) - x*(7*A + 25*I*B)/(8*a**3) + (7*A + 25*I*B)*tan(c + d*x)/(8*a**3*d) - (I*A - 3*B)*log(cos(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_52():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**3
    F = -I*B*log(cos(c + d*x))/(a**3*d) + (A + 7*I*B)/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - B)*tan(c + d*x)**3/(6*d*(I*a*tan(c + d*x) + a)**3) + (A + 3*I*B)*tan(c + d*x)**2/(8*a*d*(I*a*tan(c + d*x) + a)**2) + x*(I*A - 7*B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_53():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = (I*A - B)*tan(c + d*x)**2/(6*d*(I*a*tan(c + d*x) + a)**3) + (I*A + 17*B)/(24*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - 7*B)/(24*a*d*(I*a*tan(c + d*x) + a)**2) - x*(A - I*B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_54():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = (A - I*B)/(8*d*(I*a**3*tan(c + d*x) + a**3)) - (A + I*B)/(6*d*(I*a*tan(c + d*x) + a)**3) + (A + 3*I*B)/(8*a*d*(I*a*tan(c + d*x) + a)**2) - x*(I*A + B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_55():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**3
    F = (I*A - B)/(6*d*(I*a*tan(c + d*x) + a)**3) + (I*A + B)/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A + B)/(8*a*d*(I*a*tan(c + d*x) + a)**2) + x*(A - I*B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_56():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = A*log(sin(c + d*x))/(a**3*d) + (A + I*B)/(6*d*(I*a*tan(c + d*x) + a)**3) + (7*A + I*B)/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (3*A + I*B)/(8*a*d*(I*a*tan(c + d*x) + a)**2) - x*(7*I*A - B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_57():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = (A + I*B)*cot(c + d*x)/(6*d*(I*a*tan(c + d*x) + a)**3) + (3*A + I*B)*cot(c + d*x)/(2*d*(I*a**3*tan(c + d*x) + a**3)) + (11*A + 5*I*B)*cot(c + d*x)/(24*a*d*(I*a*tan(c + d*x) + a)**2) - x*(25*A + 7*I*B)/(8*a**3) - (25*A + 7*I*B)*cot(c + d*x)/(8*a**3*d) - (3*I*A - B)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_58():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**3
    F = (A + I*B)*cot(c + d*x)**2/(6*d*(I*a*tan(c + d*x) + a)**3) + (55*A + 25*I*B)*cot(c + d*x)**2/(24*d*(I*a**3*tan(c + d*x) + a**3)) + (13*A + 7*I*B)*cot(c + d*x)**2/(24*a*d*(I*a*tan(c + d*x) + a)**2) + x*(55*I*A - 25*B)/(8*a**3) - (7*A + 3*I*B)*log(sin(c + d*x))/(a**3*d) - (7*A + 3*I*B)*cot(c + d*x)**2/(2*a**3*d) + (55*I*A - 25*B)*cot(c + d*x)/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_59():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**4
    F = -B*log(cos(c + d*x))/(a**4*d) + (I*A - B)*tan(c + d*x)**4/(8*d*(I*a*tan(c + d*x) + a)**4) + (A + 3*I*B)*tan(c + d*x)**3/(12*a*d*(I*a*tan(c + d*x) + a)**3) + x*(A + 15*I*B)/(16*a**4) - (I*A - 15*B)/(16*a**4*d*(I*tan(c + d*x) + 1)) - (I*A - 7*B)*tan(c + d*x)**2/(16*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_60():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**4
    F = (I*A - B)*tan(c + d*x)**3/(8*d*(I*a*tan(c + d*x) + a)**4) + (A + 5*I*B)*tan(c + d*x)**2/(24*a*d*(I*a*tan(c + d*x) + a)**3) + x*(I*A + B)/(16*a**4) - (A - 13*I*B)/(48*a**4*d*(I*tan(c + d*x) + 1)**2) + (5*A - 29*I*B)/(48*a**4*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_61():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = -B/(6*a*d*(I*a*tan(c + d*x) + a)**3) + (I*A - B)*tan(c + d*x)**2/(8*d*(I*a*tan(c + d*x) + a)**4) - x*(A - I*B)/(16*a**4) - (I*A + B)/(16*a**4*d*(I*tan(c + d*x) + 1)) + (I*A + 5*B)/(16*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_62():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = (A - I*B)/(16*d*(I*a**4*tan(c + d*x) + a**4)) + (A - I*B)/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) - (A + I*B)/(8*d*(I*a*tan(c + d*x) + a)**4) + (A + 3*I*B)/(12*a*d*(I*a*tan(c + d*x) + a)**3) - x*(I*A + B)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_63():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**4
    F = (I*A - B)/(8*d*(I*a*tan(c + d*x) + a)**4) + (I*A + B)/(16*d*(I*a**4*tan(c + d*x) + a**4)) + (I*A + B)/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + (I*A + B)/(12*a*d*(I*a*tan(c + d*x) + a)**3) + x*(A - I*B)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_64():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = A*log(sin(c + d*x))/(a**4*d) + (A + I*B)/(8*d*(I*a*tan(c + d*x) + a)**4) + (3*A + I*B)/(12*a*d*(I*a*tan(c + d*x) + a)**3) - x*(15*I*A - B)/(16*a**4) + (7*A + I*B)/(16*a**4*d*(I*tan(c + d*x) + 1)**2) + (15*A + I*B)/(16*a**4*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_65():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = (A + I*B)*cot(c + d*x)/(8*d*(I*a*tan(c + d*x) + a)**4) + (7*A + 3*I*B)*cot(c + d*x)/(24*a*d*(I*a*tan(c + d*x) + a)**3) - x*(65*A + 15*I*B)/(16*a**4) + (4*A + I*B)*cot(c + d*x)/(2*a**4*d*(I*tan(c + d*x) + 1)) + (31*A + 9*I*B)*cot(c + d*x)/(48*a**4*d*(I*tan(c + d*x) + 1)**2) - (65*A + 15*I*B)*cot(c + d*x)/(16*a**4*d) - (4*I*A - B)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_66():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**4
    F = (A + I*B)*cot(c + d*x)**2/(8*d*(I*a*tan(c + d*x) + a)**4) + (2*A + I*B)*cot(c + d*x)**2/(6*a*d*(I*a*tan(c + d*x) + a)**3) + x*(175*I*A - 65*B)/(16*a**4) - (11*A + 4*I*B)*log(sin(c + d*x))/(a**4*d) - (11*A + 4*I*B)*cot(c + d*x)**2/(2*a**4*d) + (43*A + 17*I*B)*cot(c + d*x)**2/(48*a**4*d*(I*tan(c + d*x) + 1)**2) + (175*A + 65*I*B)*cot(c + d*x)**2/(48*a**4*d*(I*tan(c + d*x) + 1)) + (175*I*A - 65*B)*cot(c + d*x)/(16*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_67():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3
    F = 2*B*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(7*d) + sqrt(2)*sqrt(a)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + (14*A - 2*I*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(35*d) - (56*A - 8*I*B)*sqrt(I*a*tan(c + d*x) + a)/(35*d) - (14*A - 62*I*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_68():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2
    F = 2*B*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(5*d) - 8*B*sqrt(I*a*tan(c + d*x) + a)/(5*d) + sqrt(2)*sqrt(a)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - (10*I*A + 2*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_69():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)
    F = 2*A*sqrt(I*a*tan(c + d*x) + a)/d - 2*I*B*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d) - sqrt(2)*sqrt(a)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_70():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)
    F = 2*B*sqrt(I*a*tan(c + d*x) + a)/d - sqrt(2)*sqrt(a)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_71():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)
    F = -2*A*sqrt(a)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + sqrt(2)*sqrt(a)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_72():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2
    F = -A*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/d + sqrt(2)*sqrt(a)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - sqrt(a)*(I*A + 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_73():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3
    F = -A*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*d) - sqrt(2)*sqrt(a)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + sqrt(a)*(7*A - 4*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - (I*A + 4*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_74():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**4
    F = -A*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3/(3*d) - sqrt(2)*sqrt(a)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + sqrt(a)*(9*I*A + 14*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(8*d) + (7*A - 2*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(8*d) - (I*A + 6*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_75():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**2
    F = 2*I*B*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(7*d) + 2*sqrt(2)*a**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*a*(7*I*A + 8*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(35*d) - 8*a*(7*I*A + 8*B)*sqrt(I*a*tan(c + d*x) + a)/(35*d) - (84*I*A + 76*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_76():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)
    F = 2*A*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*I*B*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d) - 2*sqrt(2)*a**(sympy.S(3)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*a*(A - I*B)*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_77():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*sqrt(2)*a**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*a*(I*A + B)*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_78():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)
    F = -2*A*a**(sympy.S(3)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 2*I*B*a*sqrt(I*a*tan(c + d*x) + a)/d + 2*sqrt(2)*a**(sympy.S(3)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_79():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**2
    F = -A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/d + 2*sqrt(2)*a**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**(sympy.S(3)/2)*(3*I*A + 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_80():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3
    F = -A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*d) - 2*sqrt(2)*a**(sympy.S(3)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + a**(sympy.S(3)/2)*(11*A - 12*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - a*(5*I*A + 4*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_81():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**4
    F = -A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3/(3*d) - 2*sqrt(2)*a**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + a**(sympy.S(3)/2)*(23*I*A + 22*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(8*d) + a*(9*A - 10*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(8*d) - a*(7*I*A + 6*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_82():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**2
    F = 2*I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**3/(9*d) + 4*sqrt(2)*a**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*a**2*(3*A - 4*I*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(21*d) + 2*a**2*(45*I*A + 46*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(105*d) - 8*a**2*(45*I*A + 46*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d) - 8*a*(60*I*A + 59*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(315*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_83():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)
    F = 2*A*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d) - 2*I*B*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d) - 4*sqrt(2)*a**(sympy.S(5)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 4*a**2*(A - I*B)*sqrt(I*a*tan(c + d*x) + a)/d + 2*a*(A - I*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_84():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d) - 4*sqrt(2)*a**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 4*a**2*(I*A + B)*sqrt(I*a*tan(c + d*x) + a)/d + 2*a*(I*A + B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_85():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)
    F = -2*A*a**(sympy.S(5)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 2*I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 4*sqrt(2)*a**(sympy.S(5)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*a**2*(A - 2*I*B)*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_86():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**2
    F = -A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)/d + 4*sqrt(2)*a**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**(sympy.S(5)/2)*(5*I*A + 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + a**2*(I*A - 2*B)*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_87():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**3
    F = -A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**2/(2*d) - 4*sqrt(2)*a**(sympy.S(5)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + a**(sympy.S(5)/2)*(23*A - 20*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - a**2*(7*I*A + 4*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_88():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**4
    F = -A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(3*d) - 4*sqrt(2)*a**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + a**(sympy.S(5)/2)*(45*I*A + 46*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(8*d) + a**2*(19*A - 18*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(8*d) - a**2*(3*I*A + 2*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_89():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5
    F = -A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**4/(4*d) + 4*sqrt(2)*a**(sympy.S(5)/2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 3*a**(sympy.S(5)/2)*(121*A - 120*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(64*d) + a**2*(107*A - 104*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(96*d) - a**2*(11*I*A + 8*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3/(24*d) + a**2*(149*I*A + 152*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(64*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_90():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**3/(d*sqrt(I*a*tan(c + d*x) + a)) - (5*A + 7*I*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(5*a*d) + (20*A + 28*I*B)*sqrt(I*a*tan(c + d*x) + a)/(5*a*d) - (25*A + 23*I*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(15*a**2*d) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_91():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**2/(d*sqrt(I*a*tan(c + d*x) + a)) - (4*I*A - 4*B)*sqrt(I*a*tan(c + d*x) + a)/(a*d) + (3*I*A - 5*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_92():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = -2*I*B*sqrt(I*a*tan(c + d*x) + a)/(a*d) - (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_93():
    f = (A + B*tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = (I*A - B)/(d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_94():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = -2*A*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) + (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_95():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)/(d*sqrt(I*a*tan(c + d*x) + a)) - (2*A + I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(a*d) + (I*A - 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_96():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**2/(d*sqrt(I*a*tan(c + d*x) + a)) - (3*A + 2*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*a*d) + (7*I*A - 8*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a*d) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d) + (11*A + 4*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_97():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*tan(c + d*x)**3/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (3*A + 5*I*B)*tan(c + d*x)**2/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - (6*A + 10*I*B)*sqrt(I*a*tan(c + d*x) + a)/(a**2*d) + (11*A + 21*I*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(6*a**3*d) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_98():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*tan(c + d*x)**2/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (5*I*A - 11*B)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) + (I*A - 7*B)*sqrt(I*a*tan(c + d*x) + a)/(3*a**2*d) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_99():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A + I*B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (A + 3*I*B)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_100():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (I*A + B)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_101():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*A*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) + (A + I*B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (3*A + I*B)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_102():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (A + I*B)*cot(c + d*x)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (13*A + 7*I*B)*cot(c + d*x)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - (7*A + 3*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(2*a**2*d) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d) + (3*I*A - 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_103():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (A + I*B)*cot(c + d*x)**2/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (17*A + 11*I*B)*cot(c + d*x)**2/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - (22*A + 13*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(6*a**2*d) + (21*I*A - 14*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a**2*d) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d) + (23*A + 12*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_104():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*tan(c + d*x)**4/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (11*A + 21*I*B)*tan(c + d*x)**3/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - (39*I*A - 89*B)*tan(c + d*x)**2/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (39*I*A - 89*B)*sqrt(I*a*tan(c + d*x) + a)/(5*a**3*d) - (151*I*A - 361*B)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(60*a**4*d) - sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_105():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*tan(c + d*x)**3/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (7*A + 17*I*B)*tan(c + d*x)**2/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (41*A + 151*I*B)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (13*A + 83*I*B)*sqrt(I*a*tan(c + d*x) + a)/(30*a**3*d) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_106():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*tan(c + d*x)**2/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (3*I*A - 13*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - (I*A - 31*B)/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_107():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A + I*B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (A + 3*I*B)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (A - I*B)/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_108():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (I*A + B)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (I*A + B)/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_109():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*A*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) + (A + I*B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (3*A + I*B)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (7*A + I*B)/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_110():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (A + I*B)*cot(c + d*x)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (19*A + 9*I*B)*cot(c + d*x)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (41*A + 15*I*B)*cot(c + d*x)/(12*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - (21*A + 7*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a**3*d) + sqrt(2)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d) + (5*I*A - 2*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_111():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (A + I*B)*cot(c + d*x)**2/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (23*A + 13*I*B)*cot(c + d*x)**2/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (337*A + 167*I*B)*cot(c + d*x)**2/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - (85*A + 41*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(12*a**3*d) + (42*I*A - 21*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a**3*d) - sqrt(2)*(A - I*B)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d) + (43*A + 20*I*B)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_112():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*I*B*a*tan(c + d*x)**(sympy.S(7)/2)/(7*d) + 2*a*(A - I*B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(I*A + B)*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - 2*a*(I*A + B)*sqrt(tan(c + d*x))/d - 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_113():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*I*B*a*tan(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*(A - I*B)*sqrt(tan(c + d*x))/d + 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d + 2*a*(I*A + B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_114():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))
    F = 2*I*B*a*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(I*A + B)*sqrt(tan(c + d*x))/d + 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_115():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/sqrt(tan(c + d*x))
    F = 2*I*B*a*sqrt(tan(c + d*x))/d - 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_116():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a/(d*sqrt(tan(c + d*x))) - 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_117():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a*(I*A + B)/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_118():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + 2*a*(A - I*B)/(d*sqrt(tan(c + d*x))) + 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a*(I*A + B)/(3*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_119():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*I*B*(I*a**2*tan(c + d*x) + a**2)*tan(c + d*x)**(sympy.S(7)/2)/(9*d) + 4*a**2*(A - I*B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*a**2*(9*A - 11*I*B)*tan(c + d*x)**(sympy.S(7)/2)/(63*d) + 4*a**2*(I*A + B)*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - 4*a**2*(I*A + B)*sqrt(tan(c + d*x))/d - 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_120():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*I*B*(I*a**2*tan(c + d*x) + a**2)*tan(c + d*x)**(sympy.S(5)/2)/(7*d) + 4*a**2*(A - I*B)*sqrt(tan(c + d*x))/d + 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a**2*(7*A - 9*I*B)*tan(c + d*x)**(sympy.S(5)/2)/(35*d) + 4*a**2*(I*A + B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_121():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*sqrt(tan(c + d*x))
    F = 2*I*B*(I*a**2*tan(c + d*x) + a**2)*tan(c + d*x)**(sympy.S(3)/2)/(5*d) - 2*a**2*(5*A - 7*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(15*d) + 4*a**2*(I*A + B)*sqrt(tan(c + d*x))/d + 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_122():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/sqrt(tan(c + d*x))
    F = 2*I*B*(I*a**2*tan(c + d*x) + a**2)*sqrt(tan(c + d*x))/(3*d) - 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a**2*(3*A - 5*I*B)*sqrt(tan(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_123():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*(I*a**2*tan(c + d*x) + a**2)/(d*sqrt(tan(c + d*x))) + 2*a**2*(I*A - B)*sqrt(tan(c + d*x))/d - 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_124():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*(I*a**2*tan(c + d*x) + a**2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a**2*(5*I*A + 3*B)/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_125():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*(I*a**2*tan(c + d*x) + a**2)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(A - I*B)/(d*sqrt(tan(c + d*x))) + 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 2*a**2*(7*I*A + 5*B)/(15*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_126():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*(I*a**2*tan(c + d*x) + a**2)/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d + 4*a**2*(A - I*B)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(I*A + B)/(d*sqrt(tan(c + d*x))) - 2*a**2*(9*I*A + 7*B)/(35*d*tan(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_127():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*I*B*a*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(5)/2)/(9*d) + 8*a**3*(A - I*B)*sqrt(tan(c + d*x))/d + 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 16*a**3*(18*A - 19*I*B)*tan(c + d*x)**(sympy.S(5)/2)/(315*d) + 8*a**3*(I*A + B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - (18*A - 26*I*B)*(I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**(sympy.S(5)/2)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_128():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*sqrt(tan(c + d*x))
    F = 2*I*B*a*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(3)/2)/(7*d) - 8*a**3*(21*A - 23*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(105*d) + 8*a**3*(I*A + B)*sqrt(tan(c + d*x))/d + 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - (14*A - 22*I*B)*(I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_129():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/sqrt(tan(c + d*x))
    F = 2*I*B*a*(I*a*tan(c + d*x) + a)**2*sqrt(tan(c + d*x))/(5*d) - 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - 16*a**3*(5*A - 6*I*B)*sqrt(tan(c + d*x))/(15*d) - (10*A - 18*I*B)*(I*a**3*tan(c + d*x) + a**3)*sqrt(tan(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_130():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**2/(d*sqrt(tan(c + d*x))) - 16*B*a**3*sqrt(tan(c + d*x))/(3*d) - 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d + (6*I*A - 2*B)*(I*a**3*tan(c + d*x) + a**3)*sqrt(tan(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_131():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/tan(c + d*x)**(sympy.S(5)/2)
    F = -16*A*a**3*sqrt(tan(c + d*x))/(3*d) - 2*A*a*(I*a*tan(c + d*x) + a)**2/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - (14*I*A + 6*B)*(I*a**3*tan(c + d*x) + a**3)/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_132():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**2/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + 16*a**3*(6*A - 5*I*B)/(15*d*sqrt(tan(c + d*x))) + 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d - (18*I*A + 10*B)*(I*a**3*tan(c + d*x) + a**3)/(15*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_133():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**2/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atanh((-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x)))/d + 8*a**3*(23*A - 21*I*B)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) + 8*a**3*(I*A + B)/(d*sqrt(tan(c + d*x))) - (22*I*A + 14*B)*(I*a**3*tan(c + d*x) + a**3)/(35*d*tan(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_134():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**(sympy.S(5)/2)/(2*d*(I*a*tan(c + d*x) + a)) - (3*A + 7*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(6*a*d) - (5*I*A - 5*B)*sqrt(tan(c + d*x))/(2*a*d) - sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(1 + 4*I) - B*(6 + I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d) - sqrt(2)*(A*(3 - 5*I) + B*(5 + 7*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d) + sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(4 + I) + B*(1 + 6*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(4 + I) + B*(1 + 6*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_135():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(2*d*(I*a*tan(c + d*x) + a)) - (A + 5*I*B)*sqrt(tan(c + d*x))/(2*a*d) + sqrt(2)*(A*(1 - 3*I) + B*(3 + 5*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(8*a*d) - sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(1 + 2*I) - B*(4 + I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(2 + I) + B*(1 + 4*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(2 + I) + B*(1 + 4*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_136():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)
    F = (I*A - B)*sqrt(tan(c + d*x))/(2*d*(I*a*tan(c + d*x) + a)) + sqrt(2)*(sympy.S(1)/4 - I/4)*(A + B*(2 - I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/4 - I/4)*(A + B*(2 - I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/8 + I/8)*(A - B*(2 + I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/8 + I/8)*(A - B*(2 + I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_137():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x)))
    F = (A + I*B)*sqrt(tan(c + d*x))/(2*d*(I*a*tan(c + d*x) + a)) + sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(2 + I) + B)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(2 + I) + B)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d) - sqrt(2)*(A*(3 + I) - B*(1 + I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d) + sqrt(2)*(A*(3 + I) - B*(1 + I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_138():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(2*d*(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - (5*A + I*B)/(2*a*d*sqrt(tan(c + d*x))) + sqrt(2)*(A*(-5 - 3*I) + B*(3 - I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(8*a*d) - sqrt(2)*(sympy.S(1)/8 - I/8)*(A*(4 + I) + B*(1 + 2*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a*d) + sqrt(2)*(A*(5 - 3*I) + B*(3 + I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d) - sqrt(2)*(A*(5 + 3*I) - B*(3 - I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_139():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(2*d*(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) - (7*A + 3*I*B)/(6*a*d*tan(c + d*x)**(sympy.S(3)/2)) + (5*I*A - 5*B)/(2*a*d*sqrt(tan(c + d*x))) + sqrt(2)*(A*(-7 - 5*I) + B*(5 - 3*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d) - sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(6 + I) + B*(1 + 4*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d) - sqrt(2)*(A*(7 - 5*I) + B*(5 + 3*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(8*a*d) + sqrt(2)*(A*(7 + 5*I) - B*(5 - 3*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_140():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**2
    F = (I*A - B)*tan(c + d*x)**(sympy.S(5)/2)/(4*d*(I*a*tan(c + d*x) + a)**2) + (3*A + 7*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(8*a**2*d*(I*tan(c + d*x) + 1)) + (5*I*A - 25*B)*sqrt(tan(c + d*x))/(8*a**2*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(7 + 2*I) + B*(2 + 23*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(7 + 2*I) + B*(2 + 23*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(A*(9 + 5*I) - B*(25 - 21*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**2*d) - sqrt(2)*(A*(9 + 5*I) - B*(25 - 21*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_141():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**2
    F = (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(4*d*(I*a*tan(c + d*x) + a)**2) + (A + 5*I*B)*sqrt(tan(c + d*x))/(8*a**2*d*(I*tan(c + d*x) + 1)) + sqrt(2)*(A*(1 - 3*I) - B*(9 - 5*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d) - sqrt(2)*(A*(1 - 3*I) - B*(9 - 5*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d) - sqrt(2)*(A*(1 + 3*I) + B*(9 + 5*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**2*d) - sqrt(2)*(A*(1 + 3*I) + B*(9 + 5*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_142():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**2
    F = (I*A - B)*sqrt(tan(c + d*x))/(4*d*(I*a*tan(c + d*x) + a)**2) + (I*A + 3*B)*sqrt(tan(c + d*x))/(8*a**2*d*(I*tan(c + d*x) + 1)) - sqrt(2)*(A*(-1 + 3*I) + B*(1 + 3*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**2*d) - sqrt(2)*(A*(-1 + 3*I) + B*(1 + 3*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**2*d) + sqrt(2)*(A*(1 + 3*I) + B*(1 - 3*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d) - sqrt(2)*(A*(1 + 3*I) + B*(1 - 3*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_143():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*sqrt(tan(c + d*x)))
    F = (A + I*B)*sqrt(tan(c + d*x))/(4*d*(I*a*tan(c + d*x) + a)**2) + (5*A + I*B)*sqrt(tan(c + d*x))/(8*a**2*d*(I*tan(c + d*x) + 1)) + sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(-7 + 2*I) + B*(2 + I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(-2 + 7*I) + B*(1 + 2*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d) + sqrt(2)*(A*(9 - 5*I) + B*(1 - 3*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**2*d) + sqrt(2)*(A*(9 + 5*I) - B*(1 + 3*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_144():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(4*d*(I*a*tan(c + d*x) + a)**2*sqrt(tan(c + d*x))) + (7*A + 3*I*B)/(8*a**2*d*(I*tan(c + d*x) + 1)*sqrt(tan(c + d*x))) - (25*A + 5*I*B)/(8*a**2*d*sqrt(tan(c + d*x))) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(2 + 23*I) - B*(7 + 2*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(23 + 2*I) + B*(2 + 7*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(23 + 2*I) + B*(2 + 7*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(A*(25 + 21*I) - B*(9 - 5*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_145():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(4*d*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(3)/2)) + (9*A + 5*I*B)/(8*a**2*d*(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)) - (49*A + 21*I*B)/(24*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + (45*I*A - 25*B)/(8*a**2*d*sqrt(tan(c + d*x))) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(2 + 47*I) - B*(23 + 2*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(47 + 2*I) + B*(2 + 23*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(47 + 2*I) + B*(2 + 23*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d) + sqrt(2)*(A*(49 + 45*I) - B*(25 - 21*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_146():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**3
    F = (I*A - B)*tan(c + d*x)**(sympy.S(9)/2)/(6*d*(I*a*tan(c + d*x) + a)**3) - (6*I*A - 15*B)*tan(c + d*x)**(sympy.S(5)/2)/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (A + 2*I*B)*tan(c + d*x)**(sympy.S(7)/2)/(4*a*d*(I*a*tan(c + d*x) + a)**2) + (28*A + 77*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(24*a**3*d) + (30*I*A - 75*B)*sqrt(tan(c + d*x))/(8*a**3*d) - sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(1 + 29*I) - B*(76 + I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(A*(28 - 30*I) + B*(75 + 77*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(29 + I) + B*(1 + 76*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(29 + I) + B*(1 + 76*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_147():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**3
    F = (I*A - B)*tan(c + d*x)**(sympy.S(7)/2)/(6*d*(I*a*tan(c + d*x) + a)**3) - (7*I*A - 28*B)*tan(c + d*x)**(sympy.S(3)/2)/(24*d*(I*a**3*tan(c + d*x) + a**3)) + (2*A + 5*I*B)*tan(c + d*x)**(sympy.S(5)/2)/(12*a*d*(I*a*tan(c + d*x) + a)**2) + (5*A + 30*I*B)*sqrt(tan(c + d*x))/(8*a**3*d) + sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + 6*I) - B*(29 + I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(A*(5 - 7*I) + B*(28 + 30*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**3*d) + sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(6 + I) + B*(1 + 29*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(6 + I) + B*(1 + 29*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_148():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**3
    F = 5*B*sqrt(tan(c + d*x))/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - B)*tan(c + d*x)**(sympy.S(5)/2)/(6*d*(I*a*tan(c + d*x) + a)**3) + (A + 4*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(12*a*d*(I*a*tan(c + d*x) + a)**2) - sqrt(2)*(2*A + B*(5 - 7*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**3*d) - sqrt(2)*(2*A + B*(5 - 7*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**3*d) - sqrt(2)*(2*A - B*(5 + 7*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(2*A - B*(5 + 7*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_149():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**3
    F = I*B*sqrt(tan(c + d*x))/(4*a*d*(I*a*tan(c + d*x) + a)**2) + (A - 2*I*B)*sqrt(tan(c + d*x))/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(6*d*(I*a*tan(c + d*x) + a)**3) - sqrt(2)*(A*(-1 + I) + 2*B)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(A*(-1 + I) + 2*B)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) - sqrt(2)*(A*(1 + I) + 2*B)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**3*d) - sqrt(2)*(A*(1 + I) + 2*B)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_150():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**3
    F = B*sqrt(tan(c + d*x))/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (I*A - B)*sqrt(tan(c + d*x))/(6*d*(I*a*tan(c + d*x) + a)**3) + (I*A + 2*B)*sqrt(tan(c + d*x))/(12*a*d*(I*a*tan(c + d*x) + a)**2) + sqrt(2)*(2*I*A + B*(1 - I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) - sqrt(2)*(2*I*A + B*(1 - I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + I) + B)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + I) + B)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_151():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*sqrt(tan(c + d*x)))
    F = 5*A*sqrt(tan(c + d*x))/(8*d*(I*a**3*tan(c + d*x) + a**3)) + (A + I*B)*sqrt(tan(c + d*x))/(6*d*(I*a*tan(c + d*x) + a)**3) + (4*A + I*B)*sqrt(tan(c + d*x))/(12*a*d*(I*a*tan(c + d*x) + a)**2) + sqrt(2)*(A*(7 - 5*I) - 2*I*B)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**3*d) + sqrt(2)*(A*(7 - 5*I) - 2*I*B)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(32*a**3*d) - sqrt(2)*(A*(7 + 5*I) - 2*I*B)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(A*(7 + 5*I) - 2*I*B)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_152():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(6*d*(I*a*tan(c + d*x) + a)**3*sqrt(tan(c + d*x))) + (28*A + 7*I*B)/(24*d*(I*a**3*tan(c + d*x) + a**3)*sqrt(tan(c + d*x))) + (5*A + 2*I*B)/(12*a*d*(I*a*tan(c + d*x) + a)**2*sqrt(tan(c + d*x))) - (30*A + 5*I*B)/(8*a**3*d*sqrt(tan(c + d*x))) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(1 + 29*I) - B*(6 + I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(29 + I) + B*(1 + 6*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(29 + I) + B*(1 + 6*I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(A*(30 + 28*I) - B*(7 - 5*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_153():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(6*d*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**(sympy.S(3)/2)) + (15*A + 6*I*B)/(8*d*(I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**(sympy.S(3)/2)) + (2*A + I*B)/(4*a*d*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(3)/2)) - (77*A + 28*I*B)/(24*a**3*d*tan(c + d*x)**(sympy.S(3)/2)) + (75*I*A - 30*B)/(8*a**3*d*sqrt(tan(c + d*x))) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(1 + 76*I) - B*(29 + I))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(76 + I) + B*(1 + 29*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(76 + I) + B*(1 + 29*I))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**3*d) + sqrt(2)*(A*(77 + 75*I) - B*(30 - 28*I))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(64*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_154():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)
    F = B*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(2*d) + sqrt(a)*(1 + I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (-1)**(sympy.S(3)/4)*sqrt(a)*(4*I*A + 7*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) + (4*A - I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_155():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))
    F = B*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/d - sqrt(a)*(1 + I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*sqrt(a)*(2*A - I*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_156():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/sqrt(tan(c + d*x))
    F = -2*(-1)**(sympy.S(3)/4)*B*sqrt(a)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - sqrt(a)*(1 + I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_157():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) + sqrt(a)*(1 + I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_158():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a)*(1 + I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (2*I*A + 6*B)*sqrt(I*a*tan(c + d*x) + a)/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_159():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a)*(-1 - I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (26*A - 10*I*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*sqrt(tan(c + d*x))) - (2*I*A + 10*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_160():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)/(7*d*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(a)*(1 - I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (62*A - 14*I*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) - (2*I*A + 14*B)*sqrt(I*a*tan(c + d*x) + a)/(35*d*tan(c + d*x)**(sympy.S(5)/2)) + (86*I*A + 182*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_161():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = I*B*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)/(3*d) + a**(sympy.S(3)/2)*(2 + 2*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*(22*I*A + 23*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(8*d) + a*(10*A - 9*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(8*d) + a*(6*I*A + 7*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_162():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))
    F = I*B*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(2*d) - a**(sympy.S(3)/2)*(2 + 2*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*(12*A - 11*I*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) + a*(4*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_163():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/sqrt(tan(c + d*x))
    F = I*B*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/d - a**(sympy.S(3)/2)*(2 + 2*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*(2*I*A + 3*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_164():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) + 2*(-1)**(sympy.S(1)/4)*B*a**(sympy.S(3)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(3)/2)*(2 + 2*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_165():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(3)/2)*(2 + 2*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a*(4*I*A + 3*B)*sqrt(I*a*tan(c + d*x) + a)/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_166():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(3)/2)*(-2 - 2*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a*(9*A - 10*I*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*sqrt(tan(c + d*x))) - 2*a*(6*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_167():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)/(7*d*tan(c + d*x)**(sympy.S(7)/2)) + a**(sympy.S(3)/2)*(2 - 2*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a*(19*A - 21*I*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a*(8*I*A + 7*B)*sqrt(I*a*tan(c + d*x) + a)/(35*d*tan(c + d*x)**(sympy.S(5)/2)) + 4*a*(67*I*A + 63*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_168():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)/(9*d*tan(c + d*x)**(sympy.S(9)/2)) + a**(sympy.S(3)/2)*(2 + 2*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a*(11*A - 12*I*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(5)/2)) - 4*a*(193*A - 201*I*B)*sqrt(I*a*tan(c + d*x) + a)/(315*d*sqrt(tan(c + d*x))) - 2*a*(10*I*A + 9*B)*sqrt(I*a*tan(c + d*x) + a)/(63*d*tan(c + d*x)**(sympy.S(7)/2)) + 4*a*(61*I*A + 57*B)*sqrt(I*a*tan(c + d*x) + a)/(315*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_169():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2)/(4*d) + a**(sympy.S(5)/2)*(4 + 4*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 3*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(120*I*A + 121*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(64*d) - a**2*(8*A - 11*I*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)/(24*d) + a**2*(152*A - 149*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(64*d) + a**2*(104*I*A + 107*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(96*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_170():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x))
    F = I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(46*A - 45*I*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(8*d) - a**2*(2*A - 3*I*B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(4*d) + a**2*(18*I*A + 19*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_171():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/sqrt(tan(c + d*x))
    F = I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(2*d) + a**(sympy.S(5)/2)*(4 - 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(20*I*A + 23*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) - a**2*(4*A - 7*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_172():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(d*sqrt(tan(c + d*x))) + a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(2*A - 5*I*B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**2*(2*I*A - B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_173():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*(-1)**(sympy.S(3)/4)*B*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 + 4*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*(2*I*A + B)*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_174():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(5)/2)*(-4 - 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(38*A - 35*I*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*sqrt(tan(c + d*x))) - 2*a**2*(8*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)/(15*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_175():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(7*d*tan(c + d*x)**(sympy.S(7)/2)) + a**(sympy.S(5)/2)*(4 - 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(80*A - 77*I*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(10*I*A + 7*B)*sqrt(I*a*tan(c + d*x) + a)/(35*d*tan(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(130*I*A + 133*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_176():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(9*d*tan(c + d*x)**(sympy.S(9)/2)) + a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(46*A - 45*I*B)*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(5)/2)) - 8*a**2*(197*A - 195*I*B)*sqrt(I*a*tan(c + d*x) + a)/(315*d*sqrt(tan(c + d*x))) - 2*a**2*(4*I*A + 3*B)*sqrt(I*a*tan(c + d*x) + a)/(21*d*tan(c + d*x)**(sympy.S(7)/2)) + 8*a**2*(59*I*A + 60*B)*sqrt(I*a*tan(c + d*x) + a)/(315*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_177():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(13)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(11*d*tan(c + d*x)**(sympy.S(11)/2)) + a**(sympy.S(5)/2)*(4 + 4*I)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(212*A - 209*I*B)*sqrt(I*a*tan(c + d*x) + a)/(693*d*tan(c + d*x)**(sympy.S(7)/2)) - 8*a**2*(655*A - 649*I*B)*sqrt(I*a*tan(c + d*x) + a)/(3465*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(14*I*A + 11*B)*sqrt(I*a*tan(c + d*x) + a)/(99*d*tan(c + d*x)**(sympy.S(9)/2)) + 4*a**2*(250*I*A + 253*B)*sqrt(I*a*tan(c + d*x) + a)/(1155*d*tan(c + d*x)**(sympy.S(5)/2)) - 8*a**2*(2155*I*A + 2167*B)*sqrt(I*a*tan(c + d*x) + a)/(3465*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_178():
    f = (B*tan(c + d*x) + 3*B*b/(2*a))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = 2*(-1)**(sympy.S(3)/4)*B*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + B*a**(sympy.S(3)/2)*(2 + 2*I)*(2*a + 3*I*b)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*B*a*(a + 3*I*b)*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) - B*b*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_179():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) - (A + 2*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(a*d) - (sympy.S.Half - I/2)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) + (-1)**(sympy.S(3)/4)*(2*I*A - B)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_180():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = -2*(-1)**(sympy.S(1)/4)*B*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) + (I*A - B)*sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S.Half + I/2)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_181():
    f = (A + B*tan(c + d*x))/(sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x)))
    F = (A + I*B)*sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S.Half - I/2)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_182():
    f = (A + B*tan(c + d*x))/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - (3*A + I*B)*sqrt(I*a*tan(c + d*x) + a)/(a*d*sqrt(tan(c + d*x))) + (sympy.S.Half + I/2)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_183():
    f = (A + B*tan(c + d*x))/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) - (5*A + 3*I*B)*sqrt(I*a*tan(c + d*x) + a)/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) + (7*I*A - 9*B)*sqrt(I*a*tan(c + d*x) + a)/(3*a*d*sqrt(tan(c + d*x))) + (sympy.S.Half + I/2)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_184():
    f = (A + B*tan(c + d*x))/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/2))
    F = (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)) - (7*A + 5*I*B)*sqrt(I*a*tan(c + d*x) + a)/(5*a*d*tan(c + d*x)**(sympy.S(5)/2)) + (61*A + 35*I*B)*sqrt(I*a*tan(c + d*x) + a)/(15*a*d*sqrt(tan(c + d*x))) + (23*I*A - 25*B)*sqrt(I*a*tan(c + d*x) + a)/(15*a*d*tan(c + d*x)**(sympy.S(3)/2)) + (sympy.S(-1)/2 - I/2)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_185():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*(-1)**(sympy.S(3)/4)*B*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (A + 3*I*B)*sqrt(tan(c + d*x))/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S(1)/4 - I/4)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_186():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*sqrt(tan(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (I*A + 5*B)*sqrt(tan(c + d*x))/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S(1)/4 + I/4)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_187():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x)))
    F = (A + I*B)*sqrt(tan(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (7*A + I*B)*sqrt(tan(c + d*x))/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/4 - I/4)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_188():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + (11*A + 5*I*B)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - (25*A + 7*I*B)*sqrt(I*a*tan(c + d*x) + a)/(6*a**2*d*sqrt(tan(c + d*x))) + (sympy.S(1)/4 + I/4)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_189():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + (5*A + 3*I*B)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) - (21*A + 11*I*B)*sqrt(I*a*tan(c + d*x) + a)/(6*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + (39*I*A - 25*B)*sqrt(I*a*tan(c + d*x) + a)/(6*a**2*d*sqrt(tan(c + d*x))) + (sympy.S(1)/4 + I/4)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_190():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*(-1)**(sympy.S(1)/4)*B*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (I*A - B)*tan(c + d*x)**(sympy.S(5)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (A + 3*I*B)*tan(c + d*x)**(sympy.S(3)/2)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - (I*A - 7*B)*sqrt(tan(c + d*x))/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/8 + I/8)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_191():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*tan(c + d*x)**(sympy.S(3)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (A + 11*I*B)*sqrt(tan(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (13*A - 37*I*B)*sqrt(tan(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S(1)/8 - I/8)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_192():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*sqrt(tan(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (3*I*A + 7*B)*sqrt(tan(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - (3*I*A - 13*B)*sqrt(tan(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S(1)/8 + I/8)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_193():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x)))
    F = (A + I*B)*sqrt(tan(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (13*A + 3*I*B)*sqrt(tan(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (67*A - 3*I*B)*sqrt(tan(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/8 - I/8)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_194():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = (A + I*B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x))) + (17*A + 7*I*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + (151*A + 41*I*B)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - (317*A + 67*I*B)*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*sqrt(tan(c + d*x))) + (sympy.S(1)/8 + I/8)*(A - I*B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_195():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = (A + I*B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)) + (21*A + 11*I*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + (89*A + 39*I*B)/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) - (361*A + 151*I*B)*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*tan(c + d*x)**(sympy.S(3)/2)) + (707*I*A - 317*B)*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*sqrt(tan(c + d*x))) + (sympy.S(1)/8 + I/8)*(I*A + B)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_196():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*B*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d - 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*x*(A - I*B)/4 + 3*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(I*A + B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_197():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)**2
    F = 3*B*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)**2/(8*d) - 9*B*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(8*d) + 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x*(A - I*B)/4 - 3*2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(2)/3)*sqrt(3)*a**(sympy.S(2)/3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) - (12*I*A + 3*B)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/3)/(20*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_198():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)
    F = 3*A*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(2*d) - 3*I*B*(I*a*tan(c + d*x) + a)**(sympy.S(5)/3)/(5*a*d) + 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x*(I*A + B)/4 + 3*2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(A - I*B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(A - I*B)*log(cos(c + d*x))/(4*d) + 2**(sympy.S(2)/3)*sqrt(3)*a**(sympy.S(2)/3)*(A - I*B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_199():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*B*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(2*d) - 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x*(A - I*B)/4 + 3*2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) + 2**(sympy.S(2)/3)*sqrt(3)*a**(sympy.S(2)/3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_200():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*cot(c + d*x)
    F = 3*A*a**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - A*a**(sympy.S(2)/3)*log(tan(c + d*x))/(2*d) + sqrt(3)*A*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d - 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x*(I*A + B)/4 - 3*2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(A - I*B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(A - I*B)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(2)/3)*sqrt(3)*a**(sympy.S(2)/3)*(A - I*B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_201():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*cot(c + d*x)**2
    F = -A*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)*cot(c + d*x)/d + 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x*(A - I*B)/4 - 3*2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(2)/3)*sqrt(3)*a**(sympy.S(2)/3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) + a**(sympy.S(2)/3)*(2*I*A + 3*B)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - a**(sympy.S(2)/3)*(2*I*A + 3*B)*log(tan(c + d*x))/(6*d) + sqrt(3)*a**(sympy.S(2)/3)*(2*I*A + 3*B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_202():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = (3*I*A - 3*B)/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*x*(A - I*B)/(8*a**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*(I*A + B)*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*(3*I*A + 3*B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_203():
    f = (A + B*tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)
    F = (3*I*A - 3*B)/(4*d*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*x*(A - I*B)/(8*a**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*(I*A + B)*log(cos(c + d*x))/(8*a**(sympy.S(2)/3)*d) - 2**(sympy.S(1)/3)*sqrt(3)*(I*A + B)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(2)/3)*d) + 2**(sympy.S(1)/3)*(3*I*A + 3*B)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(2)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_204():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**4*tan(c + d*x)**m
    F = I*B*a*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**(m + 1)/(d*(m + 4)) + 8*a**4*(A - I*B)*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), I*tan(c + d*x))/(d*(m + 1)) - 2*a**4*(A*(2*m**3 + 19*m**2 + 60*m + 64) - I*B*(2*m**3 + 19*m**2 + 60*m + 67))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 2)*(m + 3)*(m + 4)) - (A*(m + 4) - I*B*(m + 7))*(I*a**2*tan(c + d*x) + a**2)**2*tan(c + d*x)**(m + 1)/(d*(m + 3)*(m + 4)) - (2*A*(m + 4)**2 - 2*I*B*(m**2 + 8*m + 19))*(I*a**4*tan(c + d*x) + a**4)*tan(c + d*x)**(m + 1)/(d*(m + 2)*(m + 3)*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_205():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*tan(c + d*x)**m
    F = I*B*a*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(m + 1)/(d*(m + 3)) + 4*a**3*(A - I*B)*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), I*tan(c + d*x))/(d*(m + 1)) - a**3*(A*(2*m**2 + 11*m + 15) - I*B*(2*m**2 + 11*m + 17))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 2)*(m + 3)) - (A*(m + 3) - I*B*(m + 5))*(I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**(m + 1)/(d*(m + 2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_206():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**m
    F = I*B*(I*a**2*tan(c + d*x) + a**2)*tan(c + d*x)**(m + 1)/(d*(m + 2)) + 2*a**2*(A - I*B)*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), I*tan(c + d*x))/(d*(m + 1)) + I*a**2*(B + (m + 2)*(I*A + B))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_207():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*tan(c + d*x)**m
    F = I*B*a*tan(c + d*x)**(m + 1)/(d*(m + 1)) + a*(A - I*B)*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), I*tan(c + d*x))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_208():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*tan(c + d*x)**(m + 1)/(2*d*(I*a*tan(c + d*x) + a)) + m*(I*A - B)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(2*a*d*(m + 2)) + (A*(1 - m) - I*B*(m + 1))*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(2*a*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_209():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**2
    F = (A + I*B)*tan(c + d*x)**(m + 1)/(4*d*(I*a*tan(c + d*x) + a)**2) + m*(I*A*(2 - m) + B*m)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(4*a**2*d*(m + 2)) + (1 - m)*(A*(1 - m) - I*B*(m + 1))*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(4*a**2*d*(m + 1)) + (A*(2 - m) - I*B*m)*tan(c + d*x)**(m + 1)/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_210():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**3
    F = (2 - m)*(A*(5 - 2*m) - I*(2*B*m + B))*tan(c + d*x)**(m + 1)/(24*d*(I*a**3*tan(c + d*x) + a**3)) + (A + I*B)*tan(c + d*x)**(m + 1)/(6*d*(I*a*tan(c + d*x) + a)**3) + (A*(7 - 2*m) + I*B*(1 - 2*m))*tan(c + d*x)**(m + 1)/(24*a*d*(I*a*tan(c + d*x) + a)**2) + m*(2 - m)*(I*A*(5 - 2*m) + 2*B*m + B)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(24*a**3*d*(m + 2)) - (1 - m)*(-A*(2*m**2 - 7*m + 3) + I*B*(-2*m**2 + m + 3))*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(24*a**3*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_211():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**4
    F = (A + I*B)*tan(c + d*x)**(m + 1)/(8*d*(I*a*tan(c + d*x) + a)**4) + (A*(5 - m) + I*B*(1 - m))*tan(c + d*x)**(m + 1)/(24*a*d*(I*a*tan(c + d*x) + a)**3) + m*(2 - m)*(I*A*(m**2 - 6*m + 8) + B*(-m**2 + 2*m + 2))*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(48*a**4*d*(m + 2)) - (2 - m)*(-A*(m**2 - 6*m + 8) + I*B*(-m**2 + 2*m + 2))*tan(c + d*x)**(m + 1)/(48*a**4*d*(I*tan(c + d*x) + 1)) - (-A*(m**2 - 7*m + 13) + I*B*(-m**2 + 3*m + 1))*tan(c + d*x)**(m + 1)/(48*a**4*d*(I*tan(c + d*x) + 1)**2) - (-A*(m**2 - 4*m + 1) + I*B*(1 - m**2))*(m**2 - 4*m + 3)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(48*a**4*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_212():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**m
    F = 2*I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(m + 1)/(d*(2*m + 5)) + 4*a**3*(A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(d*(m + 1)*sqrt(I*a*tan(c + d*x) + a)) + 2*a**2*(-A*(2*m + 5) + 2*I*B*(m + 4))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(m + 1)/(d*(2*m + 3)*(2*m + 5)) + 2*a**2*(I*A*(8*m**2 + 34*m + 35) + 2*B*(4*m**2 + 17*m + 19))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(d*(-I*tan(c + d*x))**m*(2*m + 3)*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_213():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**m
    F = 2*I*B*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(m + 1)/(d*(2*m + 3)) + 2*a**2*(A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(d*(m + 1)*sqrt(I*a*tan(c + d*x) + a)) + 2*a*(B + (2*m + 3)*(I*A + B))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(d*(-I*tan(c + d*x))**m*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_214():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m
    F = 2*B*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(d*(-I*tan(c + d*x))**m) + a*(A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(d*(m + 1)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_215():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/sqrt(I*a*tan(c + d*x) + a)
    F = (A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(2*d*(m + 1)*sqrt(I*a*tan(c + d*x) + a)) + (A + I*B)*tan(c + d*x)**(m + 1)/(d*sqrt(I*a*tan(c + d*x) + a)) + (2*m + 1)*(I*A - B)*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(a*d*(-I*tan(c + d*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_216():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (A + I*B)*tan(c + d*x)**(m + 1)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(4*a*d*(m + 1)*sqrt(I*a*tan(c + d*x) + a)) + (A*(5 - 4*m) - I*(4*B*m + B))*tan(c + d*x)**(m + 1)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) + (2*m + 1)*sqrt(I*a*tan(c + d*x) + a)*(I*A*(5 - 4*m) + 4*B*m + B)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(6*a**2*d*(-I*tan(c + d*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_217():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (A + I*B)*tan(c + d*x)**(m + 1)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (A*(11 - 4*m) + I*B*(1 - 4*m))*tan(c + d*x)**(m + 1)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (A - I*B)*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -I*tan(c + d*x), I*tan(c + d*x))/(8*a**2*d*(m + 1)*sqrt(I*a*tan(c + d*x) + a)) - (-A*(16*m**2 - 52*m + 37) + I*B*(-16*m**2 + 12*m + 13))*tan(c + d*x)**(m + 1)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (2*m + 1)*(I*A*(16*m**2 - 52*m + 37) + B*(-16*m**2 + 12*m + 13))*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), I*tan(c + d*x) + 1)/(60*a**3*d*(-I*tan(c + d*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_218():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**m
    F = I*B*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(m + 1)*hyper((m + 1, 1 - n), (m + 2,), -I*tan(c + d*x))/(d*(m + 1)*(I*tan(c + d*x) + 1)**n) + (A - I*B)*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, 1 - n, m + 2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(m + 1)*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_219():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**3
    F = B*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**3/(d*(n + 3)) - (-A*(n + 3) + I*B*n)*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**2/(d*(n + 2)*(n + 3)) + (A - I*B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n) + (-2*A*(n + 3) + 2*I*B*n)*(I*a*tan(c + d*x) + a)**n/(d*n*(n + 2)*(n + 3)) - (I*a*tan(c + d*x) + a)**(n + 1)*(A*n*(n + 3) - I*B*(n**2 + 3*n + 6))/(a*d*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_220():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**2
    F = B*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**2/(d*(n + 2)) - 2*B*(I*a*tan(c + d*x) + a)**n/(d*n*(n + 2)) + (I*A + B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n) - (I*A*(n + 2) + B*n)*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_221():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)
    F = A*(I*a*tan(c + d*x) + a)**n/(d*n) - I*B*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(n + 1)) - (A - I*B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_222():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n
    F = B*(I*a*tan(c + d*x) + a)**n/(d*n) - (I*A + B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_223():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)
    F = -A*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x) + 1)/(d*n) + (A - I*B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_224():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**2
    F = -A*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)/d + (I*A + B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n) - (I*A*n + B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x) + 1)/(d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_225():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**3
    F = -A*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**2/(2*d) - (I*A*n + 2*B)*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)/(2*d) - (A - I*B)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n) - (-A*(n**2 - n + 2) + 2*I*B*n)*(I*a*tan(c + d*x) + a)**n*hyper((1, n), (n + 1,), I*tan(c + d*x) + 1)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_226():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(5)/2)/(d*(2*n + 5)) + (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) - (-2*A*(2*n + 5) + 4*I*B*n)*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 3)*(2*n + 5)) - (I*a*tan(c + d*x) + a)**n*(4*I*A*n*(2*n + 5) + 2*B*(4*n**2 + 10*n + 15))*sqrt(tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)*(2*n + 5)) - (2*I*A*(8*n**3 + 32*n**2 + 36*n + 15) + 8*B*n*(2*n**2 + 8*n + 9))*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)*(2*n + 5)*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_227():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 3)) - (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) - (-2*A*(2*n + 3) + 4*I*B*n)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)) + (I*a*tan(c + d*x) + a)**n*(4*A*n*(2*n + 3) - 2*I*B*(4*n**2 + 6*n + 3))*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_228():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))
    F = 2*B*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))/(d*(2*n + 1)) - (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) + (2*I*A*(2*n + 1) + 4*B*n)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_229():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/sqrt(tan(c + d*x))
    F = 2*I*B*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) + (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_230():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*I*A*(1 - 2*n)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) - 2*A*(I*a*tan(c + d*x) + a)**n/(d*sqrt(tan(c + d*x))) + (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_231():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*(I*a*tan(c + d*x) + a)**n/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - (2 - 4*n)*(-2*A*n + 3*I*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(3*d*(I*tan(c + d*x) + 1)**n) - (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n) - (4*I*A*n + 6*B)*(I*a*tan(c + d*x) + a)**n/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_232():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*tan(c + d*x)**2
    F = B*b*tan(c + d*x)**3/(3*d) - x*(A*a - B*b) + (A*a - B*b)*tan(c + d*x)/d + (A*b + B*a)*log(cos(c + d*x))/d + (A*b + B*a)*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_233():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*tan(c + d*x)
    F = B*b*tan(c + d*x)**2/(2*d) - x*(A*b + B*a) - (A*a - B*b)*log(cos(c + d*x))/d + (A*b + B*a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_234():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))
    F = B*b*tan(c + d*x)/d + x*(A*a - B*b) - (A*b + B*a)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_235():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)
    F = A*a*log(sin(c + d*x))/d - B*b*log(cos(c + d*x))/d + x*(A*b + B*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_236():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**2
    F = -A*a*cot(c + d*x)/d - x*(A*a - B*b) + (A*b + B*a)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_237():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**3
    F = -A*a*cot(c + d*x)**2/(2*d) - x*(A*b + B*a) - (A*a - B*b)*log(sin(c + d*x))/d - (A*b + B*a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_238():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**4
    F = -A*a*cot(c + d*x)**3/(3*d) + x*(A*a - B*b) + (A*a - B*b)*cot(c + d*x)/d - (A*b + B*a)*log(sin(c + d*x))/d - (A*b + B*a)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_239():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**5
    F = -A*a*cot(c + d*x)**4/(4*d) + x*(A*b + B*a) + (A*a - B*b)*log(sin(c + d*x))/d + (A*a - B*b)*cot(c + d*x)**2/(2*d) - (A*b + B*a)*cot(c + d*x)**3/(3*d) + (A*b + B*a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_240():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*tan(c + d*x)**2
    F = -B*(a + b*tan(c + d*x))**2/(2*d) + B*(a + b*tan(c + d*x))**3*tan(c + d*x)/(4*b*d) - b*(A*b + B*a)*tan(c + d*x)/d - x*(A*a**2 - A*b**2 - 2*B*a*b) + (2*A*a*b + B*a**2 - B*b**2)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**3*(4*A*b - B*a)/(12*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_241():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*tan(c + d*x)
    F = A*(a + b*tan(c + d*x))**2/(2*d) + B*(a + b*tan(c + d*x))**3/(3*b*d) + b*(A*a - B*b)*tan(c + d*x)/d - x*(2*A*a*b + B*a**2 - B*b**2) - (A*a**2 - A*b**2 - 2*B*a*b)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_242():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2
    F = B*(a + b*tan(c + d*x))**2/(2*d) + b*(A*b + B*a)*tan(c + d*x)/d + x*(A*a**2 - A*b**2 - 2*B*a*b) - (2*A*a*b + B*a**2 - B*b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_243():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)
    F = A*a**2*log(sin(c + d*x))/d + B*b**2*tan(c + d*x)/d - b*(A*b + 2*B*a)*log(cos(c + d*x))/d + x*(2*A*a*b + B*a**2 - B*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_244():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**2
    F = -A*a**2*cot(c + d*x)/d - B*b**2*log(cos(c + d*x))/d + a*(2*A*b + B*a)*log(sin(c + d*x))/d - x*(A*a**2 - A*b**2 - 2*B*a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_245():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**3
    F = -A*a**2*cot(c + d*x)**2/(2*d) - a*(2*A*b + B*a)*cot(c + d*x)/d + x*(B*b**2 - a*(2*A*b + B*a)) - (A*a**2 - A*b**2 - 2*B*a*b)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_246():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**4
    F = -A*a**2*cot(c + d*x)**3/(3*d) - a*(2*A*b + B*a)*cot(c + d*x)**2/(2*d) + x*(A*a**2 - A*b**2 - 2*B*a*b) + (B*b**2 - a*(2*A*b + B*a))*log(sin(c + d*x))/d + (A*a**2 - A*b**2 - 2*B*a*b)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_247():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**5
    F = -A*a**2*cot(c + d*x)**4/(4*d) - a*(2*A*b + B*a)*cot(c + d*x)**3/(3*d) + x*(2*A*a*b + B*a**2 - B*b**2) - (B*b**2 - a*(2*A*b + B*a))*cot(c + d*x)/d + (A*a**2 - A*b**2 - 2*B*a*b)*log(sin(c + d*x))/d + (A*a**2 - A*b**2 - 2*B*a*b)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_248():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*tan(c + d*x)**2
    F = -B*(a + b*tan(c + d*x))**3/(3*d) + B*(a + b*tan(c + d*x))**4*tan(c + d*x)/(5*b*d) - b*(2*A*a*b + B*a**2 - B*b**2)*tan(c + d*x)/d - x*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3) - (a + b*tan(c + d*x))**2*(A*b + B*a)/(2*d) + (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**4*(5*A*b - B*a)/(20*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_249():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*tan(c + d*x)
    F = A*(a + b*tan(c + d*x))**3/(3*d) + B*(a + b*tan(c + d*x))**4/(4*b*d) + b*(A*a**2 - A*b**2 - 2*B*a*b)*tan(c + d*x)/d - x*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2) + (a + b*tan(c + d*x))**2*(A*a - B*b)/(2*d) - (A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_250():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3
    F = B*(a + b*tan(c + d*x))**3/(3*d) + b*(2*A*a*b + B*a**2 - B*b**2)*tan(c + d*x)/d + x*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3) + (a + b*tan(c + d*x))**2*(A*b + B*a)/(2*d) - (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_251():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)
    F = A*a**3*log(sin(c + d*x))/d + B*b*(a + b*tan(c + d*x))**2/(2*d) + b**2*(A*b + 2*B*a)*tan(c + d*x)/d - b*(3*A*a*b + 3*B*a**2 - B*b**2)*log(cos(c + d*x))/d + x*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_252():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**2
    F = -A*a*(a + b*tan(c + d*x))**2*cot(c + d*x)/d + a**2*(3*A*b + B*a)*log(sin(c + d*x))/d + b**2*(A*a + B*b)*tan(c + d*x)/d - b**2*(A*b + 3*B*a)*log(cos(c + d*x))/d - x*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_253():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**3
    F = -A*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**2/(2*d) - B*b**3*log(cos(c + d*x))/d - a**2*(2*A*b + B*a)*cot(c + d*x)/d - a*(A*a**2 - 3*A*b**2 - 3*B*a*b)*log(sin(c + d*x))/d - x*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_254():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**4
    F = -A*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**3/(3*d) - a**2*(5*A*b + 3*B*a)*cot(c + d*x)**2/(6*d) + a*(3*A*a**2 - 8*A*b**2 - 9*B*a*b)*cot(c + d*x)/(3*d) + x*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3) - (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_255():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**5
    F = -A*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**4/(4*d) - a**2*(3*A*b + 2*B*a)*cot(c + d*x)**3/(6*d) + a*(2*A*a**2 - 5*A*b**2 - 6*B*a*b)*cot(c + d*x)**2/(4*d) + x*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2) + (A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)*log(sin(c + d*x))/d + (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_256():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**6
    F = -A*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**5/(5*d) - a**2*(7*A*b + 5*B*a)*cot(c + d*x)**4/(20*d) + a*(5*A*a**2 - 12*A*b**2 - 15*B*a*b)*cot(c + d*x)**3/(15*d) - x*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3) - (A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)*cot(c + d*x)/d + (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*log(sin(c + d*x))/d + (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_257():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*tan(c + d*x)**2
    F = -B*(a + b*tan(c + d*x))**4/(4*d) + B*(a + b*tan(c + d*x))**5*tan(c + d*x)/(6*b*d) - b*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*tan(c + d*x)/d - x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3) - (a + b*tan(c + d*x))**3*(A*b + B*a)/(3*d) - (a + b*tan(c + d*x))**2*(2*A*a*b + B*a**2 - B*b**2)/(2*d) + (4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**5*(6*A*b - B*a)/(30*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_258():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*tan(c + d*x)
    F = A*(a + b*tan(c + d*x))**4/(4*d) + B*(a + b*tan(c + d*x))**5/(5*b*d) + b*(A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)*tan(c + d*x)/d - x*(4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4) + (a + b*tan(c + d*x))**3*(A*a - B*b)/(3*d) + (a + b*tan(c + d*x))**2*(A*a**2 - A*b**2 - 2*B*a*b)/(2*d) - (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_259():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4
    F = B*(a + b*tan(c + d*x))**4/(4*d) + b*(3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*tan(c + d*x)/d + x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3) + (a + b*tan(c + d*x))**3*(A*b + B*a)/(3*d) + (a + b*tan(c + d*x))**2*(2*A*a*b + B*a**2 - B*b**2)/(2*d) - (4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_260():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)
    F = A*a**4*log(sin(c + d*x))/d + B*b*(a + b*tan(c + d*x))**3/(3*d) + b**2*(3*A*a*b + 3*B*a**2 - B*b**2)*tan(c + d*x)/d + b*(a + b*tan(c + d*x))**2*(A*b + 2*B*a)/(2*d) - b*(6*A*a**2*b - A*b**3 + 4*B*a**3 - 4*B*a*b**2)*log(cos(c + d*x))/d + x*(4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_261():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**2
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)/d + a**3*(4*A*b + B*a)*log(sin(c + d*x))/d + b**2*(A*a**2 + A*b**2 + 3*B*a*b)*tan(c + d*x)/d - b**2*(4*A*a*b + 6*B*a**2 - B*b**2)*log(cos(c + d*x))/d + b*(a + b*tan(c + d*x))**2*(2*A*a + B*b)/(2*d) - x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_262():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**3
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)**2/(2*d) - a**2*(A*a**2 - 6*A*b**2 - 4*B*a*b)*log(sin(c + d*x))/d - a*(a + b*tan(c + d*x))**2*(5*A*b + 2*B*a)*cot(c + d*x)/(2*d) - b**3*(A*b + 4*B*a)*log(cos(c + d*x))/d + b**2*(3*A*a*b + B*a**2 + B*b**2)*tan(c + d*x)/d - x*(4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_263():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**4
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)**3/(3*d) - B*b**4*log(cos(c + d*x))/d + a**2*(A*a**2 - 3*A*b**2 - 3*B*a*b)*cot(c + d*x)/d - a*(a + b*tan(c + d*x))**2*(2*A*b + B*a)*cot(c + d*x)**2/(2*d) - a*(4*A*a**2*b - 4*A*b**3 + B*a**3 - 6*B*a*b**2)*log(sin(c + d*x))/d + x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_264():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**5
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)**4/(4*d) + a**2*(6*A*a**2 - 13*A*b**2 - 16*B*a*b)*cot(c + d*x)**2/(12*d) - a*(a + b*tan(c + d*x))**2*(7*A*b + 4*B*a)*cot(c + d*x)**3/(12*d) + a*(24*A*a**2*b - 19*A*b**3 + 6*B*a**3 - 34*B*a*b**2)*cot(c + d*x)/(6*d) + x*(4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4) + (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_265():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**6
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)**5/(5*d) + a**2*(10*A*a**2 - 18*A*b**2 - 25*B*a*b)*cot(c + d*x)**3/(30*d) - a*(a + b*tan(c + d*x))**2*(8*A*b + 5*B*a)*cot(c + d*x)**4/(20*d) + a*(40*A*a**2*b - 28*A*b**3 + 10*B*a**3 - 55*B*a*b**2)*cot(c + d*x)**2/(20*d) - x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3) - (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*cot(c + d*x)/d + (4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_266():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*cot(c + d*x)**7
    F = -A*a*(a + b*tan(c + d*x))**3*cot(c + d*x)**6/(6*d) + a**2*(5*A*a**2 - 8*A*b**2 - 12*B*a*b)*cot(c + d*x)**4/(20*d) - a*(a + b*tan(c + d*x))**2*(3*A*b + 2*B*a)*cot(c + d*x)**5/(10*d) + a*(20*A*a**2*b - 13*A*b**3 + 5*B*a**3 - 27*B*a*b**2)*cot(c + d*x)**3/(15*d) - x*(4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4) - (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*log(sin(c + d*x))/d - (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*cot(c + d*x)**2/(2*d) - (4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_267():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))
    F = B*tan(c + d*x)**2/(2*b*d) - a**3*(A*b - B*a)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)) - x*(A*b - B*a)/(a**2 + b**2) + (A*a + B*b)*log(cos(c + d*x))/(d*(a**2 + b**2)) + (A*b - B*a)*tan(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_268():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))
    F = B*tan(c + d*x)/(b*d) + a**2*(A*b - B*a)*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)) - x*(A*a + B*b)/(a**2 + b**2) - (A*b - B*a)*log(cos(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_269():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))
    F = -B*log(cos(c + d*x))/(b*d) - a*(A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(b*d*(a**2 + b**2)) + x*(A*b - B*a)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_270():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))
    F = x*(A*a + B*b)/(a**2 + b**2) + (A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_271():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))
    F = A*log(sin(c + d*x))/(a*d) - x*(A*b - B*a)/(a**2 + b**2) - b*(A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_272():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))
    F = -A*cot(c + d*x)/(a*d) - x*(A*a + B*b)/(a**2 + b**2) + b**2*(A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)) - (A*b - B*a)*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_273():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))
    F = -A*cot(c + d*x)**2/(2*a*d) + x*(A*b - B*a)/(a**2 + b**2) + (A*b - B*a)*cot(c + d*x)/(a**2*d) - b**3*(A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)) - (A*a**2 - A*b**2 + B*a*b)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_274():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**4/(a + b*tan(c + d*x))
    F = -A*cot(c + d*x)**3/(3*a*d) + x*(A*a + B*b)/(a**2 + b**2) + (A*b - B*a)*cot(c + d*x)**2/(2*a**2*d) + (A*a**2 - A*b**2 + B*a*b)*cot(c + d*x)/(a**3*d) + b**4*(A*b - B*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)) + (a**2 - b**2)*(A*b - B*a)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_275():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = a**2*(A*a**2*b + 3*A*b**3 - 2*B*a**3 - 4*B*a*b**2)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**2) + a*(A*b - B*a)*tan(c + d*x)**2/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - x*(2*A*a*b - B*a**2 + B*b**2)/(a**2 + b**2)**2 + (A*a**2 - A*b**2 + 2*B*a*b)*log(cos(c + d*x))/(d*(a**2 + b**2)**2) - (A*a*b - 2*B*a**2 - B*b**2)*tan(c + d*x)/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_276():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -a**2*(A*b - B*a)/(b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - a*(2*A*b**3 - B*a*(a**2 + 3*b**2))*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)**2) - x*(A*a**2 - A*b**2 + 2*B*a*b)/(a**2 + b**2)**2 - (2*A*a*b - B*a**2 + B*b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_277():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**2
    F = a*(A*b - B*a)/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + x*(2*A*a*b - B*a**2 + B*b**2)/(a**2 + b**2)**2 - (A*a**2 - A*b**2 + 2*B*a*b)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_278():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**2
    F = x*(A*a**2 - A*b**2 + 2*B*a*b)/(a**2 + b**2)**2 + (2*A*a*b - B*a**2 + B*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) - (A*b - B*a)/(d*(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_279():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**2
    F = A*log(sin(c + d*x))/(a**2*d) - x*(2*A*a*b - B*a**2 + B*b**2)/(a**2 + b**2)**2 + b*(A*b - B*a)/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - b*(3*A*a**2*b + A*b**3 - 2*B*a**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_280():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -A*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))) - x*(A*a**2 - A*b**2 + 2*B*a*b)/(a**2 + b**2)**2 - b*(A*a**2 + 2*A*b**2 - B*a*b)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + b**2*(4*A*a**2*b + 2*A*b**3 - 3*B*a**3 - B*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**2) - (2*A*b - B*a)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_281():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = -A*cot(c + d*x)**2/(2*a*d*(a + b*tan(c + d*x))) + x*(2*A*a*b - B*a**2 + B*b**2)/(a**2 + b**2)**2 + (3*A*b - 2*B*a)*cot(c + d*x)/(2*a**2*d*(a + b*tan(c + d*x))) + b*(2*A*a**2*b + 3*A*b**3 - B*a**3 - 2*B*a*b**2)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - b**3*(5*A*a**2*b + 3*A*b**3 - 4*B*a**3 - 2*B*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**2) - (A*a**2 - 3*A*b**2 + 2*B*a*b)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_282():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = a**2*(A*a**4*b + 3*A*a**2*b**3 + 6*A*b**5 - 3*B*a**5 - 9*B*a**3*b**2 - 10*B*a*b**4)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**3) + a*(A*b - B*a)*tan(c + d*x)**3/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(A*a**2*b + 5*A*b**3 - 3*B*a**3 - 7*B*a*b**2)*tan(c + d*x)**2/(2*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + x*(A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(a**2 + b**2)**3 + (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**3) - (A*a**3*b + 3*A*a*b**3 - 3*B*a**4 - 6*B*a**2*b**2 - B*b**4)*tan(c + d*x)/(b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_283():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = -a**2*(2*A*b**3 - B*a*(a**2 + 3*b**2))/(b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + a*(A*b - B*a)*tan(c + d*x)**2/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(A*a**2*b**3 - 3*A*b**5 + B*a**5 + 3*B*a**3*b**2 + 6*B*a*b**4)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**3) - x*(3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(a**2 + b**2)**3 + (A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)*log(cos(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_284():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -a**2*(A*b - B*a)/(2*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(2*A*b**3 - B*a*(a**2 + 3*b**2))/(b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - x*(A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(a**2 + b**2)**3 - (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_285():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**3
    F = a*(A*b - B*a)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + x*(3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(a**2 + b**2)**3 - (A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + (A*a**2 - A*b**2 + 2*B*a*b)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_286():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**3
    F = x*(A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(a**2 + b**2)**3 + (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) - (2*A*a*b - B*a**2 + B*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (A*b - B*a)/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_287():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**3
    F = A*log(sin(c + d*x))/(a**3*d) - x*(3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(a**2 + b**2)**3 + b*(A*b - B*a)/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b*(3*A*a**2*b + A*b**3 - 2*B*a**3)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b*(6*A*a**4*b + 3*A*a**2*b**3 + A*b**5 - 3*B*a**5 + B*a**3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_288():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -A*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**2) - x*(A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(a**2 + b**2)**3 - b*(2*A*a**2 + 3*A*b**2 - B*a*b)/(2*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - b*(A*a**4 + 6*A*a**2*b**2 + 3*A*b**4 - 3*B*a**3*b - B*a*b**3)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + b**2*(10*A*a**4*b + 9*A*a**2*b**3 + 3*A*b**5 - 6*B*a**5 - 3*B*a**3*b**2 - B*a*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**3) - (3*A*b - B*a)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_289():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = -A*cot(c + d*x)**2/(2*a*d*(a + b*tan(c + d*x))**2) + x*(3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(a**2 + b**2)**3 + (2*A*b - B*a)*cot(c + d*x)/(a**2*d*(a + b*tan(c + d*x))**2) + b*(5*A*a**2*b + 6*A*b**3 - 2*B*a**3 - 3*B*a*b**2)/(2*a**3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b*(3*A*a**4*b + 11*A*a**2*b**3 + 6*A*b**5 - B*a**5 - 6*B*a**3*b**2 - 3*B*a*b**4)/(a**4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b**3*(15*A*a**4*b + 17*A*a**2*b**3 + 6*A*b**5 - 10*B*a**5 - 9*B*a**3*b**2 - 3*B*a*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**5*d*(a**2 + b**2)**3) - (A*a**2 - 6*A*b**2 + 3*B*a*b)*log(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_290():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(a + b*tan(c + d*x))**4
    F = a**2*(A*a**2*b**3 - 3*A*b**5 + B*a**5 + 3*B*a**3*b**2 + 6*B*a*b**4)/(b**4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + a*(A*b - B*a)*tan(c + d*x)**3/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + a*(2*A*b**3 - B*a*(a**2 + 3*b**2))*tan(c + d*x)**2/(2*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + a*(4*A*a**2*b**5 - 4*A*b**7 + B*a**7 + 4*B*a**5*b**2 + 5*B*a**3*b**4 + 10*B*a*b**6)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**4) + x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)/(a**2 + b**2)**4 + (4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)*log(cos(c + d*x))/(d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_291():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**4
    F = a**2*(A*a**2*b - 5*A*b**3 + 2*B*a**3 + 8*B*a*b**2)/(6*b**3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + a*(A*b - B*a)*tan(c + d*x)**2/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - a*(A*a**4*b + 5*A*a**2*b**3 - 8*A*b**5 + 2*B*a**5 + 7*B*a**3*b**2 + 17*B*a*b**4)/(3*b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - x*(4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)/(a**2 + b**2)**4 + (A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_292():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -a**2*(A*b - B*a)/(3*b**2*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + a*(2*A*b**3 - B*a*(a**2 + 3*b**2))/(2*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)/(a**2 + b**2)**4 - (4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_293():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**4
    F = a*(A*b - B*a)/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + x*(4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)/(a**2 + b**2)**4 - (A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + (A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + (A*a**2 - A*b**2 + 2*B*a*b)/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_294():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**4
    F = x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)/(a**2 + b**2)**4 + (4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) - (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - (2*A*a*b - B*a**2 + B*b**2)/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - (A*b - B*a)/(d*(a + b*tan(c + d*x))**3*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_295():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**4
    F = A*log(sin(c + d*x))/(a**4*d) - x*(4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)/(a**2 + b**2)**4 + b*(A*b - B*a)/(3*a*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + b*(3*A*a**2*b + A*b**3 - 2*B*a**3)/(2*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + b*(6*A*a**4*b + 3*A*a**2*b**3 + A*b**5 - 3*B*a**5 + B*a**3*b**2)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b*(10*A*a**6*b + 5*A*a**4*b**3 + 4*A*a**2*b**5 + A*b**7 - 4*B*a**7 + 4*B*a**5*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_296():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -A*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**3) - x*(A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)/(a**2 + b**2)**4 - b*(3*A*a**2 + 4*A*b**2 - B*a*b)/(3*a**2*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - b*(2*A*a**4 + 8*A*a**2*b**2 + 4*A*b**4 - 3*B*a**3*b - B*a*b**3)/(2*a**3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - b*(A*a**6 + 13*A*a**4*b**2 + 12*A*a**2*b**4 + 4*A*b**6 - 6*B*a**5*b - 3*B*a**3*b**3 - B*a*b**5)/(a**4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + b**2*(20*A*a**6*b + 24*A*a**4*b**3 + 16*A*a**2*b**5 + 4*A*b**7 - 10*B*a**7 - 5*B*a**5*b**2 - 4*B*a**3*b**4 - B*a*b**6)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**5*d*(a**2 + b**2)**4) - (4*A*b - B*a)*log(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_297():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**4
    F = -A*cot(c + d*x)**2/(2*a*d*(a + b*tan(c + d*x))**3) + x*(4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)/(a**2 + b**2)**4 + (5*A*b - 2*B*a)*cot(c + d*x)/(2*a**2*d*(a + b*tan(c + d*x))**3) + b*(9*A*a**2*b + 10*A*b**3 - 3*B*a**3 - 4*B*a*b**2)/(3*a**3*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + b*(7*A*a**4*b + 19*A*a**2*b**3 + 10*A*b**5 - 2*B*a**5 - 8*B*a**3*b**2 - 4*B*a*b**4)/(2*a**4*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + b*(4*A*a**6*b + 27*A*a**4*b**3 + 29*A*a**2*b**5 + 10*A*b**7 - B*a**7 - 13*B*a**5*b**2 - 12*B*a**3*b**4 - 4*B*a*b**6)/(a**5*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b**3*(35*A*a**6*b + 56*A*a**4*b**3 + 39*A*a**2*b**5 + 10*A*b**7 - 20*B*a**7 - 24*B*a**5*b**2 - 16*B*a**3*b**4 - 4*B*a*b**6)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**6*d*(a**2 + b**2)**4) - (A*a**2 - 10*A*b**2 + 4*B*a*b)*log(sin(c + d*x))/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_298():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))
    F = B*log(cos(c + d*x))/d + B*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_299():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))
    F = -B*x + B*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_300():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))
    F = -B*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_301():
    f = (B*a + B*b*tan(c + d*x))/(a + b*tan(c + d*x))
    F = B*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_302():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))
    F = B*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_303():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))
    F = -B*x - B*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_304():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))
    F = -B*log(sin(c + d*x))/d - B*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_305():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**4/(a + b*tan(c + d*x))
    F = B*x - B*cot(c + d*x)**3/(3*d) + B*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_306():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = B*a**4*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)) + B*a*x/(a**2 + b**2) - B*a*tan(c + d*x)/(b**2*d) + B*b*log(cos(c + d*x))/(d*(a**2 + b**2)) + B*tan(c + d*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_307():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = -B*a**3*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)) + B*a*log(cos(c + d*x))/(d*(a**2 + b**2)) - B*b*x/(a**2 + b**2) + B*tan(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_308():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = B*a**3*x/(b**2*(a**2 + b**2)) + B*a**2*log(a*cos(c + d*x) + b*sin(c + d*x))/(b*d*(a**2 + b**2)) - B*a*x/b**2 - B*log(cos(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_309():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**2
    F = -B*a*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)) + B*b*x/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_310():
    f = (B*a + B*b*tan(c + d*x))/(a + b*tan(c + d*x))**2
    F = B*a*x/(a**2 + b**2) + B*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_311():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**2
    F = -B*b*x/(a**2 + b**2) - B*b**2*log(a*cos(c + d*x) + b*sin(c + d*x))/(a*d*(a**2 + b**2)) + B*log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_312():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -B*a*x/(a**2 + b**2) - B*cot(c + d*x)/(a*d) + B*b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)) - B*b*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_313():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = B*b*x/(a**2 + b**2) - B*cot(c + d*x)**2/(2*a*d) + B*b*cot(c + d*x)/(a**2*d) - B*b**4*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)) - B*(a**2 - b**2)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_314():
    f = (tan(c + d*x) + 3)/(2 - tan(c + d*x))
    F = x - log(-sin(c + d*x) + 2*cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_315():
    f = (B*tan(c + d*x) + B*b/a)/(a + b*tan(c + d*x))
    F = 2*B*b*x/(a**2 + b**2) - B*(a - b**2/a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_316():
    f = (a + b*tan(c + d*x))/(a*tan(c + d*x) + b)**2
    F = -a*x*(a**2 - 3*b**2)/(a**2 + b**2)**2 + b*(3*a**2 - b**2)*log(a*sin(c + d*x) + b*cos(c + d*x))/(d*(a**2 + b**2)**2) - (a**2 - b**2)/(d*(a**2 + b**2)*(a*tan(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_317():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**3
    F = -2*A*sqrt(a + b*tan(c + d*x))/d + 2*B*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**2/(7*b*d) + (A - I*B)*sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (A + I*B)*sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(14*A*b - 8*B*a)*tan(c + d*x)/(35*b**2*d) - (a + b*tan(c + d*x))**(sympy.S(3)/2)*(28*A*a*b - 16*B*a**2 + 70*B*b**2)/(105*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_318():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**2
    F = -2*B*sqrt(a + b*tan(c + d*x))/d + 2*B*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*b*d) + sqrt(a - I*b)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - sqrt(a + I*b)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(10*A*b - 4*B*a)/(15*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_319():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*tan(c + d*x)
    F = 2*A*sqrt(a + b*tan(c + d*x))/d + 2*B*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*b*d) - (A - I*B)*sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_320():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))
    F = 2*B*sqrt(a + b*tan(c + d*x))/d - sqrt(a - I*b)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + sqrt(a + I*b)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_321():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)
    F = -2*A*sqrt(a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + (A - I*B)*sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (A + I*B)*sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_322():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2
    F = -A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/d + sqrt(a - I*b)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - sqrt(a + I*b)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - (A*b + 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_323():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**3
    F = -A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - (A - I*B)*sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*(A*b + 4*B*a)*cot(c + d*x)/(4*a*d) + (8*A*a**2 + A*b**2 - 4*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_324():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**4
    F = -A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**3/(3*d) - sqrt(a - I*b)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + sqrt(a + I*b)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*(A*b + 6*B*a)*cot(c + d*x)**2/(12*a*d) + sqrt(a + b*tan(c + d*x))*(8*A*a**2 + A*b**2 - 2*B*a*b)*cot(c + d*x)/(8*a**2*d) + (8*A*a**2*b - A*b**3 + 16*B*a**3 + 2*B*a*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_325():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**2
    F = -2*B*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*B*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*b*d) + (a - I*b)**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(3)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 2*B*a)/d + (a + b*tan(c + d*x))**(sympy.S(5)/2)*(14*A*b - 4*B*a)/(35*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_326():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)
    F = 2*A*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*B*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*b*d) - (A - I*B)*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + sqrt(a + b*tan(c + d*x))*(2*A*a - 2*B*b)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_327():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*B*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) - (a - I*b)**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(3)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + sqrt(a + b*tan(c + d*x))*(2*A*b + 2*B*a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_328():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)
    F = -2*A*a**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + 2*B*b*sqrt(a + b*tan(c + d*x))/d + (A - I*B)*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (A + I*B)*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_329():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**2
    F = -A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/d - sqrt(a)*(3*A*b + 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + (a - I*b)**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(3)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_330():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**3
    F = -A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - (A - I*B)*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*(5*A*b + 4*B*a)*cot(c + d*x)/(4*d) + (8*A*a**2 - 3*A*b**2 - 12*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_331():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**4
    F = -A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**3/(3*d) - (a - I*b)**(sympy.S(3)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(3)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*(7*A*b + 6*B*a)*cot(c + d*x)**2/(12*d) + sqrt(a + b*tan(c + d*x))*(8*A*a**2 - A*b**2 - 10*B*a*b)*cot(c + d*x)/(8*a*d) + (24*A*a**2*b + A*b**3 + 16*B*a**3 - 6*B*a*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(8*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_332():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**2
    F = -2*B*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) + 2*B*(a + b*tan(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(9*b*d) + (a - I*b)**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - (a + b*tan(c + d*x))**(sympy.S(3)/2)*(2*A*b + 2*B*a)/(3*d) - sqrt(a + b*tan(c + d*x))*(4*A*a*b + 2*B*a**2 - 2*B*b**2)/d + (a + b*tan(c + d*x))**(sympy.S(7)/2)*(18*A*b - 4*B*a)/(63*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_333():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)
    F = 2*A*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) + 2*B*(a + b*tan(c + d*x))**(sympy.S(7)/2)/(7*b*d) - (A - I*B)*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(2*A*a - 2*B*b)/(3*d) + sqrt(a + b*tan(c + d*x))*(2*A*a**2 - 2*A*b**2 - 4*B*a*b)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_334():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*B*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) - (a - I*b)**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(5)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(2*A*b + 2*B*a)/(3*d) + sqrt(a + b*tan(c + d*x))*(4*A*a*b + 2*B*a**2 - 2*B*b**2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_335():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)
    F = -2*A*a**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + 2*B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*b*sqrt(a + b*tan(c + d*x))*(A*b + 2*B*a)/d + (A - I*B)*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (A + I*B)*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_336():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**2
    F = -A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)/d - a**(sympy.S(3)/2)*(5*A*b + 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + b*sqrt(a + b*tan(c + d*x))*(A*a + 2*B*b)/d + (a - I*b)**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_337():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**3
    F = -A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**2/(2*d) + sqrt(a)*(8*A*a**2 - 15*A*b**2 - 20*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*d) - a*sqrt(a + b*tan(c + d*x))*(7*A*b + 4*B*a)*cot(c + d*x)/(4*d) - (A - I*B)*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (A + I*B)*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_338():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**4
    F = -A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**3/(3*d) - a*sqrt(a + b*tan(c + d*x))*(3*A*b + 2*B*a)*cot(c + d*x)**2/(4*d) - (a - I*b)**(sympy.S(5)/2)*(I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(5)/2)*(I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + sqrt(a + b*tan(c + d*x))*(8*A*a**2 - 11*A*b**2 - 18*B*a*b)*cot(c + d*x)/(8*d) + (40*A*a**2*b - 5*A*b**3 + 16*B*a**3 - 30*B*a*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_339():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**5
    F = -A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**4/(4*d) - a*sqrt(a + b*tan(c + d*x))*(11*A*b + 8*B*a)*cot(c + d*x)**3/(24*d) + (A - I*B)*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (A + I*B)*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + sqrt(a + b*tan(c + d*x))*(48*A*a**2 - 59*A*b**2 - 104*B*a*b)*cot(c + d*x)**2/(96*d) + sqrt(a + b*tan(c + d*x))*(144*A*a**2*b - 5*A*b**3 + 64*B*a**3 - 88*B*a*b**2)*cot(c + d*x)/(64*a*d) - (128*A*a**4 - 240*A*a**2*b**2 - 5*A*b**4 - 320*B*a**3*b + 40*B*a*b**3)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(64*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_340():
    f = (-a + b*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) - 2*b*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)/d + (a - I*b)**(sympy.S(5)/2)*(I*a - b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*(I*a + b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_341():
    f = (-a + b*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) - sqrt(2)*b*(a**2 + b**2)*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*(a**2 + b**2)*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*(a**2 + b**2)*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + sqrt(2)*b*(a**2 + b**2)*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_342():
    f = (-a + b*tan(c + d*x))*sqrt(a + b*tan(c + d*x))
    F = 2*b*sqrt(a + b*tan(c + d*x))/d + sqrt(2)*b*sqrt(a**2 + b**2)*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*sqrt(a**2 + b**2)*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*sqrt(a**2 + b**2)*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + sqrt(2)*b*sqrt(a**2 + b**2)*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_343():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/sqrt(a + b*tan(c + d*x))
    F = 2*B*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**2/(5*b*d) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + sqrt(a + b*tan(c + d*x))*(10*A*b - 8*B*a)*tan(c + d*x)/(15*b**2*d) - sqrt(a + b*tan(c + d*x))*(20*A*a*b - 16*B*a**2 + 30*B*b**2)/(15*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_344():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/sqrt(a + b*tan(c + d*x))
    F = 2*B*sqrt(a + b*tan(c + d*x))*tan(c + d*x)/(3*b*d) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + sqrt(a + b*tan(c + d*x))*(6*A*b - 4*B*a)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_345():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/sqrt(a + b*tan(c + d*x))
    F = 2*B*sqrt(a + b*tan(c + d*x))/(b*d) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_346():
    f = (A + B*tan(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) - (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_347():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/sqrt(a + b*tan(c + d*x))
    F = -2*A*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(sqrt(a)*d) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_348():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/sqrt(a + b*tan(c + d*x))
    F = -A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(a*d) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + (A*b - 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_349():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/sqrt(a + b*tan(c + d*x))
    F = -A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*a*d) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + sqrt(a + b*tan(c + d*x))*(3*A*b - 4*B*a)*cot(c + d*x)/(4*a**2*d) + (8*A*a**2 - 3*A*b**2 + 4*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_350():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(A*b - B*a)*tan(c + d*x)**2/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - sqrt(a + b*tan(c + d*x))*(6*A*a*b - 8*B*a**2 - 2*B*b**2)*tan(c + d*x)/(3*b**2*d*(a**2 + b**2)) + sqrt(a + b*tan(c + d*x))*(12*A*a**2*b + 6*A*b**3 - 16*B*a**3 - 10*B*a*b**2)/(3*b**3*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_351():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(a + b*tan(c + d*x))/(b**2*d) - 2*a**2*(A*b - B*a)/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_352():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(A*b - B*a)/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_353():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -(2*A*b - 2*B*a)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_354():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*A*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + 2*b*(A*b - B*a)/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_355():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -A*cot(c + d*x)/(a*d*sqrt(a + b*tan(c + d*x))) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) - b*(A*a**2 + 3*A*b**2 - 2*B*a*b)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + (3*A*b - 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_356():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -A*cot(c + d*x)**2/(2*a*d*sqrt(a + b*tan(c + d*x))) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + (5*A*b - 4*B*a)*cot(c + d*x)/(4*a**2*d*sqrt(a + b*tan(c + d*x))) + b*(7*A*a**2*b + 15*A*b**3 - 4*B*a**3 - 12*B*a*b**2)/(4*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + (8*A*a**2 - 15*A*b**2 + 12*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_357():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**4/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*tan(c + d*x)**3/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*a*(A*a**2*b + 3*A*b**3 - 2*B*a**3 - 4*B*a*b**2)*tan(c + d*x)**2/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) - sqrt(a + b*tan(c + d*x))*(8*A*a**3*b + 20*A*a*b**3 - 16*B*a**4 - 30*B*a**2*b**2 - 2*B*b**4)*tan(c + d*x)/(3*b**3*d*(a**2 + b**2)**2) + sqrt(a + b*tan(c + d*x))*(16*A*a**4*b + 34*A*a**2*b**3 + 6*A*b**5 - 32*B*a**5 - 60*B*a**3*b**2 - 16*B*a*b**4)/(3*b**4*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_358():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(A*a**2*b + 7*A*b**3 - 4*B*a**3 - 10*B*a*b**2)/(3*b**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + 2*a*(A*b - B*a)*tan(c + d*x)**2/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - sqrt(a + b*tan(c + d*x))*(2*A*a*b - 8*B*a**2 - 6*B*b**2)/(3*b**3*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_359():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(A*b - B*a)/(3*b**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*a*(2*A*b**3 - B*a*(a**2 + 3*b**2))/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_360():
    f = (A + B*tan(c + d*x))*tan(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (2*A*a**2 - 2*A*b**2 + 4*B*a*b)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_361():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -(4*A*a*b - 2*B*a**2 + 2*B*b**2)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (2*A*b - 2*B*a)/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_362():
    f = (A + B*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*A*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d) + (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) + (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + 2*b*(A*b - B*a)/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*b*(3*A*a**2*b + A*b**3 - 2*B*a**3)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_363():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -A*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) - (I*A - B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (I*A + B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) - b*(3*A*a**2 + 5*A*b**2 - 2*B*a*b)/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - b*(A*a**4 + 10*A*a**2*b**2 + 5*A*b**4 - 6*B*a**3*b - 2*B*a*b**3)/(a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (5*A*b - 2*B*a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_364():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -A*cot(c + d*x)**2/(2*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) - (A - I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) - (A + I*B)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (7*A*b - 4*B*a)*cot(c + d*x)/(4*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + b*(27*A*a**2*b + 35*A*b**3 - 12*B*a**3 - 20*B*a*b**2)/(12*a**3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + b*(11*A*a**4*b + 62*A*a**2*b**3 + 35*A*b**5 - 4*B*a**5 - 40*B*a**3*b**2 - 20*B*a*b**4)/(4*a**4*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (8*A*a**2 - 35*A*b**2 + 20*B*a*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_365():
    f = (B*a + B*b*tan(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = sqrt(2)*B*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*B*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*B*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(2)*B*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_366():
    f = (B*a + B*b*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -sqrt(2)*B*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*B*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*B*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*B*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_367():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) - 2*B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_368():
    f = (B*a + B*b*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*B*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + I*B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - I*B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_369():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + 2*B*b**2/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*B*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_370():
    f = (-a + b*tan(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = -(I*a + b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + (I*a - b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_371():
    f = (-a + b*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 4*a*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - (I*a + b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + (I*a - b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_372():
    f = (-a + b*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 4*a*b/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + 2*b*(3*a**2 - b**2)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (I*a + b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + (I*a - b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_373():
    f = (I*tan(c + d*x) + 1)/sqrt(a + b*tan(c + d*x))
    F = -2*I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_374():
    f = (-I*tan(c + d*x) + 1)/sqrt(a + b*tan(c + d*x))
    F = 2*I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_375():
    f = (tan(x) + 3)/sqrt(3*tan(x) + 4)
    F = -sqrt(2)*atan(sqrt(2)*(1 - 3*tan(x))/(2*sqrt(3*tan(x) + 4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_376():
    f = (1 - 3*tan(x))/sqrt(3*tan(x) + 4)
    F = sqrt(2)*atanh(sqrt(2)*(tan(x) + 3)/(2*sqrt(3*tan(x) + 4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_377():
    f = (4 - 3*tan(a + b*x))/sqrt(3*tan(a + b*x) + 4)
    F = -9*sqrt(2)*atan(sqrt(2)*(1 - 3*tan(a + b*x))/(2*sqrt(3*tan(a + b*x) + 4)))/(10*b) + 13*sqrt(2)*atanh(sqrt(2)*(tan(a + b*x) + 3)/(2*sqrt(3*tan(a + b*x) + 4)))/(10*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_378():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*tan(c + d*x)**(sympy.S(7)/2)/(7*d) + (2*A*a - 2*B*b)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*A*b + 2*B*a)*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - (2*A*b + 2*B*a)*sqrt(tan(c + d*x))/d - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_379():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*tan(c + d*x)**(sympy.S(5)/2)/(5*d) + (2*A*a - 2*B*b)*sqrt(tan(c + d*x))/d + (2*A*b + 2*B*a)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(2)*(a*(A - B) - b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A - B) - b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_380():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*sqrt(tan(c + d*x))
    F = 2*B*b*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*A*b + 2*B*a)*sqrt(tan(c + d*x))/d + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_381():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))/sqrt(tan(c + d*x))
    F = 2*B*b*sqrt(tan(c + d*x))/d - sqrt(2)*(a*(A - B) - b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A - B) - b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_382():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a/(d*sqrt(tan(c + d*x))) - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_383():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - (2*A*b + 2*B*a)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a*(A - B) - b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A - B) - b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_384():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + (2*A*a - 2*B*b)/(d*sqrt(tan(c + d*x))) - (2*A*b + 2*B*a)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_385():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*b*(9*A*b + 11*B*a)*tan(c + d*x)**(sympy.S(7)/2)/(63*d) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + (4*A*a*b + 2*B*a**2 - 2*B*b**2)*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - (4*A*a*b + 2*B*a**2 - 2*B*b**2)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_386():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*b*(7*A*b + 9*B*a)*tan(c + d*x)**(sympy.S(5)/2)/(35*d) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*sqrt(tan(c + d*x))/d + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + (4*A*a*b + 2*B*a**2 - 2*B*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_387():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*sqrt(tan(c + d*x))
    F = 2*B*b*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*b*(5*A*b + 7*B*a)*tan(c + d*x)**(sympy.S(3)/2)/(15*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + (4*A*a*b + 2*B*a**2 - 2*B*b**2)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_388():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2/sqrt(tan(c + d*x))
    F = 2*B*b*(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(3*d) + 2*b*(3*A*b + 5*B*a)*sqrt(tan(c + d*x))/(3*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_389():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a**2/(d*sqrt(tan(c + d*x))) + 2*B*b**2*sqrt(tan(c + d*x))/d - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_390():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a**2/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a*(2*A*b + B*a)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_391():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a**2/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*a*(2*A*b + B*a)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_392():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(5)/2)/(9*d) + 2*b**2*(9*A*b + 13*B*a)*tan(c + d*x)**(sympy.S(7)/2)/(63*d) + 2*b*(27*A*a*b + 22*B*a**2 - 9*B*b**2)*tan(c + d*x)**(sympy.S(5)/2)/(45*d) + (2*A*a**3 - 6*A*a*b**2 - 6*B*a**2*b + 2*B*b**3)*sqrt(tan(c + d*x))/d + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_393():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*sqrt(tan(c + d*x))
    F = 2*B*b*(a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*b**2*(7*A*b + 11*B*a)*tan(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*b*(21*A*a*b + 18*B*a**2 - 7*B*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(21*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_394():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3/sqrt(tan(c + d*x))
    F = 2*B*b*(a + b*tan(c + d*x))**2*sqrt(tan(c + d*x))/(5*d) + 2*b**2*(5*A*b + 9*B*a)*tan(c + d*x)**(sympy.S(3)/2)/(15*d) + 2*b*(15*A*a*b + 14*B*a**2 - 5*B*b**2)*sqrt(tan(c + d*x))/(5*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_395():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**2/(d*sqrt(tan(c + d*x))) + 2*b**2*(3*A*a + B*b)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b*(2*A*a**2 + A*b**2 + 3*B*a*b)*sqrt(tan(c + d*x))/d - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_396():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**2/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(7*A*b + 3*B*a)/(3*d*sqrt(tan(c + d*x))) + 2*b**2*(A*a + 3*B*b)*sqrt(tan(c + d*x))/(3*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_397():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**2/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*a**2*(9*A*b + 5*B*a)/(15*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*a*(5*A*a**2 - 14*A*b**2 - 15*B*a*b)/(5*d*sqrt(tan(c + d*x))) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_398():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = 2*B*tan(c + d*x)**(sympy.S(3)/2)/(3*b*d) - 2*a**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + (2*A*b - 2*B*a)*sqrt(tan(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_399():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = 2*B*sqrt(tan(c + d*x))/(b*d) + 2*a**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_400():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))
    F = -2*sqrt(a)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(b)*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_401():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = -sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + 2*sqrt(b)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(a)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_402():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2))
    F = -2*A/(a*d*sqrt(tan(c + d*x))) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - 2*b**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_403():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2))
    F = -2*A/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + (2*A*b - 2*B*a)/(a**2*d*sqrt(tan(c + d*x))) + 2*b**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_404():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**2
    F = a**(sympy.S(3)/2)*(A*a**2*b + 5*A*b**3 - 3*B*a**3 - 7*B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)**2) + a*(A*b - B*a)*tan(c + d*x)**(sympy.S(3)/2)/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - (A*a*b - 3*B*a**2 - 2*B*b**2)*sqrt(tan(c + d*x))/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_405():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**2
    F = sqrt(a)*(A*a**2*b - 3*A*b**3 + B*a**3 + 5*B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)**2) + a*(A*b - B*a)*sqrt(tan(c + d*x))/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_406():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**2
    F = sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - (A*b - B*a)*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) - (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(a)*sqrt(b)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_407():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*sqrt(tan(c + d*x)))
    F = -sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + b*(A*b - B*a)*sqrt(tan(c + d*x))/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + sqrt(b)*(5*A*a**2*b + A*b**3 - 3*B*a**3 + B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_408():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b*(A*b - B*a)/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(tan(c + d*x))) - (2*A*a**2 + 3*A*b**2 - B*a*b)/(a**2*d*(a**2 + b**2)*sqrt(tan(c + d*x))) - b**(sympy.S(3)/2)*(7*A*a**2*b + 3*A*b**3 - 5*B*a**3 - B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_409():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(5)/2))
    F = sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + b*(A*b - B*a)/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)*tan(c + d*x)**(sympy.S(3)/2)) - (2*A*a**2 + 5*A*b**2 - 3*B*a*b)/(3*a**2*d*(a**2 + b**2)*tan(c + d*x)**(sympy.S(3)/2)) + (4*A*a**2*b + 5*A*b**3 - 2*B*a**3 - 3*B*a*b**2)/(a**3*d*(a**2 + b**2)*sqrt(tan(c + d*x))) + b**(sympy.S(5)/2)*(9*A*a**2*b + 5*A*b**3 - 7*B*a**3 - 3*B*a*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(7)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_410():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))**3
    F = a**(sympy.S(3)/2)*(3*A*a**4*b + 6*A*a**2*b**3 + 35*A*b**5 - 15*B*a**5 - 46*B*a**3*b**2 - 63*B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(7)/2)*d*(a**2 + b**2)**3) + a*(A*b - B*a)*tan(c + d*x)**(sympy.S(5)/2)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(A*a**2*b + 9*A*b**3 - 5*B*a**3 - 13*B*a*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(4*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - (3*A*a**3*b + 11*A*a*b**3 - 15*B*a**4 - 31*B*a**2*b**2 - 8*B*b**4)*sqrt(tan(c + d*x))/(4*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_411():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**3
    F = sqrt(a)*(A*a**4*b + 18*A*a**2*b**3 - 15*A*b**5 + 3*B*a**5 + 6*B*a**3*b**2 + 35*B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(5)/2)*d*(a**2 + b**2)**3) + a*(A*b - B*a)*tan(c + d*x)**(sympy.S(3)/2)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a*(A*a**2*b - 7*A*b**3 + 3*B*a**3 + 11*B*a*b**2)*sqrt(tan(c + d*x))/(4*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_412():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**3
    F = a*(A*b - B*a)*sqrt(tan(c + d*x))/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + (3*A*a**2*b - 5*A*b**3 + B*a**3 + 9*B*a*b**2)*sqrt(tan(c + d*x))/(4*b*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (3*A*a**4*b - 26*A*a**2*b**3 + 3*A*b**5 + B*a**5 + 18*B*a**3*b**2 - 15*B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*sqrt(a)*b**(sympy.S(3)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_413():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**3
    F = sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - (A*b - B*a)*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2)) - (7*A*a**2*b - A*b**3 - 3*B*a**3 + 5*B*a*b**2)*sqrt(tan(c + d*x))/(4*a*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (15*A*a**4*b - 18*A*a**2*b**3 - A*b**5 - 3*B*a**5 + 26*B*a**3*b**2 - 3*B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_414():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*sqrt(tan(c + d*x)))
    F = -sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b*(A*b - B*a)*sqrt(tan(c + d*x))/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b*(11*A*a**2*b + 3*A*b**3 - 7*B*a**3 + B*a*b**2)*sqrt(tan(c + d*x))/(4*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + sqrt(b)*(35*A*a**4*b + 6*A*a**2*b**3 + 3*A*b**5 - 15*B*a**5 + 18*B*a**3*b**2 + B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(5)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_415():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b*(A*b - B*a)/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)*sqrt(tan(c + d*x))) + b*(13*A*a**2*b + 5*A*b**3 - 9*B*a**3 - B*a*b**2)/(4*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(tan(c + d*x))) - (8*A*a**4 + 31*A*a**2*b**2 + 15*A*b**4 - 11*B*a**3*b - 3*B*a*b**3)/(4*a**3*d*(a**2 + b**2)**2*sqrt(tan(c + d*x))) - b**(sympy.S(3)/2)*(63*A*a**4*b + 46*A*a**2*b**3 + 15*A*b**5 - 35*B*a**5 - 6*B*a**3*b**2 - 3*B*a*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(7)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_416():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + 2*B*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_417():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + 2*B*sqrt(tan(c + d*x))/d - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_418():
    f = (B*a + B*b*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_419():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_420():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 2*B/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_421():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 2*B/(3*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_422():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**2
    F = -2*B*a**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) - sqrt(2)*B*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*B*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + 2*B*sqrt(tan(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_423():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**2
    F = 2*B*a**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(b)*d*(a**2 + b**2)) - sqrt(2)*B*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*B*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*B*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_424():
    f = (B*a + B*b*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**2
    F = -2*B*sqrt(a)*sqrt(b)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(d*(a**2 + b**2)) + sqrt(2)*B*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*B*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_425():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**2*sqrt(tan(c + d*x)))
    F = sqrt(2)*B*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*B*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*B*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + 2*B*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(a)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_426():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*B*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*B*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*B*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - 2*B/(a*d*sqrt(tan(c + d*x))) - 2*B*b**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_427():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(2*b*d) + (I*A - B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(4*A*b - B*a)*sqrt(tan(c + d*x))/(4*b*d) + (4*A*a*b - B*a**2 - 8*B*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_428():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))
    F = B*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/d - (A - I*B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (2*A*b + B*a)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_429():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/sqrt(tan(c + d*x))
    F = 2*B*sqrt(b)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A - B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_430():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) + (A - I*B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_431():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + (I*A - B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 6*B*a)/(3*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_432():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - (A - I*B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 10*B*a)/(15*a*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(30*A*a**2 + 4*A*b**2 - 10*B*a*b)/(15*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_433():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - (I*A - B)*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 14*B*a)/(35*a*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 + 8*A*b**2 - 14*B*a*b)/(105*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(70*A*a**2*b - 16*A*b**3 + 210*B*a**3 + 28*B*a*b**2)/(105*a**3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_434():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x))/(3*b*d) + (A - I*B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(6*A*b - B*a)*sqrt(tan(c + d*x))/(12*b*d) + sqrt(a + b*tan(c + d*x))*(6*A*a*b - B*a**2 - 8*B*b**2)*sqrt(tan(c + d*x))/(8*b*d) + (6*A*a**2*b - 16*A*b**3 - B*a**3 - 24*B*a*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_435():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))
    F = B*b*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(2*d) + (a + I*b)**2*(I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*(4*A*b + 5*B*a)*sqrt(tan(c + d*x))/(4*d) + (I*A + B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (12*A*a*b + 3*B*a**2 - 8*B*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_436():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/sqrt(tan(c + d*x))
    F = B*b*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/d + sqrt(b)*(2*A*b + 3*B*a)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A - I*B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_437():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) + 2*B*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (a + I*b)**2*(I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_438():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + (A - I*B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(8*A*b + 6*B*a)/(3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_439():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + (a + I*b)**2*(I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - sqrt(a + b*tan(c + d*x))*(12*A*b + 10*B*a)/(15*d*tan(c + d*x)**(sympy.S(3)/2)) + (I*A + B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 6*A*b**2 - 40*B*a*b)/(15*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_440():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - (A - I*B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(16*A*b + 14*B*a)/(35*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 - 6*A*b**2 - 84*B*a*b)/(105*a*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(280*A*a**2*b + 12*A*b**3 + 210*B*a**3 - 42*B*a*b**2)/(105*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_441():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))/(9*d*tan(c + d*x)**(sympy.S(9)/2)) - sqrt(a + b*tan(c + d*x))*(20*A*b + 18*B*a)/(63*d*tan(c + d*x)**(sympy.S(7)/2)) + (I*A - B)*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(42*A*a**2 - 2*A*b**2 - 48*B*a*b)/(105*a*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(252*A*a**2*b + 8*A*b**3 + 210*B*a**3 - 18*B*a*b**2)/(315*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*tan(c + d*x))*(630*A*a**4 - 126*A*a**2*b**2 + 16*A*b**4 - 840*B*a**3*b - 36*B*a*b**3)/(315*a**3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_442():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(7)/2)*sqrt(tan(c + d*x))/(4*b*d) - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(5)/2)*(8*A*b - B*a)*sqrt(tan(c + d*x))/(24*b*d) + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(40*A*a*b - 5*B*a**2 - 48*B*b**2)*sqrt(tan(c + d*x))/(96*b*d) + sqrt(a + b*tan(c + d*x))*(40*A*a**2*b - 64*A*b**3 - 5*B*a**3 - 112*B*a*b**2)*sqrt(tan(c + d*x))/(64*b*d) + (40*A*a**3*b - 320*A*a*b**3 - 5*B*a**4 - 240*B*a**2*b**2 + 128*B*b**4)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(64*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_443():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x))
    F = B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + (A - I*B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(2*A*b + 3*B*a)*sqrt(tan(c + d*x))/(4*d) + sqrt(a + b*tan(c + d*x))*(14*A*a*b + 5*B*a**2 - 8*B*b**2)*sqrt(tan(c + d*x))/(8*d) + (30*A*a**2*b - 16*A*b**3 + 5*B*a**3 - 40*B*a*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_444():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/sqrt(tan(c + d*x))
    F = B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(2*d) + sqrt(b)*(20*A*a*b + 15*B*a**2 - 8*B*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*d) + b*sqrt(a + b*tan(c + d*x))*(4*A*b + 7*B*a)*sqrt(tan(c + d*x))/(4*d) + (I*A - B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_445():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(d*sqrt(tan(c + d*x))) + b**(sympy.S(3)/2)*(2*A*b + 5*B*a)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b*sqrt(a + b*tan(c + d*x))*(2*A*a + B*b)*sqrt(tan(c + d*x))/d - (A - I*B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_446():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*B*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*a*sqrt(a + b*tan(c + d*x))*(2*A*b + B*a)/(d*sqrt(tan(c + d*x))) - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_447():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*a*sqrt(a + b*tan(c + d*x))*(8*A*b + 5*B*a)/(15*d*tan(c + d*x)**(sympy.S(3)/2)) + (A - I*B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 46*A*b**2 - 70*B*a*b)/(15*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_448():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - 2*a*sqrt(a + b*tan(c + d*x))*(10*A*b + 7*B*a)/(35*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 - 90*A*b**2 - 154*B*a*b)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) + (I*A - B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(490*A*a**2*b - 30*A*b**3 + 210*B*a**3 - 322*B*a*b**2)/(105*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_449():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(9*d*tan(c + d*x)**(sympy.S(9)/2)) - 2*a*sqrt(a + b*tan(c + d*x))*(4*A*b + 3*B*a)/(21*d*tan(c + d*x)**(sympy.S(7)/2)) - (A - I*B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(42*A*a**2 - 50*A*b**2 - 90*B*a*b)/(105*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(462*A*a**2*b - 10*A*b**3 + 210*B*a**3 - 270*B*a*b**2)/(315*a*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*tan(c + d*x))*(630*A*a**4 - 966*A*a**2*b**2 - 20*A*b**4 - 1470*B*a**3*b + 90*B*a*b**3)/(315*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_450():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(13)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(11*d*tan(c + d*x)**(sympy.S(11)/2)) - 2*a*sqrt(a + b*tan(c + d*x))*(14*A*b + 11*B*a)/(99*d*tan(c + d*x)**(sympy.S(9)/2)) + sqrt(a + b*tan(c + d*x))*(198*A*a**2 - 226*A*b**2 - 418*B*a*b)/(693*d*tan(c + d*x)**(sympy.S(7)/2)) - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(990*A*a**2*b - 10*A*b**3 + 462*B*a**3 - 550*B*a*b**2)/(1155*a*d*tan(c + d*x)**(sympy.S(5)/2)) - sqrt(a + b*tan(c + d*x))*(2310*A*a**4 - 2970*A*a**2*b**2 - 40*A*b**4 - 5082*B*a**3*b + 110*B*a*b**3)/(3465*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*tan(c + d*x))*(16170*A*a**4*b - 990*A*a**2*b**3 + 80*A*b**5 + 6930*B*a**5 - 10626*B*a**3*b**2 - 220*B*a*b**4)/(3465*a**3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_451():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*(B*tan(c + d*x) + 3*B*b/(2*a))/tan(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(d*tan(c + d*x)**(sympy.S(3)/2)) - B*sqrt(a + b*tan(c + d*x))*(2*a**2 + 6*b**2)/(d*sqrt(tan(c + d*x))) + B*(2*a - 3*I*b)*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(2*a*d) - B*(2*a + 3*I*b)*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_452():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*tan(c + d*x))
    F = B*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(b*d) - (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + (2*A*b - B*a)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_453():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = 2*B*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d) + (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_454():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_455():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2))
    F = -2*A*sqrt(a + b*tan(c + d*x))/(a*d*sqrt(tan(c + d*x))) - (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_456():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2))
    F = -2*A*sqrt(a + b*tan(c + d*x))/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) - (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*(4*A*b - 6*B*a)/(3*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_457():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2))
    F = -2*A*sqrt(a + b*tan(c + d*x))/(5*a*d*tan(c + d*x)**(sympy.S(5)/2)) + (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + sqrt(a + b*tan(c + d*x))*(8*A*b - 10*B*a)/(15*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 16*A*b**2 + 20*B*a*b)/(15*a**3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_458():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*B*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d) + 2*a*(A*b - B*a)*sqrt(tan(c + d*x))/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_459():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (2*A*b - 2*B*a)*sqrt(tan(c + d*x))/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_460():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x)))
    F = (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + 2*b*(A*b - B*a)*sqrt(tan(c + d*x))/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_461():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = -2*A/(a*d*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))) - (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2*b*(A*a**2 + 2*A*b**2 - B*a*b)*sqrt(tan(c + d*x))/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_462():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = -2*A/(3*a*d*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)) - (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + (8*A*b - 6*B*a)/(3*a**2*d*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))) + 2*b*(5*A*a**2*b + 8*A*b**3 - 3*B*a**3 - 6*B*a*b**2)*sqrt(tan(c + d*x))/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_463():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*B*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(5)/2)*d) + 2*a*(A*b - B*a)*tan(c + d*x)**(sympy.S(3)/2)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*a*(2*A*b**3 - B*a*(a**2 + 3*b**2))*sqrt(tan(c + d*x))/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_464():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*sqrt(tan(c + d*x))/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*A*a**2*b - 8*A*b**3 + 2*B*a**3 + 14*B*a*b**2)*sqrt(tan(c + d*x))/(3*b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_465():
    f = (A + B*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -(I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - (2*A*b - 2*B*a)*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) - (10*A*a**2*b - 2*A*b**3 - 4*B*a**3 + 8*B*a*b**2)*sqrt(tan(c + d*x))/(3*a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_466():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x)))
    F = -(A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + 2*b*(A*b - B*a)*sqrt(tan(c + d*x))/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*b*(8*A*a**2*b + 2*A*b**3 - 5*B*a**3 + B*a*b**2)*sqrt(tan(c + d*x))/(3*a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_467():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = -2*A/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + (I*A - B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - (I*A + B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - 2*b*(3*A*a**2 + 4*A*b**2 - B*a*b)*sqrt(tan(c + d*x))/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 2*b*(3*A*a**4 + 17*A*a**2*b**2 + 8*A*b**4 - 8*B*a**3*b - 2*B*a*b**3)*sqrt(tan(c + d*x))/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_468():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = -2*A/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + (A - I*B)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + (A + I*B)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*A*b - 2*B*a)/(a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + 2*b*(7*A*a**2*b + 8*A*b**3 - 3*B*a**3 - 4*B*a*b**2)*sqrt(tan(c + d*x))/(3*a**3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*b*(8*A*a**4*b + 30*A*a**2*b**3 + 16*A*b**5 - 3*B*a**5 - 17*B*a**3*b**2 - 8*B*a*b**4)*sqrt(tan(c + d*x))/(3*a**4*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_469():
    f = (B*a + B*b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -B*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - B*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + 2*B*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_470():
    f = (B*a + B*b*tan(c + d*x))*sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -I*B*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + I*B*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_471():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x)))
    F = B*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + B*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_472():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = I*B*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - I*B*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*B*sqrt(a + b*tan(c + d*x))/(a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_473():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(2)/3)
    F = 3*B*(a + b*tan(c + d*x))**(sympy.S(2)/3)/(2*d) - x*(A - I*B)*(a - I*b)**(sympy.S(2)/3)/4 - x*(A + I*B)*(a + I*b)**(sympy.S(2)/3)/4 + 3*(a - I*b)**(sympy.S(2)/3)*(I*A + B)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) + (a - I*b)**(sympy.S(2)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) + sqrt(3)*(a - I*b)**(sympy.S(2)/3)*(I*A + B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d) - 3*(a + I*b)**(sympy.S(2)/3)*(I*A - B)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) - (a + I*b)**(sympy.S(2)/3)*(I*A - B)*log(cos(c + d*x))/(4*d) - sqrt(3)*(a + I*b)**(sympy.S(2)/3)*(I*A - B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_474():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(1)/3)
    F = 3*B*(a + b*tan(c + d*x))**(sympy.S(1)/3)/d - x*(A - I*B)*(a - I*b)**(sympy.S(1)/3)/4 - x*(A + I*B)*(a + I*b)**(sympy.S(1)/3)/4 + 3*(a - I*b)**(sympy.S(1)/3)*(I*A + B)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) + (a - I*b)**(sympy.S(1)/3)*(I*A + B)*log(cos(c + d*x))/(4*d) - sqrt(3)*(a - I*b)**(sympy.S(1)/3)*(I*A + B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d) - 3*(a + I*b)**(sympy.S(1)/3)*(I*A - B)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) - (a + I*b)**(sympy.S(1)/3)*(I*A - B)*log(cos(c + d*x))/(4*d) + sqrt(3)*(a + I*b)**(sympy.S(1)/3)*(I*A - B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_475():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(1)/3)
    F = -x*(A - I*B)/(4*(a - I*b)**(sympy.S(1)/3)) - x*(A + I*B)/(4*(a + I*b)**(sympy.S(1)/3)) - (I*A - B)*log(cos(c + d*x))/(4*d*(a + I*b)**(sympy.S(1)/3)) - sqrt(3)*(I*A - B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d*(a + I*b)**(sympy.S(1)/3)) - (3*I*A - 3*B)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a + I*b)**(sympy.S(1)/3)) + (I*A + B)*log(cos(c + d*x))/(4*d*(a - I*b)**(sympy.S(1)/3)) + sqrt(3)*(I*A + B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d*(a - I*b)**(sympy.S(1)/3)) + (3*I*A + 3*B)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a - I*b)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_476():
    f = (A + B*tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(2)/3)
    F = -x*(A - I*B)/(4*(a - I*b)**(sympy.S(2)/3)) - x*(A + I*B)/(4*(a + I*b)**(sympy.S(2)/3)) - (I*A - B)*log(cos(c + d*x))/(4*d*(a + I*b)**(sympy.S(2)/3)) + sqrt(3)*(I*A - B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d*(a + I*b)**(sympy.S(2)/3)) - (3*I*A - 3*B)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a + I*b)**(sympy.S(2)/3)) + (I*A + B)*log(cos(c + d*x))/(4*d*(a - I*b)**(sympy.S(2)/3)) - sqrt(3)*(I*A + B)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d*(a - I*b)**(sympy.S(2)/3)) + (3*I*A + 3*B)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a - I*b)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_477():
    f = (-tan(e + f*x) + I)/(c + d*tan(e + f*x))**(sympy.S(1)/3)
    F = -I*x/(2*(c - I*d)**(sympy.S(1)/3)) - 3*log((c - I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(2*f*(c - I*d)**(sympy.S(1)/3)) - log(cos(e + f*x))/(2*f*(c - I*d)**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - I*d)**(sympy.S(1)/3))/3)/(f*(c - I*d)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_478():
    f = (-c*tan(e + f*x) + d)/(c + d*tan(e + f*x))**(sympy.S(2)/3)
    F = -I*x*(c - I*d)**(sympy.S(1)/3)/4 + I*x*(c + I*d)**(sympy.S(1)/3)/4 - 3*(c - I*d)**(sympy.S(1)/3)*log((c - I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c - I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c - I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - I*d)**(sympy.S(1)/3))/3)/(2*f) - 3*(c + I*d)**(sympy.S(1)/3)*log((c + I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c + I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c + I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + I*d)**(sympy.S(1)/3))/3)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_479():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**4*tan(c + d*x)**m
    F = B*b*(a + b*tan(c + d*x))**3*tan(c + d*x)**(m + 1)/(d*(m + 4)) + b**2*(2*A*a*b*(m + 4)**2 + B*a**2*(m**2 + 9*m + 26) - B*b**2*(m**2 + 7*m + 12))*tan(c + d*x)**(m + 2)/(d*(m + 2)*(m + 3)*(m + 4)) + b*(a + b*tan(c + d*x))**2*(A*b*(m + 4) + B*a*(m + 7))*tan(c + d*x)**(m + 1)/(d*(m + 3)*(m + 4)) - b*(-A*a**2*b*(5*m**2 + 37*m + 68) + A*b**3*(m**2 + 7*m + 12) - 2*B*a**3*(m**2 + 8*m + 19) + 4*B*a*b**2*(m**2 + 7*m + 12))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 3)*(m + 4)) + (4*A*a**3*b - 4*A*a*b**3 + B*a**4 - 6*B*a**2*b**2 + B*b**4)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(m + 2)) + (A*a**4 - 6*A*a**2*b**2 + A*b**4 - 4*B*a**3*b + 4*B*a*b**3)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_480():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*tan(c + d*x)**m
    F = B*b*(a + b*tan(c + d*x))**2*tan(c + d*x)**(m + 1)/(d*(m + 3)) + b**2*(A*b*(m + 3) + B*a*(m + 5))*tan(c + d*x)**(m + 2)/(d*(m + 2)*(m + 3)) + b*(3*A*a*b*(m + 3) + 2*B*a**2*(m + 4) - B*b**2*(m + 3))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 3)) + (3*A*a**2*b - A*b**3 + B*a**3 - 3*B*a*b**2)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(m + 2)) + (A*a**3 - 3*A*a*b**2 - 3*B*a**2*b + B*b**3)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_481():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*tan(c + d*x)**m
    F = B*b*(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)/(d*(m + 2)) + b*(A*b*(m + 2) + B*a*(m + 3))*tan(c + d*x)**(m + 1)/(d*(m + 1)*(m + 2)) + (2*A*a*b + B*a**2 - B*b**2)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(m + 2)) + (A*a**2 - A*b**2 - 2*B*a*b)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_482():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*tan(c + d*x)**m
    F = B*b*tan(c + d*x)**(m + 1)/(d*(m + 1)) + (A*b + B*a)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(m + 2)) + (A*a - B*b)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_483():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))
    F = -(A*b - B*a)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)*(m + 2)) + (A*a + B*b)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)*(m + 1)) + b*(A*b - B*a)*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -b*tan(c + d*x)/a)/(a*d*(a**2 + b**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_484():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))**2
    F = -(2*A*a*b - B*a**2 + B*b**2)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**2*(m + 2)) + (A*a**2 - A*b**2 + 2*B*a*b)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**2*(m + 1)) + b*(A*b - B*a)*tan(c + d*x)**(m + 1)/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + b*(A*a**2*b*(2 - m) - A*b**3*m + B*a*b**2*(m + 1) - a**3*(-B*m + B))*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -b*tan(c + d*x)/a)/(a**2*d*(a**2 + b**2)**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_485():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))**3
    F = -(3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**3*(m + 2)) + (A*a**3 - 3*A*a*b**2 + 3*B*a**2*b - B*b**3)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**3*(m + 1)) + b*(A*b - B*a)*tan(c + d*x)**(m + 1)/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b*(A*a**2*b*(5 - m) + A*b**3*(1 - m) - B*a**3*(3 - m) + B*a*b**2*(m + 1))*tan(c + d*x)**(m + 1)/(2*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b*(-A*a**4*b*(m**2 - 5*m + 6) + 2*A*a**2*b**3*(-m**2 + 3*m + 1) + A*b**5*m*(1 - m) + B*a**5*(m**2 - 3*m + 2) - 2*B*a**3*b**2*(-m**2 + m + 3) + B*a*b**4*m*(m + 1))*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -b*tan(c + d*x)/a)/(2*a**3*d*(a**2 + b**2)**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_486():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))**4
    F = -(4*A*a**3*b - 4*A*a*b**3 - B*a**4 + 6*B*a**2*b**2 - B*b**4)*tan(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**4*(m + 2)) + (A*a**4 - 6*A*a**2*b**2 + A*b**4 + 4*B*a**3*b - 4*B*a*b**3)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(a**2 + b**2)**4*(m + 1)) + b*(A*b - B*a)*tan(c + d*x)**(m + 1)/(3*a*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + b*(A*a**2*b*(8 - m) + A*b**3*(2 - m) - B*a**3*(5 - m) + B*a*b**2*(m + 1))*tan(c + d*x)**(m + 1)/(6*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + b*(A*a**4*b*(m**2 - 9*m + 26) + 2*A*a**2*b**3*(m**2 - 6*m + 2) + A*b**5*(m**2 - 3*m + 2) - B*a**5*(m**2 - 6*m + 11) + 2*B*a**3*b**2*(-m**2 + 3*m + 7) + B*a*b**4*(1 - m**2))*tan(c + d*x)**(m + 1)/(6*a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b*(-A*a**6*b*(-m**3 + 9*m**2 - 26*m + 24) + 3*A*a**4*b**3*(m**3 - 7*m**2 + 10*m + 8) + 3*A*a**2*b**5*m*(m**2 - 5*m + 2) + A*b**7*m*(m**2 - 3*m + 2) + B*a**7*(-m**3 + 6*m**2 - 11*m + 6) - 3*B*a**5*b**2*(m**3 - 4*m**2 - m + 12) + 3*B*a**3*b**4*(-m**3 + 2*m**2 + 5*m + 2) + B*a*b**6*m*(1 - m**2))*tan(c + d*x)**(m + 1)*hyper((1, m + 1), (m + 2,), -b*tan(c + d*x)/a)/(6*a**4*d*(a**2 + b**2)**4*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_487():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**m
    F = a**2*(A - I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-5)/2, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1)) + a**2*(A + I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-5)/2, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_488():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**m
    F = a*(A - I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1)) + a*(A + I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_489():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**m
    F = (A - I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1)) + (A + I*B)*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_490():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/sqrt(a + b*tan(c + d*x))
    F = sqrt(1 + b*tan(c + d*x)/a)*(A - I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*(m + 1)) + sqrt(1 + b*tan(c + d*x)/a)*(A + I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_491():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = sqrt(1 + b*tan(c + d*x)/a)*(A - I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a*d*sqrt(a + b*tan(c + d*x))*(m + 1)) + sqrt(1 + b*tan(c + d*x)/a)*(A + I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a*d*sqrt(a + b*tan(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_492():
    f = (A + B*tan(c + d*x))*tan(c + d*x)**m/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = sqrt(1 + b*tan(c + d*x)/a)*(A - I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(5)/2, m + 2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a**2*d*sqrt(a + b*tan(c + d*x))*(m + 1)) + sqrt(1 + b*tan(c + d*x)/a)*(A + I*B)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(5)/2, m + 2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a**2*d*sqrt(a + b*tan(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_493():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)**m
    F = (A - I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, -n, m + 2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*d*(1 + b*tan(c + d*x)/a)**n*(m + 1)) + (A + I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, -n, m + 2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*d*(1 + b*tan(c + d*x)/a)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_494():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)**4
    F = B*(a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)**3/(b*d*(n + 4)) + (A - I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b)) - (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(n + 1)*(2*I*a - 2*b)) - (a + b*tan(c + d*x))**(n + 1)*(-A*b*(n + 4) + 3*B*a)*tan(c + d*x)**2/(b**2*d*(n + 3)*(n + 4)) - (a + b*tan(c + d*x))**(n + 1)*(B*b**2*(n + 3)*(n + 4) - 2*a*(-A*b*(n + 4) + 3*B*a))*tan(c + d*x)/(b**3*d*(n + 2)*(n + 3)*(n + 4)) - (a + b*tan(c + d*x))**(n + 1)*(A*b**3*(n + 2)*(n + 3)*(n + 4) - a*(B*b**2*(n + 3)*(n + 4) - 2*a*(-A*b*(n + 4) + 3*B*a)))/(b**4*d*(n + 1)*(n + 2)*(n + 3)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_495():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)**3
    F = B*(a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)**2/(b*d*(n + 3)) + (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*(I*A + B)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b)) - (a + b*tan(c + d*x))**(n + 1)*(-A*b*(n + 3) + 2*B*a)*tan(c + d*x)/(b**2*d*(n + 2)*(n + 3)) + (a + b*tan(c + d*x))**(n + 1)*(-A*a*b*(n + 3) + 2*B*a**2 - B*b**2*(n**2 + 5*n + 6))/(b**3*d*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_496():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)**2
    F = B*(a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)/(b*d*(n + 2)) + (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(n + 1)*(2*I*a - 2*b)) + (a + b*tan(c + d*x))**(n + 1)*(I*A + B)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*(-A*b*(n + 2) + B*a)/(b**2*d*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_497():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)
    F = B*(a + b*tan(c + d*x))**(n + 1)/(b*d*(n + 1)) - (A - I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1)) - (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_498():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n
    F = (A - I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b)) + (a + b*tan(c + d*x))**(n + 1)*(I*A - B)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_499():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*cot(c + d*x)
    F = -A*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(a*d*(n + 1)) + (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*(I*A + B)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_500():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*cot(c + d*x)**2
    F = -A*(a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)/(a*d) - (A - I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b)) + (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(n + 1)*(2*I*a - 2*b)) - (a + b*tan(c + d*x))**(n + 1)*(A*b*n + B*a)*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(a**2*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_501():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*cot(c + d*x)**3
    F = -A*(a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)**2/(2*a*d) - (A + I*B)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*(I*A + B)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(n + 1)*(2*I*a + 2*b)) - (a + b*tan(c + d*x))**(n + 1)*(-A*b*(1 - n) + 2*B*a)*cot(c + d*x)/(2*a**2*d) + (a + b*tan(c + d*x))**(n + 1)*(2*A*a**2 + A*b**2*n*(1 - n) - 2*B*a*b*n)*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(2*a**3*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_502():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*(A - I*B)*sqrt(cot(c + d*x))/d + 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a*(I*A + B)*cot(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_503():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*a*(I*A + B)*sqrt(cot(c + d*x))/d - 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_504():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(cot(c + d*x))/d - 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_505():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))
    F = 2*I*B*a/(d*sqrt(cot(c + d*x))) + 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_506():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/sqrt(cot(c + d*x))
    F = 2*I*B*a/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*(-1)**(sympy.S(3)/4)*a*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 2*a*(I*A + B)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_507():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)/cot(c + d*x)**(sympy.S(3)/2)
    F = 2*I*B*a/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*a*(A - I*B)/(d*sqrt(cot(c + d*x))) - 2*(-1)**(sympy.S(3)/4)*a*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 2*a*(I*A + B)/(3*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_508():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*(a**2*cot(c + d*x) + I*a**2)*cot(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a**2*(A - I*B)*sqrt(cot(c + d*x))/d + 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a**2*(7*I*A + 5*B)*cot(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_509():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*(a**2*cot(c + d*x) + I*a**2)*sqrt(cot(c + d*x))/(3*d) - 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a**2*(5*I*A + 3*B)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_510():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(3)/2)
    F = 2*I*B*(a**2*cot(c + d*x) + I*a**2)/(d*sqrt(cot(c + d*x))) - 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a**2*(A + I*B)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_511():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2*sqrt(cot(c + d*x))
    F = 2*I*B*(a**2*cot(c + d*x) + I*a**2)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(3*A - 5*I*B)/(3*d*sqrt(cot(c + d*x))) + 4*(-1)**(sympy.S(3)/4)*a**2*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_512():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**2/sqrt(cot(c + d*x))
    F = 2*I*B*(a**2*cot(c + d*x) + I*a**2)/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 4*(-1)**(sympy.S(3)/4)*a**2*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a**2*(5*A - 7*I*B)/(15*d*cot(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(I*A + B)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_513():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(a*cot(c + d*x) + I*a)**2*cot(c + d*x)**(sympy.S(3)/2)/(7*d) + 8*a**3*(23*A - 21*I*B)*cot(c + d*x)**(sympy.S(3)/2)/(105*d) + 8*a**3*(I*A + B)*sqrt(cot(c + d*x))/d + 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - (22*I*A + 14*B)*(a**3*cot(c + d*x) + I*a**3)*cot(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_514():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a*cot(c + d*x) + I*a)**2*sqrt(cot(c + d*x))/(5*d) + 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 16*a**3*(6*A - 5*I*B)*sqrt(cot(c + d*x))/(15*d) - (18*I*A + 10*B)*(a**3*cot(c + d*x) + I*a**3)*sqrt(cot(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_515():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(5)/2)
    F = -16*I*A*a**3*sqrt(cot(c + d*x))/(3*d) + 2*I*B*a*(a*cot(c + d*x) + I*a)**2/(d*sqrt(cot(c + d*x))) - 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - (2*A + 6*I*B)*(a**3*cot(c + d*x) + I*a**3)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_516():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(3)/2)
    F = -16*I*B*a**3*sqrt(cot(c + d*x))/(3*d) + 2*I*B*a*(a*cot(c + d*x) + I*a)**2/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - (6*A - 14*I*B)*(a**3*cot(c + d*x) + I*a**3)/(3*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_517():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3*sqrt(cot(c + d*x))
    F = 2*I*B*a*(a*cot(c + d*x) + I*a)**2/(5*d*cot(c + d*x)**(sympy.S(5)/2)) - 16*a**3*(5*A - 6*I*B)/(15*d*sqrt(cot(c + d*x))) + 8*(-1)**(sympy.S(3)/4)*a**3*(I*A + B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - (10*A - 18*I*B)*(a**3*cot(c + d*x) + I*a**3)/(15*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_518():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**3/sqrt(cot(c + d*x))
    F = 2*I*B*a*(a*cot(c + d*x) + I*a)**2/(7*d*cot(c + d*x)**(sympy.S(7)/2)) + 8*(-1)**(sympy.S(3)/4)*a**3*(A - I*B)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 8*a**3*(21*A - 23*I*B)/(105*d*cot(c + d*x)**(sympy.S(3)/2)) + 8*a**3*(I*A + B)/(d*sqrt(cot(c + d*x))) - (14*A - 22*I*B)*(a**3*cot(c + d*x) + I*a**3)/(35*d*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_519():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**(sympy.S(5)/2)/(2*d*(a*cot(c + d*x) + I*a)) - (7*A + 3*I*B)*cot(c + d*x)**(sympy.S(3)/2)/(6*a*d) + (5*I*A - 5*B)*sqrt(cot(c + d*x))/(2*a*d) + sqrt(2)*(A*(-7 - 5*I) + B*(5 - 3*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d) + sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(6 + I) + B*(1 + 4*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) + sqrt(2)*(A*(7 - 5*I) + B*(5 + 3*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(8*a*d) + sqrt(2)*(A*(7 + 5*I) - B*(5 - 3*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_520():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cot(c + d*x) + I*a)) - (5*A + I*B)*sqrt(cot(c + d*x))/(2*a*d) - sqrt(2)*(A*(-5 - 3*I) + B*(3 - I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(8*a*d) - sqrt(2)*(sympy.S(1)/8 - I/8)*(A*(4 + I) + B*(1 + 2*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(A*(5 - 3*I) + B*(3 + I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d) + sqrt(2)*(A*(5 + 3*I) - B*(3 - I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_521():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)
    F = (A + I*B)*sqrt(cot(c + d*x))/(2*d*(a*cot(c + d*x) + I*a)) - sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(2 + I) + B)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/4 - I/4)*(A*(2 + I) + B)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d) - sqrt(2)*(A*(3 + I) - B*(1 + I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d) + sqrt(2)*(A*(3 + I) - B*(1 + I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_522():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x)))
    F = (I*A - B)*sqrt(cot(c + d*x))/(2*d*(a*cot(c + d*x) + I*a)) - sqrt(2)*(sympy.S(1)/4 - I/4)*(A + B*(2 - I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/4 - I/4)*(A + B*(2 - I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/8 + I/8)*(A - B*(2 + I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/8 + I/8)*(A - B*(2 + I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_523():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2))
    F = (I*A - B)/(2*d*(a*cot(c + d*x) + I*a)*sqrt(cot(c + d*x))) - (A + 5*I*B)/(2*a*d*sqrt(cot(c + d*x))) - sqrt(2)*(A*(1 - 3*I) + B*(3 + 5*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(8*a*d) + sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(1 + 2*I) - B*(4 + I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(2 + I) + B*(1 + 4*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(2 + I) + B*(1 + 4*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_524():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2))
    F = (I*A - B)/(2*d*(a*cot(c + d*x) + I*a)*cot(c + d*x)**(sympy.S(3)/2)) - (3*A + 7*I*B)/(6*a*d*cot(c + d*x)**(sympy.S(3)/2)) - (5*I*A - 5*B)/(2*a*d*sqrt(cot(c + d*x))) + sqrt(2)*(sympy.S(1)/8 + I/8)*(A*(1 + 4*I) - B*(6 + I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(A*(3 - 5*I) + B*(5 + 7*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(16*a*d) - sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(4 + I) + B*(1 + 6*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/4 + I/4)*(A*(4 + I) + B*(1 + 6*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_525():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**2
    F = (A + I*B)*cot(c + d*x)**(sympy.S(5)/2)/(4*d*(a*cot(c + d*x) + I*a)**2) + (7*A + 3*I*B)*cot(c + d*x)**(sympy.S(3)/2)/(8*a**2*d*(cot(c + d*x) + I)) - (25*A + 5*I*B)*sqrt(cot(c + d*x))/(8*a**2*d) + sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(2 + 23*I) - B*(7 + 2*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(23 + 2*I) + B*(2 + 7*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(23 + 2*I) + B*(2 + 7*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(A*(25 + 21*I) - B*(9 - 5*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_526():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**2
    F = (A + I*B)*cot(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cot(c + d*x) + I*a)**2) + (5*A + I*B)*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I)) + sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(-7 + 2*I) + B*(2 + I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(-2 + 7*I) + B*(1 + 2*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) - sqrt(2)*(A*(9 - 5*I) + B*(1 - 3*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(32*a**2*d) + sqrt(2)*(A*(9 + 5*I) - B*(1 + 3*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_527():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*sqrt(cot(c + d*x)))
    F = (A + I*B)*sqrt(cot(c + d*x))/(4*d*(a*cot(c + d*x) + I*a)**2) + (3*I*A + B)*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I)) + sqrt(2)*(A*(-1 + 3*I) + B*(1 + 3*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(32*a**2*d) + sqrt(2)*(A*(-1 + 3*I) + B*(1 + 3*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**2*d) + sqrt(2)*(A*(1 + 3*I) + B*(1 - 3*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**2*d) - sqrt(2)*(A*(1 + 3*I) + B*(1 - 3*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_528():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(3)/2))
    F = (I*A - B)*sqrt(cot(c + d*x))/(4*d*(a*cot(c + d*x) + I*a)**2) + (A + 5*I*B)*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I)) + sqrt(2)*(A*(1 - 3*I) - B*(9 - 5*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**2*d) + sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(1 + 2*I) + B*(2 - 7*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(A*(1 + 3*I) + B*(9 + 5*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**2*d) + sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(2 + I) + B*(7 - 2*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_529():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(5)/2))
    F = (I*A - B)/(4*d*(a*cot(c + d*x) + I*a)**2*sqrt(cot(c + d*x))) + (3*A + 7*I*B)/(8*a**2*d*(cot(c + d*x) + I)*sqrt(cot(c + d*x))) + (5*I*A - 25*B)/(8*a**2*d*sqrt(cot(c + d*x))) + sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(2 + 7*I) - B*(23 + 2*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(7 + 2*I) + B*(2 + 23*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(7 + 2*I) + B*(2 + 23*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(A*(9 + 5*I) - B*(25 - 21*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_530():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**3
    F = (A + I*B)*cot(c + d*x)**(sympy.S(7)/2)/(6*d*(a*cot(c + d*x) + I*a)**3) + (28*A + 7*I*B)*cot(c + d*x)**(sympy.S(3)/2)/(24*d*(a**3*cot(c + d*x) + I*a**3)) + (5*A + 2*I*B)*cot(c + d*x)**(sympy.S(5)/2)/(12*a*d*(a*cot(c + d*x) + I*a)**2) - (30*A + 5*I*B)*sqrt(cot(c + d*x))/(8*a**3*d) + sqrt(2)*(sympy.S(1)/16 - I/16)*(A*(1 + 29*I) - B*(6 + I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(29 + I) + B*(1 + 6*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) + sqrt(2)*(sympy.S(1)/32 - I/32)*(A*(29 + I) + B*(1 + 6*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) + sqrt(2)*(A*(30 + 28*I) - B*(7 - 5*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_531():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**3
    F = 5*A*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + (A + I*B)*cot(c + d*x)**(sympy.S(5)/2)/(6*d*(a*cot(c + d*x) + I*a)**3) + (4*A + I*B)*cot(c + d*x)**(sympy.S(3)/2)/(12*a*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*(A*(-7 + 5*I) + 2*I*B)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(32*a**3*d) + sqrt(2)*(A*(-7 + 5*I) + 2*I*B)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**3*d) - sqrt(2)*(A*(7 + 5*I) - 2*I*B)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(A*(7 + 5*I) - 2*I*B)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_532():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*sqrt(cot(c + d*x)))
    F = A*sqrt(cot(c + d*x))/(4*a*d*(a*cot(c + d*x) + I*a)**2) + (A + I*B)*cot(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cot(c + d*x) + I*a)**3) + (2*I*A + B)*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + sqrt(2)*(2*I*A + B*(1 - I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) - sqrt(2)*(2*I*A + B*(1 - I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + I) + B)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**3*d) + sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + I) + B)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_533():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(3)/2))
    F = A*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + (A + I*B)*sqrt(cot(c + d*x))/(6*d*(a*cot(c + d*x) + I*a)**3) + (2*I*A + B)*sqrt(cot(c + d*x))/(12*a*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*(A*(-1 + I) + 2*B)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(A*(-1 + I) + 2*B)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(A*(1 + I) + 2*B)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(32*a**3*d) + sqrt(2)*(A*(1 + I) + 2*B)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_534():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(5)/2))
    F = 5*B*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + (I*A - B)*sqrt(cot(c + d*x))/(6*d*(a*cot(c + d*x) + I*a)**3) + (A + 4*I*B)*sqrt(cot(c + d*x))/(12*a*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*(2*A + B*(5 - 7*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(32*a**3*d) + sqrt(2)*(2*A + B*(5 - 7*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**3*d) - sqrt(2)*(2*A - B*(5 + 7*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d) + sqrt(2)*(2*A - B*(5 + 7*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(64*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_535():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(7)/2))
    F = (I*A - B)/(6*d*(a*cot(c + d*x) + I*a)**3*sqrt(cot(c + d*x))) - (7*I*A - 28*B)/(24*d*(a**3*cot(c + d*x) + I*a**3)*sqrt(cot(c + d*x))) + (2*A + 5*I*B)/(12*a*d*(a*cot(c + d*x) + I*a)**2*sqrt(cot(c + d*x))) + (5*A + 30*I*B)/(8*a**3*d*sqrt(cot(c + d*x))) - sqrt(2)*(sympy.S(1)/16 + I/16)*(A*(1 + 6*I) - B*(29 + I))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**3*d) + sqrt(2)*(A*(5 - 7*I) + B*(28 + 30*I))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(32*a**3*d) + sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(6 + I) + B*(1 + 29*I))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(1)/32 + I/32)*(A*(6 + I) + B*(1 + 29*I))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_536():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + sqrt(a)*(-1 - I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (26*A - 10*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(15*d) - (2*I*A + 10*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_537():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(a)*(1 + I)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (2*I*A + 6*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_538():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d + sqrt(a)*(1 + I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_539():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))
    F = -2*(-1)**(sympy.S(3)/4)*B*sqrt(a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + sqrt(a)*(1 - I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_540():
    f = (A + B*tan(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/sqrt(cot(c + d*x))
    F = B*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(cot(c + d*x))) - sqrt(a)*(1 + I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*sqrt(a)*(2*A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_541():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)/(7*d) + a**(sympy.S(3)/2)*(2 - 2*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a*(19*A - 21*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(105*d) - 2*a*(8*I*A + 7*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(35*d) + 4*a*(67*I*A + 63*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_542():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + a**(sympy.S(3)/2)*(-2 - 2*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a*(9*A - 10*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(15*d) - 2*a*(6*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_543():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + a**(sympy.S(3)/2)*(2 + 2*I)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a*(4*I*A + 3*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_544():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d + 2*(-1)**(sympy.S(1)/4)*B*a**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(3)/2)*(2 + 2*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_545():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))
    F = I*B*a*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(cot(c + d*x))) + a**(sympy.S(3)/2)*(2 - 2*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*(2*I*A + 3*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_546():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cot(c + d*x))
    F = I*B*a*sqrt(I*a*tan(c + d*x) + a)/(2*d*cot(c + d*x)**(sympy.S(3)/2)) - a**(sympy.S(3)/2)*(2 + 2*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*(12*A - 11*I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) + a*(4*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)/(4*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_547():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(9)/2)/(9*d) + a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(46*A - 45*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(105*d) - 8*a**2*(197*A - 195*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(315*d) - 2*a**2*(4*I*A + 3*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)/(21*d) + 8*a**2*(59*I*A + 60*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(315*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_548():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)/(7*d) + a**(sympy.S(5)/2)*(4 - 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(80*A - 77*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(105*d) - 2*a**2*(10*I*A + 7*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(35*d) + 4*a**2*(130*I*A + 133*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_549():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + a**(sympy.S(5)/2)*(-4 - 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 2*a**2*(38*A - 35*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(15*d) - 2*a**2*(8*I*A + 5*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_550():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*(-1)**(sympy.S(3)/4)*B*a**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 + 4*I)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*(2*I*A + B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_551():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))/d + a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(2*A - 5*I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**2*(2*I*A - B)*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_552():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))
    F = I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(2*d*sqrt(cot(c + d*x))) + a**(sympy.S(5)/2)*(4 - 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(20*I*A + 23*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) - a**2*(4*A - 7*I*B)*sqrt(I*a*tan(c + d*x) + a)/(4*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_553():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cot(c + d*x))
    F = I*B*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - a**(sympy.S(5)/2)*(4 + 4*I)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - (-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*(46*A - 45*I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(8*d) - a**2*(2*A - 3*I*B)*sqrt(I*a*tan(c + d*x) + a)/(4*d*cot(c + d*x)**(sympy.S(3)/2)) + a**2*(18*I*A + 19*B)*sqrt(I*a*tan(c + d*x) + a)/(8*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_554():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = (A + I*B)*cot(c + d*x)**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) - (5*A + 3*I*B)*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) + (7*I*A - 9*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(3*a*d) + (sympy.S.Half + I/2)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_555():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = (A + I*B)*sqrt(cot(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) - (3*A + I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(a*d) + (sympy.S.Half + I/2)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_556():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = (A + I*B)/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S.Half - I/2)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_557():
    f = (A + B*tan(c + d*x))/(sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x)))
    F = -2*(-1)**(sympy.S(1)/4)*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) + (I*A - B)/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - (sympy.S.Half + I/2)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_558():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (A + I*B)*sqrt(cot(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (11*A + 5*I*B)*sqrt(cot(c + d*x))/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - (25*A + 7*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(6*a**2*d) + (sympy.S(1)/4 + I/4)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_559():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = (A + I*B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + (7*A + I*B)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/4 - I/4)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_560():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x)))
    F = (I*A - B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + (I*A + 5*B)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - (sympy.S(1)/4 + I/4)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_561():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*(-1)**(sympy.S(3)/4)*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (I*A - B)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + (A + 3*I*B)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/4 + I/4)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_562():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (A + I*B)*sqrt(cot(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + (17*A + 7*I*B)*sqrt(cot(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + (151*A + 41*I*B)*sqrt(cot(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - (317*A + 67*I*B)*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(60*a**3*d) + (sympy.S(1)/8 + I/8)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_563():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = (A + I*B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))) + (13*A + 3*I*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + (67*A - 3*I*B)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/8 - I/8)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_564():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x)))
    F = (I*A - B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))) + (3*I*A + 7*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) - (3*I*A - 13*B)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - (sympy.S(1)/8 + I/8)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_565():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = (I*A - B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2)) + (A + 11*I*B)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + (13*A - 37*I*B)/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/8 + I/8)*(I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_566():
    f = (A + B*tan(c + d*x))/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = 2*(-1)**(sympy.S(1)/4)*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (I*A - B)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)) + (A + 3*I*B)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) - (I*A - 7*B)/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/8 + I/8)*(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_567():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**m
    F = I*B*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(m - 1)*hyper((1 - m, 1 - n), (2 - m,), -I*tan(c + d*x))/(d*(1 - m)*(I*tan(c + d*x) + 1)**n) + (A - I*B)*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(m - 1)*appellf1(1 - m, 1, 1 - n, 2 - m, I*tan(c + d*x), -I*tan(c + d*x))/(d*(1 - m)*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_568():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - (2 - 4*n)*(-2*A*n + 3*I*B)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(3*d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) - (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) - (4*I*A*n + 6*B)*(I*a*tan(c + d*x) + a)**n*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_569():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*I*A*(1 - 2*n)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) - 2*A*(I*a*tan(c + d*x) + a)**n*sqrt(cot(c + d*x))/d + (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_570():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n*sqrt(cot(c + d*x))
    F = 2*I*B*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) + (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_571():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/sqrt(cot(c + d*x))
    F = 2*B*(I*a*tan(c + d*x) + a)**n/(d*(2*n + 1)*sqrt(cot(c + d*x))) - (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) + (2*I*A*(2*n + 1) + 4*B*n)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_572():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/cot(c + d*x)**(sympy.S(3)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**n/(d*(2*n + 3)*cot(c + d*x)**(sympy.S(3)/2)) - (2*A - 2*I*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) - (-2*A*(2*n + 3) + 4*I*B*n)*(I*a*tan(c + d*x) + a)**n/(d*(2*n + 1)*(2*n + 3)*sqrt(cot(c + d*x))) + (I*a*tan(c + d*x) + a)**n*(4*A*n*(2*n + 3) - 2*I*B*(4*n**2 + 6*n + 3))*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_573():
    f = (A + B*tan(c + d*x))*(I*a*tan(c + d*x) + a)**n/cot(c + d*x)**(sympy.S(5)/2)
    F = 2*B*(I*a*tan(c + d*x) + a)**n/(d*(2*n + 5)*cot(c + d*x)**(sympy.S(5)/2)) + (2*I*A + 2*B)*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x))) - (-2*A*(2*n + 5) + 4*I*B*n)*(I*a*tan(c + d*x) + a)**n/(d*(2*n + 3)*(2*n + 5)*cot(c + d*x)**(sympy.S(3)/2)) - (I*a*tan(c + d*x) + a)**n*(4*I*A*n*(2*n + 5) + 2*B*(4*n**2 + 10*n + 15))/(d*(2*n + 1)*(2*n + 3)*(2*n + 5)*sqrt(cot(c + d*x))) - (2*I*A*(8*n**3 + 32*n**2 + 36*n + 15) + 8*B*n*(2*n**2 + 8*n + 9))*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), -I*tan(c + d*x))/(d*(2*n + 1)*(2*n + 3)*(2*n + 5)*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_574():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - (2*A*b + 2*B*a)*sqrt(cot(c + d*x))/d + sqrt(2)*(a*(A - B) - b*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A - B) - b*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_575():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(cot(c + d*x))/d + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_576():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))*sqrt(cot(c + d*x))
    F = 2*B*b/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a*(A - B) - b*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a*(A - B) - b*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_577():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))/sqrt(cot(c + d*x))
    F = 2*B*b/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + (2*A*b + 2*B*a)/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a*(A - B) - b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_578():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a*cot(c + d*x) + b)*cot(c + d*x)**(sympy.S(3)/2)/(5*d) - 2*a*(7*A*b + 5*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*sqrt(cot(c + d*x))/d - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_579():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(a*cot(c + d*x) + b)*sqrt(cot(c + d*x))/(3*d) - 2*a*(5*A*b + 3*B*a)*sqrt(cot(c + d*x))/(3*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_580():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a**2*sqrt(cot(c + d*x))/d + 2*B*b**2/(d*sqrt(cot(c + d*x))) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_581():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2*sqrt(cot(c + d*x))
    F = 2*B*b**2/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*b*(A*b + 2*B*a)/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_582():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**2/sqrt(cot(c + d*x))
    F = 2*B*b**2/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*b*(A*b + 2*B*a)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2*(A - B) - 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2*(A + B) + 2*a*b*(A - B) - b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + (4*A*a*b + 2*B*a**2 - 2*B*b**2)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_583():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(a*cot(c + d*x) + b)**2*cot(c + d*x)**(sympy.S(3)/2)/(7*d) - 2*a**2*(11*A*b + 7*B*a)*cot(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*a*(7*A*a**2 - 18*A*b**2 - 21*B*a*b)*cot(c + d*x)**(sympy.S(3)/2)/(21*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_584():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a*cot(c + d*x) + b)**2*sqrt(cot(c + d*x))/(5*d) - 2*a**2*(9*A*b + 5*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + 2*a*(5*A*a**2 - 14*A*b**2 - 15*B*a*b)*sqrt(cot(c + d*x))/(5*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_585():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*(a*cot(c + d*x) + b)**2/(d*sqrt(cot(c + d*x))) - 2*a**2*(A*a + 3*B*b)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*a*(3*A*a*b + B*a**2 + 2*B*b**2)*sqrt(cot(c + d*x))/d + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_586():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a*cot(c + d*x) + b)**2/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(3*A*a + B*b)*sqrt(cot(c + d*x))/(3*d) + 2*b**2*(3*A*b + 7*B*a)/(3*d*sqrt(cot(c + d*x))) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_587():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3*sqrt(cot(c + d*x))
    F = 2*B*b*(a*cot(c + d*x) + b)**2/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*b**2*(5*A*b + 9*B*a)/(15*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*b*(15*A*a*b + 14*B*a**2 - 5*B*b**2)/(5*d*sqrt(cot(c + d*x))) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_588():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**3/sqrt(cot(c + d*x))
    F = 2*B*b*(a*cot(c + d*x) + b)**2/(7*d*cot(c + d*x)**(sympy.S(7)/2)) + 2*b**2*(7*A*b + 11*B*a)/(35*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*b*(21*A*a*b + 18*B*a**2 - 7*B*b**2)/(21*d*cot(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**3*(A - B) - 3*a**2*b*(A + B) - 3*a*b**2*(A - B) + b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**3*(A + B) + 3*a**2*b*(A - B) - 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_589():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = -2*A*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) + sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + (2*A*b - 2*B*a)*sqrt(cot(c + d*x))/(a**2*d) - 2*b**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(5)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_590():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = -2*A*sqrt(cot(c + d*x))/(a*d) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + 2*b**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(3)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_591():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))
    F = -sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - 2*sqrt(b)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(a)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_592():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*sqrt(cot(c + d*x)))
    F = 2*sqrt(a)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(b)*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_593():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*B/(b*d*sqrt(cot(c + d*x))) - 2*a**(sympy.S(3)/2)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a*(A - B) + b*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_594():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2))
    F = 2*B/(3*b*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*a**(sympy.S(5)/2)*(A*b - B*a)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a*(A - B) + b*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(-a*(A + B) + b*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(-a*(A + B) + b*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + (2*A*b - 2*B*a)/(b**2*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_595():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**2
    F = sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b*(A*b - B*a)*cot(c + d*x)**(sympy.S(3)/2)/(a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - (2*A*a**2 + 3*A*b**2 - B*a*b)*sqrt(cot(c + d*x))/(a**2*d*(a**2 + b**2)) + b**(sympy.S(3)/2)*(7*A*a**2*b + 3*A*b**3 - 5*B*a**3 - B*a*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(5)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_596():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**2
    F = -sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + b*(A*b - B*a)*sqrt(cot(c + d*x))/(a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - sqrt(b)*(5*A*a**2*b + A*b**3 - 3*B*a**3 + B*a*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(3)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_597():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*sqrt(cot(c + d*x)))
    F = -(A*b - B*a)*sqrt(cot(c + d*x))/(d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + (3*A*a**2*b - A*b**3 - B*a**3 + 3*B*a*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(a)*sqrt(b)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_598():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(a)*(A*a**2*b - 3*A*b**3 + B*a**3 + 5*B*a*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)**2) + a*(A*b - B*a)*sqrt(cot(c + d*x))/(b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_599():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(5)/2))
    F = -a**(sympy.S(3)/2)*(A*a**2*b + 5*A*b**3 - 3*B*a**3 - 7*B*a*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)**2) + a*(A*b - B*a)/(b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)*sqrt(cot(c + d*x))) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2*(A - B) + 2*a*b*(A + B) - b**2*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(-a**2*(A + B) + 2*a*b*(A - B) + b**2*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - (A*a*b - 3*B*a**2 - 2*B*b**2)/(b**2*d*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_600():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**3
    F = sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b*(A*b - B*a)*cot(c + d*x)**(sympy.S(5)/2)/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) + b*(13*A*a**2*b + 5*A*b**3 - 9*B*a**3 - B*a*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(4*a**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - (8*A*a**4 + 31*A*a**2*b**2 + 15*A*b**4 - 11*B*a**3*b - 3*B*a*b**3)*sqrt(cot(c + d*x))/(4*a**3*d*(a**2 + b**2)**2) + b**(sympy.S(3)/2)*(63*A*a**4*b + 46*A*a**2*b**3 + 15*A*b**5 - 35*B*a**5 - 6*B*a**3*b**2 - 3*B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(7)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_601():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**3
    F = -sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b*(A*b - B*a)*cot(c + d*x)**(sympy.S(3)/2)/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) + b*(11*A*a**2*b + 3*A*b**3 - 7*B*a**3 + B*a*b**2)*sqrt(cot(c + d*x))/(4*a**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - sqrt(b)*(35*A*a**4*b + 6*A*a**2*b**3 + 3*A*b**5 - 15*B*a**5 + 18*B*a**3*b**2 + B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(5)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_602():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*sqrt(cot(c + d*x)))
    F = -sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b*(A*b - B*a)*sqrt(cot(c + d*x))/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) - (9*A*a**2*b + A*b**3 - 5*B*a**3 + 3*B*a*b**2)*sqrt(cot(c + d*x))/(4*a*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) + (15*A*a**4*b - 18*A*a**2*b**3 - A*b**5 - 3*B*a**5 + 26*B*a**3*b**2 - 3*B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(3)/2)*sqrt(b)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_603():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(3)/2))
    F = -(A*b - B*a)*sqrt(cot(c + d*x))/(d*(2*a**2 + 2*b**2)*(a*cot(c + d*x) + b)**2) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + (5*A*a**2*b - 3*A*b**3 - B*a**3 + 7*B*a*b**2)*sqrt(cot(c + d*x))/(4*b*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - (3*A*a**4*b - 26*A*a**2*b**3 + 3*A*b**5 + B*a**5 + 18*B*a**3*b**2 - 15*B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*sqrt(a)*b**(sympy.S(3)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_604():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(a)*(A*a**4*b + 18*A*a**2*b**3 - 15*A*b**5 + 3*B*a**5 + 6*B*a**3*b**2 + 35*B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*b**(sympy.S(5)/2)*d*(a**2 + b**2)**3) + a*(A*b - B*a)*sqrt(cot(c + d*x))/(2*b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) - a*(A*a**2*b - 7*A*b**3 + 3*B*a**3 + 11*B*a*b**2)*sqrt(cot(c + d*x))/(4*b**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_605():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(7)/2))
    F = -a**(sympy.S(3)/2)*(3*A*a**4*b + 6*A*a**2*b**3 + 35*A*b**5 - 15*B*a**5 - 46*B*a**3*b**2 - 63*B*a*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*b**(sympy.S(7)/2)*d*(a**2 + b**2)**3) + a*(A*b - B*a)/(2*b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2*sqrt(cot(c + d*x))) + a*(A*a**2*b + 9*A*b**3 - 5*B*a**3 - 13*B*a*b**2)/(4*b**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)*sqrt(cot(c + d*x))) - sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a**3*(A - B) + 3*a**2*b*(A + B) - 3*a*b**2*(A - B) - b**3*(A + B))*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(-a**3*(A + B) + 3*a**2*b*(A - B) + 3*a*b**2*(A + B) - b**3*(A - B))*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - (3*A*a**3*b + 11*A*a*b**3 - 15*B*a**4 - 31*B*a**2*b**2 - 8*B*b**4)/(4*b**3*d*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_606():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - 2*B*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_607():
    f = (B*a + B*b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - 2*B*sqrt(cot(c + d*x))/d + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_608():
    f = (B*a + B*b*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_609():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*sqrt(cot(c + d*x)))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_610():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2))
    F = sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + 2*B/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_611():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(2)*B*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*B*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*B*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + 2*B/(3*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_612():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - (I*A - B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 14*B*a)*cot(c + d*x)**(sympy.S(5)/2)/(35*a*d) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 + 8*A*b**2 - 14*B*a*b)*cot(c + d*x)**(sympy.S(3)/2)/(105*a**2*d) + sqrt(a + b*tan(c + d*x))*(70*A*a**2*b - 16*A*b**3 + 210*B*a**3 + 28*B*a*b**2)*sqrt(cot(c + d*x))/(105*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_613():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - (A - I*B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 10*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*a*d) + sqrt(a + b*tan(c + d*x))*(30*A*a**2 + 4*A*b**2 - 10*B*a*b)*sqrt(cot(c + d*x))/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_614():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + (I*A - B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(2*A*b + 6*B*a)*sqrt(cot(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_615():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/d + (A - I*B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_616():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))
    F = 2*B*sqrt(b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A - B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_617():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/sqrt(cot(c + d*x))
    F = B*sqrt(a + b*tan(c + d*x))/(d*sqrt(cot(c + d*x))) - (A - I*B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (2*A*b + B*a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_618():
    f = (A + B*tan(c + d*x))*sqrt(a + b*tan(c + d*x))/cot(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(2*b*d*sqrt(cot(c + d*x))) + (I*A - B)*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(4*A*b - B*a)/(4*b*d*sqrt(cot(c + d*x))) + (4*A*a*b - B*a**2 - 8*B*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_619():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(9)/2)/(9*d) - sqrt(a + b*tan(c + d*x))*(20*A*b + 18*B*a)*cot(c + d*x)**(sympy.S(7)/2)/(63*d) + (I*A - B)*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(42*A*a**2 - 2*A*b**2 - 48*B*a*b)*cot(c + d*x)**(sympy.S(5)/2)/(105*a*d) + sqrt(a + b*tan(c + d*x))*(252*A*a**2*b + 8*A*b**3 + 210*B*a**3 - 18*B*a*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(315*a**2*d) - sqrt(a + b*tan(c + d*x))*(630*A*a**4 - 126*A*a**2*b**2 + 16*A*b**4 - 840*B*a**3*b - 36*B*a*b**3)*sqrt(cot(c + d*x))/(315*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_620():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - (A - I*B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(16*A*b + 14*B*a)*cot(c + d*x)**(sympy.S(5)/2)/(35*d) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 - 6*A*b**2 - 84*B*a*b)*cot(c + d*x)**(sympy.S(3)/2)/(105*a*d) + sqrt(a + b*tan(c + d*x))*(280*A*a**2*b + 12*A*b**3 + 210*B*a**3 - 42*B*a*b**2)*sqrt(cot(c + d*x))/(105*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_621():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + (a + I*b)**2*(I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - sqrt(a + b*tan(c + d*x))*(12*A*b + 10*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + (I*A + B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 6*A*b**2 - 40*B*a*b)*sqrt(cot(c + d*x))/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_622():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + (A - I*B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(a + b*tan(c + d*x))*(8*A*b + 6*B*a)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_623():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/d + 2*B*b**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (a + I*b)**2*(I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_624():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x))
    F = B*b*sqrt(a + b*tan(c + d*x))/(d*sqrt(cot(c + d*x))) + sqrt(b)*(2*A*b + 3*B*a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A - I*B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_625():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/sqrt(cot(c + d*x))
    F = B*b*sqrt(a + b*tan(c + d*x))/(2*d*cot(c + d*x)**(sympy.S(3)/2)) + (a + I*b)**2*(I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*(4*A*b + 5*B*a)/(4*d*sqrt(cot(c + d*x))) + (I*A + B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (12*A*a*b + 3*B*a**2 - 8*B*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_626():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(3)/2)/cot(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(3*b*d*sqrt(cot(c + d*x))) + (A - I*B)*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(6*A*b - B*a)/(12*b*d*sqrt(cot(c + d*x))) + sqrt(a + b*tan(c + d*x))*(6*A*a*b - B*a**2 - 8*B*b**2)/(8*b*d*sqrt(cot(c + d*x))) + (6*A*a**2*b - 16*A*b**3 - B*a**3 - 24*B*a*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_627():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(13)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(11)/2)/(11*d) - 2*a*sqrt(a + b*tan(c + d*x))*(14*A*b + 11*B*a)*cot(c + d*x)**(sympy.S(9)/2)/(99*d) + sqrt(a + b*tan(c + d*x))*(198*A*a**2 - 226*A*b**2 - 418*B*a*b)*cot(c + d*x)**(sympy.S(7)/2)/(693*d) - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(990*A*a**2*b - 10*A*b**3 + 462*B*a**3 - 550*B*a*b**2)*cot(c + d*x)**(sympy.S(5)/2)/(1155*a*d) - sqrt(a + b*tan(c + d*x))*(2310*A*a**4 - 2970*A*a**2*b**2 - 40*A*b**4 - 5082*B*a**3*b + 110*B*a*b**3)*cot(c + d*x)**(sympy.S(3)/2)/(3465*a**2*d) - sqrt(a + b*tan(c + d*x))*(16170*A*a**4*b - 990*A*a**2*b**3 + 80*A*b**5 + 6930*B*a**5 - 10626*B*a**3*b**2 - 220*B*a*b**4)*sqrt(cot(c + d*x))/(3465*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_628():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(11)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(9)/2)/(9*d) - 2*a*sqrt(a + b*tan(c + d*x))*(4*A*b + 3*B*a)*cot(c + d*x)**(sympy.S(7)/2)/(21*d) - (A - I*B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(42*A*a**2 - 50*A*b**2 - 90*B*a*b)*cot(c + d*x)**(sympy.S(5)/2)/(105*d) + sqrt(a + b*tan(c + d*x))*(462*A*a**2*b - 10*A*b**3 + 210*B*a**3 - 270*B*a*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(315*a*d) - sqrt(a + b*tan(c + d*x))*(630*A*a**4 - 966*A*a**2*b**2 - 20*A*b**4 - 1470*B*a**3*b + 90*B*a*b**3)*sqrt(cot(c + d*x))/(315*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_629():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - 2*a*sqrt(a + b*tan(c + d*x))*(10*A*b + 7*B*a)*cot(c + d*x)**(sympy.S(5)/2)/(35*d) + sqrt(a + b*tan(c + d*x))*(70*A*a**2 - 90*A*b**2 - 154*B*a*b)*cot(c + d*x)**(sympy.S(3)/2)/(105*d) + (I*A - B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(490*A*a**2*b - 30*A*b**3 + 210*B*a**3 - 322*B*a*b**2)*sqrt(cot(c + d*x))/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_630():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 2*a*sqrt(a + b*tan(c + d*x))*(8*A*b + 5*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + (A - I*B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 46*A*b**2 - 70*B*a*b)*sqrt(cot(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_631():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*B*b**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*a*sqrt(a + b*tan(c + d*x))*(2*A*b + B*a)*sqrt(cot(c + d*x))/d - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_632():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*A*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x))/d + b**(sympy.S(3)/2)*(2*A*b + 5*B*a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b*sqrt(a + b*tan(c + d*x))*(2*A*a + B*b)/(d*sqrt(cot(c + d*x))) - (A - I*B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (A + I*B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_633():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(cot(c + d*x))
    F = B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(2*d*sqrt(cot(c + d*x))) + sqrt(b)*(20*A*a*b + 15*B*a**2 - 8*B*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*d) + b*sqrt(a + b*tan(c + d*x))*(4*A*b + 7*B*a)/(4*d*sqrt(cot(c + d*x))) + (I*A - B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*A + B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_634():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/sqrt(cot(c + d*x))
    F = B*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + (A - I*B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (A + I*B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(2*A*b + 3*B*a)/(4*d*sqrt(cot(c + d*x))) + sqrt(a + b*tan(c + d*x))*(14*A*a*b + 5*B*a**2 - 8*B*b**2)/(8*d*sqrt(cot(c + d*x))) + (30*A*a**2*b - 16*A*b**3 + 5*B*a**3 - 40*B*a*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_635():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**(sympy.S(5)/2)/cot(c + d*x)**(sympy.S(3)/2)
    F = B*(a + b*tan(c + d*x))**(sympy.S(7)/2)/(4*b*d*sqrt(cot(c + d*x))) - (I*A - B)*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*A + B)*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(5)/2)*(8*A*b - B*a)/(24*b*d*sqrt(cot(c + d*x))) + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(40*A*a*b - 5*B*a**2 - 48*B*b**2)/(96*b*d*sqrt(cot(c + d*x))) + sqrt(a + b*tan(c + d*x))*(40*A*a**2*b - 64*A*b**3 - 5*B*a**3 - 112*B*a*b**2)/(64*b*d*sqrt(cot(c + d*x))) + (40*A*a**3*b - 320*A*a*b**3 - 5*B*a**4 - 240*B*a**2*b**2 + 128*B*b**4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(64*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_636():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/sqrt(a + b*tan(c + d*x))
    F = -2*A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*a*d) + (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + sqrt(a + b*tan(c + d*x))*(8*A*b - 10*B*a)*cot(c + d*x)**(sympy.S(3)/2)/(15*a**2*d) + sqrt(a + b*tan(c + d*x))*(30*A*a**2 - 16*A*b**2 + 20*B*a*b)*sqrt(cot(c + d*x))/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_637():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*tan(c + d*x))
    F = -2*A*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) - (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*(4*A*b - 6*B*a)*sqrt(cot(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_638():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*tan(c + d*x))
    F = -2*A*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(a*d) - (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_639():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_640():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x)))
    F = 2*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d) + (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_641():
    f = (A + B*tan(c + d*x))/(sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2))
    F = B*sqrt(a + b*tan(c + d*x))/(b*d*sqrt(cot(c + d*x))) - (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + (2*A*b - B*a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_642():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*A*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d*sqrt(a + b*tan(c + d*x))) - (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + (8*A*b - 6*B*a)*sqrt(cot(c + d*x))/(3*a**2*d*sqrt(a + b*tan(c + d*x))) + 2*b*(5*A*a**2*b + 8*A*b**3 - 3*B*a**3 - 6*B*a*b**2)/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_643():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*A*sqrt(cot(c + d*x))/(a*d*sqrt(a + b*tan(c + d*x))) - (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2*b*(A*a**2 + 2*A*b**2 - B*a*b)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_644():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + 2*b*(A*b - B*a)/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_645():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x)))
    F = (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (2*A*b - 2*B*a)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_646():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d) + 2*a*(A*b - B*a)/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x))) - (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_647():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*A*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*A*b - 2*B*a)*sqrt(cot(c + d*x))/(a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + 2*b*(7*A*a**2*b + 8*A*b**3 - 3*B*a**3 - 4*B*a*b**2)/(3*a**3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + 2*b*(8*A*a**4*b + 30*A*a**2*b**3 + 16*A*b**5 - 3*B*a**5 - 17*B*a**3*b**2 - 8*B*a*b**4)/(3*a**4*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_648():
    f = (A + B*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*A*sqrt(cot(c + d*x))/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - 2*b*(3*A*a**2 + 4*A*b**2 - B*a*b)/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) - 2*b*(3*A*a**4 + 17*A*a**2*b**2 + 8*A*b**4 - 8*B*a**3*b - 2*B*a*b**3)/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_649():
    f = (A + B*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -(A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + 2*b*(A*b - B*a)/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + 2*b*(8*A*a**2*b + 2*A*b**3 - 5*B*a**3 + B*a*b**2)/(3*a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_650():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(cot(c + d*x)))
    F = -(I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - (2*A*b - 2*B*a)/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*sqrt(cot(c + d*x))) - (10*A*a**2*b - 2*A*b**3 - 4*B*a**3 + 8*B*a*b**2)/(3*a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_651():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*a*(A*b - B*a)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + (A - I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + (A + I*B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*A*a**2*b - 8*A*b**3 + 2*B*a**3 + 14*B*a*b**2)/(3*b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_652():
    f = (A + B*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = 2*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(5)/2)*d) + 2*a*(A*b - B*a)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*cot(c + d*x)**(sympy.S(3)/2)) + 2*a*(2*A*b**3 - B*a*(a**2 + 3*b**2))/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x))) + (I*A - B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - (I*A + B)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_653():
    f = (B*a + B*b*tan(c + d*x))*sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_654():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x)))
    F = -I*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + I*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_655():
    f = (B*a + B*b*tan(c + d*x))/((a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = -B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + 2*B*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_656():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*cot(c + d*x)**m
    F = (A - I*B)*(a + b*tan(c + d*x))**n*cot(c + d*x)**(m - 1)*appellf1(1 - m, 1, -n, 2 - m, I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*d*(1 - m)*(1 + b*tan(c + d*x)/a)**n) + (A + I*B)*(a + b*tan(c + d*x))**n*cot(c + d*x)**(m - 1)*appellf1(1 - m, 1, -n, 2 - m, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*d*(1 - m)*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_657():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*cot(c + d*x)**(sympy.S(3)/2)
    F = -(A - I*B)*(a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n) - (A + I*B)*(a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_658():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))
    F = (A - I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(cot(c + d*x))) + (A + I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_659():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n/sqrt(cot(c + d*x))
    F = (A - I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(3)/2)) + (A + I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_660():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n/cot(c + d*x)**(sympy.S(3)/2)
    F = (A - I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(5)/2)) + (A + I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_661():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)
    F = (A - I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n) + (A + I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_662():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))
    F = (A - I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n) + (A + I*B)*(a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_663():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n/sqrt(tan(c + d*x))
    F = (A - I*B)*(a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n) + (A + I*B)*(a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_664():
    f = (A + B*tan(c + d*x))*(a + b*tan(c + d*x))**n/tan(c + d*x)**(sympy.S(3)/2)
    F = -(A - I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(tan(c + d*x))) - (A + I*B)*(a + b*tan(c + d*x))**n*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_665():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**n
    F = -B*a*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1)) + a*(I*A + B)*(-I*c*tan(e + f*x) + c)**n/(f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_666():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**4
    F = -B*a*c**4*(-I*tan(e + f*x) + 1)**5/(5*f) + a*c**4*(I*A + B)*(-I*tan(e + f*x) + 1)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_667():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**3
    F = -B*a*c**3*(-I*tan(e + f*x) + 1)**4/(4*f) + a*c**3*(I*A + B)*(-I*tan(e + f*x) + 1)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_668():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**2
    F = A*a*c**2*tan(e + f*x)/f - I*B*a*c**2*tan(e + f*x)**3/(3*f) - a*c**2*(I*A - B)*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_669():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)
    F = A*a*c*tan(e + f*x)/f + B*a*c*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_670():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)
    F = I*B*a*tan(e + f*x)/f + a*x*(A - I*B) - a*(I*A + B)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_671():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)
    F = I*B*a*x/c + B*a*log(cos(e + f*x))/(c*f) + a*(A - I*B)/(c*f*(tan(e + f*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_672():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**2
    F = a*(A + B*tan(e + f*x))**2/(c**2*f*(2*I*A + 2*B)*(-I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_673():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**3
    F = -B*a/(2*c**3*f*(tan(e + f*x) + I)**2) - a*(A - I*B)/(3*c**3*f*(tan(e + f*x) + I)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_674():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**4
    F = -I*B*a/(3*c**4*f*(tan(e + f*x) + I)**3) - a*(I*A + B)/(4*c**4*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_675():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**5
    F = B*a/(4*c**5*f*(tan(e + f*x) + I)**4) + a*(A - I*B)/(5*c**5*f*(tan(e + f*x) + I)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_676():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**n
    F = B*a**2*(-I*c*tan(e + f*x) + c)**(n + 2)/(c**2*f*(n + 2)) + 2*a**2*(I*A + B)*(-I*c*tan(e + f*x) + c)**n/(f*n) - a**2*(I*A + 3*B)*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_677():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**5
    F = B*a**2*c**5*(-I*tan(e + f*x) + 1)**7/(7*f) + 2*a**2*c**5*(I*A + B)*(-I*tan(e + f*x) + 1)**5/(5*f) - a**2*c**5*(I*A + 3*B)*(-I*tan(e + f*x) + 1)**6/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_678():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**4
    F = B*a**2*c**4*(-I*tan(e + f*x) + 1)**6/(6*f) + a**2*c**4*(I*A + B)*(-I*tan(e + f*x) + 1)**4/(2*f) - a**2*c**4*(I*A + 3*B)*(-I*tan(e + f*x) + 1)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_679():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**3
    F = B*a**2*c**3*(-I*tan(e + f*x) + 1)**5/(5*f) + 2*a**2*c**3*(I*A + B)*(-I*tan(e + f*x) + 1)**3/(3*f) - a**2*c**3*(I*A + 3*B)*(-I*tan(e + f*x) + 1)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_680():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**2
    F = A*a**2*c**2*tan(e + f*x)**3/(3*f) + A*a**2*c**2*tan(e + f*x)/f + B*a**2*c**2*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_681():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)
    F = A*a**2*c*tan(e + f*x)/f + I*B*a**2*c*tan(e + f*x)**3/(3*f) + a**2*c*(I*A + B)*tan(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_682():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2
    F = B*(I*a*tan(e + f*x) + a)**2/(2*f) + 2*a**2*x*(A - I*B) - a**2*(A - I*B)*tan(e + f*x)/f - 2*a**2*(I*A + B)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_683():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)
    F = -I*B*a**2*tan(e + f*x)/(c*f) - a**2*x*(A - 3*I*B)/c + 2*a**2*(A - I*B)/(c*f*(tan(e + f*x) + I)) + a**2*(I*A + 3*B)*log(cos(e + f*x))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_684():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**2
    F = -I*B*a**2*x/c**2 - B*a**2*log(cos(e + f*x))/(c**2*f) - a**2*(A - 3*I*B)/(c**2*f*(tan(e + f*x) + I)) + a**2*(I*A + B)/(c**2*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_685():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**3
    F = -I*B*a**2/(c**3*f*(tan(e + f*x) + I)) - 2*a**2*(A - I*B)/(3*c**3*f*(tan(e + f*x) + I)**3) - a**2*(I*A + 3*B)/(2*c**3*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_686():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**4
    F = B*a**2/(2*c**4*f*(tan(e + f*x) + I)**2) + a**2*(A - 3*I*B)/(3*c**4*f*(tan(e + f*x) + I)**3) - a**2*(I*A + B)/(2*c**4*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_687():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**5
    F = I*B*a**2/(3*c**5*f*(tan(e + f*x) + I)**3) + 2*a**2*(A - I*B)/(5*c**5*f*(tan(e + f*x) + I)**5) + a**2*(I*A + 3*B)/(4*c**5*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_688():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**6
    F = -B*a**2/(4*c**6*f*(tan(e + f*x) + I)**4) - a**2*(A - 3*I*B)/(5*c**6*f*(tan(e + f*x) + I)**5) + a**2*(I*A + B)/(3*c**6*f*(tan(e + f*x) + I)**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_689():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**n
    F = -B*a**3*(-I*c*tan(e + f*x) + c)**(n + 3)/(c**3*f*(n + 3)) + 4*a**3*(I*A + B)*(-I*c*tan(e + f*x) + c)**n/(f*n) - 4*a**3*(I*A + 2*B)*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1)) + a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(n + 2)/(c**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_690():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**6
    F = -B*a**3*c**6*(-I*tan(e + f*x) + 1)**9/(9*f) + 2*a**3*c**6*(I*A + B)*(-I*tan(e + f*x) + 1)**6/(3*f) - 4*a**3*c**6*(I*A + 2*B)*(-I*tan(e + f*x) + 1)**7/(7*f) + a**3*c**6*(I*A + 5*B)*(-I*tan(e + f*x) + 1)**8/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_691():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**5
    F = -B*a**3*c**5*(-I*tan(e + f*x) + 1)**8/(8*f) + 4*a**3*c**5*(I*A + B)*(-I*tan(e + f*x) + 1)**5/(5*f) - 2*a**3*c**5*(I*A + 2*B)*(-I*tan(e + f*x) + 1)**6/(3*f) + a**3*c**5*(I*A + 5*B)*(-I*tan(e + f*x) + 1)**7/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_692():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**4
    F = -B*a**3*c**4*(-I*tan(e + f*x) + 1)**7/(7*f) + a**3*c**4*(I*A + B)*(-I*tan(e + f*x) + 1)**4/f - 4*a**3*c**4*(I*A + 2*B)*(-I*tan(e + f*x) + 1)**5/(5*f) + a**3*c**4*(I*A + 5*B)*(-I*tan(e + f*x) + 1)**6/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_693():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**3
    F = A*a**3*c**3*tan(e + f*x)**5/(5*f) + 2*A*a**3*c**3*tan(e + f*x)**3/(3*f) + A*a**3*c**3*tan(e + f*x)/f + B*a**3*c**3*sec(e + f*x)**6/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_694():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**2
    F = B*a**3*c**2*(I*tan(e + f*x) + 1)**5/(5*f) + a**3*c**2*(I*A - 3*B)*(I*tan(e + f*x) + 1)**4/(4*f) - 2*a**3*c**2*(I*A - B)*(I*tan(e + f*x) + 1)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_695():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)
    F = -B*a**3*c*(I*tan(e + f*x) + 1)**4/(4*f) - a**3*c*(I*A - B)*(I*tan(e + f*x) + 1)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_696():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3
    F = B*(I*a*tan(e + f*x) + a)**3/(3*f) + 4*a**3*x*(A - I*B) - 2*a**3*(A - I*B)*tan(e + f*x)/f - 4*a**3*(I*A + B)*log(cos(e + f*x))/f + a*(I*A + B)*(I*a*tan(e + f*x) + a)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_697():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)
    F = B*a**3*tan(e + f*x)**2/(2*c*f) - 4*a**3*x*(A - 2*I*B)/c + a**3*(A - 4*I*B)*tan(e + f*x)/(c*f) + 4*a**3*(A - I*B)/(c*f*(tan(e + f*x) + I)) + 4*a**3*(I*A + 2*B)*log(cos(e + f*x))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_698():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**2
    F = I*B*a**3*tan(e + f*x)/(c**2*f) + a**3*x*(A - 5*I*B)/c**2 - 4*a**3*(A - 2*I*B)/(c**2*f*(tan(e + f*x) + I)) + 2*a**3*(I*A + B)/(c**2*f*(tan(e + f*x) + I)**2) - a**3*(I*A + 5*B)*log(cos(e + f*x))/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_699():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**3
    F = I*B*a**3*x/c**3 + B*a**3*log(cos(e + f*x))/(c**3*f) - 4*I*B*a**3/(c**3*f*(tan(e + f*x) + I)) - 2*B*a**3/(c**3*f*(tan(e + f*x) + I)**2) - a**3*(I*A + B)*(I*tan(e + f*x) + 1)**3/(6*c**3*f*(-I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_700():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**4
    F = -a**3*(I*A - 7*B)*(I*tan(e + f*x) + 1)**3/(48*c**4*f*(-I*tan(e + f*x) + 1)**3) - a**3*(I*A + B)*(I*tan(e + f*x) + 1)**3/(8*c**4*f*(-I*tan(e + f*x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_701():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**5
    F = -B*a**3/(2*c**5*f*(tan(e + f*x) + I)**2) - a**3*(A - 5*I*B)/(3*c**5*f*(tan(e + f*x) + I)**3) + 4*a**3*(A - I*B)/(5*c**5*f*(tan(e + f*x) + I)**5) + a**3*(I*A + 2*B)/(c**5*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_702():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**6
    F = -I*B*a**3/(3*c**6*f*(tan(e + f*x) + I)**3) - 4*a**3*(A - 2*I*B)/(5*c**6*f*(tan(e + f*x) + I)**5) + 2*a**3*(I*A + B)/(3*c**6*f*(tan(e + f*x) + I)**6) - a**3*(I*A + 5*B)/(4*c**6*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_703():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**7
    F = B*a**3/(4*c**7*f*(tan(e + f*x) + I)**4) + a**3*(A - 5*I*B)/(5*c**7*f*(tan(e + f*x) + I)**5) - 4*a**3*(A - I*B)/(7*c**7*f*(tan(e + f*x) + I)**7) - 2*a**3*(I*A + 2*B)/(3*c**7*f*(tan(e + f*x) + I)**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_704():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**8
    F = I*B*a**3/(5*c**8*f*(tan(e + f*x) + I)**5) + 4*a**3*(A - 2*I*B)/(7*c**8*f*(tan(e + f*x) + I)**7) - a**3*(I*A + B)/(2*c**8*f*(tan(e + f*x) + I)**8) + a**3*(I*A + 5*B)/(6*c**8*f*(tan(e + f*x) + I)**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_705():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**n/(2*a*f*(I*tan(e + f*x) + 1)) + (I*A*(1 - n) + B*(n + 1))*(-I*c*tan(e + f*x) + c)**n*hyper((1, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(4*a*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_706():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)
    F = -I*B*c**4*tan(e + f*x)**3/(3*a*f) + c**4*x*(-12*A - 20*I*B)/a + c**4*(5*A + 12*I*B)*tan(e + f*x)/(a*f) - c**4*(8*A + 8*I*B)/(a*f*(-tan(e + f*x) + I)) - c**4*(I*A - 5*B)*tan(e + f*x)**2/(2*a*f) - c**4*(12*I*A - 20*B)*log(cos(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_707():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)
    F = B*c**3*tan(e + f*x)**2/(2*a*f) + c**3*x*(-4*A - 8*I*B)/a + c**3*(A + 4*I*B)*tan(e + f*x)/(a*f) - c**3*(4*A + 4*I*B)/(a*f*(-tan(e + f*x) + I)) - c**3*(4*I*A - 8*B)*log(cos(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_708():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)
    F = I*B*c**2*tan(e + f*x)/(a*f) - c**2*x*(A + 3*I*B)/a - c**2*(2*A + 2*I*B)/(a*f*(-tan(e + f*x) + I)) - c**2*(I*A - 3*B)*log(cos(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_709():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)
    F = -I*B*c*x/a + B*c*log(cos(e + f*x))/(a*f) - c*(A + I*B)/(a*f*(-tan(e + f*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_710():
    f = (A + B*tan(e + f*x))/(I*a*tan(e + f*x) + a)
    F = (I*A - B)/(2*f*(I*a*tan(e + f*x) + a)) + x*(A - I*B)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_711():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c))
    F = A*x/(2*a*c) - (-A*tan(e + f*x) + B)*cos(e + f*x)**2/(2*a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_712():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**2)
    F = A/(4*a*c**2*f*(tan(e + f*x) + I)) + x*(3*A + I*B)/(8*a*c**2) - (A + I*B)/(8*a*c**2*f*(-tan(e + f*x) + I)) + (I*A + B)/(8*a*c**2*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_713():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**3)
    F = I*A/(8*a*c**3*f*(tan(e + f*x) + I)**2) + x*(2*A + I*B)/(8*a*c**3) - (A - I*B)/(12*a*c**3*f*(tan(e + f*x) + I)**3) - (A + I*B)/(16*a*c**3*f*(-tan(e + f*x) + I)) + (3*A + I*B)/(16*a*c**3*f*(tan(e + f*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_714():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**4)
    F = -A/(12*a*c**4*f*(tan(e + f*x) + I)**3) + x*(5*A + 3*I*B)/(32*a*c**4) - (A + I*B)/(32*a*c**4*f*(-tan(e + f*x) + I)) + (2*A + I*B)/(16*a*c**4*f*(tan(e + f*x) + I)) - (I*A + B)/(16*a*c**4*f*(tan(e + f*x) + I)**4) + (3*I*A - B)/(32*a*c**4*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_715():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)**2
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**n/(4*a**2*f*(I*tan(e + f*x) + 1)**2) + (I*A*(2 - n) + B*(n + 2))*(-I*c*tan(e + f*x) + c)**n*hyper((2, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(16*a**2*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_716():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**5/(I*a*tan(e + f*x) + a)**2
    F = I*B*c**5*tan(e + f*x)**3/(3*a**2*f) + c**5*x*(24*A + 56*I*B)/a**2 - c**5*(7*A + 24*I*B)*tan(e + f*x)/(a**2*f) + c**5*(32*A + 48*I*B)/(a**2*f*(-tan(e + f*x) + I)) + c**5*(I*A - 7*B)*tan(e + f*x)**2/(2*a**2*f) - c**5*(8*I*A - 8*B)/(a**2*f*(-tan(e + f*x) + I)**2) + c**5*(24*I*A - 56*B)*log(cos(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_717():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**2
    F = -B*c**4*tan(e + f*x)**2/(2*a**2*f) + c**4*x*(6*A + 18*I*B)/a**2 - c**4*(A + 6*I*B)*tan(e + f*x)/(a**2*f) + c**4*(12*A + 20*I*B)/(a**2*f*(-tan(e + f*x) + I)) - c**4*(4*I*A - 4*B)/(a**2*f*(-tan(e + f*x) + I)**2) + c**4*(6*I*A - 18*B)*log(cos(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_718():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**2
    F = -I*B*c**3*tan(e + f*x)/(a**2*f) + c**3*x*(A + 5*I*B)/a**2 + c**3*(4*A + 8*I*B)/(a**2*f*(-tan(e + f*x) + I)) + c**3*(I*A - 5*B)*log(cos(e + f*x))/(a**2*f) - c**3*(2*I*A - 2*B)/(a**2*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_719():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)**2
    F = I*B*c**2*x/a**2 - B*c**2*log(cos(e + f*x))/(a**2*f) + c**2*(A + 3*I*B)/(a**2*f*(-tan(e + f*x) + I)) - c**2*(I*A - B)/(a**2*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_720():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**2
    F = -c*(A + B*tan(e + f*x))**2/(2*a**2*f*(I*A - B)*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_721():
    f = (A + B*tan(e + f*x))/(I*a*tan(e + f*x) + a)**2
    F = (I*A - B)/(4*f*(I*a*tan(e + f*x) + a)**2) + (I*A + B)/(4*f*(I*a**2*tan(e + f*x) + a**2)) + x*(A - I*B)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_722():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c))
    F = -A/(4*a**2*c*f*(-tan(e + f*x) + I)) + x*(3*A - I*B)/(8*a**2*c) + (A - I*B)/(8*a**2*c*f*(tan(e + f*x) + I)) - (I*A - B)/(8*a**2*c*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_723():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**2)
    F = 3*A*x/(8*a**2*c**2) + 3*A*sin(e + f*x)*cos(e + f*x)/(8*a**2*c**2*f) - (-A*tan(e + f*x) + B)*cos(e + f*x)**4/(4*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_724():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**3)
    F = 3*A/(16*a**2*c**3*f*(tan(e + f*x) + I)) + x*(5*A + I*B)/(16*a**2*c**3) - (A - I*B)/(24*a**2*c**3*f*(tan(e + f*x) + I)**3) - (2*A + I*B)/(16*a**2*c**3*f*(-tan(e + f*x) + I)) - (I*A - B)/(32*a**2*c**3*f*(-tan(e + f*x) + I)**2) + (3*I*A + B)/(32*a**2*c**3*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_725():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**4)
    F = 3*I*A/(32*a**2*c**4*f*(tan(e + f*x) + I)**2) + x*(15*A + 5*I*B)/(64*a**2*c**4) - (3*A - I*B)/(48*a**2*c**4*f*(tan(e + f*x) + I)**3) + (5*A + I*B)/(32*a**2*c**4*f*(tan(e + f*x) + I)) - (5*A + 3*I*B)/(64*a**2*c**4*f*(-tan(e + f*x) + I)) - (I*A - B)/(64*a**2*c**4*f*(-tan(e + f*x) + I)**2) - (I*A + B)/(32*a**2*c**4*f*(tan(e + f*x) + I)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_726():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**5)
    F = -A/(16*a**2*c**5*f*(tan(e + f*x) + I)**3) + x*(21*A + 9*I*B)/(128*a**2*c**5) + (A - I*B)/(40*a**2*c**5*f*(tan(e + f*x) + I)**5) - (3*A + 2*I*B)/(64*a**2*c**5*f*(-tan(e + f*x) + I)) + (15*A + 5*I*B)/(128*a**2*c**5*f*(tan(e + f*x) + I)) - (I*A - B)/(128*a**2*c**5*f*(-tan(e + f*x) + I)**2) - (3*I*A + B)/(64*a**2*c**5*f*(tan(e + f*x) + I)**4) + (5*I*A - B)/(64*a**2*c**5*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_727():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)**3
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**n/(6*a**3*f*(I*tan(e + f*x) + 1)**3) + (I*A*(3 - n) + B*(n + 3))*(-I*c*tan(e + f*x) + c)**n*hyper((3, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(48*a**3*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_728():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**5/(I*a*tan(e + f*x) + a)**3
    F = B*c**5*tan(e + f*x)**2/(2*a**3*f) + c**5*x*(-8*A - 32*I*B)/a**3 + c**5*(A + 8*I*B)*tan(e + f*x)/(a**3*f) + c**5*(16*A + 16*I*B)/(3*a**3*f*(-tan(e + f*x) + I)**3) - c**5*(24*A + 56*I*B)/(a**3*f*(-tan(e + f*x) + I)) - c**5*(8*I*A - 32*B)*log(cos(e + f*x))/(a**3*f) + c**5*(16*I*A - 24*B)/(a**3*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_729():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**3
    F = I*B*c**4*tan(e + f*x)/(a**3*f) - c**4*x*(A + 7*I*B)/a**3 - c**4*(6*A + 18*I*B)/(a**3*f*(-tan(e + f*x) + I)) + c**4*(8*A + 8*I*B)/(3*a**3*f*(-tan(e + f*x) + I)**3) - c**4*(I*A - 7*B)*log(cos(e + f*x))/(a**3*f) + c**4*(6*I*A - 10*B)/(a**3*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_730():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**3
    F = -I*B*c**3*x/a**3 + B*c**3*log(cos(e + f*x))/(a**3*f) - 4*I*B*c**3/(a**3*f*(-tan(e + f*x) + I)) - 2*B*c**3/(a**3*f*(-tan(e + f*x) + I)**2) + c**3*(I*A - B)*(-I*tan(e + f*x) + 1)**3/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_731():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)**3
    F = -I*B*c**2/(a**3*f*(-tan(e + f*x) + I)) + c**2*(2*A + 2*I*B)/(3*a**3*f*(-tan(e + f*x) + I)**3) + c**2*(I*A - 3*B)/(2*a**3*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_732():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**3
    F = -B*c/(2*a**3*f*(-tan(e + f*x) + I)**2) + c*(A + I*B)/(3*a**3*f*(-tan(e + f*x) + I)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_733():
    f = (A + B*tan(e + f*x))/(I*a*tan(e + f*x) + a)**3
    F = (I*A - B)/(6*f*(I*a*tan(e + f*x) + a)**3) + (I*A + B)/(8*f*(I*a**3*tan(e + f*x) + a**3)) + (I*A + B)/(8*a*f*(I*a*tan(e + f*x) + a)**2) + x*(A - I*B)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_734():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c))
    F = -I*A/(8*a**3*c*f*(-tan(e + f*x) + I)**2) + x*(2*A - I*B)/(8*a**3*c) + (A - I*B)/(16*a**3*c*f*(tan(e + f*x) + I)) + (A + I*B)/(12*a**3*c*f*(-tan(e + f*x) + I)**3) - (3*A - I*B)/(16*a**3*c*f*(-tan(e + f*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_735():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**2)
    F = -3*A/(16*a**3*c**2*f*(-tan(e + f*x) + I)) + x*(5*A - I*B)/(16*a**3*c**2) + (A + I*B)/(24*a**3*c**2*f*(-tan(e + f*x) + I)**3) + (2*A - I*B)/(16*a**3*c**2*f*(tan(e + f*x) + I)) + (I*A + B)/(32*a**3*c**2*f*(tan(e + f*x) + I)**2) - (3*I*A - B)/(32*a**3*c**2*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_736():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**3)
    F = 5*A*x/(16*a**3*c**3) + 5*A*sin(e + f*x)*cos(e + f*x)**3/(24*a**3*c**3*f) + 5*A*sin(e + f*x)*cos(e + f*x)/(16*a**3*c**3*f) - (-A*tan(e + f*x) + B)*cos(e + f*x)**6/(6*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_737():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**4)
    F = 5*A/(32*a**3*c**4*f*(tan(e + f*x) + I)) + x*(35*A + 5*I*B)/(128*a**3*c**4) + (A + I*B)/(96*a**3*c**4*f*(-tan(e + f*x) + I)**3) - (2*A - I*B)/(48*a**3*c**4*f*(tan(e + f*x) + I)**3) - (15*A + 5*I*B)/(128*a**3*c**4*f*(-tan(e + f*x) + I)) - (I*A + B)/(64*a**3*c**4*f*(tan(e + f*x) + I)**4) - (5*I*A - 3*B)/(128*a**3*c**4*f*(-tan(e + f*x) + I)**2) + (5*I*A + B)/(64*a**3*c**4*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_738():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**5)
    F = 5*I*A/(64*a**3*c**5*f*(tan(e + f*x) + I)**2) + x*(28*A + 7*I*B)/(128*a**3*c**5) + (A - I*B)/(80*a**3*c**5*f*(tan(e + f*x) + I)**5) + (A + I*B)/(192*a**3*c**5*f*(-tan(e + f*x) + I)**3) - (5*A - I*B)/(96*a**3*c**5*f*(tan(e + f*x) + I)**3) - (21*A + 9*I*B)/(256*a**3*c**5*f*(-tan(e + f*x) + I)) + (35*A + 5*I*B)/(256*a**3*c**5*f*(tan(e + f*x) + I)) - (2*I*A + B)/(64*a**3*c**5*f*(tan(e + f*x) + I)**4) - (3*I*A - 2*B)/(128*a**3*c**5*f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_739():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**6)
    F = -5*A/(96*a**3*c**6*f*(tan(e + f*x) + I)**3) + x*(21*A + 7*I*B)/(128*a**3*c**6) + (A + I*B)/(384*a**3*c**6*f*(-tan(e + f*x) + I)**3) + (2*A - I*B)/(80*a**3*c**6*f*(tan(e + f*x) + I)**5) - (14*A + 7*I*B)/(256*a**3*c**6*f*(-tan(e + f*x) + I)) + (28*A + 7*I*B)/(256*a**3*c**6*f*(tan(e + f*x) + I)) + (I*A + B)/(96*a**3*c**6*f*(tan(e + f*x) + I)**6) - (5*I*A + B)/(128*a**3*c**6*f*(tan(e + f*x) + I)**4) - (7*I*A - 5*B)/(512*a**3*c**6*f*(-tan(e + f*x) + I)**2) + (35*I*A - 5*B)/(512*a**3*c**6*f*(tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_740():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c*f) + 2*a*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_741():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c*f) + 2*a*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_742():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c*f) + 2*a*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_743():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)
    F = -2*B*a*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c*f) + 2*a*(I*A + B)*sqrt(-I*c*tan(e + f*x) + c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_744():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/sqrt(-I*c*tan(e + f*x) + c)
    F = -2*B*a*sqrt(-I*c*tan(e + f*x) + c)/(c*f) - 2*a*(I*A + B)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_745():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 2*B*a/(c*f*sqrt(-I*c*tan(e + f*x) + c)) - 2*a*(I*A + B)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_746():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*B*a/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*a*(I*A + B)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_747():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = 2*B*a/(5*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 2*a*(I*A + B)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_748():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = 2*B*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)/(11*c**2*f) + 4*a**2*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*f) - 2*a**2*(I*A + 3*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_749():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*B*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c**2*f) + 4*a**2*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f) - 2*a**2*(I*A + 3*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_750():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 2*B*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c**2*f) + 4*a**2*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f) - 2*a**2*(I*A + 3*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_751():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2*sqrt(-I*c*tan(e + f*x) + c)
    F = 2*B*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c**2*f) + 4*a**2*(I*A + B)*sqrt(-I*c*tan(e + f*x) + c)/f - 2*a**2*(I*A + 3*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_752():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/sqrt(-I*c*tan(e + f*x) + c)
    F = 2*B*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c**2*f) - 4*a**2*(I*A + B)/(f*sqrt(-I*c*tan(e + f*x) + c)) - 2*a**2*(I*A + 3*B)*sqrt(-I*c*tan(e + f*x) + c)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_753():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 2*B*a**2*sqrt(-I*c*tan(e + f*x) + c)/(c**2*f) - 4*a**2*(I*A + B)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**2*(I*A + 3*B)/(c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_754():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a**2/(c**2*f*sqrt(-I*c*tan(e + f*x) + c)) - 4*a**2*(I*A + B)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**2*(I*A + 3*B)/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_755():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a**2/(3*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 4*a**2*(I*A + B)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) + 2*a**2*(I*A + 3*B)/(5*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_756():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)/(13*c**3*f) + 8*a**3*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*f) - 8*a**3*(I*A + 2*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c*f) + 2*a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)/(11*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_757():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)/(11*c**3*f) + 8*a**3*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f) - 8*a**3*(I*A + 2*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c*f) + 2*a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_758():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c**3*f) + 8*a**3*(I*A + B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f) - 8*a**3*(I*A + 2*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c*f) + 2*a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_759():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3*sqrt(-I*c*tan(e + f*x) + c)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c**3*f) + 8*a**3*(I*A + B)*sqrt(-I*c*tan(e + f*x) + c)/f - 8*a**3*(I*A + 2*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c*f) + 2*a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_760():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/sqrt(-I*c*tan(e + f*x) + c)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c**3*f) - 8*a**3*(I*A + B)/(f*sqrt(-I*c*tan(e + f*x) + c)) - 8*a**3*(I*A + 2*B)*sqrt(-I*c*tan(e + f*x) + c)/(c*f) + 2*a**3*(I*A + 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_761():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c**3*f) - 8*a**3*(I*A + B)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*a**3*(I*A + 2*B)/(c*f*sqrt(-I*c*tan(e + f*x) + c)) + 2*a**3*(I*A + 5*B)*sqrt(-I*c*tan(e + f*x) + c)/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_762():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a**3*sqrt(-I*c*tan(e + f*x) + c)/(c**3*f) - 8*a**3*(I*A + B)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 8*a**3*(I*A + 2*B)/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*a**3*(I*A + 5*B)/(c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_763():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = 2*B*a**3/(c**3*f*sqrt(-I*c*tan(e + f*x) + c)) - 8*a**3*(I*A + B)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) + 8*a**3*(I*A + 2*B)/(5*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 2*a**3*(I*A + 5*B)/(3*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_764():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)
    F = -2*sqrt(2)*c**(sympy.S(7)/2)*(5*I*A - 9*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(a*f) + c**3*(10*I*A - 18*B)*sqrt(-I*c*tan(e + f*x) + c)/(a*f) + c**2*(5*I*A - 9*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*a*f) + c*(5*I*A - 9*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(10*a*f) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(2*a*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_765():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)
    F = -sqrt(2)*c**(sympy.S(5)/2)*(3*I*A - 7*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(a*f) + c**2*(3*I*A - 7*B)*sqrt(-I*c*tan(e + f*x) + c)/(a*f) + c*(3*I*A - 7*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(6*a*f) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(2*a*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_766():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)
    F = -sqrt(2)*c**(sympy.S(3)/2)*(I*A - 5*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(2*a*f) + c*(I*A - 5*B)*sqrt(-I*c*tan(e + f*x) + c)/(2*a*f) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(2*a*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_767():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)
    F = sqrt(2)*sqrt(c)*(I*A + 3*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(4*a*f) + (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(2*a*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_768():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    F = (I*A - B)/(2*a*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) - (3*I*A + B)/(4*a*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(3*I*A + B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(8*a*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_769():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = (I*A - B)/(2*a*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (5*I*A - B)/(12*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (5*I*A - B)/(8*a*c*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(5*I*A - B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(16*a*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_770():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = (I*A - B)/(2*a*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (7*I*A - 3*B)/(20*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (7*I*A - 3*B)/(24*a*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (7*I*A - 3*B)/(16*a*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(7*I*A - 3*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_771():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(I*a*tan(e + f*x) + a)**2
    F = sqrt(2)*c**(sympy.S(9)/2)*(35*I*A - 91*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(2*a**2*f) - c**4*(35*I*A - 91*B)*sqrt(-I*c*tan(e + f*x) + c)/(2*a**2*f) - c**3*(35*I*A - 91*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(12*a**2*f) - c**2*(35*I*A - 91*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(40*a**2*f) - c*(5*I*A - 13*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(8*a**2*f*(I*tan(e + f*x) + 1)) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_772():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**2
    F = sqrt(2)*c**(sympy.S(7)/2)*(15*I*A - 55*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(4*a**2*f) - c**3*(15*I*A - 55*B)*sqrt(-I*c*tan(e + f*x) + c)/(4*a**2*f) - c**2*(15*I*A - 55*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(24*a**2*f) - c*(3*I*A - 11*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(8*a**2*f*(I*tan(e + f*x) + 1)) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_773():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**2
    F = sqrt(2)*c**(sympy.S(5)/2)*(3*I*A - 27*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(8*a**2*f) - c**2*(3*I*A - 27*B)*sqrt(-I*c*tan(e + f*x) + c)/(8*a**2*f) - c*(I*A - 9*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(8*a**2*f*(I*tan(e + f*x) + 1)) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_774():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**2
    F = -sqrt(2)*c**(sympy.S(3)/2)*(I*A + 7*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(16*a**2*f) + c*(I*A + 7*B)*sqrt(-I*c*tan(e + f*x) + c)/(8*a**2*f*(I*tan(e + f*x) + 1)) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_775():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**2
    F = sqrt(2)*sqrt(c)*(3*I*A + 5*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a**2*f) + (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(4*a**2*f*(I*tan(e + f*x) + 1)**2) + (3*I*A + 5*B)*sqrt(-I*c*tan(e + f*x) + c)/(16*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_776():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*sqrt(-I*c*tan(e + f*x) + c))
    F = (I*A - B)/(4*a**2*f*(I*tan(e + f*x) + 1)**2*sqrt(-I*c*tan(e + f*x) + c)) + (5*I*A + 3*B)/(16*a**2*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) - (15*I*A + 9*B)/(32*a**2*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(15*I*A + 9*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(64*a**2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_777():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = (I*A - B)/(4*a**2*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + (7*I*A + B)/(16*a**2*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (35*I*A + 5*B)/(96*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (35*I*A + 5*B)/(64*a**2*c*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(35*I*A + 5*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(128*a**2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_778():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = (I*A - B)/(4*a**2*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + (9*I*A - B)/(16*a**2*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (63*I*A - 7*B)/(160*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (63*I*A - 7*B)/(192*a**2*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (63*I*A - 7*B)/(128*a**2*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(63*I*A - 7*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(256*a**2*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_779():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(I*a*tan(e + f*x) + a)**3
    F = -sqrt(2)*c**(sympy.S(9)/2)*(35*I*A - 175*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(8*a**3*f) + c**4*(35*I*A - 175*B)*sqrt(-I*c*tan(e + f*x) + c)/(8*a**3*f) + c**3*(35*I*A - 175*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(48*a**3*f) + c**2*(7*I*A - 35*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(16*a**3*f*(I*tan(e + f*x) + 1)) - c*(I*A - 5*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(8*a**3*f*(I*tan(e + f*x) + 1)**2) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_780():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**3
    F = sqrt(2)*c**(sympy.S(7)/2)*(-5*I*A + 65*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(16*a**3*f) + c**3*(5*I*A - 65*B)*sqrt(-I*c*tan(e + f*x) + c)/(16*a**3*f) + c**2*(5*I*A - 65*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(48*a**3*f*(I*tan(e + f*x) + 1)) - c*(I*A - 13*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(24*a**3*f*(I*tan(e + f*x) + 1)**2) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_781():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**3
    F = sqrt(2)*c**(sympy.S(5)/2)*(I*A + 11*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a**3*f) - c**2*(I*A + 11*B)*sqrt(-I*c*tan(e + f*x) + c)/(16*a**3*f*(I*tan(e + f*x) + 1)) + c*(I*A + 11*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(24*a**3*f*(I*tan(e + f*x) + 1)**2) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_782():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**3
    F = -sqrt(2)*c**(sympy.S(3)/2)*(I*A + 3*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(64*a**3*f) - c*(I*A + 3*B)*sqrt(-I*c*tan(e + f*x) + c)/(32*a**3*f*(I*tan(e + f*x) + 1)) + c*(I*A + 3*B)*sqrt(-I*c*tan(e + f*x) + c)/(8*a**3*f*(I*tan(e + f*x) + 1)**2) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_783():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**3
    F = sqrt(2)*sqrt(c)*(5*I*A + 7*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(128*a**3*f) + (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(6*a**3*f*(I*tan(e + f*x) + 1)**3) + (5*I*A + 7*B)*sqrt(-I*c*tan(e + f*x) + c)/(64*a**3*f*(I*tan(e + f*x) + 1)) + (5*I*A + 7*B)*sqrt(-I*c*tan(e + f*x) + c)/(48*a**3*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_784():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*sqrt(-I*c*tan(e + f*x) + c))
    F = (I*A - B)/(6*a**3*f*(I*tan(e + f*x) + 1)**3*sqrt(-I*c*tan(e + f*x) + c)) + (7*I*A + 5*B)/(48*a**3*f*(I*tan(e + f*x) + 1)**2*sqrt(-I*c*tan(e + f*x) + c)) - (35*I*A + 25*B)/(128*a**3*f*sqrt(-I*c*tan(e + f*x) + c)) + (35*I*A + 25*B)/(192*a**3*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(35*I*A + 25*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(256*a**3*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_785():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = (I*A - B)/(6*a**3*f*(I*tan(e + f*x) + 1)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + (3*I*A + B)/(16*a**3*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + (21*I*A + 7*B)/(64*a**3*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (105*I*A + 35*B)/(384*a**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (105*I*A + 35*B)/(256*a**3*c*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(105*I*A + 35*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(512*a**3*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_786():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = (I*A - B)/(6*a**3*f*(I*tan(e + f*x) + 1)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + (11*I*A + B)/(48*a**3*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + (33*I*A + 3*B)/(64*a**3*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (231*I*A + 21*B)/(640*a**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (77*I*A + 7*B)/(256*a**3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (231*I*A + 21*B)/(512*a**3*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + sqrt(2)*(231*I*A + 21*B)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(1024*a**3*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_787():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = B*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(4*f) - 5*sqrt(a)*c**(sympy.S(7)/2)*(4*I*A - 3*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) - c**3*(20*I*A - 15*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(8*f) - c**2*(20*I*A - 15*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(24*f) - c*(4*I*A - 3*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_788():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = B*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(3*f) - sqrt(a)*c**(sympy.S(5)/2)*(3*I*A - 2*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f - c**2*(3*I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f) - c*(3*I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_789():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = B*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(2*f) - sqrt(a)*c**(sympy.S(3)/2)*(2*I*A - B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f - c*(2*I*A - B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_790():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)
    F = -2*I*A*sqrt(a)*sqrt(c)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + B*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_791():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/sqrt(-I*c*tan(e + f*x) + c)
    F = 2*B*sqrt(a)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - (I*A + B)*sqrt(I*a*tan(e + f*x) + a)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_792():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -(I*A + B)*sqrt(I*a*tan(e + f*x) + a)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_793():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -(I*A + B)*sqrt(I*a*tan(e + f*x) + a)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (2*I*A - 3*B)*sqrt(I*a*tan(e + f*x) + a)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (2*I*A - 3*B)*sqrt(I*a*tan(e + f*x) + a)/(15*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_794():
    f = (A + B*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -(I*A + B)*sqrt(I*a*tan(e + f*x) + a)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (3*I*A - 4*B)*sqrt(I*a*tan(e + f*x) + a)/(35*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (6*I*A - 8*B)*sqrt(I*a*tan(e + f*x) + a)/(105*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (6*I*A - 8*B)*sqrt(I*a*tan(e + f*x) + a)/(105*c**3*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_795():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(5*f) - a**(sympy.S(3)/2)*c**(sympy.S(7)/2)*(5*I*A - 2*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + a*c**3*(5*A + 2*I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) - c**2*(5*I*A - 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(12*f) - c*(5*I*A - 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_796():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(4*f) - a**(sympy.S(3)/2)*c**(sympy.S(5)/2)*(4*I*A - B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + a*c**2*(4*A + I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) - c*(4*I*A - B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_797():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -I*A*a**(sympy.S(3)/2)*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + A*a*c*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(2*f) + B*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_798():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(2*f) - a**(sympy.S(3)/2)*sqrt(c)*(2*I*A + B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + a*(2*I*A + B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_799():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 2*a**(sympy.S(3)/2)*(I*A + 2*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - a*(I*A + 2*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c*f) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_800():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + 2*B*a*sqrt(I*a*tan(e + f*x) + a)/(c*f*sqrt(-I*c*tan(e + f*x) + c)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_801():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (I*A - 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_802():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (2*I*A - 5*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(35*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (2*I*A - 5*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(105*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_803():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(9*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (I*A - 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(21*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (2*I*A - 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(105*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (2*I*A - 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(315*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_804():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(11*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (4*I*A - 7*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(99*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (4*I*A - 7*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(231*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (8*I*A - 14*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(1155*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (8*I*A - 14*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3465*c**4*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_805():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(6*f) - a**(sympy.S(5)/2)*c**(sympy.S(7)/2)*(6*I*A - B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(8*f) + a**2*c**3*(6*A + I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(16*f) + a*c**2*(6*A + I*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(24*f) - c*(6*I*A - B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_806():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -3*I*A*a**(sympy.S(5)/2)*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + 3*A*a**2*c**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) + A*a*c*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(4*f) + B*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_807():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(4*f) - a**(sympy.S(5)/2)*c**(sympy.S(3)/2)*(4*I*A + B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + a**2*c*(4*A - I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) + a*(4*I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_808():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)/(3*f) - a**(sympy.S(5)/2)*sqrt(c)*(3*I*A + 2*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + a**2*(3*I*A + 2*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f) + a*(3*I*A + 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_809():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 3*a**(sympy.S(5)/2)*(2*I*A + 3*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - 3*a**2*(2*I*A + 3*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c*f) - a*(2*I*A + 3*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(2*c*f) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_810():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(5)/2)*(I*A + 4*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + a**2*(I*A + 4*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c**2*f) + 2*a*(I*A + 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_811():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*B*a**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(5)/2)*f) - 2*B*a**2*sqrt(I*a*tan(e + f*x) + a)/(c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 2*B*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_812():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (I*A - 6*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(35*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_813():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(9*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (2*I*A - 7*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(63*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (2*I*A - 7*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(315*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_814():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(11*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (3*I*A - 8*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(99*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (6*I*A - 16*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(693*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (6*I*A - 16*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3465*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_815():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(13*f*(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)) - (4*I*A - 9*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(143*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (4*I*A - 9*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(429*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (8*I*A - 18*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3003*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)) - (8*I*A - 18*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(15015*c**4*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_816():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(8*f) - 5*a**(sympy.S(7)/2)*c**(sympy.S(9)/2)*(8*I*A - B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(64*f) + 5*a**3*c**4*(8*A + I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(128*f) + 5*a**2*c**3*(8*A + I*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(192*f) + a*c**2*(8*A + I*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(48*f) - c*(8*I*A - B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(56*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_817():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -5*I*A*a**(sympy.S(7)/2)*c**(sympy.S(7)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(8*f) + 5*A*a**3*c**3*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(16*f) + 5*A*a**2*c**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(24*f) + A*a*c*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(6*f) + B*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_818():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(6*f) - a**(sympy.S(7)/2)*c**(sympy.S(5)/2)*(6*I*A + B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(8*f) + a**3*c**2*(6*A - I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(16*f) + a**2*c*(6*A - I*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(24*f) + a*(6*I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_819():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(5*f) - a**(sympy.S(7)/2)*c**(sympy.S(3)/2)*(5*I*A + 2*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + a**3*c*(5*A - 2*I*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) + a**2*(5*I*A + 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(12*f) + a*(5*I*A + 2*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_820():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-I*c*tan(e + f*x) + c)
    F = B*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-I*c*tan(e + f*x) + c)/(4*f) - 5*a**(sympy.S(7)/2)*sqrt(c)*(4*I*A + 3*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + 5*a**3*(4*I*A + 3*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(8*f) + 5*a**2*(4*I*A + 3*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(24*f) + a*(4*I*A + 3*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_821():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 5*a**(sympy.S(7)/2)*(3*I*A + 4*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - 5*a**3*(3*I*A + 4*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c*f) - 5*a**2*(3*I*A + 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(6*c*f) - a*(3*I*A + 4*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)/(3*c*f) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_822():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -5*a**(sympy.S(7)/2)*(2*I*A + 5*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + 5*a**3*(2*I*A + 5*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c**2*f) + 5*a**2*(2*I*A + 5*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(6*c**2*f) + 2*a*(2*I*A + 5*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_823():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(7)/2)*(I*A + 6*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(5)/2)*f) - a**3*(I*A + 6*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c**3*f) - 2*a**2*(I*A + 6*B)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 2*a*(I*A + 6*B)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_824():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a**(sympy.S(7)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(7)/2)*f) + 2*B*a**3*sqrt(I*a*tan(e + f*x) + a)/(c**3*f*sqrt(-I*c*tan(e + f*x) + c)) - 2*B*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 2*B*a*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(5*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(7*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_825():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(9*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (I*A - 8*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(63*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_826():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(11*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (2*I*A - 9*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(99*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (2*I*A - 9*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(693*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_827():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(13*f*(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)) - (3*I*A - 10*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(143*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (6*I*A - 20*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(1287*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (6*I*A - 20*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(9009*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_828():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(15)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(15*f*(-I*c*tan(e + f*x) + c)**(sympy.S(15)/2)) - (4*I*A - 11*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(195*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)) - (4*I*A - 11*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(715*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (8*I*A - 22*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(6435*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (8*I*A - 22*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(45045*c**4*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_829():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(17)/2)
    F = -(I*A + B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(17*f*(-I*c*tan(e + f*x) + c)**(sympy.S(17)/2)) - (5*I*A - 12*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(255*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(15)/2)) - (20*I*A - 48*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(3315*c**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(13)/2)) - (20*I*A - 48*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(12155*c**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(11)/2)) - (40*I*A - 96*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(109395*c**4*f*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)) - (40*I*A - 96*B)*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(765765*c**5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_830():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(f*sqrt(I*a*tan(e + f*x) + a)) + c**2*(6*I*A - 9*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*a*f) + c*(2*I*A - 3*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(2*a*f) + c**(sympy.S(5)/2)*(6*I*A - 9*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_831():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(f*sqrt(I*a*tan(e + f*x) + a)) + c*(I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(a*f) + c**(sympy.S(3)/2)*(2*I*A - 4*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_832():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/sqrt(I*a*tan(e + f*x) + a)
    F = -2*B*sqrt(c)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(a)*f) + (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_833():
    f = (A + B*tan(e + f*x))/(sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    F = I*A*sqrt(-I*c*tan(e + f*x) + c)/(c*f*sqrt(I*a*tan(e + f*x) + a)) - (I*A + B)/(f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_834():
    f = (A + B*tan(e + f*x))/(sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = (I*A - B)/(f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (2*I*A - B)*sqrt(I*a*tan(e + f*x) + a)/(3*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (2*I*A - B)*sqrt(I*a*tan(e + f*x) + a)/(3*a*c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_835():
    f = (A + B*tan(e + f*x))/(sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = (I*A - B)/(f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (3*I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)/(5*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (6*I*A - 4*B)*sqrt(I*a*tan(e + f*x) + a)/(15*a*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (6*I*A - 4*B)*sqrt(I*a*tan(e + f*x) + a)/(15*a*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_836():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) - c*(4*I*A - 10*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(3*a*f*sqrt(I*a*tan(e + f*x) + a)) - c**3*(10*I*A - 25*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*a**2*f) - c**2*(10*I*A - 25*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(6*a**2*f) + c**(sympy.S(7)/2)*(-10*I*A + 25*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_837():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) - c*(2*I*A - 8*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*a*f*sqrt(I*a*tan(e + f*x) + a)) - c**2*(I*A - 4*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(a**2*f) + c**(sympy.S(5)/2)*(-2*I*A + 8*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_838():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*B*c*sqrt(-I*c*tan(e + f*x) + c)/(a*f*sqrt(I*a*tan(e + f*x) + a)) + 2*B*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(3)/2)*f) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_839():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (I*A + 2*B)*sqrt(-I*c*tan(e + f*x) + c)/(3*a*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_840():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c))
    F = -(I*A + B)/(f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + (2*I*A + B)*sqrt(-I*c*tan(e + f*x) + c)/(3*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (2*I*A + B)*sqrt(-I*c*tan(e + f*x) + c)/(3*a*c*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_841():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = I*A/(3*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 2*A*tan(e + f*x)/(3*a*c*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)) + (-I*A - B)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_842():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = (I*A - B)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + (4*I*A - B)/(3*a*f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (4*I*A - B)*sqrt(I*a*tan(e + f*x) + a)/(5*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - (8*I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)/(15*a**2*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - (8*I*A - 2*B)*sqrt(I*a*tan(e + f*x) + a)/(15*a**2*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_843():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) - c*(4*I*A - 14*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + c**2*(28*I*A - 98*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(15*a**2*f*sqrt(I*a*tan(e + f*x) + a)) + c**4*(14*I*A - 49*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*a**3*f) + c**3*(14*I*A - 49*B)*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(6*a**3*f) + c**(sympy.S(9)/2)*(14*I*A - 49*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_844():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) - c*(2*I*A - 12*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + c**2*(2*I*A - 12*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*a**2*f*sqrt(I*a*tan(e + f*x) + a)) + c**3*(I*A - 6*B)*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(a**3*f) + c**(sympy.S(7)/2)*(2*I*A - 12*B)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_845():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*B*c*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) - 2*B*c**2*sqrt(-I*c*tan(e + f*x) + c)/(a**2*f*sqrt(I*a*tan(e + f*x) + a)) - 2*B*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(5)/2)*f) + (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_846():
    f = (A + B*tan(e + f*x))*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (I*A + 4*B)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_847():
    f = (A + B*tan(e + f*x))*sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*sqrt(-I*c*tan(e + f*x) + c)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (2*I*A + 3*B)*sqrt(-I*c*tan(e + f*x) + c)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (2*I*A + 3*B)*sqrt(-I*c*tan(e + f*x) + c)/(15*a**2*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_848():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c))
    F = -(I*A + B)/(f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)) + (3*I*A + 2*B)*sqrt(-I*c*tan(e + f*x) + c)/(5*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (6*I*A + 4*B)*sqrt(-I*c*tan(e + f*x) + c)/(15*a*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (6*I*A + 4*B)*sqrt(-I*c*tan(e + f*x) + c)/(15*a**2*c*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_849():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = (-I*A - B)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + (4*I*A + B)/(15*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)) + (4*I*A + B)/(15*a*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + (8*A - 2*I*B)*tan(e + f*x)/(15*a**2*c*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_850():
    f = (A + B*tan(e + f*x))/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = I*A/(5*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 4*A*tan(e + f*x)/(15*a*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*A*tan(e + f*x)/(15*a**2*c**2*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)) + (-I*A - B)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_851():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**n
    F = -2**(n - 1)*(I*A*(m + n) + B*(m - n))*(I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**n*hyper((m, -n), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(f*m*n*(-I*tan(e + f*x) + 1)**n) + (I*A + B)*(I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**n/(2*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_852():
    f = (A + B*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(m + 1)*(-I*c*tan(e + f*x) + c)**(-m - 1)
    F = 2**m*B*a*(I*a*tan(e + f*x) + a)**m*hyper((-m, -m), (1 - m,), -I*tan(e + f*x)/2 + sympy.S.Half)/(c*f*m*(I*tan(e + f*x) + 1)**m*(-I*c*tan(e + f*x) + c)**m) - (I*A + B)*(I*a*tan(e + f*x) + a)**(m + 1)*(-I*c*tan(e + f*x) + c)**(-m - 1)/(2*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_853():
    f = ((n - 2)*tan(e + f*x) - I*(n + 2))*(-I*c*tan(e + f*x) + c)**n/(tan(e + f*x) - I)**2
    F = (-I*c*tan(e + f*x) + c)**n/(f*(-tan(e + f*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_854():
    f = (A + B*tan(e + f*x))*(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**2
    F = (c + I*d)*(I*A - B)/(4*f*(I*a*tan(e + f*x) + a)**2) + x*(A - I*B)*(c - I*d)/(4*a**2) + (A*(I*c + d) + B*(c + 3*I*d))/(4*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_3_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_855():
    f = (A + B*tan(e + f*x))*(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (c + I*d)*(I*A - B)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (A*(I*c + d) + B*(c + 3*I*d))/(2*a*f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*(c - I*d)*(I*A + B)*atanh(sqrt(2)*sqrt(I*a*tan(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F

