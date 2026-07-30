"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.2.1 (a+b tan)^m (c+d tan)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1():
    f = (I*a*tan(c + d*x) + a)*tan(c + d*x)**5
    F = -I*a*x - a*log(cos(c + d*x))/d + I*a*tan(c + d*x)**5/(5*d) + a*tan(c + d*x)**4/(4*d) - I*a*tan(c + d*x)**3/(3*d) - a*tan(c + d*x)**2/(2*d) + I*a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_2():
    f = (I*a*tan(c + d*x) + a)*tan(c + d*x)**4
    F = a*x - I*a*log(cos(c + d*x))/d + I*a*tan(c + d*x)**4/(4*d) + a*tan(c + d*x)**3/(3*d) - I*a*tan(c + d*x)**2/(2*d) - a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_3():
    f = (I*a*tan(c + d*x) + a)*tan(c + d*x)**3
    F = I*a*x + a*log(cos(c + d*x))/d + I*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)**2/(2*d) - I*a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_4():
    f = (I*a*tan(c + d*x) + a)*tan(c + d*x)**2
    F = -a*x + I*a*log(cos(c + d*x))/d + I*a*tan(c + d*x)**2/(2*d) + a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_5():
    f = (I*a*tan(c + d*x) + a)*tan(c + d*x)
    F = -I*a*x - a*log(cos(c + d*x))/d + I*a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_6():
    f = I*a*tan(c + d*x) + a
    F = a*x - I*a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_7():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)
    F = I*a*x + a*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_8():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**2
    F = -a*x + I*a*log(sin(c + d*x))/d - a*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_9():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**3
    F = -I*a*x - a*log(sin(c + d*x))/d - a*cot(c + d*x)**2/(2*d) - I*a*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_10():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**4
    F = a*x - I*a*log(sin(c + d*x))/d - a*cot(c + d*x)**3/(3*d) - I*a*cot(c + d*x)**2/(2*d) + a*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_11():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**5
    F = I*a*x + a*log(sin(c + d*x))/d - a*cot(c + d*x)**4/(4*d) - I*a*cot(c + d*x)**3/(3*d) + a*cot(c + d*x)**2/(2*d) + I*a*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_12():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**6
    F = -a*x + I*a*log(sin(c + d*x))/d - a*cot(c + d*x)**5/(5*d) - I*a*cot(c + d*x)**4/(4*d) + a*cot(c + d*x)**3/(3*d) + I*a*cot(c + d*x)**2/(2*d) - a*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_13():
    f = (I*a*tan(c + d*x) + a)**2*tan(c + d*x)**4
    F = 2*a**2*x - 2*I*a**2*log(cos(c + d*x))/d - a**2*tan(c + d*x)**5/(5*d) + I*a**2*tan(c + d*x)**4/(2*d) + 2*a**2*tan(c + d*x)**3/(3*d) - I*a**2*tan(c + d*x)**2/d - 2*a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_14():
    f = (I*a*tan(c + d*x) + a)**2*tan(c + d*x)**3
    F = 2*I*a**2*x + 2*a**2*log(cos(c + d*x))/d - a**2*tan(c + d*x)**4/(4*d) + 2*I*a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)**2/d - 2*I*a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_15():
    f = (I*a*tan(c + d*x) + a)**2*tan(c + d*x)**2
    F = -2*a**2*x + 2*I*a**2*log(cos(c + d*x))/d + a**2*tan(c + d*x)/d - I*(I*a*tan(c + d*x) + a)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_16():
    f = (I*a*tan(c + d*x) + a)**2*tan(c + d*x)
    F = -2*I*a**2*x - 2*a**2*log(cos(c + d*x))/d + I*a**2*tan(c + d*x)/d + (I*a*tan(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_17():
    f = (I*a*tan(c + d*x) + a)**2
    F = 2*a**2*x - 2*I*a**2*log(cos(c + d*x))/d - a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_18():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)
    F = 2*I*a**2*x + a**2*log(sin(c + d*x))/d + a**2*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_19():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**2
    F = -2*a**2*x + 2*I*a**2*log(sin(c + d*x))/d - a**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_20():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**3
    F = -2*I*a**2*x - 2*a**2*log(sin(c + d*x))/d - a**2*cot(c + d*x)**2/(2*d) - 2*I*a**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_21():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**4
    F = 2*a**2*x - 2*I*a**2*log(sin(c + d*x))/d - a**2*cot(c + d*x)**3/(3*d) - I*a**2*cot(c + d*x)**2/d + 2*a**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_22():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**5
    F = 2*I*a**2*x + 2*a**2*log(sin(c + d*x))/d - a**2*cot(c + d*x)**4/(4*d) - 2*I*a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)**2/d + 2*I*a**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_23():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**6
    F = -2*a**2*x + 2*I*a**2*log(sin(c + d*x))/d - a**2*cot(c + d*x)**5/(5*d) - I*a**2*cot(c + d*x)**4/(2*d) + 2*a**2*cot(c + d*x)**3/(3*d) + I*a**2*cot(c + d*x)**2/d - 2*a**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_24():
    f = (I*a*tan(c + d*x) + a)**3*tan(c + d*x)**3
    F = 4*I*a**3*x + 4*a**3*log(cos(c + d*x))/d - 11*a**3*tan(c + d*x)**4/(20*d) + 4*I*a**3*tan(c + d*x)**3/(3*d) + 2*a**3*tan(c + d*x)**2/d - 4*I*a**3*tan(c + d*x)/d - (I*a**3*tan(c + d*x) + a**3)*tan(c + d*x)**4/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_25():
    f = (I*a*tan(c + d*x) + a)**3*tan(c + d*x)**2
    F = -4*a**3*x + 4*I*a**3*log(cos(c + d*x))/d + 2*a**3*tan(c + d*x)/d - I*a*(I*a*tan(c + d*x) + a)**2/(2*d) - I*(I*a*tan(c + d*x) + a)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_26():
    f = (I*a*tan(c + d*x) + a)**3*tan(c + d*x)
    F = -4*I*a**3*x - 4*a**3*log(cos(c + d*x))/d + 2*I*a**3*tan(c + d*x)/d + a*(I*a*tan(c + d*x) + a)**2/(2*d) + (I*a*tan(c + d*x) + a)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_27():
    f = (I*a*tan(c + d*x) + a)**3
    F = 4*a**3*x - 4*I*a**3*log(cos(c + d*x))/d - 2*a**3*tan(c + d*x)/d + I*a*(I*a*tan(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_28():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)
    F = 4*I*a**3*x + a**3*log(sin(c + d*x))/d + 3*a**3*log(cos(c + d*x))/d - (I*a**3*tan(c + d*x) + a**3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_29():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**2
    F = -4*a**3*x + 3*I*a**3*log(sin(c + d*x))/d + I*a**3*log(cos(c + d*x))/d - (I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_30():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3
    F = -4*I*a**3*x - 4*a**3*log(sin(c + d*x))/d - 2*I*a**3*cot(c + d*x)/d - a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_31():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**4
    F = 4*a**3*x - 4*I*a**3*log(sin(c + d*x))/d + 2*a**3*cot(c + d*x)/d - I*a*(I*a*tan(c + d*x) + a)**2*cot(c + d*x)**2/(2*d) - (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_32():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**5
    F = 4*I*a**3*x + 4*a**3*log(sin(c + d*x))/d - 3*I*a**3*cot(c + d*x)**3/(4*d) + 2*a**3*cot(c + d*x)**2/d + 4*I*a**3*cot(c + d*x)/d - (I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_33():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**6
    F = -4*a**3*x + 4*I*a**3*log(sin(c + d*x))/d - 11*I*a**3*cot(c + d*x)**4/(20*d) + 4*a**3*cot(c + d*x)**3/(3*d) + 2*I*a**3*cot(c + d*x)**2/d - 4*a**3*cot(c + d*x)/d - (I*a**3*tan(c + d*x) + a**3)*cot(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_34():
    f = (I*a*tan(c + d*x) + a)**4*tan(c + d*x)**3
    F = 8*I*a**4*x + 8*a**4*log(cos(c + d*x))/d - 67*a**4*tan(c + d*x)**4/(60*d) + 8*I*a**4*tan(c + d*x)**3/(3*d) + 4*a**4*tan(c + d*x)**2/d - 8*I*a**4*tan(c + d*x)/d - (I*a**2*tan(c + d*x) + a**2)**2*tan(c + d*x)**4/(6*d) - 7*(I*a**4*tan(c + d*x) + a**4)*tan(c + d*x)**4/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_35():
    f = (I*a*tan(c + d*x) + a)**4*tan(c + d*x)**2
    F = -8*a**4*x + 8*I*a**4*log(cos(c + d*x))/d + 4*a**4*tan(c + d*x)/d - I*a*(I*a*tan(c + d*x) + a)**3/(3*d) - I*(I*a**2*tan(c + d*x) + a**2)**2/d - I*(I*a*tan(c + d*x) + a)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_36():
    f = (I*a*tan(c + d*x) + a)**4*tan(c + d*x)
    F = -8*I*a**4*x - 8*a**4*log(cos(c + d*x))/d + 4*I*a**4*tan(c + d*x)/d + a*(I*a*tan(c + d*x) + a)**3/(3*d) + (I*a*tan(c + d*x) + a)**4/(4*d) + (I*a**2*tan(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_37():
    f = (I*a*tan(c + d*x) + a)**4
    F = 8*a**4*x - 8*I*a**4*log(cos(c + d*x))/d - 4*a**4*tan(c + d*x)/d + I*a*(I*a*tan(c + d*x) + a)**3/(3*d) + I*(I*a**2*tan(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_38():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)
    F = 8*I*a**4*x + a**4*log(sin(c + d*x))/d + 7*a**4*log(cos(c + d*x))/d - (I*a**2*tan(c + d*x) + a**2)**2/(2*d) - (3*I*a**4*tan(c + d*x) + 3*a**4)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_39():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**2
    F = -8*a**4*x + 4*I*a**4*log(sin(c + d*x))/d + 4*I*a**4*log(cos(c + d*x))/d - (I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_40():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**3
    F = -8*I*a**4*x - 7*a**4*log(sin(c + d*x))/d - a**4*log(cos(c + d*x))/d - (I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**2/(2*d) - 3*I*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_41():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**4
    F = 8*a**4*x - 8*I*a**4*log(sin(c + d*x))/d + 4*a**4*cot(c + d*x)/d - a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3/(3*d) - I*(I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_42():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**5
    F = 8*I*a**4*x + 8*a**4*log(sin(c + d*x))/d + 4*I*a**4*cot(c + d*x)/d - I*a*(I*a*tan(c + d*x) + a)**3*cot(c + d*x)**3/(3*d) - (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**4/(4*d) + (I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_43():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**6
    F = -8*a**4*x + 8*I*a**4*log(sin(c + d*x))/d + 23*a**4*cot(c + d*x)**3/(15*d) + 4*I*a**4*cot(c + d*x)**2/d - 8*a**4*cot(c + d*x)/d - (I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**5/(5*d) - 3*I*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)**4/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_44():
    f = (I*a*tan(c + d*x) + a)**4*cot(c + d*x)**7
    F = -8*I*a**4*x - 8*a**4*log(sin(c + d*x))/d + 67*a**4*cot(c + d*x)**4/(60*d) + 8*I*a**4*cot(c + d*x)**3/(3*d) - 4*a**4*cot(c + d*x)**2/d - 8*I*a**4*cot(c + d*x)/d - (I*a**2*tan(c + d*x) + a**2)**2*cot(c + d*x)**6/(6*d) - 7*I*(I*a**4*tan(c + d*x) + a**4)*cot(c + d*x)**5/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_45():
    f = tan(c + d*x)**6/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**5/(2*d*(I*a*tan(c + d*x) + a)) + 5*x/(2*a) + 3*I*log(cos(c + d*x))/(a*d) - 3*I*tan(c + d*x)**4/(4*a*d) + 5*tan(c + d*x)**3/(6*a*d) + 3*I*tan(c + d*x)**2/(2*a*d) - 5*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_46():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**4/(2*d*(I*a*tan(c + d*x) + a)) - 5*I*x/(2*a) + 2*log(cos(c + d*x))/(a*d) - 5*I*tan(c + d*x)**3/(6*a*d) + tan(c + d*x)**2/(a*d) + 5*I*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_47():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**3/(2*d*(I*a*tan(c + d*x) + a)) - 3*x/(2*a) - 2*I*log(cos(c + d*x))/(a*d) - I*tan(c + d*x)**2/(a*d) + 3*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_48():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**2/(2*d*(I*a*tan(c + d*x) + a)) + 3*I*x/(2*a) - log(cos(c + d*x))/(a*d) - 3*I*tan(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_49():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = -I/(2*d*(I*a*tan(c + d*x) + a)) + x/(2*a) + I*log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_50():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)
    F = -1/(2*d*(I*a*tan(c + d*x) + a)) - I*x/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_51():
    f = 1/(I*a*tan(c + d*x) + a)
    F = I/(2*d*(I*a*tan(c + d*x) + a)) + x/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_52():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)
    F = 1/(2*d*(I*a*tan(c + d*x) + a)) - I*x/(2*a) + log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_53():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)/(2*d*(I*a*tan(c + d*x) + a)) - 3*x/(2*a) - I*log(sin(c + d*x))/(a*d) - 3*cot(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_54():
    f = cot(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)**2/(2*d*(I*a*tan(c + d*x) + a)) + 3*I*x/(2*a) - 2*log(sin(c + d*x))/(a*d) - cot(c + d*x)**2/(a*d) + 3*I*cot(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_55():
    f = cot(c + d*x)**4/(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)**3/(2*d*(I*a*tan(c + d*x) + a)) + 5*x/(2*a) + 2*I*log(sin(c + d*x))/(a*d) - 5*cot(c + d*x)**3/(6*a*d) + I*cot(c + d*x)**2/(a*d) + 5*cot(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_56():
    f = tan(c + d*x)**6/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**5/(4*d*(I*a*tan(c + d*x) + a)**2) - 25*x/(4*a**2) - 6*I*log(cos(c + d*x))/(a**2*d) - 25*tan(c + d*x)**3/(12*a**2*d) - 3*I*tan(c + d*x)**2/(a**2*d) + 25*tan(c + d*x)/(4*a**2*d) + 3*I*tan(c + d*x)**4/(2*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_57():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**4/(4*d*(I*a*tan(c + d*x) + a)**2) + 15*I*x/(4*a**2) - 4*log(cos(c + d*x))/(a**2*d) - 2*tan(c + d*x)**2/(a**2*d) - 15*I*tan(c + d*x)/(4*a**2*d) + 5*I*tan(c + d*x)**3/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_58():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**3/(4*d*(I*a*tan(c + d*x) + a)**2) + 9*x/(4*a**2) + 2*I*log(cos(c + d*x))/(a**2*d) - 9*tan(c + d*x)/(4*a**2*d) + I*tan(c + d*x)**2/(a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_59():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**2/(4*d*(I*a*tan(c + d*x) + a)**2) - 3*I*x/(4*a**2) + log(cos(c + d*x))/(a**2*d) - 3/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_60():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = -I/(4*d*(I*a*tan(c + d*x) + a)**2) - x/(4*a**2) + 3*I/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_61():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = 1/(4*d*(I*a**2*tan(c + d*x) + a**2)) - 1/(4*d*(I*a*tan(c + d*x) + a)**2) - I*x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_62():
    f = (I*a*tan(c + d*x) + a)**(-2)
    F = I/(4*d*(I*a**2*tan(c + d*x) + a**2)) + I/(4*d*(I*a*tan(c + d*x) + a)**2) + x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_63():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = 1/(4*d*(I*a*tan(c + d*x) + a)**2) - 3*I*x/(4*a**2) + log(sin(c + d*x))/(a**2*d) + 3/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_64():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = cot(c + d*x)/(4*d*(I*a*tan(c + d*x) + a)**2) - 9*x/(4*a**2) - 2*I*log(sin(c + d*x))/(a**2*d) - 9*cot(c + d*x)/(4*a**2*d) + cot(c + d*x)/(a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_65():
    f = cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = cot(c + d*x)**2/(4*d*(I*a*tan(c + d*x) + a)**2) + 15*I*x/(4*a**2) - 4*log(sin(c + d*x))/(a**2*d) - 2*cot(c + d*x)**2/(a**2*d) + 15*I*cot(c + d*x)/(4*a**2*d) + 5*cot(c + d*x)**2/(4*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_66():
    f = tan(c + d*x)**6/(I*a*tan(c + d*x) + a)**3
    F = 55*tan(c + d*x)**3/(24*d*(I*a**3*tan(c + d*x) + a**3)) - tan(c + d*x)**5/(6*d*(I*a*tan(c + d*x) + a)**3) + 13*I*tan(c + d*x)**4/(24*a*d*(I*a*tan(c + d*x) + a)**2) + 55*x/(8*a**3) + 7*I*log(cos(c + d*x))/(a**3*d) + 7*I*tan(c + d*x)**2/(2*a**3*d) - 55*tan(c + d*x)/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_67():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)**3
    F = 3*tan(c + d*x)**2/(2*d*(I*a**3*tan(c + d*x) + a**3)) - tan(c + d*x)**4/(6*d*(I*a*tan(c + d*x) + a)**3) + 11*I*tan(c + d*x)**3/(24*a*d*(I*a*tan(c + d*x) + a)**2) - 25*I*x/(8*a**3) + 3*log(cos(c + d*x))/(a**3*d) + 25*I*tan(c + d*x)/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_68():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**3
    F = 7*I/(8*d*(I*a**3*tan(c + d*x) + a**3)) - tan(c + d*x)**3/(6*d*(I*a*tan(c + d*x) + a)**3) + 3*I*tan(c + d*x)**2/(8*a*d*(I*a*tan(c + d*x) + a)**2) - 7*x/(8*a**3) - I*log(cos(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_69():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**3
    F = I*tan(c + d*x)**3/(6*d*(I*a*tan(c + d*x) + a)**3) - 1/(8*a*d*(I*a*tan(c + d*x) + a)**2) + I*x/(8*a**3) + 3/(8*a**3*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_70():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = -I/(8*d*(I*a**3*tan(c + d*x) + a**3)) - I/(6*d*(I*a*tan(c + d*x) + a)**3) + 3*I/(8*a*d*(I*a*tan(c + d*x) + a)**2) - x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_71():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = 1/(8*d*(I*a**3*tan(c + d*x) + a**3)) - 1/(6*d*(I*a*tan(c + d*x) + a)**3) + 1/(8*a*d*(I*a*tan(c + d*x) + a)**2) - I*x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_72():
    f = (I*a*tan(c + d*x) + a)**(-3)
    F = I/(8*d*(I*a**3*tan(c + d*x) + a**3)) + I/(6*d*(I*a*tan(c + d*x) + a)**3) + I/(8*a*d*(I*a*tan(c + d*x) + a)**2) + x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_73():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = 7/(8*d*(I*a**3*tan(c + d*x) + a**3)) + 1/(6*d*(I*a*tan(c + d*x) + a)**3) + 3/(8*a*d*(I*a*tan(c + d*x) + a)**2) - 7*I*x/(8*a**3) + log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_74():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = 3*cot(c + d*x)/(2*d*(I*a**3*tan(c + d*x) + a**3)) + cot(c + d*x)/(6*d*(I*a*tan(c + d*x) + a)**3) + 11*cot(c + d*x)/(24*a*d*(I*a*tan(c + d*x) + a)**2) - 25*x/(8*a**3) - 3*I*log(sin(c + d*x))/(a**3*d) - 25*cot(c + d*x)/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_75():
    f = tan(c + d*x)**6/(I*a*tan(c + d*x) + a)**4
    F = -tan(c + d*x)**5/(8*d*(I*a*tan(c + d*x) + a)**4) + 7*I*tan(c + d*x)**4/(24*a*d*(I*a*tan(c + d*x) + a)**3) - 65*x/(16*a**4) - 4*I*log(cos(c + d*x))/(a**4*d) + 65*tan(c + d*x)/(16*a**4*d) - 2*I*tan(c + d*x)**2/(a**4*d*(I*tan(c + d*x) + 1)) + 31*tan(c + d*x)**3/(48*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_76():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)**4
    F = -tan(c + d*x)**4/(8*d*(I*a*tan(c + d*x) + a)**4) + I*tan(c + d*x)**3/(4*a*d*(I*a*tan(c + d*x) + a)**3) + 15*I*x/(16*a**4) - log(cos(c + d*x))/(a**4*d) + 15/(16*a**4*d*(I*tan(c + d*x) + 1)) + 7*tan(c + d*x)**2/(16*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_77():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**4
    F = I/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + I*tan(c + d*x)**4/(8*d*(I*a*tan(c + d*x) + a)**4) + tan(c + d*x)**3/(12*a*d*(I*a*tan(c + d*x) + a)**3) + x/(16*a**4) - 3*I/(16*a**4*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_78():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**4
    F = -1/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + tan(c + d*x)**4/(8*d*(I*a*tan(c + d*x) + a)**4) + I*tan(c + d*x)**3/(12*a*d*(I*a*tan(c + d*x) + a)**3) + I*x/(16*a**4) + 3/(16*a**4*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_79():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = -I/(16*d*(I*a**4*tan(c + d*x) + a**4)) - I/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) - I/(8*d*(I*a*tan(c + d*x) + a)**4) + I/(4*a*d*(I*a*tan(c + d*x) + a)**3) - x/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_80():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = 1/(16*d*(I*a**4*tan(c + d*x) + a**4)) + 1/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) - 1/(8*d*(I*a*tan(c + d*x) + a)**4) + 1/(12*a*d*(I*a*tan(c + d*x) + a)**3) - I*x/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_81():
    f = (I*a*tan(c + d*x) + a)**(-4)
    F = I/(16*d*(I*a**4*tan(c + d*x) + a**4)) + I/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + I/(8*d*(I*a*tan(c + d*x) + a)**4) + I/(12*a*d*(I*a*tan(c + d*x) + a)**3) + x/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_82():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = 1/(8*d*(I*a*tan(c + d*x) + a)**4) + 1/(4*a*d*(I*a*tan(c + d*x) + a)**3) - 15*I*x/(16*a**4) + log(sin(c + d*x))/(a**4*d) + 15/(16*a**4*d*(I*tan(c + d*x) + 1)) + 7/(16*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_83():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = cot(c + d*x)/(8*d*(I*a*tan(c + d*x) + a)**4) + 7*cot(c + d*x)/(24*a*d*(I*a*tan(c + d*x) + a)**3) - 65*x/(16*a**4) - 4*I*log(sin(c + d*x))/(a**4*d) - 65*cot(c + d*x)/(16*a**4*d) + 2*cot(c + d*x)/(a**4*d*(I*tan(c + d*x) + 1)) + 31*cot(c + d*x)/(48*a**4*d*(I*tan(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_84():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**4
    F = -sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(7*d) - 2*I*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(35*d) + 8*I*sqrt(I*a*tan(c + d*x) + a)/(35*d) + 62*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_85():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3
    F = sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(5*d) - 8*sqrt(I*a*tan(c + d*x) + a)/(5*d) - 2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_86():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2
    F = sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_87():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)
    F = -sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_88():
    f = sqrt(I*a*tan(c + d*x) + a)
    F = -sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_89():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_90():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2
    F = -I*sqrt(a)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_91():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3
    F = 7*sqrt(a)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*d) - I*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_92():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**3
    F = 2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*a**2*tan(c + d*x)**4/(7*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*a**2*tan(c + d*x)**3/(7*d*sqrt(I*a*tan(c + d*x) + a)) + 16*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(35*d) - 64*a*sqrt(I*a*tan(c + d*x) + a)/(35*d) - 76*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_93():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**2
    F = 2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*I*a*sqrt(I*a*tan(c + d*x) + a)/d - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_94():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)
    F = -2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*a*sqrt(I*a*tan(c + d*x) + a)/d + 2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_95():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 2*I*a*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_96():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_97():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**2
    F = -3*I*a**(sympy.S(3)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**2*cot(c + d*x)/(d*sqrt(I*a*tan(c + d*x) + a)) - I*a**2/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_98():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3
    F = 11*a**(sympy.S(3)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - 2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**2*cot(c + d*x)**2/(2*d*sqrt(I*a*tan(c + d*x) + a)) - I*a**2*cot(c + d*x)/(2*d*sqrt(I*a*tan(c + d*x) + a)) - 5*I*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_99():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**3
    F = 4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**4/(9*d) + 38*I*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(63*d) + 92*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(105*d) - 368*a**2*sqrt(I*a*tan(c + d*x) + a)/(105*d) - 472*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(315*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_100():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**2
    F = 4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 4*I*a**2*sqrt(I*a*tan(c + d*x) + a)/d - 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_101():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)
    F = -4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 4*a**2*sqrt(I*a*tan(c + d*x) + a)/d + 2*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_102():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 4*I*a**2*sqrt(I*a*tan(c + d*x) + a)/d + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_103():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 2*a**2*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_104():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**2
    F = -5*I*a**(sympy.S(5)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/d + 4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_105():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**3
    F = 23*a**(sympy.S(5)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*d) - 4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*d) - 9*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_106():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**4
    F = 45*I*a**(sympy.S(5)/2)*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(8*d) - 4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**3/(3*d) - 13*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(12*d) + 19*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_107():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -8*sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d + 8*I*a**3*sqrt(I*a*tan(c + d*x) + a)/d + 4*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_108():
    f = tan(c + d*x)**5/sqrt(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**4/(d*sqrt(I*a*tan(c + d*x) + a)) - 9*I*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**3/(7*a*d) + 47*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(35*a*d) - 188*sqrt(I*a*tan(c + d*x) + a)/(35*a*d) + 223*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(105*a**2*d) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_109():
    f = tan(c + d*x)**4/sqrt(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**3/(d*sqrt(I*a*tan(c + d*x) + a)) - 7*I*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(5*a*d) + 28*I*sqrt(I*a*tan(c + d*x) + a)/(5*a*d) - 23*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(15*a**2*d) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_110():
    f = tan(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**2/(d*sqrt(I*a*tan(c + d*x) + a)) + 4*sqrt(I*a*tan(c + d*x) + a)/(a*d) - 5*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_111():
    f = tan(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = -I/(d*sqrt(I*a*tan(c + d*x) + a)) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(a*d) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_112():
    f = tan(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = -1/(d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_113():
    f = 1/sqrt(I*a*tan(c + d*x) + a)
    F = I/(d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_114():
    f = cot(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = 1/(d*sqrt(I*a*tan(c + d*x) + a)) - 2*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_115():
    f = cot(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)/(d*sqrt(I*a*tan(c + d*x) + a)) - 2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(a*d) + I*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_116():
    f = cot(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)**2/(d*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(2*a*d) + 7*I*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a*d) + 11*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_117():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)**4/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 19*I*tan(c + d*x)**3/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 39*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**2/(10*a**2*d) + 78*sqrt(I*a*tan(c + d*x) + a)/(5*a**2*d) - 151*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(30*a**3*d) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_118():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)**3/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*I*tan(c + d*x)**2/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - 10*I*sqrt(I*a*tan(c + d*x) + a)/(a**2*d) + 7*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(2*a**3*d) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_119():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)**2/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 11/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*sqrt(I*a*tan(c + d*x) + a)/(3*a**2*d) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_120():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -I/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 3*I/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_121():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 1/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_122():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-3)/2)
    F = I/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + I/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_123():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 3/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - 2*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_124():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = cot(c + d*x)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 13*cot(c + d*x)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(2*a**2*d) + 3*I*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_125():
    f = cot(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = cot(c + d*x)**2/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 17*cot(c + d*x)**2/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 11*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**2/(3*a**2*d) + 21*I*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a**2*d) + 23*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(4*a**(sympy.S(3)/2)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_126():
    f = tan(c + d*x)**5/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)**4/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 7*I*tan(c + d*x)**3/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 89*tan(c + d*x)**2/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 89*sqrt(I*a*tan(c + d*x) + a)/(5*a**3*d) + 361*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(60*a**4*d) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_127():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)**3/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 17*I*tan(c + d*x)**2/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 151*I/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 83*I*sqrt(I*a*tan(c + d*x) + a)/(30*a**3*d) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_128():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)**2/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - 13/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 31/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_129():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -I/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I/(2*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - I/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_130():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 1/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 1/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_131():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-5)/2)
    F = I/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + I/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_132():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 1/(2*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 7/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 2*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_133():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = cot(c + d*x)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 19*cot(c + d*x)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 41*cot(c + d*x)/(12*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 21*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)/(4*a**3*d) + 5*I*atanh(sqrt(I*a*tan(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) + sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_134():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-7)/2)
    F = I/(7*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + I/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I/(12*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + I/(8*a**3*d*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(16*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_135():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)
    F = 2*(-1)**(sympy.S(1)/4)*a*d**(sympy.S(5)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f - 2*I*a*d**2*sqrt(d*tan(e + f*x))/f + 2*a*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*I*a*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_136():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)
    F = 2*(-1)**(sympy.S(3)/4)*a*d**(sympy.S(3)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 2*a*d*sqrt(d*tan(e + f*x))/f + 2*I*a*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_137():
    f = sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)
    F = -2*(-1)**(sympy.S(1)/4)*a*sqrt(d)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 2*I*a*sqrt(d*tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_138():
    f = (I*a*tan(e + f*x) + a)/sqrt(d*tan(e + f*x))
    F = -2*(-1)**(sympy.S(3)/4)*a*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_139():
    f = (I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a/(d*f*sqrt(d*tan(e + f*x))) + 2*(-1)**(sympy.S(1)/4)*a*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_140():
    f = (I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) - 2*I*a/(d**2*f*sqrt(d*tan(e + f*x))) + 2*(-1)**(sympy.S(3)/4)*a*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_141():
    f = (I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -2*a/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2)) - 2*I*a/(3*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 2*a/(d**3*f*sqrt(d*tan(e + f*x))) - 2*(-1)**(sympy.S(1)/4)*a*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_142():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(-I*a*tan(e + f*x) + a)
    F = -2*(-1)**(sympy.S(1)/4)*a*d**(sympy.S(5)/2)*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 2*I*a*d**2*sqrt(d*tan(e + f*x))/f + 2*a*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) - 2*I*a*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_143():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(-I*a*tan(e + f*x) + a)
    F = 2*(-1)**(sympy.S(3)/4)*a*d**(sympy.S(3)/2)*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 2*a*d*sqrt(d*tan(e + f*x))/f - 2*I*a*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_144():
    f = sqrt(d*tan(e + f*x))*(-I*a*tan(e + f*x) + a)
    F = 2*(-1)**(sympy.S(1)/4)*a*sqrt(d)*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f - 2*I*a*sqrt(d*tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_145():
    f = (-I*a*tan(e + f*x) + a)/sqrt(d*tan(e + f*x))
    F = -2*(-1)**(sympy.S(3)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_146():
    f = (-I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a/(d*f*sqrt(d*tan(e + f*x))) - 2*(-1)**(sympy.S(1)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_147():
    f = (-I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 2*I*a/(d**2*f*sqrt(d*tan(e + f*x))) + 2*(-1)**(sympy.S(3)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_148():
    f = (-I*a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -2*a/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2)) + 2*I*a/(3*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 2*a/(d**3*f*sqrt(d*tan(e + f*x))) + 2*(-1)**(sympy.S(1)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_149():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**2
    F = 4*(-1)**(sympy.S(1)/4)*a**2*d**(sympy.S(5)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f - 4*I*a**2*d**2*sqrt(d*tan(e + f*x))/f + 4*a**2*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 4*I*a**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) - 2*a**2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_150():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2
    F = 4*(-1)**(sympy.S(3)/4)*a**2*d**(sympy.S(3)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 4*a**2*d*sqrt(d*tan(e + f*x))/f + 4*I*a**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) - 2*a**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_151():
    f = sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2
    F = -4*(-1)**(sympy.S(1)/4)*a**2*sqrt(d)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 4*I*a**2*sqrt(d*tan(e + f*x))/f - 2*a**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_152():
    f = (I*a*tan(e + f*x) + a)**2/sqrt(d*tan(e + f*x))
    F = -2*a**2*sqrt(d*tan(e + f*x))/(d*f) - 4*(-1)**(sympy.S(3)/4)*a**2*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_153():
    f = (I*a*tan(e + f*x) + a)**2/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a**2/(d*f*sqrt(d*tan(e + f*x))) + 4*(-1)**(sympy.S(1)/4)*a**2*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_154():
    f = (I*a*tan(e + f*x) + a)**2/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a**2/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) - 4*I*a**2/(d**2*f*sqrt(d*tan(e + f*x))) + 4*(-1)**(sympy.S(3)/4)*a**2*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_155():
    f = (I*a*tan(e + f*x) + a)**2/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -2*a**2/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2)) - 4*I*a**2/(3*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 4*a**2/(d**3*f*sqrt(d*tan(e + f*x))) - 4*(-1)**(sympy.S(1)/4)*a**2*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_156():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**3
    F = 8*(-1)**(sympy.S(1)/4)*a**3*d**(sympy.S(5)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f - 8*I*a**3*d**2*sqrt(d*tan(e + f*x))/f + 8*a**3*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 8*I*a**3*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) - 40*a**3*(d*tan(e + f*x))**(sympy.S(7)/2)/(63*d*f) - 2*(d*tan(e + f*x))**(sympy.S(7)/2)*(I*a**3*tan(e + f*x) + a**3)/(9*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_157():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**3
    F = 8*(-1)**(sympy.S(3)/4)*a**3*d**(sympy.S(3)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 8*a**3*d*sqrt(d*tan(e + f*x))/f + 8*I*a**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) - 32*a**3*(d*tan(e + f*x))**(sympy.S(5)/2)/(35*d*f) - 2*(d*tan(e + f*x))**(sympy.S(5)/2)*(I*a**3*tan(e + f*x) + a**3)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_158():
    f = sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3
    F = -8*(-1)**(sympy.S(1)/4)*a**3*sqrt(d)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 8*I*a**3*sqrt(d*tan(e + f*x))/f - 8*a**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(5*d*f) - 2*(d*tan(e + f*x))**(sympy.S(3)/2)*(I*a**3*tan(e + f*x) + a**3)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_159():
    f = (I*a*tan(e + f*x) + a)**3/sqrt(d*tan(e + f*x))
    F = -16*a**3*sqrt(d*tan(e + f*x))/(3*d*f) - 8*(-1)**(sympy.S(3)/4)*a**3*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f) - 2*sqrt(d*tan(e + f*x))*(I*a**3*tan(e + f*x) + a**3)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_160():
    f = (I*a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 8*(-1)**(sympy.S(1)/4)*a**3*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f) - (2*I*a**3*tan(e + f*x) + 2*a**3)/(d*f*sqrt(d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_161():
    f = (I*a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -16*I*a**3/(3*d**2*f*sqrt(d*tan(e + f*x))) + 8*(-1)**(sympy.S(3)/4)*a**3*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f) - (2*I*a**3*tan(e + f*x) + 2*a**3)/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_162():
    f = (I*a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -8*I*a**3/(5*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 8*a**3/(d**3*f*sqrt(d*tan(e + f*x))) - 8*(-1)**(sympy.S(1)/4)*a**3*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(7)/2)*f) - (2*I*a**3*tan(e + f*x) + 2*a**3)/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_163():
    f = (I*a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(9)/2)
    F = -32*I*a**3/(35*d**2*f*(d*tan(e + f*x))**(sympy.S(5)/2)) + 8*a**3/(3*d**3*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 8*I*a**3/(d**4*f*sqrt(d*tan(e + f*x))) - 8*(-1)**(sympy.S(3)/4)*a**3*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(9)/2)*f) - (2*I*a**3*tan(e + f*x) + 2*a**3)/(7*d*f*(d*tan(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_164():
    f = (d*tan(e + f*x))**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)
    F = -d*(d*tan(e + f*x))**(sympy.S(5)/2)/(2*f*(I*a*tan(e + f*x) + a)) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/8 + 7*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/8 + 7*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/4 - 7*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/4 - 7*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f) + 5*d**3*sqrt(d*tan(e + f*x))/(2*a*f) - 7*I*d**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(6*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_165():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)
    F = -d*(d*tan(e + f*x))**(sympy.S(3)/2)/(2*f*(I*a*tan(e + f*x) + a)) + sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(3)/8 - 5*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) - sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(3)/8 - 5*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) - sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(3)/4 + 5*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f) + sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(3)/4 + 5*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f) - 5*I*d**2*sqrt(d*tan(e + f*x))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_166():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)
    F = -d*sqrt(d*tan(e + f*x))/(2*f*(I*a*tan(e + f*x) + a)) - sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/8 + 3*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) + sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/8 + 3*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*f) - sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/4 - 3*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f) + sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/4 - 3*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_167():
    f = sqrt(d*tan(e + f*x))/(I*a*tan(e + f*x) + a)
    F = I*sqrt(d*tan(e + f*x))/(2*f*(I*a*tan(e + f*x) + a)) - (-1)**(sympy.S(1)/4)*sqrt(d)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_168():
    f = 1/(sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a))
    F = sqrt(d*tan(e + f*x))/(2*d*f*(I*a*tan(e + f*x) + a)) - sqrt(2)*(sympy.S(3)/8 + I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*sqrt(d)*f) + sqrt(2)*(sympy.S(3)/8 + I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*sqrt(d)*f) - sqrt(2)*(sympy.S(3)/4 - I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*sqrt(d)*f) + sqrt(2)*(sympy.S(3)/4 - I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_169():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a))
    F = 1/(2*d*f*sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)) - 5/(2*a*d*f*sqrt(d*tan(e + f*x))) - sqrt(2)*(sympy.S(5)/8 - 3*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(5)/8 - 3*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(5)/4 + 3*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*d**(sympy.S(3)/2)*f) - sqrt(2)*(sympy.S(5)/4 + 3*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_170():
    f = 1/((d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a))
    F = 1/(2*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)) - 7/(6*a*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 5*I/(2*a*d**2*f*sqrt(d*tan(e + f*x))) + sqrt(2)*(sympy.S(7)/8 + 5*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*d**(sympy.S(5)/2)*f) - sqrt(2)*(sympy.S(7)/8 + 5*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a*d**(sympy.S(5)/2)*f) + sqrt(2)*(sympy.S(7)/4 - 5*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*d**(sympy.S(5)/2)*f) - sqrt(2)*(sympy.S(7)/4 - 5*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_171():
    f = (d*tan(e + f*x))**(sympy.S(9)/2)/(I*a*tan(e + f*x) + a)**2
    F = -d*(d*tan(e + f*x))**(sympy.S(7)/2)/(4*f*(I*a*tan(e + f*x) + a)**2) + sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(49)/32 - 45*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(49)/32 - 45*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(49)/16 + 45*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(49)/16 + 45*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - 45*I*d**4*sqrt(d*tan(e + f*x))/(8*a**2*f) - 49*d**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(24*a**2*f) + 9*I*d**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_172():
    f = (d*tan(e + f*x))**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**2
    F = -d*(d*tan(e + f*x))**(sympy.S(5)/2)/(4*f*(I*a*tan(e + f*x) + a)**2) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(25)/32 + 21*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(25)/32 + 21*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(25)/16 - 21*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(25)/16 - 21*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - 25*d**3*sqrt(d*tan(e + f*x))/(8*a**2*f) + 7*I*d**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_173():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**2
    F = -d*(d*tan(e + f*x))**(sympy.S(3)/2)/(4*f*(I*a*tan(e + f*x) + a)**2) - sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(9)/32 - 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) + sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(9)/32 - 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) + sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(9)/16 + 5*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - sqrt(2)*d**(sympy.S(5)/2)*(sympy.S(9)/16 + 5*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + 5*I*d**2*sqrt(d*tan(e + f*x))/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_174():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**2
    F = -d*sqrt(d*tan(e + f*x))/(4*f*(I*a*tan(e + f*x) + a)**2) + sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/32 - 3*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/32 - 3*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) + sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/16 + 3*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - sqrt(2)*d**(sympy.S(3)/2)*(sympy.S(1)/16 + 3*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + 3*d*sqrt(d*tan(e + f*x))/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_175():
    f = sqrt(d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**2
    F = I*sqrt(d*tan(e + f*x))/(4*f*(I*a*tan(e + f*x) + a)**2) + sqrt(2)*sqrt(d)*(sympy.S(1)/32 + 3*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*sqrt(d)*(sympy.S(1)/32 + 3*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*f) - sqrt(2)*sqrt(d)*(sympy.S(1)/16 - 3*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + sqrt(2)*sqrt(d)*(sympy.S(1)/16 - 3*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + I*sqrt(d*tan(e + f*x))/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_176():
    f = 1/(sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2)
    F = sqrt(d*tan(e + f*x))/(4*d*f*(I*a*tan(e + f*x) + a)**2) + 5*sqrt(d*tan(e + f*x))/(8*a**2*d*f*(I*tan(e + f*x) + 1)) - sqrt(2)*(sympy.S(9)/32 + 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*sqrt(d)*f) + sqrt(2)*(sympy.S(9)/32 + 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*sqrt(d)*f) - sqrt(2)*(sympy.S(9)/16 - 5*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*sqrt(d)*f) + sqrt(2)*(sympy.S(9)/16 - 5*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_177():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2)
    F = 1/(4*d*f*sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2) - 25/(8*a**2*d*f*sqrt(d*tan(e + f*x))) + 7/(8*a**2*d*f*sqrt(d*tan(e + f*x))*(I*tan(e + f*x) + 1)) - sqrt(2)*(sympy.S(25)/32 - 21*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(25)/32 - 21*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(25)/16 + 21*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(3)/2)*f) - sqrt(2)*(sympy.S(25)/16 + 21*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_178():
    f = 1/((d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**2)
    F = 1/(4*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2) - 49/(24*a**2*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 9/(8*a**2*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(I*tan(e + f*x) + 1)) + 45*I/(8*a**2*d**2*f*sqrt(d*tan(e + f*x))) + sqrt(2)*(sympy.S(49)/32 + 45*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*d**(sympy.S(5)/2)*f) - sqrt(2)*(sympy.S(49)/32 + 45*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**2*d**(sympy.S(5)/2)*f) + sqrt(2)*(sympy.S(49)/16 - 45*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(5)/2)*f) - sqrt(2)*(sympy.S(49)/16 - 45*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_179():
    f = (d*tan(e + f*x))**(sympy.S(9)/2)/(I*a*tan(e + f*x) + a)**3
    F = 7*d**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(6*f*(I*a**3*tan(e + f*x) + a**3)) - d*(d*tan(e + f*x))**(sympy.S(7)/2)/(6*f*(I*a*tan(e + f*x) + a)**3) + 5*I*d**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(12*a*f*(I*a*tan(e + f*x) + a)**2) - sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(7)/8 - 15*I/16)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*f) + sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(7)/8 - 15*I/16)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*f) + sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(7)/4 + 15*I/8)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*f) - sqrt(2)*d**(sympy.S(9)/2)*(sympy.S(7)/4 + 15*I/8)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*f) + 15*I*d**4*sqrt(d*tan(e + f*x))/(4*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_180():
    f = (d*tan(e + f*x))**(sympy.S(7)/2)/(I*a*tan(e + f*x) + a)**3
    F = 5*d**3*sqrt(d*tan(e + f*x))/(8*f*(I*a**3*tan(e + f*x) + a**3)) - d*(d*tan(e + f*x))**(sympy.S(5)/2)/(6*f*(I*a*tan(e + f*x) + a)**3) + I*d**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*a*f*(I*a*tan(e + f*x) + a)**2) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/32 + 7*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*f) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/32 + 7*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*f) + sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/16 - 7*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*f) - sqrt(2)*d**(sympy.S(7)/2)*(sympy.S(5)/16 - 7*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_181():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**3
    F = -I*d**2*sqrt(d*tan(e + f*x))/(4*f*(I*a**3*tan(e + f*x) + a**3)) - d*(d*tan(e + f*x))**(sympy.S(3)/2)/(6*f*(I*a*tan(e + f*x) + a)**3) + I*d**2*sqrt(d*tan(e + f*x))/(4*a*f*(I*a*tan(e + f*x) + a)**2) - sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(32*a**3*f) + sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(32*a**3*f) + sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(16*a**3*f) - sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(16*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_182():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**3
    F = d*sqrt(d*tan(e + f*x))/(8*f*(I*a**3*tan(e + f*x) + a**3)) - d*sqrt(d*tan(e + f*x))/(6*f*(I*a*tan(e + f*x) + a)**3) + d*sqrt(d*tan(e + f*x))/(6*a*f*(I*a*tan(e + f*x) + a)**2) + (-1)**(sympy.S(3)/4)*d**(sympy.S(3)/2)*atanh((-1)**(sympy.S(1)/4)*sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_183():
    f = sqrt(d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**3
    F = I*sqrt(d*tan(e + f*x))/(6*f*(I*a*tan(e + f*x) + a)**3) + I*sqrt(d*tan(e + f*x))/(12*a*f*(I*a*tan(e + f*x) + a)**2) + sqrt(2)*I*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(32*a**3*f) - sqrt(2)*I*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(32*a**3*f) + sqrt(2)*I*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(16*a**3*f) - sqrt(2)*I*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(16*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_184():
    f = 1/(sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3)
    F = 5*sqrt(d*tan(e + f*x))/(8*d*f*(I*a**3*tan(e + f*x) + a**3)) + sqrt(d*tan(e + f*x))/(6*d*f*(I*a*tan(e + f*x) + a)**3) + sqrt(d*tan(e + f*x))/(3*a*d*f*(I*a*tan(e + f*x) + a)**2) - sqrt(2)*(sympy.S(7)/32 + 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*sqrt(d)*f) + sqrt(2)*(sympy.S(7)/32 + 5*I/32)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*sqrt(d)*f) - sqrt(2)*(sympy.S(7)/16 - 5*I/16)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*sqrt(d)*f) + sqrt(2)*(sympy.S(7)/16 - 5*I/16)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_185():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**3)
    F = 7/(6*d*f*sqrt(d*tan(e + f*x))*(I*a**3*tan(e + f*x) + a**3)) + 1/(6*d*f*sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3) + 5/(12*a*d*f*sqrt(d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2) - 15/(4*a**3*d*f*sqrt(d*tan(e + f*x))) - sqrt(2)*(sympy.S(15)/16 - 7*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(15)/16 - 7*I/8)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*a**3*d**(sympy.S(3)/2)*f) + sqrt(2)*(sympy.S(15)/8 + 7*I/4)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*d**(sympy.S(3)/2)*f) - sqrt(2)*(sympy.S(15)/8 + 7*I/4)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**3*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_186():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)
    F = 7*(-1)**(sympy.S(3)/4)*sqrt(a)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) + sqrt(a)*(1 + I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(2*d) - I*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_187():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)
    F = -(-1)**(sympy.S(1)/4)*sqrt(a)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - sqrt(a)*(1 - I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_188():
    f = sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))
    F = -2*(-1)**(sympy.S(3)/4)*sqrt(a)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - sqrt(a)*(1 + I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_189():
    f = sqrt(I*a*tan(c + d*x) + a)/sqrt(tan(c + d*x))
    F = sqrt(a)*(1 - I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_190():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a)*(1 + I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_191():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(5)/2)
    F = sqrt(a)*(-1 + I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*I*sqrt(I*a*tan(c + d*x) + a)/(3*d*sqrt(tan(c + d*x))) - 2*sqrt(I*a*tan(c + d*x) + a)/(3*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_192():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(7)/2)
    F = sqrt(a)*(-1 - I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 26*sqrt(I*a*tan(c + d*x) + a)/(15*d*sqrt(tan(c + d*x))) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(15*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*sqrt(I*a*tan(c + d*x) + a)/(5*d*tan(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_193():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2)
    F = 23*(-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(8*d) + a**(sympy.S(3)/2)*(2 + 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*tan(c + d*x)**(sympy.S(7)/2)/(3*d*sqrt(I*a*tan(c + d*x) + a)) + I*a**2*tan(c + d*x)**(sympy.S(5)/2)/(3*d*sqrt(I*a*tan(c + d*x) + a)) + 7*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(12*d) - 9*I*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_194():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = -11*(-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) - a**(sympy.S(3)/2)*(2 - 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*tan(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(I*a*tan(c + d*x) + a)) + I*a**2*tan(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(I*a*tan(c + d*x) + a)) + 5*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_195():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))
    F = -3*(-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**(sympy.S(3)/2)*(2 + 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*tan(c + d*x)**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) + I*a**2*sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_196():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/sqrt(tan(c + d*x))
    F = 2*(-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(3)/2)*(2 - 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_197():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*(2 + 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_198():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = a**(sympy.S(3)/2)*(-2 + 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*I*a*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) - 2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_199():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = a**(sympy.S(3)/2)*(-2 - 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*I*a**2/(5*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2/(5*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)) + 12*a*sqrt(I*a*tan(c + d*x) + a)/(5*d*sqrt(tan(c + d*x))) - 4*I*a*sqrt(I*a*tan(c + d*x) + a)/(5*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_200():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = a**(sympy.S(3)/2)*(2 - 2*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*I*a**2/(7*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)) - 2*a**2/(7*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/2)) + 268*I*a*sqrt(I*a*tan(c + d*x) + a)/(105*d*sqrt(tan(c + d*x))) + 76*a*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(3)/2)) - 16*I*a*sqrt(I*a*tan(c + d*x) + a)/(35*d*tan(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_201():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)
    F = 363*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(64*d) + a**(sympy.S(5)/2)*(4 + 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/2)/(4*d) + 17*I*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)/(24*d) + 107*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(96*d) - 149*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(64*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_202():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = -45*(-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(8*d) - a**(sympy.S(5)/2)*(4 - 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)/(3*d) + 13*I*a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(12*d) + 19*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_203():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x))
    F = -23*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) - a**(sympy.S(5)/2)*(4 + 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(2*d) + 9*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_204():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/sqrt(tan(c + d*x))
    F = 5*(-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 - 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_205():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = 2*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 + 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_206():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = a**(sympy.S(5)/2)*(-4 + 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 4*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) - 2*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_207():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = a**(sympy.S(5)/2)*(-4 - 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a**2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(tan(c + d*x))) - 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d*tan(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_208():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = a**(sympy.S(5)/2)*(4 - 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 104*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(21*d*sqrt(tan(c + d*x))) + 32*a**2*sqrt(I*a*tan(c + d*x) + a)/(21*d*tan(c + d*x)**(sympy.S(3)/2)) - 6*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(7*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*a**2*sqrt(I*a*tan(c + d*x) + a)/(7*d*tan(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_209():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = a**(sympy.S(5)/2)*(4 + 4*I)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 1576*a**2*sqrt(I*a*tan(c + d*x) + a)/(315*d*sqrt(tan(c + d*x))) + 472*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(315*d*tan(c + d*x)**(sympy.S(3)/2)) + 92*a**2*sqrt(I*a*tan(c + d*x) + a)/(105*d*tan(c + d*x)**(sympy.S(5)/2)) - 38*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(63*d*tan(c + d*x)**(sympy.S(7)/2)) - 2*a**2*sqrt(I*a*tan(c + d*x) + a)/(9*d*tan(c + d*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_210():
    f = tan(c + d*x)**(sympy.S(7)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**(sympy.S(5)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) - 3*I*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)/(2*a*d) + 7*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*a*d) + 11*(-1)**(sympy.S(1)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*sqrt(a)*d) + (sympy.S.Half - I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_211():
    f = tan(c + d*x)**(sympy.S(5)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) - 2*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(a*d) - (-1)**(sympy.S(3)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) + (sympy.S.Half + I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_212():
    f = tan(c + d*x)**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = -sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) - 2*(-1)**(sympy.S(1)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) - (sympy.S.Half - I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_213():
    f = sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = I*sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(-1)/2 - I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_214():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x)))
    F = sqrt(tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S.Half - I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_215():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2))
    F = 1/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - 3*sqrt(I*a*tan(c + d*x) + a)/(a*d*sqrt(tan(c + d*x))) + (sympy.S.Half + I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_216():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2))
    F = 1/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) + 7*I*sqrt(I*a*tan(c + d*x) + a)/(3*a*d*sqrt(tan(c + d*x))) - 5*sqrt(I*a*tan(c + d*x) + a)/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) - (sympy.S.Half - I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_217():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/2))
    F = 1/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/2)) + 61*sqrt(I*a*tan(c + d*x) + a)/(15*a*d*sqrt(tan(c + d*x))) + 23*I*sqrt(I*a*tan(c + d*x) + a)/(15*a*d*tan(c + d*x)**(sympy.S(3)/2)) - 7*sqrt(I*a*tan(c + d*x) + a)/(5*a*d*tan(c + d*x)**(sympy.S(5)/2)) - (sympy.S.Half + I/2)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_218():
    f = tan(c + d*x)**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)**(sympy.S(5)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 13*I*tan(c + d*x)**(sympy.S(3)/2)/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(2*a**2*d) - 3*(-1)**(sympy.S(1)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (sympy.S(1)/4 - I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_219():
    f = tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 3*I*sqrt(tan(c + d*x))/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) + 2*(-1)**(sympy.S(3)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (sympy.S(1)/4 + I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_220():
    f = tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = I*tan(c + d*x)**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(tan(c + d*x))/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(-1)/4 + I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_221():
    f = sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = tan(c + d*x)**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + I*sqrt(tan(c + d*x))/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - (sympy.S(1)/4 + I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_222():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x)))
    F = sqrt(tan(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 7*sqrt(tan(c + d*x))/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/4 - I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_223():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = 1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + 11/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - 25*sqrt(I*a*tan(c + d*x) + a)/(6*a**2*d*sqrt(tan(c + d*x))) + (sympy.S(1)/4 + I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_224():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = 1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 5/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) + 13*I*sqrt(I*a*tan(c + d*x) + a)/(2*a**2*d*sqrt(tan(c + d*x))) - 7*sqrt(I*a*tan(c + d*x) + a)/(2*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + (sympy.S(-1)/4 + I/4)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_225():
    f = tan(c + d*x)**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)**(sympy.S(7)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 19*I*tan(c + d*x)**(sympy.S(5)/2)/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 41*tan(c + d*x)**(sympy.S(3)/2)/(12*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 21*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))/(4*a**3*d) + 5*(-1)**(sympy.S(3)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - (sympy.S(1)/8 + I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_226():
    f = tan(c + d*x)**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)**(sympy.S(5)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I*tan(c + d*x)**(sympy.S(3)/2)/(2*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 7*sqrt(tan(c + d*x))/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 2*(-1)**(sympy.S(1)/4)*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (sympy.S(1)/8 - I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_227():
    f = tan(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*tan(c + d*x)**(sympy.S(5)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + tan(c + d*x)**(sympy.S(3)/2)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - I*sqrt(tan(c + d*x))/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/8 + I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_228():
    f = tan(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = tan(c + d*x)**(sympy.S(5)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I*tan(c + d*x)**(sympy.S(3)/2)/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(tan(c + d*x))/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(-1)/8 + I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_229():
    f = sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*sqrt(tan(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + I*sqrt(tan(c + d*x))/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - I*sqrt(tan(c + d*x))/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(-1)/8 - I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_230():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x)))
    F = sqrt(tan(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 13*sqrt(tan(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 67*sqrt(tan(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + (sympy.S(1)/8 - I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_231():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = 1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(tan(c + d*x))) + 17/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + 151/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(tan(c + d*x))) - 317*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*sqrt(tan(c + d*x))) + (sympy.S(1)/8 + I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_232():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = 1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 7/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 89/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(3)/2)) + 707*I*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*sqrt(tan(c + d*x))) - 361*sqrt(I*a*tan(c + d*x) + a)/(60*a**3*d*tan(c + d*x)**(sympy.S(3)/2)) + (sympy.S(-1)/8 + I/8)*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_233():
    f = tan(c + d*x)**(sympy.S(10)/3)/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**(sympy.S(7)/3)/(2*d*(I*a*tan(c + d*x) + a)) - 5*I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(6*a*d) + 7*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - 7*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) + 5*I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(12*a*d) - 5*I*tan(c + d*x)**(sympy.S(4)/3)/(4*a*d) + 7*tan(c + d*x)**(sympy.S(1)/3)/(2*a*d) - 5*sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(6*a*d) - 7*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) - 7*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) - 7*atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_234():
    f = tan(c + d*x)**(sympy.S(8)/3)/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**(sympy.S(5)/3)/(2*d*(I*a*tan(c + d*x) + a)) + 2*I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) + 5*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - 5*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) - 2*I*tan(c + d*x)**(sympy.S(2)/3)/(a*d) - 2*sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(3*a*d) + 5*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) + 5*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) + 5*atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_235():
    f = tan(c + d*x)**(sympy.S(4)/3)/(I*a*tan(c + d*x) + a)
    F = -tan(c + d*x)**(sympy.S(1)/3)/(2*d*(I*a*tan(c + d*x) + a)) + I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) - sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) + sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(6*a*d) + sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(3*a*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) + atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_236():
    f = tan(c + d*x)**(sympy.S(2)/3)/(I*a*tan(c + d*x) + a)
    F = I*tan(c + d*x)**(sympy.S(2)/3)/(2*d*(I*a*tan(c + d*x) + a)) - I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(6*a*d) + sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) + I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(12*a*d) + sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(6*a*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) + atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_237():
    f = 1/((I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3))
    F = tan(c + d*x)**(sympy.S(2)/3)/(2*d*(I*a*tan(c + d*x) + a)) + log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) - sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) + sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(6*a*d) - sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(3*a*d) - I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) - I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) - I*atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_238():
    f = 1/((I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/3))
    F = 1/(2*d*(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(2)/3)) + 2*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) + 5*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - 5*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(3*a*d) + 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(3*a*d) - 5*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) - 5*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) - 5*I*atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d) - 2/(a*d*tan(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_239():
    f = 1/((I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/3))
    F = 1/(2*d*(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(4)/3)) - 5*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(6*a*d) + 7*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) - 7*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(24*a*d) + 5*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(12*a*d) + 5*sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(6*a*d) + 7*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(12*a*d) + 7*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(12*a*d) + 7*I*atan(tan(c + d*x)**(sympy.S(1)/3))/(6*a*d) + 7*I/(2*a*d*tan(c + d*x)**(sympy.S(1)/3)) - 5/(4*a*d*tan(c + d*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_240():
    f = tan(c + d*x)**(sympy.S(14)/3)/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**(sympy.S(11)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) + 14*I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + 121*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 121*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 7*I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) - 121*tan(c + d*x)**(sympy.S(5)/3)/(60*a**2*d) - 14*I*tan(c + d*x)**(sympy.S(2)/3)/(3*a**2*d) - 14*sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) + 121*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) + 121*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) + 121*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + 7*I*tan(c + d*x)**(sympy.S(8)/3)/(6*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_241():
    f = tan(c + d*x)**(sympy.S(10)/3)/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**(sympy.S(7)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) + 5*I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) - 49*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + 49*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 5*I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(18*a**2*d) - 49*tan(c + d*x)**(sympy.S(1)/3)/(12*a**2*d) + 5*sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) + 49*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) + 49*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) + 49*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + 5*I*tan(c + d*x)**(sympy.S(4)/3)/(6*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_242():
    f = tan(c + d*x)**(sympy.S(8)/3)/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**(sympy.S(5)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) - 2*I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) - 25*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + 25*sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + 2*sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) - 25*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) - 25*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) - 25*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + 2*I*tan(c + d*x)**(sympy.S(2)/3)/(3*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_243():
    f = tan(c + d*x)**(sympy.S(4)/3)/(I*a*tan(c + d*x) + a)**2
    F = -tan(c + d*x)**(sympy.S(1)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) + I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(18*a**2*d) + sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) - atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) - atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) - atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + tan(c + d*x)**(sympy.S(1)/3)/(3*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_244():
    f = tan(c + d*x)**(sympy.S(2)/3)/(I*a*tan(c + d*x) + a)**2
    F = tan(c + d*x)**(sympy.S(5)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) - I*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - sqrt(3)*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + I*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(18*a**2*d) + sqrt(3)*I*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) + atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) + atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + I*tan(c + d*x)**(sympy.S(2)/3)/(3*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_245():
    f = 1/((I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(1)/3))
    F = tan(c + d*x)**(sympy.S(2)/3)/(4*d*(I*a*tan(c + d*x) + a)**2) + 2*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) - 7*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + 7*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) - 7*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) - 7*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) - 7*I*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + 7*tan(c + d*x)**(sympy.S(2)/3)/(12*a**2*d*(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_246():
    f = 1/((I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(5)/3))
    F = 1/(4*d*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(2)/3)) + 8*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + 55*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 55*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 4*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(9*a**2*d) + 8*sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(9*a**2*d) - 55*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) - 55*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) - 55*I*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) - 8/(3*a**2*d*tan(c + d*x)**(sympy.S(2)/3)) + 11/(12*a**2*d*(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_247():
    f = 1/((I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(7)/3))
    F = 1/(4*d*(I*a*tan(c + d*x) + a)**2*tan(c + d*x)**(sympy.S(4)/3)) - 25*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(18*a**2*d) + 91*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) - 91*sqrt(3)*I*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(144*a**2*d) + 25*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(36*a**2*d) + 25*sqrt(3)*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(18*a**2*d) + 91*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(72*a**2*d) + 91*I*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(72*a**2*d) + 91*I*atan(tan(c + d*x)**(sympy.S(1)/3))/(36*a**2*d) + 91*I/(12*a**2*d*tan(c + d*x)**(sympy.S(1)/3)) - 25/(12*a**2*d*tan(c + d*x)**(sympy.S(4)/3)) + 13/(12*a**2*d*(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_248():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(4)/3)
    F = 3*a*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, sympy.S.Half, 1, sympy.S(10)/3, -I*tan(c + d*x), I*tan(c + d*x))/(7*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_249():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(2)/3)
    F = 3*a*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, sympy.S.Half, 1, sympy.S(8)/3, -I*tan(c + d*x), I*tan(c + d*x))/(5*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_250():
    f = sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3)
    F = 3*a*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -I*tan(c + d*x), I*tan(c + d*x))/(4*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_251():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(1)/3)
    F = 3*a*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -I*tan(c + d*x), I*tan(c + d*x))/(2*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_252():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(2)/3)
    F = 3*a*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -I*tan(c + d*x), I*tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_253():
    f = sqrt(I*a*tan(c + d*x) + a)/tan(c + d*x)**(sympy.S(4)/3)
    F = -3*a*sqrt(I*tan(c + d*x) + 1)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -I*tan(c + d*x), I*tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_254():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(4)/3)
    F = 3*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, sympy.S(-1)/2, 1, sympy.S(10)/3, -I*tan(c + d*x), I*tan(c + d*x))/(7*d*sqrt(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_255():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(2)/3)
    F = 3*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, sympy.S(-1)/2, 1, sympy.S(8)/3, -I*tan(c + d*x), I*tan(c + d*x))/(5*d*sqrt(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_256():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(1)/3)
    F = 3*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, sympy.S(-1)/2, 1, sympy.S(7)/3, -I*tan(c + d*x), I*tan(c + d*x))/(4*d*sqrt(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_257():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(1)/3)
    F = 3*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, sympy.S(-1)/2, 1, sympy.S(5)/3, -I*tan(c + d*x), I*tan(c + d*x))/(2*d*sqrt(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_258():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(2)/3)
    F = 3*a*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, sympy.S(-1)/2, 1, sympy.S(4)/3, -I*tan(c + d*x), I*tan(c + d*x))/(d*sqrt(I*tan(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_259():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(4)/3)
    F = -3*a*sqrt(I*a*tan(c + d*x) + a)*appellf1(sympy.S(-1)/3, sympy.S(-1)/2, 1, sympy.S(2)/3, -I*tan(c + d*x), I*tan(c + d*x))/(d*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_260():
    f = tan(c + d*x)**(sympy.S(4)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, 1, sympy.S(3)/2, sympy.S(10)/3, I*tan(c + d*x), -I*tan(c + d*x))/(7*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_261():
    f = tan(c + d*x)**(sympy.S(2)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, 1, sympy.S(3)/2, sympy.S(8)/3, I*tan(c + d*x), -I*tan(c + d*x))/(5*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_262():
    f = tan(c + d*x)**(sympy.S(1)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, 1, sympy.S(3)/2, sympy.S(7)/3, I*tan(c + d*x), -I*tan(c + d*x))/(4*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_263():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3))
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, 1, sympy.S(3)/2, sympy.S(5)/3, I*tan(c + d*x), -I*tan(c + d*x))/(2*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_264():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(2)/3))
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, 1, sympy.S(3)/2, sympy.S(4)/3, I*tan(c + d*x), -I*tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_265():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(4)/3))
    F = -3*sqrt(I*tan(c + d*x) + 1)*appellf1(sympy.S(-1)/3, 1, sympy.S(3)/2, sympy.S(2)/3, I*tan(c + d*x), -I*tan(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_266():
    f = tan(c + d*x)**(sympy.S(4)/3)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, 1, sympy.S(5)/2, sympy.S(10)/3, I*tan(c + d*x), -I*tan(c + d*x))/(7*a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_267():
    f = tan(c + d*x)**(sympy.S(2)/3)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, 1, sympy.S(5)/2, sympy.S(8)/3, I*tan(c + d*x), -I*tan(c + d*x))/(5*a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_268():
    f = tan(c + d*x)**(sympy.S(1)/3)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, 1, sympy.S(5)/2, sympy.S(7)/3, I*tan(c + d*x), -I*tan(c + d*x))/(4*a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_269():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(1)/3))
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, 1, sympy.S(5)/2, sympy.S(5)/3, I*tan(c + d*x), -I*tan(c + d*x))/(2*a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_270():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(2)/3))
    F = 3*sqrt(I*tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, 1, sympy.S(5)/2, sympy.S(4)/3, I*tan(c + d*x), -I*tan(c + d*x))/(a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_271():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(4)/3))
    F = -3*sqrt(I*tan(c + d*x) + 1)*appellf1(sympy.S(-1)/3, 1, sympy.S(5)/2, sympy.S(2)/3, I*tan(c + d*x), -I*tan(c + d*x))/(a*d*sqrt(I*a*tan(c + d*x) + a)*tan(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_272():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)**3
    F = -2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*x/4 - 3*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) + 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) + 3*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)**2/(7*d) - 18*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/(7*d) - 3*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)/(28*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_273():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)**2
    F = 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*x/4 - 3*2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) + 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) - 3*I*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_274():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)
    F = 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*x/4 + 3*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) + 3*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_275():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*x/4 + 3*2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) - 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_276():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)
    F = -2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*x/4 + 3*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 3*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) - a**(sympy.S(1)/3)*log(tan(c + d*x))/(2*d) + 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) - sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_277():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)**2
    F = 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*x/4 + I*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 3*2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) - 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) - I*a**(sympy.S(1)/3)*log(tan(c + d*x))/(6*d) + 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) - sqrt(3)*I*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*d) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_278():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)**3
    F = 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*x/4 - 4*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*d) + 3*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d) + 4*a**(sympy.S(1)/3)*log(tan(c + d*x))/(9*d) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d) + 8*sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*d) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)**2/(2*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_279():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(2)/3)
    F = -2**(sympy.S(2)/3)*a**(sympy.S(2)/3)*x/4 + 3*2**(sympy.S(2)/3)*I*a**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(4*d) + 2**(sympy.S(2)/3)*I*a**(sympy.S(2)/3)*log(cos(c + d*x))/(4*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*a**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_280():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*tan(c + d*x)**3
    F = -2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*x/2 - 3*2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) + 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d - 3*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d + 3*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*tan(c + d*x)**2/(10*d) - 9*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)/(20*d) - 6*(I*a*tan(c + d*x) + a)**(sympy.S(7)/3)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_281():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*tan(c + d*x)**2
    F = 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*x/2 - 3*2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) + 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d - 3*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d - 3*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/3)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_282():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*tan(c + d*x)
    F = 2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*x/2 + 3*2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) + 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d + 3*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d + 3*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_283():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = -2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*x/2 + 3*2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) + 2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) - 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d + 3*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_284():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*cot(c + d*x)
    F = -2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*x/2 + 3*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 3*2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) - a**(sympy.S(4)/3)*log(tan(c + d*x))/(2*d) + 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d - sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_285():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*cot(c + d*x)**2
    F = 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*x/2 + 2*I*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/d - 3*2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) - 2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) - 2*I*a**(sympy.S(4)/3)*log(tan(c + d*x))/(3*d) + 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d - 4*sqrt(3)*I*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*d) + I*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)/d - (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_286():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*cot(c + d*x)**3
    F = 2**(sympy.S(1)/3)*I*a**(sympy.S(4)/3)*x/2 - 11*a**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(6*d) + 3*2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) + 2**(sympy.S(1)/3)*a**(sympy.S(4)/3)*log(cos(c + d*x))/(2*d) + 11*a**(sympy.S(4)/3)*log(tan(c + d*x))/(18*d) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d + 11*sqrt(3)*a**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*d) - 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)*cot(c + d*x)/(3*d) - (I*a*tan(c + d*x) + a)**(sympy.S(4)/3)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_287():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/3)
    F = -2**(sympy.S(2)/3)*a**(sympy.S(5)/3)*x/2 + 3*2**(sympy.S(2)/3)*I*a**(sympy.S(5)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*d) + 2**(sympy.S(2)/3)*I*a**(sympy.S(5)/3)*log(cos(c + d*x))/(2*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*a**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/d + 3*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_288():
    f = tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = (I*tan(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(4)/3, m + 2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(m + 1)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_289():
    f = sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = 2*(I*tan(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, sympy.S(4)/3, sympy.S(5)/2, I*tan(c + d*x), -I*tan(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_290():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*tan(c + d*x)**3/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 15*I*tan(c + d*x)**2/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 45*I*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(8*a*d) - 39*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/3)/(20*a**2*d) - 2**(sympy.S(2)/3)*x/(8*a**(sympy.S(1)/3)) + 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_291():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*tan(c + d*x)**2/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 21/(10*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 3*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(10*a*d) - 2**(sympy.S(2)/3)*I*x/(8*a**(sympy.S(1)/3)) - 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_292():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = -3*I/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 3*I*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(2*a*d) + 2**(sympy.S(2)/3)*x/(8*a**(sympy.S(1)/3)) - 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_293():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = -3/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*I*x/(8*a**(sympy.S(1)/3)) + 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_294():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-1)/3)
    F = 3*I/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*x/(8*a**(sympy.S(1)/3)) + 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_295():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = 3/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*I*x/(8*a**(sympy.S(1)/3)) + 3*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)*d) - 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) - log(tan(c + d*x))/(2*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_296():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)
    F = -cot(c + d*x)/(d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 5*I/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*x/(8*a**(sympy.S(1)/3)) - I*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)*d) - 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(8*a**(sympy.S(1)/3)*d) + I*log(tan(c + d*x))/(6*a**(sympy.S(1)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(1)/3)*d) - sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_297():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-2)/3)
    F = 3*I/(4*d*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*x/(8*a**(sympy.S(2)/3)) + 3*2**(sympy.S(1)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(8*a**(sympy.S(2)/3)*d) + 2**(sympy.S(1)/3)*I*log(cos(c + d*x))/(8*a**(sympy.S(2)/3)*d) - 2**(sympy.S(1)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(4*a**(sympy.S(2)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_298():
    f = tan(c + d*x)**m/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = (I*tan(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(7)/3, m + 2, I*tan(c + d*x), -I*tan(c + d*x))/(a*d*(m + 1)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_299():
    f = sqrt(tan(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = 2*(I*tan(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, sympy.S(7)/3, sympy.S(5)/2, I*tan(c + d*x), -I*tan(c + d*x))/(3*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_300():
    f = tan(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = 3*tan(c + d*x)**3/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) - 39*I*tan(c + d*x)**2/(40*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) - 51*I/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 87*I*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)/(40*a**2*d) - 2**(sympy.S(2)/3)*x/(16*a**(sympy.S(4)/3)) + 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_301():
    f = tan(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = 3*tan(c + d*x)**2/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) + 15/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) - 27/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*I*x/(16*a**(sympy.S(4)/3)) - 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_302():
    f = tan(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = -3*I/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) + 9*I/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*x/(16*a**(sympy.S(4)/3)) - 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_303():
    f = tan(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = -3/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) + 3/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*I*x/(16*a**(sympy.S(4)/3)) + 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_304():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-4)/3)
    F = 3*I/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) + 3*I/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*x/(16*a**(sympy.S(4)/3)) + 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) + 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_305():
    f = cot(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = 3/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) + 9/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*I*x/(16*a**(sympy.S(4)/3)) + 3*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)*d) - 3*2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) - log(tan(c + d*x))/(2*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d) + sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_306():
    f = cot(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)
    F = -cot(c + d*x)/(d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) - 11*I/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(4)/3)) - 19*I/(4*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(2)/3)*x/(16*a**(sympy.S(4)/3)) - 2*I*log(a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(a**(sympy.S(4)/3)*d) - 3*2**(sympy.S(2)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*I*log(cos(c + d*x))/(16*a**(sympy.S(4)/3)*d) + 2*I*log(tan(c + d*x))/(3*a**(sympy.S(4)/3)*d) - 2**(sympy.S(2)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(4)/3)*d) - 4*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_307():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(-5)/3)
    F = 3*I/(10*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/3)) + 3*I/(8*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*x/(16*a**(sympy.S(5)/3)) + 3*2**(sympy.S(1)/3)*I*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(16*a**(sympy.S(5)/3)*d) + 2**(sympy.S(1)/3)*I*log(cos(c + d*x))/(16*a**(sympy.S(5)/3)*d) - 2**(sympy.S(1)/3)*sqrt(3)*I*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(I*a*tan(c + d*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(8*a**(sympy.S(5)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_308():
    f = (e*tan(c + d*x))**m*(I*a*tan(c + d*x) + a)
    F = a*(e*tan(c + d*x))**(m + 1)*hyper((1, m + 1), (m + 2,), I*tan(c + d*x))/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_309():
    f = (e*tan(c + d*x))**m*(-I*a*tan(c + d*x) + a)
    F = a*(e*tan(c + d*x))**(m + 1)*hyper((1, m + 1), (m + 2,), -I*tan(c + d*x))/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_310():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**4
    F = 8*a**4*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), I*tan(e + f*x))/(d*f*(n + 1)) - 2*a**4*(d*tan(e + f*x))**(n + 1)*(2*n**2 + 11*n + 16)/(d*f*(n + 1)*(n + 2)*(n + 3)) - (d*tan(e + f*x))**(n + 1)*(I*a**2*tan(e + f*x) + a**2)**2/(d*f*(n + 3)) - (d*tan(e + f*x))**(n + 1)*(2*n + 8)*(I*a**4*tan(e + f*x) + a**4)/(d*f*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_311():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**3
    F = 4*a**3*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), I*tan(e + f*x))/(d*f*(n + 1)) - a**3*(d*tan(e + f*x))**(n + 1)*(2*n + 5)/(d*f*(n + 1)*(n + 2)) - (d*tan(e + f*x))**(n + 1)*(I*a**3*tan(e + f*x) + a**3)/(d*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_312():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**2
    F = 2*a**2*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), I*tan(e + f*x))/(d*f*(n + 1)) - a**2*(d*tan(e + f*x))**(n + 1)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_313():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)
    F = a*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), I*tan(e + f*x))/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_314():
    f = (d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)
    F = (d*tan(e + f*x))**(n + 1)/(2*d*f*(I*a*tan(e + f*x) + a)) + (d*tan(e + f*x))**(n + 1)*(1 - n)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(2*a*d*f*(n + 1)) + I*n*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(2*a*d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_315():
    f = (d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**2
    F = (d*tan(e + f*x))**(n + 1)/(4*d*f*(I*a*tan(e + f*x) + a)**2) + (d*tan(e + f*x))**(n + 1)*(1 - n)**2*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(4*a**2*d*f*(n + 1)) + (d*tan(e + f*x))**(n + 1)*(2 - n)/(4*a**2*d*f*(I*tan(e + f*x) + 1)) + I*n*(d*tan(e + f*x))**(n + 2)*(2 - n)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(4*a**2*d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_316():
    f = (d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**3
    F = (d*tan(e + f*x))**(n + 1)*(2 - n)*(5 - 2*n)/(24*d*f*(I*a**3*tan(e + f*x) + a**3)) + (d*tan(e + f*x))**(n + 1)/(6*d*f*(I*a*tan(e + f*x) + a)**3) + (d*tan(e + f*x))**(n + 1)*(7 - 2*n)/(24*a*d*f*(I*a*tan(e + f*x) + a)**2) + (d*tan(e + f*x))**(n + 1)*(1 - 2*n)*(1 - n)*(3 - n)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(24*a**3*d*f*(n + 1)) + I*n*(d*tan(e + f*x))**(n + 2)*(2 - n)*(5 - 2*n)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(24*a**3*d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_317():
    f = (d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**4
    F = (d*tan(e + f*x))**(n + 1)/(8*d*f*(I*a*tan(e + f*x) + a)**4) + (d*tan(e + f*x))**(n + 1)*(5 - n)/(24*a*d*f*(I*a*tan(e + f*x) + a)**3) + (d*tan(e + f*x))**(n + 1)*(1 - n)*(3 - n)*(n**2 - 4*n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(48*a**4*d*f*(n + 1)) + (d*tan(e + f*x))**(n + 1)*(2 - n)**2*(4 - n)/(48*a**4*d*f*(I*tan(e + f*x) + 1)) + (d*tan(e + f*x))**(n + 1)*(n**2 - 7*n + 13)/(48*a**4*d*f*(I*tan(e + f*x) + 1)**2) + I*n*(d*tan(e + f*x))**(n + 2)*(2 - n)**2*(4 - n)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(48*a**4*d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_318():
    f = (d*tan(e + f*x))**n*(-I*a*tan(e + f*x) + a)
    F = a*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), -I*tan(e + f*x))/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_319():
    f = (d*tan(e + f*x))**n/(-I*a*tan(e + f*x) + a)
    F = (d*tan(e + f*x))**(n + 1)/(2*d*f*(-I*a*tan(e + f*x) + a)) + (d*tan(e + f*x))**(n + 1)*(1 - n)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(2*a*d*f*(n + 1)) - I*n*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(2*a*d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_320():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = a*(d*tan(e + f*x))**(n + 1)*sqrt(I*a*tan(e + f*x) + a)*appellf1(n + 1, sympy.S(-1)/2, 1, n + 2, -I*tan(e + f*x), I*tan(e + f*x))/(d*f*(n + 1)*sqrt(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_321():
    f = (d*tan(e + f*x))**n*sqrt(I*a*tan(e + f*x) + a)
    F = a*(d*tan(e + f*x))**(n + 1)*sqrt(I*tan(e + f*x) + 1)*appellf1(n + 1, sympy.S.Half, 1, n + 2, -I*tan(e + f*x), I*tan(e + f*x))/(d*f*(n + 1)*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_322():
    f = (d*tan(e + f*x))**n/sqrt(I*a*tan(e + f*x) + a)
    F = (d*tan(e + f*x))**(n + 1)*sqrt(I*tan(e + f*x) + 1)*appellf1(n + 1, 1, sympy.S(3)/2, n + 2, I*tan(e + f*x), -I*tan(e + f*x))/(d*f*(n + 1)*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_323():
    f = (d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (d*tan(e + f*x))**(n + 1)*sqrt(I*tan(e + f*x) + 1)*appellf1(n + 1, 1, sympy.S(5)/2, n + 2, I*tan(e + f*x), -I*tan(e + f*x))/(a*d*f*(n + 1)*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_324():
    f = (d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**m
    F = (d*tan(e + f*x))**(n + 1)*(I*a*tan(e + f*x) + a)**m*appellf1(n + 1, 1, 1 - m, n + 2, I*tan(e + f*x), -I*tan(e + f*x))/(d*f*(n + 1)*(I*tan(e + f*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_325():
    f = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**4
    F = -I*m*(I*a*tan(c + d*x) + a)**m*tan(c + d*x)**2/(d*(m**2 + 5*m + 6)) + 2*I*(I*a*tan(c + d*x) + a)**m/(d*(m**2 + 5*m + 6)) + (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**3/(d*(m + 3)) - I*(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m) + I*(I*a*tan(c + d*x) + a)**(m + 1)*(m**2 + 3*m + 6)/(a*d*(m + 1)*(m + 2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_326():
    f = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**3
    F = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**2/(d*(m + 2)) + (I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m) - 2*(I*a*tan(c + d*x) + a)**m/(d*m*(m + 2)) - m*(I*a*tan(c + d*x) + a)**(m + 1)/(a*d*(m**2 + 3*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_327():
    f = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**2
    F = I*(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m) - I*(I*a*tan(c + d*x) + a)**(m + 1)/(a*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_328():
    f = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)
    F = -(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m) + (I*a*tan(c + d*x) + a)**m/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_329():
    f = (I*a*tan(c + d*x) + a)**m
    F = -I*(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_330():
    f = (I*a*tan(c + d*x) + a)**m*cot(c + d*x)
    F = (I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m) - (I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x) + 1)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_331():
    f = (I*a*tan(c + d*x) + a)**m*cot(c + d*x)**2
    F = -(I*a*tan(c + d*x) + a)**m*cot(c + d*x)/d - I*(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x) + 1)/d + I*(I*a*tan(c + d*x) + a)**m*hyper((1, m), (m + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_332():
    f = (I*a*tan(c + d*x) + a)**m*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*(I*a*tan(c + d*x) + a)**m*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, 1, 1 - m, sympy.S(7)/2, I*tan(c + d*x), -I*tan(c + d*x))/(5*d*(I*tan(c + d*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_333():
    f = (I*a*tan(c + d*x) + a)**m*sqrt(tan(c + d*x))
    F = 2*(I*a*tan(c + d*x) + a)**m*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, 1 - m, sympy.S(5)/2, I*tan(c + d*x), -I*tan(c + d*x))/(3*d*(I*tan(c + d*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_334():
    f = (I*a*tan(c + d*x) + a)**m/sqrt(tan(c + d*x))
    F = 2*(I*a*tan(c + d*x) + a)**m*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, 1 - m, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_335():
    f = (I*a*tan(c + d*x) + a)**m/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*(I*a*tan(c + d*x) + a)**m*appellf1(sympy.S(-1)/2, 1, 1 - m, sympy.S.Half, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**m*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_336():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a)
    F = sqrt(2)*a*d**(sympy.S(5)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f - 2*a*d**2*sqrt(d*tan(e + f*x))/f + 2*a*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*a*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_337():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)
    F = sqrt(2)*a*d**(sympy.S(3)/2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f + 2*a*d*sqrt(d*tan(e + f*x))/f + 2*a*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_338():
    f = sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)
    F = -sqrt(2)*a*sqrt(d)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f + 2*a*sqrt(d*tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_339():
    f = (a*tan(e + f*x) + a)/sqrt(d*tan(e + f*x))
    F = -sqrt(2)*a*atan(sqrt(2)*sqrt(d)*(1 - tan(e + f*x))/(2*sqrt(d*tan(e + f*x))))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_340():
    f = (a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a/(d*f*sqrt(d*tan(e + f*x))) + sqrt(2)*a*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_341():
    f = (a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) - 2*a/(d**2*f*sqrt(d*tan(e + f*x))) + sqrt(2)*a*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_342():
    f = (a*tan(e + f*x) + a)/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -2*a/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2)) - 2*a/(3*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 2*a/(d**3*f*sqrt(d*tan(e + f*x))) - sqrt(2)*a*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_343():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a)**2
    F = -sqrt(2)*a**2*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) + sqrt(2)*a**2*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) - sqrt(2)*a**2*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f + sqrt(2)*a**2*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f - 4*a**2*d**2*sqrt(d*tan(e + f*x))/f + 4*a**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 2*a**2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_344():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)**2
    F = -sqrt(2)*a**2*d**(sympy.S(3)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) + sqrt(2)*a**2*d**(sympy.S(3)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) + sqrt(2)*a**2*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f - sqrt(2)*a**2*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 4*a**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*a**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_345():
    f = sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)**2
    F = sqrt(2)*a**2*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) - sqrt(2)*a**2*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*f) + sqrt(2)*a**2*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f - sqrt(2)*a**2*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/f + 4*a**2*sqrt(d*tan(e + f*x))/f + 2*a**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_346():
    f = (a*tan(e + f*x) + a)**2/sqrt(d*tan(e + f*x))
    F = 2*a**2*sqrt(d*tan(e + f*x))/(d*f) + sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*sqrt(d)*f) - sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*sqrt(d)*f) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_347():
    f = (a*tan(e + f*x) + a)**2/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a**2/(d*f*sqrt(d*tan(e + f*x))) - sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*d**(sympy.S(3)/2)*f) + sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*d**(sympy.S(3)/2)*f) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_348():
    f = (a*tan(e + f*x) + a)**2/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a**2/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) - 4*a**2/(d**2*f*sqrt(d*tan(e + f*x))) - sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(2*d**(sympy.S(5)/2)*f) + sqrt(2)*a**2*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(2*d**(sympy.S(5)/2)*f) + sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f) - sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_349():
    f = (d*tan(e + f*x))**(sympy.S(7)/2)*(a*tan(e + f*x) + a)**3
    F = -2*sqrt(2)*a**3*d**(sympy.S(7)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f + 4*a**3*d**3*sqrt(d*tan(e + f*x))/f - 4*a**3*d**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) - 4*a**3*d*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 4*a**3*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*f) + 16*a**3*(d*tan(e + f*x))**(sympy.S(9)/2)/(33*d*f) + 2*(d*tan(e + f*x))**(sympy.S(9)/2)*(a**3*tan(e + f*x) + a**3)/(11*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_350():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a)**3
    F = -2*sqrt(2)*a**3*d**(sympy.S(5)/2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f - 4*a**3*d**2*sqrt(d*tan(e + f*x))/f - 4*a**3*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 4*a**3*(d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 40*a**3*(d*tan(e + f*x))**(sympy.S(7)/2)/(63*d*f) + 2*(d*tan(e + f*x))**(sympy.S(7)/2)*(a**3*tan(e + f*x) + a**3)/(9*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_351():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)**3
    F = 2*sqrt(2)*a**3*d**(sympy.S(3)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f - 4*a**3*d*sqrt(d*tan(e + f*x))/f + 4*a**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 32*a**3*(d*tan(e + f*x))**(sympy.S(5)/2)/(35*d*f) + 2*(d*tan(e + f*x))**(sympy.S(5)/2)*(a**3*tan(e + f*x) + a**3)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_352():
    f = sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)**3
    F = 2*sqrt(2)*a**3*sqrt(d)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/f + 4*a**3*sqrt(d*tan(e + f*x))/f + 8*a**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(5*d*f) + 2*(d*tan(e + f*x))**(sympy.S(3)/2)*(a**3*tan(e + f*x) + a**3)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_353():
    f = (a*tan(e + f*x) + a)**3/sqrt(d*tan(e + f*x))
    F = 16*a**3*sqrt(d*tan(e + f*x))/(3*d*f) - 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(sqrt(d)*f) + 2*sqrt(d*tan(e + f*x))*(a**3*tan(e + f*x) + a**3)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_354():
    f = (a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 4*a**3*sqrt(d*tan(e + f*x))/(d**2*f) - 2*sqrt(2)*a**3*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) - (2*a**3*tan(e + f*x) + 2*a**3)/(d*f*sqrt(d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_355():
    f = (a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(5)/2)
    F = -16*a**3/(3*d**2*f*sqrt(d*tan(e + f*x))) + 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f) - (2*a**3*tan(e + f*x) + 2*a**3)/(3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_356():
    f = (a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(7)/2)
    F = -8*a**3/(5*d**2*f*(d*tan(e + f*x))**(sympy.S(3)/2)) - 4*a**3/(d**3*f*sqrt(d*tan(e + f*x))) + 2*sqrt(2)*a**3*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(7)/2)*f) - (2*a**3*tan(e + f*x) + 2*a**3)/(5*d*f*(d*tan(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_357():
    f = (a*tan(e + f*x) + a)**3/(d*tan(e + f*x))**(sympy.S(9)/2)
    F = -32*a**3/(35*d**2*f*(d*tan(e + f*x))**(sympy.S(5)/2)) - 4*a**3/(3*d**3*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 4*a**3/(d**4*f*sqrt(d*tan(e + f*x))) - 2*sqrt(2)*a**3*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(d**(sympy.S(9)/2)*f) - (2*a**3*tan(e + f*x) + 2*a**3)/(7*d*f*(d*tan(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_358():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(a*tan(e + f*x) + a)
    F = -d**(sympy.S(5)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*f) + sqrt(2)*d**(sympy.S(5)/2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(2*a*f) + 2*d**2*sqrt(d*tan(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_359():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(a*tan(e + f*x) + a)
    F = d**(sympy.S(3)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*f) - sqrt(2)*d**(sympy.S(3)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_360():
    f = sqrt(d*tan(e + f*x))/(a*tan(e + f*x) + a)
    F = -sqrt(d)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*f) - sqrt(2)*sqrt(d)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_361():
    f = 1/(sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a))
    F = atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*sqrt(d)*f) + sqrt(2)*atanh(sqrt(2)*sqrt(d)*(tan(e + f*x) + 1)/(2*sqrt(d*tan(e + f*x))))/(2*a*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_362():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a))
    F = -2/(a*d*f*sqrt(d*tan(e + f*x))) - atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(2*a*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_363():
    f = 1/((d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a))
    F = -2/(3*a*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 2/(a*d**2*f*sqrt(d*tan(e + f*x))) + atan(sqrt(d*tan(e + f*x))/sqrt(d))/(a*d**(sympy.S(5)/2)*f) - sqrt(2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(2*a*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_364():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(a*tan(e + f*x) + a)**2
    F = -d**2*sqrt(d*tan(e + f*x))/(2*f*(a**2*tan(e + f*x) + a**2)) + sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) - sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) + 3*d**(sympy.S(5)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) + sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f) - sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_365():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(a*tan(e + f*x) + a)**2
    F = d*sqrt(d*tan(e + f*x))/(2*f*(a**2*tan(e + f*x) + a**2)) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) - d**(sympy.S(3)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_366():
    f = sqrt(d*tan(e + f*x))/(a*tan(e + f*x) + a)**2
    F = -sqrt(d*tan(e + f*x))/(2*f*(a**2*tan(e + f*x) + a**2)) - sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) + sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*f) - sqrt(d)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_367():
    f = 1/(sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)**2)
    F = sqrt(d*tan(e + f*x))/(2*d*f*(a**2*tan(e + f*x) + a**2)) - sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*sqrt(d)*f) + sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*sqrt(d)*f) + 3*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*sqrt(d)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*sqrt(d)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_368():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)**2)
    F = 1/(2*d*f*sqrt(d*tan(e + f*x))*(a**2*tan(e + f*x) + a**2)) - 5/(2*a**2*d*f*sqrt(d*tan(e + f*x))) + sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*d**(sympy.S(3)/2)*f) - sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*d**(sympy.S(3)/2)*f) - 5*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_369():
    f = 1/((d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a)**2)
    F = 1/(2*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(a**2*tan(e + f*x) + a**2)) - 7/(6*a**2*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 9/(2*a**2*d**2*f*sqrt(d*tan(e + f*x))) + sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*d**(sympy.S(5)/2)*f) - sqrt(2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(8*a**2*d**(sympy.S(5)/2)*f) + 7*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(2*a**2*d**(sympy.S(5)/2)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*d**(sympy.S(5)/2)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(4*a**2*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_370():
    f = (d*tan(e + f*x))**(sympy.S(9)/2)/(a*tan(e + f*x) + a)**3
    F = -d**2*(d*tan(e + f*x))**(sympy.S(5)/2)/(4*a*f*(a*tan(e + f*x) + a)**2) - 31*d**(sympy.S(9)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f) + sqrt(2)*d**(sympy.S(9)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*f) + 27*d**4*sqrt(d*tan(e + f*x))/(8*a**3*f) - 9*d**3*(d*tan(e + f*x))**(sympy.S(3)/2)/(8*a**3*f*(tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_371():
    f = (d*tan(e + f*x))**(sympy.S(7)/2)/(a*tan(e + f*x) + a)**3
    F = -d**2*(d*tan(e + f*x))**(sympy.S(3)/2)/(4*a*f*(a*tan(e + f*x) + a)**2) + 11*d**(sympy.S(7)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f) + sqrt(2)*d**(sympy.S(7)/2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*f) - 7*d**3*sqrt(d*tan(e + f*x))/(8*a**3*f*(tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_372():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)/(a*tan(e + f*x) + a)**3
    F = -d**2*sqrt(d*tan(e + f*x))/(4*a*f*(a*tan(e + f*x) + a)**2) + d**(sympy.S(5)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f) - sqrt(2)*d**(sympy.S(5)/2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*f) + 5*d**2*sqrt(d*tan(e + f*x))/(8*a**3*f*(tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_373():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(a*tan(e + f*x) + a)**3
    F = -d*sqrt(d*tan(e + f*x))/(8*f*(a**3*tan(e + f*x) + a**3)) + d*sqrt(d*tan(e + f*x))/(4*a*f*(a*tan(e + f*x) + a)**2) - 5*d**(sympy.S(3)/2)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_374():
    f = sqrt(d*tan(e + f*x))/(a*tan(e + f*x) + a)**3
    F = -3*sqrt(d*tan(e + f*x))/(8*f*(a**3*tan(e + f*x) + a**3)) - sqrt(d*tan(e + f*x))/(4*a*f*(a*tan(e + f*x) + a)**2) + sqrt(d)*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*f) + sqrt(2)*sqrt(d)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_375():
    f = 1/(sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)**3)
    F = sqrt(d*tan(e + f*x))/(4*a*d*f*(a*tan(e + f*x) + a)**2) + 7*sqrt(d*tan(e + f*x))/(8*a**3*d*f*(tan(e + f*x) + 1)) + 11*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*sqrt(d)*f) + sqrt(2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_376():
    f = 1/((d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)**3)
    F = 1/(4*a*d*f*sqrt(d*tan(e + f*x))*(a*tan(e + f*x) + a)**2) - 27/(8*a**3*d*f*sqrt(d*tan(e + f*x))) + 9/(8*a**3*d*f*sqrt(d*tan(e + f*x))*(tan(e + f*x) + 1)) - 31*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*d**(sympy.S(3)/2)*f) - sqrt(2)*atanh(sqrt(2)*(sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_377():
    f = 1/((d*tan(e + f*x))**(sympy.S(5)/2)*(a*tan(e + f*x) + a)**3)
    F = 1/(4*a*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(a*tan(e + f*x) + a)**2) - 55/(24*a**3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)) + 11/(8*a**3*d*f*(d*tan(e + f*x))**(sympy.S(3)/2)*(tan(e + f*x) + 1)) + 63/(8*a**3*d**2*f*sqrt(d*tan(e + f*x))) + 59*atan(sqrt(d*tan(e + f*x))/sqrt(d))/(8*a**3*d**(sympy.S(5)/2)*f) - sqrt(2)*atan(sqrt(2)*(-sqrt(d)*tan(e + f*x) + sqrt(d))/(2*sqrt(d*tan(e + f*x))))/(4*a**3*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_378():
    f = sqrt(tan(e + f*x) + 1)*tan(e + f*x)**5
    F = 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**3/(9*f) - 4*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**2/(21*f) - 26*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)/(105*f) + 52*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(315*f) + 2*sqrt(tan(e + f*x) + 1)/f - sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_379():
    f = sqrt(tan(e + f*x) + 1)*tan(e + f*x)**3
    F = 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)/(5*f) - 4*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(15*f) - 2*sqrt(tan(e + f*x) + 1)/f + sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_380():
    f = sqrt(tan(e + f*x) + 1)*tan(e + f*x)
    F = 2*sqrt(tan(e + f*x) + 1)/f - sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_381():
    f = sqrt(tan(e + f*x) + 1)*cot(e + f*x)
    F = sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - 2*atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_382():
    f = sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(2*f) - sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(4*f) - sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + 9*atanh(sqrt(tan(e + f*x) + 1))/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_383():
    f = sqrt(tan(e + f*x) + 1)*cot(e + f*x)**5
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**4/(4*f) - sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(24*f) + 53*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(96*f) + 11*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(64*f) + sqrt(sympy.S(-1)/2 + sqrt(2)/2)*atan(((2 - sqrt(2))*tan(e + f*x) - 3*sqrt(2) + 4)/(2*sqrt(-7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atanh(((sqrt(2) + 2)*tan(e + f*x) + 4 + 3*sqrt(2))/(2*sqrt(7 + 5*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - 139*atanh(sqrt(tan(e + f*x) + 1))/(64*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_384():
    f = sqrt(tan(e + f*x) + 1)*tan(e + f*x)**4
    F = 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**2/(7*f) - 8*(tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)/(35*f) - 18*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(35*f) + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_385():
    f = sqrt(tan(e + f*x) + 1)*tan(e + f*x)**2
    F = 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(3*f) - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) + sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_386():
    f = sqrt(tan(e + f*x) + 1)
    F = log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_387():
    f = sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)/f - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) + sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_388():
    f = sqrt(tan(e + f*x) + 1)*cot(e + f*x)**4
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(3*f) - sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(12*f) + 9*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(8*f) + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + 7*atanh(sqrt(tan(e + f*x) + 1))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_389():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**5
    F = 2*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)**3/(11*f) - 4*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)**2/(33*f) - 50*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)/(231*f) + 20*(tan(e + f*x) + 1)**(sympy.S(5)/2)/(231*f) + 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(3*f) + 2*sqrt(tan(e + f*x) + 1)/f + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_390():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**3
    F = 2*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)/(7*f) - 4*(tan(e + f*x) + 1)**(sympy.S(5)/2)/(35*f) - 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(3*f) - 2*sqrt(tan(e + f*x) + 1)/f - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_391():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)
    F = 2*(tan(e + f*x) + 1)**(sympy.S(3)/2)/(3*f) + 2*sqrt(tan(e + f*x) + 1)/f + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_392():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*cot(e + f*x)
    F = -log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - 2*atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_393():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*cot(e + f*x)**3
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(2*f) - 5*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(4*f) + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + 5*atanh(sqrt(tan(e + f*x) + 1))/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_394():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*cot(e + f*x)**5
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**4/(4*f) - 3*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(8*f) + 15*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(32*f) + 83*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(64*f) - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(2*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/f - 83*atanh(sqrt(tan(e + f*x) + 1))/(64*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_395():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**4
    F = 2*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)**2/(9*f) - 8*(tan(e + f*x) + 1)**(sympy.S(5)/2)*tan(e + f*x)/(63*f) - 22*(tan(e + f*x) + 1)**(sympy.S(5)/2)/(63*f) + 2*sqrt(tan(e + f*x) + 1)/f - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_396():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*tan(e + f*x)**2
    F = 2*(tan(e + f*x) + 1)**(sympy.S(5)/2)/(5*f) - 2*sqrt(tan(e + f*x) + 1)/f + sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_397():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)
    F = 2*sqrt(tan(e + f*x) + 1)/f - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_398():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*cot(e + f*x)**2
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)/f + sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - 3*atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_399():
    f = (tan(e + f*x) + 1)**(sympy.S(3)/2)*cot(e + f*x)**4
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(3*f) - 7*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(12*f) + 7*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(8*f) - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/f + 25*atanh(sqrt(tan(e + f*x) + 1))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_400():
    f = tan(e + f*x)**5/sqrt(tan(e + f*x) + 1)
    F = 2*sqrt(tan(e + f*x) + 1)*tan(e + f*x)**3/(7*f) - 12*sqrt(tan(e + f*x) + 1)*tan(e + f*x)**2/(35*f) - 22*sqrt(tan(e + f*x) + 1)*tan(e + f*x)/(105*f) + 44*sqrt(tan(e + f*x) + 1)/(105*f) - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_401():
    f = tan(e + f*x)**3/sqrt(tan(e + f*x) + 1)
    F = 2*sqrt(tan(e + f*x) + 1)*tan(e + f*x)/(3*f) - 4*sqrt(tan(e + f*x) + 1)/(3*f) + sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_402():
    f = tan(e + f*x)/sqrt(tan(e + f*x) + 1)
    F = -sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_403():
    f = cot(e + f*x)/sqrt(tan(e + f*x) + 1)
    F = sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) - 2*atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_404():
    f = cot(e + f*x)**3/sqrt(tan(e + f*x) + 1)
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(2*f) + 3*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(4*f) - sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) - sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) + 5*atanh(sqrt(tan(e + f*x) + 1))/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_405():
    f = cot(e + f*x)**5/sqrt(tan(e + f*x) + 1)
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**4/(4*f) + 7*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(24*f) + 13*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(96*f) - 13*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(64*f) + sqrt(-1 + sqrt(2))*atan(((1 - sqrt(2))*tan(e + f*x) - 2*sqrt(2) + 3)/(sqrt(-14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) + sqrt(1 + sqrt(2))*atanh(((1 + sqrt(2))*tan(e + f*x) + 2*sqrt(2) + 3)/(sqrt(14 + 10*sqrt(2))*sqrt(tan(e + f*x) + 1)))/(2*f) - 115*atanh(sqrt(tan(e + f*x) + 1))/(64*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_406():
    f = tan(e + f*x)**4/sqrt(tan(e + f*x) + 1)
    F = 2*sqrt(tan(e + f*x) + 1)*tan(e + f*x)**2/(5*f) - 8*sqrt(tan(e + f*x) + 1)*tan(e + f*x)/(15*f) - 14*sqrt(tan(e + f*x) + 1)/(15*f) - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_407():
    f = tan(e + f*x)**2/sqrt(tan(e + f*x) + 1)
    F = 2*sqrt(tan(e + f*x) + 1)/f + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) - sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_408():
    f = 1/sqrt(tan(e + f*x) + 1)
    F = -log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_409():
    f = cot(e + f*x)**2/sqrt(tan(e + f*x) + 1)
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)/f + log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) - log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) + sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) - sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) + atanh(sqrt(tan(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_410():
    f = cot(e + f*x)**4/sqrt(tan(e + f*x) + 1)
    F = -sqrt(tan(e + f*x) + 1)*cot(e + f*x)**3/(3*f) + 5*sqrt(tan(e + f*x) + 1)*cot(e + f*x)**2/(12*f) + 3*sqrt(tan(e + f*x) + 1)*cot(e + f*x)/(8*f) - log(-sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) + log(sqrt(2 + 2*sqrt(2))*sqrt(tan(e + f*x) + 1) + tan(e + f*x) + 1 + sqrt(2))/(4*f*sqrt(1 + sqrt(2))) - sqrt(1 + sqrt(2))*atan((-2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) + sqrt(1 + sqrt(2))*atan((2*sqrt(tan(e + f*x) + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/(2*f) - 3*atanh(sqrt(tan(e + f*x) + 1))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_411():
    f = (d*tan(e + f*x))**n*(a*tan(e + f*x) + a)**m
    F = (d*tan(e + f*x))**(n + 1)*(a*tan(e + f*x) + a)**m*appellf1(n + 1, 1, -m, n + 2, -I*tan(e + f*x), -tan(e + f*x))/(2*d*f*(n + 1)*(tan(e + f*x) + 1)**m) + (d*tan(e + f*x))**(n + 1)*(a*tan(e + f*x) + a)**m*appellf1(n + 1, 1, -m, n + 2, I*tan(e + f*x), -tan(e + f*x))/(2*d*f*(n + 1)*(tan(e + f*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_412():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**5
    F = -a*log(cos(c + d*x))/d + a*tan(c + d*x)**4/(4*d) - a*tan(c + d*x)**2/(2*d) - b*x + b*tan(c + d*x)**5/(5*d) - b*tan(c + d*x)**3/(3*d) + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_413():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**4
    F = a*x + a*tan(c + d*x)**3/(3*d) - a*tan(c + d*x)/d - b*log(cos(c + d*x))/d + b*tan(c + d*x)**4/(4*d) - b*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_414():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**3
    F = a*log(cos(c + d*x))/d + a*tan(c + d*x)**2/(2*d) + b*x + b*tan(c + d*x)**3/(3*d) - b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_415():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**2
    F = -a*x + a*tan(c + d*x)/d + b*log(cos(c + d*x))/d + b*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_416():
    f = (a + b*tan(c + d*x))*tan(c + d*x)
    F = -a*log(cos(c + d*x))/d - b*x + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_417():
    f = a + b*tan(c + d*x)
    F = a*x - b*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_418():
    f = (a + b*tan(c + d*x))*cot(c + d*x)
    F = a*log(sin(c + d*x))/d + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_419():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**2
    F = -a*x - a*cot(c + d*x)/d + b*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_420():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**3
    F = -a*log(sin(c + d*x))/d - a*cot(c + d*x)**2/(2*d) - b*x - b*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_421():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**4
    F = a*x - a*cot(c + d*x)**3/(3*d) + a*cot(c + d*x)/d - b*log(sin(c + d*x))/d - b*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_422():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**5
    F = a*log(sin(c + d*x))/d - a*cot(c + d*x)**4/(4*d) + a*cot(c + d*x)**2/(2*d) + b*x - b*cot(c + d*x)**3/(3*d) + b*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_423():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**6
    F = -a*x - a*cot(c + d*x)**5/(5*d) + a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + b*log(sin(c + d*x))/d - b*cot(c + d*x)**4/(4*d) + b*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_424():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)**4
    F = -2*a*b*log(cos(c + d*x))/d + a*b*tan(c + d*x)**4/(2*d) - a*b*tan(c + d*x)**2/d + b**2*tan(c + d*x)**5/(5*d) + x*(a**2 - b**2) + (a**2 - b**2)*tan(c + d*x)**3/(3*d) - (a**2 - b**2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_425():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)**3
    F = 2*a*b*x + 2*a*b*tan(c + d*x)**3/(3*d) - 2*a*b*tan(c + d*x)/d + b**2*tan(c + d*x)**4/(4*d) + (a**2 - b**2)*log(cos(c + d*x))/d + (a**2 - b**2)*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_426():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)**2
    F = 2*a*b*log(cos(c + d*x))/d - b**2*tan(c + d*x)/d - x*(a**2 - b**2) + (a + b*tan(c + d*x))**3/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_427():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)
    F = -2*a*b*x + a*b*tan(c + d*x)/d + (a + b*tan(c + d*x))**2/(2*d) - (a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_428():
    f = (a + b*tan(c + d*x))**2
    F = -2*a*b*log(cos(c + d*x))/d + b**2*tan(c + d*x)/d + x*(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_429():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)
    F = a**2*log(sin(c + d*x))/d + 2*a*b*x - b**2*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_430():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**2
    F = -a**2*cot(c + d*x)/d + 2*a*b*log(sin(c + d*x))/d - x*(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_431():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**3
    F = -a**2*cot(c + d*x)**2/(2*d) - 2*a*b*x - 2*a*b*cot(c + d*x)/d - (a**2 - b**2)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_432():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**4
    F = -a**2*cot(c + d*x)**3/(3*d) - 2*a*b*log(sin(c + d*x))/d - a*b*cot(c + d*x)**2/d + x*(a**2 - b**2) + (a**2 - b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_433():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**5
    F = -a**2*cot(c + d*x)**4/(4*d) + 2*a*b*x - 2*a*b*cot(c + d*x)**3/(3*d) + 2*a*b*cot(c + d*x)/d + (a**2 - b**2)*log(sin(c + d*x))/d + (a**2 - b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_434():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**6
    F = -a**2*cot(c + d*x)**5/(5*d) + 2*a*b*log(sin(c + d*x))/d - a*b*cot(c + d*x)**4/(2*d) + a*b*cot(c + d*x)**2/d - x*(a**2 - b**2) + (a**2 - b**2)*cot(c + d*x)**3/(3*d) - (a**2 - b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_435():
    f = (a + b*tan(c + d*x))**3*tan(c + d*x)**3
    F = -a*(a + b*tan(c + d*x))**2/(2*d) + a*(a**2 - 3*b**2)*log(cos(c + d*x))/d - a*(a + b*tan(c + d*x))**4/(20*b**2*d) + b*x*(3*a**2 - b**2) - b*(a**2 - b**2)*tan(c + d*x)/d - (a + b*tan(c + d*x))**3/(3*d) + (a + b*tan(c + d*x))**4*tan(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_436():
    f = (a + b*tan(c + d*x))**3*tan(c + d*x)**2
    F = -2*a*b**2*tan(c + d*x)/d - a*x*(a**2 - 3*b**2) - b*(a + b*tan(c + d*x))**2/(2*d) + b*(3*a**2 - b**2)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**4/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_437():
    f = (a + b*tan(c + d*x))**3*tan(c + d*x)
    F = a*(a + b*tan(c + d*x))**2/(2*d) - a*(a**2 - 3*b**2)*log(cos(c + d*x))/d - b*x*(3*a**2 - b**2) + b*(a**2 - b**2)*tan(c + d*x)/d + (a + b*tan(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_438():
    f = (a + b*tan(c + d*x))**3
    F = 2*a*b**2*tan(c + d*x)/d + a*x*(a**2 - 3*b**2) + b*(a + b*tan(c + d*x))**2/(2*d) - b*(3*a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_439():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)
    F = a**3*log(sin(c + d*x))/d - 3*a*b**2*log(cos(c + d*x))/d + b**2*(a + b*tan(c + d*x))/d + b*x*(3*a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_440():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**2
    F = 3*a**2*b*log(sin(c + d*x))/d - a**2*(a + b*tan(c + d*x))*cot(c + d*x)/d - a*x*(a**2 - 3*b**2) - b**3*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_441():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**3
    F = -5*a**2*b*cot(c + d*x)/(2*d) - a**2*(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - a*(a**2 - 3*b**2)*log(sin(c + d*x))/d - b*x*(3*a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_442():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**4
    F = -7*a**2*b*cot(c + d*x)**2/(6*d) - a**2*(a + b*tan(c + d*x))*cot(c + d*x)**3/(3*d) + a*x*(a**2 - 3*b**2) + a*(a**2 - 3*b**2)*cot(c + d*x)/d - b*(3*a**2 - b**2)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_443():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**5
    F = -3*a**2*b*cot(c + d*x)**3/(4*d) - a**2*(a + b*tan(c + d*x))*cot(c + d*x)**4/(4*d) + a*(a**2 - 3*b**2)*log(sin(c + d*x))/d + a*(a**2 - 3*b**2)*cot(c + d*x)**2/(2*d) + b*x*(3*a**2 - b**2) + b*(3*a**2 - b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_444():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**6
    F = -11*a**2*b*cot(c + d*x)**4/(20*d) - a**2*(a + b*tan(c + d*x))*cot(c + d*x)**5/(5*d) - a*x*(a**2 - 3*b**2) + a*(a**2 - 3*b**2)*cot(c + d*x)**3/(3*d) - a*(a**2 - 3*b**2)*cot(c + d*x)/d + b*(3*a**2 - b**2)*log(sin(c + d*x))/d + b*(3*a**2 - b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_445():
    f = (a + b*tan(c + d*x))**4*tan(c + d*x)**3
    F = 4*a*b*x*(a**2 - b**2) - a*b*(a**2 - 3*b**2)*tan(c + d*x)/d - a*(a + b*tan(c + d*x))**3/(3*d) - a*(a + b*tan(c + d*x))**5/(30*b**2*d) - (a + b*tan(c + d*x))**4/(4*d) - (a + b*tan(c + d*x))**2*(a**2 - b**2)/(2*d) + (a**4 - 6*a**2*b**2 + b**4)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**5*tan(c + d*x)/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_446():
    f = (a + b*tan(c + d*x))**4*tan(c + d*x)**2
    F = -a*b*(a + b*tan(c + d*x))**2/d + 4*a*b*(a**2 - b**2)*log(cos(c + d*x))/d - b**2*(3*a**2 - b**2)*tan(c + d*x)/d - b*(a + b*tan(c + d*x))**3/(3*d) - x*(a**4 - 6*a**2*b**2 + b**4) + (a + b*tan(c + d*x))**5/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_447():
    f = (a + b*tan(c + d*x))**4*tan(c + d*x)
    F = -4*a*b*x*(a**2 - b**2) + a*b*(a**2 - 3*b**2)*tan(c + d*x)/d + a*(a + b*tan(c + d*x))**3/(3*d) + (a + b*tan(c + d*x))**4/(4*d) + (a + b*tan(c + d*x))**2*(a**2 - b**2)/(2*d) - (a**4 - 6*a**2*b**2 + b**4)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_448():
    f = (a + b*tan(c + d*x))**4
    F = a*b*(a + b*tan(c + d*x))**2/d - 4*a*b*(a**2 - b**2)*log(cos(c + d*x))/d + b**2*(3*a**2 - b**2)*tan(c + d*x)/d + b*(a + b*tan(c + d*x))**3/(3*d) + x*(a**4 - 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_449():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)
    F = a**4*log(sin(c + d*x))/d + 3*a*b**3*tan(c + d*x)/d + 4*a*b*x*(a**2 - b**2) + b**2*(a + b*tan(c + d*x))**2/(2*d) - b**2*(6*a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_450():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**2
    F = 4*a**3*b*log(sin(c + d*x))/d - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)/d - 4*a*b**3*log(cos(c + d*x))/d + b**2*(a**2 + b**2)*tan(c + d*x)/d - x*(a**4 - 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_451():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**3
    F = -3*a**3*b*cot(c + d*x)/d - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)**2/(2*d) - a**2*(a**2 - 6*b**2)*log(sin(c + d*x))/d - 4*a*b*x*(a**2 - b**2) - b**4*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_452():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**4
    F = -4*a**3*b*cot(c + d*x)**2/(3*d) - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)**3/(3*d) + a**2*(3*a**2 - 17*b**2)*cot(c + d*x)/(3*d) - 4*a*b*(a**2 - b**2)*log(sin(c + d*x))/d + x*(a**4 - 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_453():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**5
    F = -5*a**3*b*cot(c + d*x)**3/(6*d) - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)**4/(4*d) + a**2*(2*a**2 - 11*b**2)*cot(c + d*x)**2/(4*d) + 4*a*b*x*(a**2 - b**2) + 4*a*b*(a**2 - b**2)*cot(c + d*x)/d + (a**4 - 6*a**2*b**2 + b**4)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_454():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**6
    F = -3*a**3*b*cot(c + d*x)**4/(5*d) - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)**5/(5*d) + a**2*(5*a**2 - 27*b**2)*cot(c + d*x)**3/(15*d) + 4*a*b*(a**2 - b**2)*log(sin(c + d*x))/d + 2*a*b*(a**2 - b**2)*cot(c + d*x)**2/d - x*(a**4 - 6*a**2*b**2 + b**4) - (a**4 - 6*a**2*b**2 + b**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_455():
    f = (a + b*tan(c + d*x))**4*cot(c + d*x)**7
    F = -7*a**3*b*cot(c + d*x)**5/(15*d) - a**2*(a + b*tan(c + d*x))**2*cot(c + d*x)**6/(6*d) + a**2*(3*a**2 - 16*b**2)*cot(c + d*x)**4/(12*d) - 4*a*b*x*(a**2 - b**2) + 4*a*b*(a**2 - b**2)*cot(c + d*x)**3/(3*d) - 4*a*b*(a**2 - b**2)*cot(c + d*x)/d - (a**4 - 6*a**2*b**2 + b**4)*log(sin(c + d*x))/d - (a**4 - 6*a**2*b**2 + b**4)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_456():
    f = tan(c + d*x)**6/(a + b*tan(c + d*x))
    F = a**6*log(a + b*tan(c + d*x))/(b**5*d*(a**2 + b**2)) - a*x/(a**2 + b**2) - a*tan(c + d*x)**3/(3*b**2*d) - a*(a**2 - b**2)*tan(c + d*x)/(b**4*d) - b*log(cos(c + d*x))/(d*(a**2 + b**2)) + tan(c + d*x)**4/(4*b*d) + (a**2 - b**2)*tan(c + d*x)**2/(2*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_457():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))
    F = -a**5*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)) - a*log(cos(c + d*x))/(d*(a**2 + b**2)) - a*tan(c + d*x)**2/(2*b**2*d) + b*x/(a**2 + b**2) + tan(c + d*x)**3/(3*b*d) + (a**2 - b**2)*tan(c + d*x)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_458():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))
    F = a**4*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)) + a*x/(a**2 + b**2) - a*tan(c + d*x)/(b**2*d) + b*log(cos(c + d*x))/(d*(a**2 + b**2)) + tan(c + d*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_459():
    f = tan(c + d*x)**3/(a + b*tan(c + d*x))
    F = -a**3*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)) + a*log(cos(c + d*x))/(d*(a**2 + b**2)) - b*x/(a**2 + b**2) + tan(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_460():
    f = tan(c + d*x)/(a + b*tan(c + d*x))
    F = -a*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)) + b*x/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_461():
    f = 1/(a + b*tan(c + d*x))
    F = a*x/(a**2 + b**2) + b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_462():
    f = cot(c + d*x)/(a + b*tan(c + d*x))
    F = -b*x/(a**2 + b**2) - b**2*log(a*cos(c + d*x) + b*sin(c + d*x))/(a*d*(a**2 + b**2)) + log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_463():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))
    F = -a*x/(a**2 + b**2) - cot(c + d*x)/(a*d) + b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)) - b*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_464():
    f = cot(c + d*x)**3/(a + b*tan(c + d*x))
    F = b*x/(a**2 + b**2) - cot(c + d*x)**2/(2*a*d) + b*cot(c + d*x)/(a**2*d) - b**4*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)) - (a**2 - b**2)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_465():
    f = cot(c + d*x)**4/(a + b*tan(c + d*x))
    F = a*x/(a**2 + b**2) - cot(c + d*x)**3/(3*a*d) + b*cot(c + d*x)**2/(2*a**2*d) + (a**2 - b**2)*cot(c + d*x)/(a**3*d) + b**5*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)) + b*(a**2 - b**2)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_466():
    f = tan(c + d*x)**6/(a + b*tan(c + d*x))**2
    F = -2*a**5*(2*a**2 + 3*b**2)*log(a + b*tan(c + d*x))/(b**5*d*(a**2 + b**2)**2) - a**2*tan(c + d*x)**4/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*a*b*log(cos(c + d*x))/(d*(a**2 + b**2)**2) - a*(2*a**2 + b**2)*tan(c + d*x)**2/(b**3*d*(a**2 + b**2)) - x*(a**2 - b**2)/(a**2 + b**2)**2 + (4*a**2 + b**2)*tan(c + d*x)**3/(3*b**2*d*(a**2 + b**2)) + (4*a**4 + 2*a**2*b**2 - b**4)*tan(c + d*x)/(b**4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_467():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))**2
    F = a**4*(3*a**2 + 5*b**2)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**2) - a**2*tan(c + d*x)**3/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + 2*a*b*x/(a**2 + b**2)**2 - a*(3*a**2 + 2*b**2)*tan(c + d*x)/(b**3*d*(a**2 + b**2)) - (a**2 - b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**2) + (3*a**2 + b**2)*tan(c + d*x)**2/(2*b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_468():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = -2*a**3*(a**2 + 2*b**2)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**2) - a**2*tan(c + d*x)**2/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + 2*a*b*log(cos(c + d*x))/(d*(a**2 + b**2)**2) + x*(a**2 - b**2)/(a**2 + b**2)**2 + (2*a**2 + b**2)*tan(c + d*x)/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_469():
    f = tan(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -a**2/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*a*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) - x*(a**2 - b**2)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_470():
    f = tan(c + d*x)/(a + b*tan(c + d*x))**2
    F = 2*a*b*x/(a**2 + b**2)**2 + a/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) - (a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_471():
    f = (a + b*tan(c + d*x))**(-2)
    F = 2*a*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) - b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) + x*(a**2 - b**2)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_472():
    f = cot(c + d*x)/(a + b*tan(c + d*x))**2
    F = -2*a*b*x/(a**2 + b**2)**2 + b**2/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - b**2*(3*a**2 + b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)**2) + log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_473():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -x*(a**2 - b**2)/(a**2 + b**2)**2 - cot(c + d*x)/(a*d*(a + b*tan(c + d*x))) - b*(a**2 + 2*b**2)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + 2*b**3*(2*a**2 + b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**2) - 2*b*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_474():
    f = cot(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = 2*a*b*x/(a**2 + b**2)**2 - cot(c + d*x)**2/(2*a*d*(a + b*tan(c + d*x))) + 3*b*cot(c + d*x)/(2*a**2*d*(a + b*tan(c + d*x))) + b**2*(2*a**2 + 3*b**2)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - b**4*(5*a**2 + 3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**2) - (a**2 - 3*b**2)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_475():
    f = tan(c + d*x)**6/(a + b*tan(c + d*x))**3
    F = a**4*(6*a**4 + 17*a**2*b**2 + 15*b**4)*log(a + b*tan(c + d*x))/(b**5*d*(a**2 + b**2)**3) - a**2*tan(c + d*x)**4/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - 2*a**2*(a**2 + 2*b**2)*tan(c + d*x)**3/(b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 - a*(6*a**4 + 11*a**2*b**2 + 3*b**4)*tan(c + d*x)/(b**4*d*(a**2 + b**2)**2) - b*(3*a**2 - b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**3) + (6*a**4 + 11*a**2*b**2 + b**4)*tan(c + d*x)**2/(2*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_476():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))**3
    F = -a**3*(3*a**4 + 9*a**2*b**2 + 10*b**4)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**3) - a**2*tan(c + d*x)**3/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a**2*(3*a**2 + 7*b**2)*tan(c + d*x)**2/(2*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - a*(a**2 - 3*b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**3) + b*x*(3*a**2 - b**2)/(a**2 + b**2)**3 + (3*a**4 + 6*a**2*b**2 + b**4)*tan(c + d*x)/(b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_477():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = a**3*(a**2 + 3*b**2)/(b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - a**2*tan(c + d*x)**2/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a**2*(a**4 + 3*a**2*b**2 + 6*b**4)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**3) + a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 + b*(3*a**2 - b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_478():
    f = tan(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = -a**2*tan(c + d*x)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a**2*(a**2 + 5*b**2)/(2*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + a*(a**2 - 3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) - b*x*(3*a**2 - b**2)/(a**2 + b**2)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_479():
    f = tan(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -a**2/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + 2*a*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 - b*(3*a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_480():
    f = tan(c + d*x)/(a + b*tan(c + d*x))**3
    F = -a*(a**2 - 3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + a/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2)) + b*x*(3*a**2 - b**2)/(a**2 + b**2)**3 + (a**2 - b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_481():
    f = (a + b*tan(c + d*x))**(-3)
    F = -2*a*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 + b*(3*a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) - b/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_482():
    f = cot(c + d*x)/(a + b*tan(c + d*x))**3
    F = -b*x*(3*a**2 - b**2)/(a**2 + b**2)**3 + b**2/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b**2*(3*a**2 + b**2)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b**2*(6*a**4 + 3*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**3) + log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_483():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 - cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**2) - b*(2*a**2 + 3*b**2)/(2*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - b*(a**4 + 6*a**2*b**2 + 3*b**4)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + b**3*(10*a**4 + 9*a**2*b**2 + 3*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**3) - 3*b*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_484():
    f = tan(c + d*x)**6/(a + b*tan(c + d*x))**4
    F = -4*a**3*(a**6 + 4*a**4*b**2 + 6*a**2*b**4 + 5*b**6)*log(a + b*tan(c + d*x))/(b**5*d*(a**2 + b**2)**4) - a**2*tan(c + d*x)**4/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - a**2*(2*a**2 + 5*b**2)*tan(c + d*x)**3/(3*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - 2*a**2*(a**4 + 3*a**2*b**2 + 4*b**4)*tan(c + d*x)**2/(b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - 4*a*b*(a**2 - b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**4) - x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4 + (4*a**6 + 12*a**4*b**2 + 13*a**2*b**4 + b**6)*tan(c + d*x)/(b**4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_485():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))**4
    F = a**3*(a**4 + 3*a**2*b**2 + 6*b**4)/(b**4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - a**2*tan(c + d*x)**3/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - a**2*(a**2 + 3*b**2)*tan(c + d*x)**2/(2*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + a**2*(a**6 + 4*a**4*b**2 + 5*a**2*b**4 + 10*b**6)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**4) + 4*a*b*x*(a**2 - b**2)/(a**2 + b**2)**4 - (a**4 - 6*a**2*b**2 + b**4)*log(cos(c + d*x))/(d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_486():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))**4
    F = a**3*(a**2 + 4*b**2)/(3*b**3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - a**2*tan(c + d*x)**2/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - a**2*(2*a**4 + 7*a**2*b**2 + 17*b**4)/(3*b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + 4*a*b*(a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_487():
    f = tan(c + d*x)**3/(a + b*tan(c + d*x))**4
    F = -a**2*tan(c + d*x)/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - a**2*(a**2 + 7*b**2)/(6*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - 4*a*b*x*(a**2 - b**2)/(a**2 + b**2)**4 - a*(a**2 - 3*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + (a**4 - 6*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_488():
    f = tan(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -a**2/(3*b*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - 4*a*b*(a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + a*b/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + b*(3*a**2 - b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_489():
    f = tan(c + d*x)/(a + b*tan(c + d*x))**4
    F = 4*a*b*x*(a**2 - b**2)/(a**2 + b**2)**4 + a*(a**2 - 3*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + a/(d*(a + b*tan(c + d*x))**3*(3*a**2 + 3*b**2)) - (a**4 - 6*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + (a**2 - b**2)/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_490():
    f = (a + b*tan(c + d*x))**(-4)
    F = 4*a*b*(a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) - a*b/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - b*(3*a**2 - b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b/(d*(a + b*tan(c + d*x))**3*(3*a**2 + 3*b**2)) + x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_491():
    f = cot(c + d*x)/(a + b*tan(c + d*x))**4
    F = -4*a*b*x*(a**2 - b**2)/(a**2 + b**2)**4 + b**2/(3*a*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) + b**2*(3*a**2 + b**2)/(2*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + b**2*(6*a**4 + 3*a**2*b**2 + b**4)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b**2*(10*a**6 + 5*a**4*b**2 + 4*a**2*b**4 + b**6)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**4) + log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_492():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4 - cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**3) - b*(3*a**2 + 4*b**2)/(3*a**2*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)) - b*(a**4 + 4*a**2*b**2 + 2*b**4)/(a**3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - b*(a**6 + 13*a**4*b**2 + 12*a**2*b**4 + 4*b**6)/(a**4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + 4*b**3*(5*a**6 + 6*a**4*b**2 + 4*a**2*b**4 + b**6)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**5*d*(a**2 + b**2)**4) - 4*b*log(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_493():
    f = 1/(5*tan(c + d*x) + 3)
    F = 3*x/34 + 5*log(5*sin(c + d*x) + 3*cos(c + d*x))/(34*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_494():
    f = (5*tan(c + d*x) + 3)**(-2)
    F = -4*x/289 + 15*log(5*sin(c + d*x) + 3*cos(c + d*x))/(578*d) - 5/(34*d*(5*tan(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_495():
    f = (5*tan(c + d*x) + 3)**(-3)
    F = -99*x/19652 + 5*log(5*sin(c + d*x) + 3*cos(c + d*x))/(19652*d) - 15/(578*d*(5*tan(c + d*x) + 3)) - 5/(68*d*(5*tan(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_496():
    f = (5*tan(c + d*x) + 3)**(-4)
    F = -161*x/334084 - 60*log(5*sin(c + d*x) + 3*cos(c + d*x))/(83521*d) - 5/(19652*d*(5*tan(c + d*x) + 3)) - 15/(1156*d*(5*tan(c + d*x) + 3)**2) - 5/(102*d*(5*tan(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_497():
    f = 1/(3*tan(c + d*x) + 5)
    F = 5*x/34 + 3*log(3*sin(c + d*x) + 5*cos(c + d*x))/(34*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_498():
    f = (3*tan(c + d*x) + 5)**(-2)
    F = 4*x/289 + 15*log(3*sin(c + d*x) + 5*cos(c + d*x))/(578*d) - 3/(34*d*(3*tan(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_499():
    f = (3*tan(c + d*x) + 5)**(-3)
    F = -5*x/19652 + 99*log(3*sin(c + d*x) + 5*cos(c + d*x))/(19652*d) - 15/(578*d*(3*tan(c + d*x) + 5)) - 3/(68*d*(3*tan(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_500():
    f = (3*tan(c + d*x) + 5)**(-4)
    F = -161*x/334084 + 60*log(3*sin(c + d*x) + 5*cos(c + d*x))/(83521*d) - 99/(19652*d*(3*tan(c + d*x) + 5)) - 15/(1156*d*(3*tan(c + d*x) + 5)**2) - 1/(34*d*(3*tan(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_501():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**4
    F = -8*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(35*b**2*d) + sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + 2*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**2/(7*b*d) + (a + b*tan(c + d*x))**(sympy.S(3)/2)*(16*a**2 - 70*b**2)/(105*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_502():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**3
    F = -4*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(15*b**2*d) + sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - 2*sqrt(a + b*tan(c + d*x))/d + 2*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_503():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**2
    F = -sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + 2*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_504():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)
    F = -sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*sqrt(a + b*tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_505():
    f = sqrt(a + b*tan(c + d*x))
    F = sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_506():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_507():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2
    F = -sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) + sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))) - sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) + sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))) - sqrt(a + b*tan(c + d*x))*cot(c + d*x)/d - b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_508():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)**3
    F = -sqrt(a - I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - sqrt(a + I*b)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(4*a*d) + (8*a**2 + b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_509():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**4
    F = -8*a*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(63*b**2*d) + 2*b*sqrt(a + b*tan(c + d*x))/d - I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**2/(9*b*d) + (a + b*tan(c + d*x))**(sympy.S(5)/2)*(16*a**2 - 126*b**2)/(315*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_510():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**3
    F = -2*a*sqrt(a + b*tan(c + d*x))/d - 4*a*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(35*b**2*d) + (a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - 2*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*(a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_511():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**2
    F = -2*b*sqrt(a + b*tan(c + d*x))/d + I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_512():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)
    F = 2*a*sqrt(a + b*tan(c + d*x))/d - (a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_513():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*b*sqrt(a + b*tan(c + d*x))/d - I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_514():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + (a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_515():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**2
    F = -3*sqrt(a)*b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d - a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/d + I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_516():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**3
    F = -a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - 5*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(4*d) - (a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + (8*a**2 - 3*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_517():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**3
    F = -2*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) - 4*a*(a + b*tan(c + d*x))**(sympy.S(7)/2)/(63*b**2*d) + (a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d - 2*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) - sqrt(a + b*tan(c + d*x))*(2*a**2 - 2*b**2)/d + 2*(a + b*tan(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(9*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_518():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**2
    F = -4*a*b*sqrt(a + b*tan(c + d*x))/d - 2*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*(a + b*tan(c + d*x))**(sympy.S(7)/2)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_519():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)
    F = 2*a*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) - (a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + 2*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) + sqrt(a + b*tan(c + d*x))*(2*a**2 - 2*b**2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_520():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 4*a*b*sqrt(a + b*tan(c + d*x))/d + 2*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) - I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_521():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d + 2*b**2*sqrt(a + b*tan(c + d*x))/d + (a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + (a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_522():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**2
    F = -5*a**(sympy.S(3)/2)*b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/d - a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/d + I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_523():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**3
    F = sqrt(a)*(8*a**2 - 15*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*d) - a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*d) - 9*a*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(4*d) - (a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d - (a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_524():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**4
    F = -a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**3/(3*d) - 13*a*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(12*d) - I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d + sqrt(a + b*tan(c + d*x))*(8*a**2 - 11*b**2)*cot(c + d*x)/(8*d) + 5*b*(8*a**2 - b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_525():
    f = (a + b*tan(c + d*x))**(sympy.S(7)/2)
    F = 4*a*b*(a + b*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*b*(a + b*tan(c + d*x))**(sympy.S(5)/2)/(5*d) + 2*b*sqrt(a + b*tan(c + d*x))*(3*a**2 - b**2)/d - I*(a - I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/d + I*(a + I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_526():
    f = tan(c + d*x)**5/sqrt(a + b*tan(c + d*x))
    F = -12*a*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**2/(35*b**2*d) - 4*a*sqrt(a + b*tan(c + d*x))*(24*a**2 - 35*b**2)/(105*b**4*d) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + 2*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**3/(7*b*d) + sqrt(a + b*tan(c + d*x))*(48*a**2 - 70*b**2)*tan(c + d*x)/(105*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_527():
    f = tan(c + d*x)**4/sqrt(a + b*tan(c + d*x))
    F = -8*a*sqrt(a + b*tan(c + d*x))*tan(c + d*x)/(15*b**2*d) - sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + 2*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**2/(5*b*d) + sqrt(a + b*tan(c + d*x))*(16*a**2 - 30*b**2)/(15*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_528():
    f = tan(c + d*x)**3/sqrt(a + b*tan(c + d*x))
    F = -4*a*sqrt(a + b*tan(c + d*x))/(3*b**2*d) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) + 2*sqrt(a + b*tan(c + d*x))*tan(c + d*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_529():
    f = tan(c + d*x)**2/sqrt(a + b*tan(c + d*x))
    F = sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + 2*sqrt(a + b*tan(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_530():
    f = tan(c + d*x)/sqrt(a + b*tan(c + d*x))
    F = -atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_531():
    f = 1/sqrt(a + b*tan(c + d*x))
    F = -sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_532():
    f = cot(c + d*x)/sqrt(a + b*tan(c + d*x))
    F = atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) - 2*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_533():
    f = cot(c + d*x)**2/sqrt(a + b*tan(c + d*x))
    F = sqrt(2)*b*log(a + b*tan(c + d*x) - sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*log(a + b*tan(c + d*x) + sqrt(2)*sqrt(a + b*tan(c + d*x))*sqrt(a + sqrt(a**2 + b**2)) + sqrt(a**2 + b**2))/(4*d*sqrt(a + sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(2)*b*atanh((-sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) + sqrt(2)*b*atanh((sqrt(2)*sqrt(a + b*tan(c + d*x)) + sqrt(a + sqrt(a**2 + b**2)))/sqrt(a - sqrt(a**2 + b**2)))/(2*d*sqrt(a - sqrt(a**2 + b**2))*sqrt(a**2 + b**2)) - sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(a*d) + b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_534():
    f = cot(c + d*x)**3/sqrt(a + b*tan(c + d*x))
    F = -atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*sqrt(a + I*b)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*sqrt(a - I*b)) - sqrt(a + b*tan(c + d*x))*cot(c + d*x)**2/(2*a*d) + 3*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)/(4*a**2*d) + (8*a**2 - 3*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_535():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)**3/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*a*sqrt(a + b*tan(c + d*x))*(8*a**2 + 3*b**2)*tan(c + d*x)/(5*b**3*d*(a**2 + b**2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(12*a**2 + 2*b**2)*tan(c + d*x)**2/(5*b**2*d*(a**2 + b**2)) + sqrt(a + b*tan(c + d*x))*(32*a**4 + 12*a**2*b**2 - 10*b**4)/(5*b**4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_536():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)**2/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*a*sqrt(a + b*tan(c + d*x))*(8*a**2 + 5*b**2)/(3*b**3*d*(a**2 + b**2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(8*a**2 + 2*b**2)*tan(c + d*x)/(3*b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_537():
    f = tan(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(4*a**2 + 2*b**2)/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_538():
    f = tan(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_539():
    f = tan(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*a/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_540():
    f = (a + b*tan(c + d*x))**(sympy.S(-3)/2)
    F = -2*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_541():
    f = cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) + 2*b**2/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - 2*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_542():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) - cot(c + d*x)/(a*d*sqrt(a + b*tan(c + d*x))) - b*(a**2 + 3*b**2)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + 3*b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_543():
    f = cot(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(3)/2)) - cot(c + d*x)**2/(2*a*d*sqrt(a + b*tan(c + d*x))) + 5*b*cot(c + d*x)/(4*a**2*d*sqrt(a + b*tan(c + d*x))) + b**2*(7*a**2 + 15*b**2)/(4*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + (8*a**2 - 15*b**2)*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_544():
    f = tan(c + d*x)**5/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)**3/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 4*a**2*(a**2 + 2*b**2)*tan(c + d*x)**2/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - 4*a*sqrt(a + b*tan(c + d*x))*(8*a**4 + 15*a**2*b**2 + 4*b**4)/(3*b**4*d*(a**2 + b**2)**2) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(16*a**4 + 30*a**2*b**2 + 2*b**4)*tan(c + d*x)/(3*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_545():
    f = tan(c + d*x)**4/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 4*a**3*(2*a**2 + 5*b**2)/(3*b**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - 2*a**2*tan(c + d*x)**2/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(8*a**2 + 6*b**2)/(3*b**3*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_546():
    f = tan(c + d*x)**3/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 4*a**2*(a**2 + 4*b**2)/(3*b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_547():
    f = tan(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 4*a*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_548():
    f = tan(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*a/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + (2*a**2 - 2*b**2)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_549():
    f = (a + b*tan(c + d*x))**(sympy.S(-5)/2)
    F = -4*a*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - 2*b/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_550():
    f = cot(c + d*x)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) + 2*b**2/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*b**2*(3*a**2 + b**2)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - 2*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_551():
    f = cot(c + d*x)**2/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(5)/2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(5)/2)) - cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) - b*(3*a**2 + 5*b**2)/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - b*(a**4 + 10*a**2*b**2 + 5*b**4)/(a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) + 5*b*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_552():
    f = (a + b*tan(c + d*x))**(sympy.S(-7)/2)
    F = -4*a*b/(3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)**2) - 2*b*(3*a**2 - b**2)/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**3) - 2*b/(d*(a + b*tan(c + d*x))**(sympy.S(5)/2)*(5*a**2 + 5*b**2)) + I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a + I*b))/(d*(a + I*b)**(sympy.S(7)/2)) - I*atanh(sqrt(a + b*tan(c + d*x))/sqrt(a - I*b))/(d*(a - I*b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_553():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)
    F = 2*a*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - 2*b*sqrt(tan(c + d*x))/d - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_554():
    f = (a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sqrt(tan(c + d*x))/d + 2*b*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_555():
    f = (a + b*tan(c + d*x))*sqrt(tan(c + d*x))
    F = 2*b*sqrt(tan(c + d*x))/d + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_556():
    f = (a + b*tan(c + d*x))/sqrt(tan(c + d*x))
    F = -sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_557():
    f = (a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*a/(d*sqrt(tan(c + d*x))) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_558():
    f = (a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*a/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*b/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_559():
    f = (a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(7)/2)
    F = 2*a/(d*sqrt(tan(c + d*x))) - 2*a/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*b/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_560():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(5)/2)
    F = 4*a*b*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - 4*a*b*sqrt(tan(c + d*x))/d + 2*b**2*tan(c + d*x)**(sympy.S(7)/2)/(7*d) + (2*a**2 - 2*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_561():
    f = (a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b**2*tan(c + d*x)**(sympy.S(5)/2)/(5*d) + (2*a**2 - 2*b**2)*sqrt(tan(c + d*x))/d + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_562():
    f = (a + b*tan(c + d*x))**2*sqrt(tan(c + d*x))
    F = 4*a*b*sqrt(tan(c + d*x))/d + 2*b**2*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_563():
    f = (a + b*tan(c + d*x))**2/sqrt(tan(c + d*x))
    F = 2*b**2*sqrt(tan(c + d*x))/d - sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_564():
    f = (a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2/(d*sqrt(tan(c + d*x))) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_565():
    f = (a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*a**2/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 4*a*b/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_566():
    f = (a + b*tan(c + d*x))**2/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*a**2/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 4*a*b/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + (2*a**2 - 2*b**2)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_567():
    f = (a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(5)/2)
    F = 40*a*b**2*tan(c + d*x)**(sympy.S(7)/2)/(63*d) + 2*a*(a**2 - 3*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b**2*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*b*(3*a**2 - b**2)*tan(c + d*x)**(sympy.S(5)/2)/(5*d) - 2*b*(3*a**2 - b**2)*sqrt(tan(c + d*x))/d - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_568():
    f = (a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(3)/2)
    F = 32*a*b**2*tan(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*a*(a**2 - 3*b**2)*sqrt(tan(c + d*x))/d + 2*b**2*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*b*(3*a**2 - b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_569():
    f = (a + b*tan(c + d*x))**3*sqrt(tan(c + d*x))
    F = 8*a*b**2*tan(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*b**2*(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*b*(3*a**2 - b**2)*sqrt(tan(c + d*x))/d + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_570():
    f = (a + b*tan(c + d*x))**3/sqrt(tan(c + d*x))
    F = 16*a*b**2*sqrt(tan(c + d*x))/(3*d) + 2*b**2*(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(3*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_571():
    f = (a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2*(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) + 2*b*(a**2 + b**2)*sqrt(tan(c + d*x))/d - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_572():
    f = (a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(5)/2)
    F = -16*a**2*b/(3*d*sqrt(tan(c + d*x))) - 2*a**2*(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_573():
    f = (a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(7)/2)
    F = -8*a**2*b/(5*d*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2*(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + 2*a*(a**2 - 3*b**2)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_574():
    f = (a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(9)/2)
    F = -32*a**2*b/(35*d*tan(c + d*x)**(sympy.S(5)/2)) - 2*a**2*(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(7)/2)) + 2*a*(a**2 - 3*b**2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*b*(3*a**2 - b**2)/(d*sqrt(tan(c + d*x))) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_575():
    f = (a + b*tan(c + d*x))**3/tan(c + d*x)**(sympy.S(11)/2)
    F = -40*a**2*b/(63*d*tan(c + d*x)**(sympy.S(7)/2)) - 2*a**2*(a + b*tan(c + d*x))/(9*d*tan(c + d*x)**(sympy.S(9)/2)) - 2*a*(a**2 - 3*b**2)/(d*sqrt(tan(c + d*x))) + 2*a*(a**2 - 3*b**2)/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + 2*b*(3*a**2 - b**2)/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_576():
    f = (a + b*tan(c + d*x))/sqrt(tan(c + d*x))
    F = -sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_577():
    f = (a + b*tan(c + d*x))/sqrt(-tan(c + d*x))
    F = -sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(-tan(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(-tan(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(-tan(c + d*x)) - tan(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(-tan(c + d*x)) - tan(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_578():
    f = (a + b*tan(c + d*x))/sqrt(e*tan(c + d*x))
    F = -sqrt(2)*(a - b)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) + sqrt(2)*(a - b)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*(a + b)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e)) + sqrt(2)*(a + b)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_579():
    f = (a + b*tan(c + d*x))/sqrt(-e*tan(c + d*x))
    F = sqrt(2)*(a - b)*atan(1 - sqrt(2)*sqrt(-e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e)) - sqrt(2)*(a - b)*atan(1 + sqrt(2)*sqrt(-e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e)) + sqrt(2)*(a + b)*log(-sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(-e*tan(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*(a + b)*log(-sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(-e*tan(c + d*x)))/(4*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_580():
    f = tan(c + d*x)**(sympy.S(9)/2)/(a + b*tan(c + d*x))
    F = -2*a**(sympy.S(9)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(7)/2)*d*(a**2 + b**2)) - 2*a*tan(c + d*x)**(sympy.S(3)/2)/(3*b**2*d) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + 2*tan(c + d*x)**(sympy.S(5)/2)/(5*b*d) + (2*a**2 - 2*b**2)*sqrt(tan(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_581():
    f = tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))
    F = 2*a**(sympy.S(7)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)) - 2*a*sqrt(tan(c + d*x))/(b**2*d) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + 2*tan(c + d*x)**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_582():
    f = tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + 2*sqrt(tan(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_583():
    f = tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(b)*d*(a**2 + b**2)) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_584():
    f = sqrt(tan(c + d*x))/(a + b*tan(c + d*x))
    F = -2*sqrt(a)*sqrt(b)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_585():
    f = 1/((a + b*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + 2*b**(sympy.S(3)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(a)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_586():
    f = 1/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - 2/(a*d*sqrt(tan(c + d*x))) - 2*b**(sympy.S(5)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_587():
    f = 1/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - 2/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*b/(a**2*d*sqrt(tan(c + d*x))) + 2*b**(sympy.S(7)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_588():
    f = 1/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2))
    F = sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - 2/(5*a*d*tan(c + d*x)**(sympy.S(5)/2)) + 2*b/(3*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + (2*a**2 - 2*b**2)/(a**3*d*sqrt(tan(c + d*x))) - 2*b**(sympy.S(9)/2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(7)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_589():
    f = tan(c + d*x)**(sympy.S(9)/2)/(a + b*tan(c + d*x))**2
    F = a**(sympy.S(7)/2)*(5*a**2 + 9*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(7)/2)*d*(a**2 + b**2)**2) - a**2*tan(c + d*x)**(sympy.S(5)/2)/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - a*(5*a**2 + 4*b**2)*sqrt(tan(c + d*x))/(b**3*d*(a**2 + b**2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + (5*a**2 + 2*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_590():
    f = tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))**2
    F = -a**(sympy.S(5)/2)*(3*a**2 + 7*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)**2) - a**2*tan(c + d*x)**(sympy.S(3)/2)/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + (3*a**2 + 2*b**2)*sqrt(tan(c + d*x))/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_591():
    f = tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**2
    F = a**(sympy.S(3)/2)*(a**2 + 5*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)**2) - a**2*sqrt(tan(c + d*x))/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_592():
    f = tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**2
    F = sqrt(a)*(a**2 - 3*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(b)*d*(a**2 + b**2)**2) + a*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_593():
    f = sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**2
    F = -b*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(b)*(3*a**2 - b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(sqrt(a)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_594():
    f = 1/((a + b*tan(c + d*x))**2*sqrt(tan(c + d*x)))
    F = sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b**2*sqrt(tan(c + d*x))/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + b**(sympy.S(3)/2)*(5*a**2 + b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_595():
    f = 1/((a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + b**2/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(tan(c + d*x))) - (2*a**2 + 3*b**2)/(a**2*d*(a**2 + b**2)*sqrt(tan(c + d*x))) - b**(sympy.S(5)/2)*(7*a**2 + 3*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(5)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_596():
    f = 1/((a + b*tan(c + d*x))**2*tan(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b**2/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)*tan(c + d*x)**(sympy.S(3)/2)) - (2*a**2 + 5*b**2)/(3*a**2*d*(a**2 + b**2)*tan(c + d*x)**(sympy.S(3)/2)) + b*(4*a**2 + 5*b**2)/(a**3*d*(a**2 + b**2)*sqrt(tan(c + d*x))) + b**(sympy.S(7)/2)*(9*a**2 + 5*b**2)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(a**(sympy.S(7)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_597():
    f = tan(c + d*x)**(sympy.S(11)/2)/(a + b*tan(c + d*x))**3
    F = a**(sympy.S(7)/2)*(35*a**4 + 102*a**2*b**2 + 99*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(9)/2)*d*(a**2 + b**2)**3) - a**2*tan(c + d*x)**(sympy.S(7)/2)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a**2*(7*a**2 + 15*b**2)*tan(c + d*x)**(sympy.S(5)/2)/(4*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - a*(35*a**4 + 67*a**2*b**2 + 24*b**4)*sqrt(tan(c + d*x))/(4*b**4*d*(a**2 + b**2)**2) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + (35*a**4 + 67*a**2*b**2 + 8*b**4)*tan(c + d*x)**(sympy.S(3)/2)/(12*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_598():
    f = tan(c + d*x)**(sympy.S(9)/2)/(a + b*tan(c + d*x))**3
    F = -a**(sympy.S(5)/2)*(15*a**4 + 46*a**2*b**2 + 63*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(7)/2)*d*(a**2 + b**2)**3) - a**2*tan(c + d*x)**(sympy.S(5)/2)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a**2*(5*a**2 + 13*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(4*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + (15*a**4 + 31*a**2*b**2 + 8*b**4)*sqrt(tan(c + d*x))/(4*b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_599():
    f = tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))**3
    F = a**(sympy.S(3)/2)*(3*a**4 + 6*a**2*b**2 + 35*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(5)/2)*d*(a**2 + b**2)**3) - a**2*tan(c + d*x)**(sympy.S(3)/2)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - a**2*(3*a**2 + 11*b**2)*sqrt(tan(c + d*x))/(4*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_600():
    f = tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**3
    F = sqrt(a)*(a**4 + 18*a**2*b**2 - 15*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*b**(sympy.S(3)/2)*d*(a**2 + b**2)**3) - a**2*sqrt(tan(c + d*x))/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(a**2 + 9*b**2)*sqrt(tan(c + d*x))/(4*b*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_601():
    f = tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**3
    F = a*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + (3*a**2 - 5*b**2)*sqrt(tan(c + d*x))/(4*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (3*a**4 - 26*a**2*b**2 + 3*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*sqrt(a)*sqrt(b)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_602():
    f = sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**3
    F = -b*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - b*(7*a**2 - b**2)*sqrt(tan(c + d*x))/(4*a*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - sqrt(b)*(15*a**4 - 18*a**2*b**2 - b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(3)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_603():
    f = 1/((a + b*tan(c + d*x))**3*sqrt(tan(c + d*x)))
    F = -sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b**2*sqrt(tan(c + d*x))/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b**2*(11*a**2 + 3*b**2)*sqrt(tan(c + d*x))/(4*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + b**(sympy.S(3)/2)*(35*a**4 + 6*a**2*b**2 + 3*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(5)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_604():
    f = 1/((a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b**2/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)*sqrt(tan(c + d*x))) + b**2*(13*a**2 + 5*b**2)/(4*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(tan(c + d*x))) - (8*a**4 + 31*a**2*b**2 + 15*b**4)/(4*a**3*d*(a**2 + b**2)**2*sqrt(tan(c + d*x))) - b**(sympy.S(5)/2)*(63*a**4 + 46*a**2*b**2 + 15*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(7)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_605():
    f = 1/((a + b*tan(c + d*x))**3*tan(c + d*x)**(sympy.S(5)/2))
    F = sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b**2/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)*tan(c + d*x)**(sympy.S(3)/2)) + b**2*(15*a**2 + 7*b**2)/(4*a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2*tan(c + d*x)**(sympy.S(3)/2)) - (8*a**4 + 67*a**2*b**2 + 35*b**4)/(12*a**3*d*(a**2 + b**2)**2*tan(c + d*x)**(sympy.S(3)/2)) + b*(24*a**4 + 67*a**2*b**2 + 35*b**4)/(4*a**4*d*(a**2 + b**2)**2*sqrt(tan(c + d*x))) + b**(sympy.S(7)/2)*(99*a**4 + 102*a**2*b**2 + 35*b**4)*atan(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a))/(4*a**(sympy.S(9)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_606():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)
    F = -a*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(4*b*d) - sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(2*b*d) - (a**2 + 8*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_607():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)
    F = a*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d) + sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/d + I*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_608():
    f = sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))
    F = 2*sqrt(b)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_609():
    f = sqrt(a + b*tan(c + d*x))/sqrt(tan(c + d*x))
    F = -I*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_610():
    f = sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*sqrt(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) - sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_611():
    f = sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*sqrt(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) + I*sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*b*sqrt(a + b*tan(c + d*x))/(3*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_612():
    f = sqrt(a + b*tan(c + d*x))/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(I*a - b)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(I*a + b)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*b*sqrt(a + b*tan(c + d*x))/(15*a*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(30*a**2 + 4*b**2)/(15*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_613():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2)
    F = -a*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(12*b*d) - a*(a**2 + 24*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*b**(sympy.S(3)/2)*d) + I*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x))/(3*b*d) - sqrt(a + b*tan(c + d*x))*(a**2 + 8*b**2)*sqrt(tan(c + d*x))/(8*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_614():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = 3*a*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(4*d) + (a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))/(2*d) + (I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (3*a**2 - 8*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_615():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))
    F = 3*a*sqrt(b)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/d - I*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_616():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/sqrt(tan(c + d*x))
    F = 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_617():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) + I*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_618():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 8*b*sqrt(a + b*tan(c + d*x))/(3*d*sqrt(tan(c + d*x))) + (I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_619():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 4*b*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(3)/2)) - I*(I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(10*a**2 - 2*b**2)/(5*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_620():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - 16*b*sqrt(a + b*tan(c + d*x))/(35*d*tan(c + d*x)**(sympy.S(5)/2)) - (I*a - b)**(sympy.S(3)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(3)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(70*a**2 - 6*b**2)/(105*a*d*tan(c + d*x)**(sympy.S(3)/2)) + 4*b*sqrt(a + b*tan(c + d*x))*(70*a**2 + 3*b**2)/(105*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_621():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)
    F = -a*(a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x))/(24*b*d) - a*sqrt(a + b*tan(c + d*x))*(5*a**2 + 112*b**2)*sqrt(tan(c + d*x))/(64*b*d) + (I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (a + b*tan(c + d*x))**(sympy.S(7)/2)*sqrt(tan(c + d*x))/(4*b*d) - (a + b*tan(c + d*x))**(sympy.S(3)/2)*(5*a**2 + 48*b**2)*sqrt(tan(c + d*x))/(96*b*d) - (5*a**4 + 240*a**2*b**2 - 128*b**4)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(64*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_622():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2)
    F = 13*a*b*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(12*d) + 5*a*(a**2 - 8*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*sqrt(b)*d) + b**2*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2)/(3*d) + sqrt(a + b*tan(c + d*x))*(11*a**2 - 8*b**2)*sqrt(tan(c + d*x))/(8*d) - I*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_623():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x))
    F = 9*a*b*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(4*d) + sqrt(b)*(15*a**2 - 8*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*d) + b**2*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(2*d) - (I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_624():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/sqrt(tan(c + d*x))
    F = 5*a*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b**2*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/d + I*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_625():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))/(d*sqrt(tan(c + d*x))) + 2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_626():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(5)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))/(3*d*tan(c + d*x)**(sympy.S(3)/2)) - 14*a*b*sqrt(a + b*tan(c + d*x))/(3*d*sqrt(tan(c + d*x))) - I*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_627():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(7)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))/(5*d*tan(c + d*x)**(sympy.S(5)/2)) - 22*a*b*sqrt(a + b*tan(c + d*x))/(15*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(30*a**2 - 46*b**2)/(15*d*sqrt(tan(c + d*x))) - (I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_628():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(9)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(7)/2)) - 6*a*b*sqrt(a + b*tan(c + d*x))/(7*d*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(14*a**2 - 18*b**2)/(21*d*tan(c + d*x)**(sympy.S(3)/2)) + I*(I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + 2*b*sqrt(a + b*tan(c + d*x))*(49*a**2 - 3*b**2)/(21*a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_629():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/tan(c + d*x)**(sympy.S(11)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))/(9*d*tan(c + d*x)**(sympy.S(9)/2)) - 38*a*b*sqrt(a + b*tan(c + d*x))/(63*d*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*tan(c + d*x))*(42*a**2 - 50*b**2)/(105*d*tan(c + d*x)**(sympy.S(5)/2)) + (I*a - b)**(sympy.S(5)/2)*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(5)/2)*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + 2*b*sqrt(a + b*tan(c + d*x))*(231*a**2 - 5*b**2)/(315*a*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*tan(c + d*x))*(630*a**4 - 966*a**2*b**2 - 20*b**4)/(315*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_630():
    f = tan(c + d*x)**(sympy.S(7)/2)/sqrt(a + b*tan(c + d*x))
    F = -3*a*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(4*b**2*d) + atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)/(2*b*d) + (3*a**2 - 8*b**2)*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_631():
    f = tan(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*tan(c + d*x))
    F = -a*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d) + I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_632():
    f = tan(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*tan(c + d*x))
    F = -atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + 2*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_633():
    f = sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = -I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_634():
    f = 1/(sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_635():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2))
    F = I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*sqrt(a + b*tan(c + d*x))/(a*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_636():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/2))
    F = -atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*sqrt(a + b*tan(c + d*x))/(3*a*d*tan(c + d*x)**(sympy.S(3)/2)) + 4*b*sqrt(a + b*tan(c + d*x))/(3*a**2*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_637():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(7)/2))
    F = -I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*sqrt(a + b*tan(c + d*x))/(5*a*d*tan(c + d*x)**(sympy.S(5)/2)) + 8*b*sqrt(a + b*tan(c + d*x))/(15*a**2*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(30*a**2 - 16*b**2)/(15*a**3*d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_638():
    f = tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)**(sympy.S(3)/2)/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - 3*a*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(5)/2)*d) + I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(3*a**2 + b**2)*sqrt(tan(c + d*x))/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_639():
    f = tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sqrt(tan(c + d*x))/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + 2*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_640():
    f = tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sqrt(tan(c + d*x))/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) - I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_641():
    f = sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*b*sqrt(tan(c + d*x))/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)) + atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_642():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x)))
    F = I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + 2*b**2*sqrt(tan(c + d*x))/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_643():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = -atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2/(a*d*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))) - 2*b*(a**2 + 2*b**2)*sqrt(tan(c + d*x))/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_644():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = -I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2/(3*a*d*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(3)/2)) + 8*b/(3*a**2*d*sqrt(a + b*tan(c + d*x))*sqrt(tan(c + d*x))) + 2*b**2*(5*a**2 + 8*b**2)*sqrt(tan(c + d*x))/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_645():
    f = tan(c + d*x)**(sympy.S(9)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)**(sympy.S(5)/2)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 2*a**2*(5*a**2 + 11*b**2)*tan(c + d*x)**(sympy.S(3)/2)/(3*b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - 5*a*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(7)/2)*d) + I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(5*a**4 + 10*a**2*b**2 + b**4)*sqrt(tan(c + d*x))/(b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_646():
    f = tan(c + d*x)**(sympy.S(7)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)**(sympy.S(3)/2)/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 2*a**2*(a**2 + 3*b**2)*sqrt(tan(c + d*x))/(b**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + 2*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_647():
    f = tan(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sqrt(tan(c + d*x))/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 2*a*(a**2 + 7*b**2)*sqrt(tan(c + d*x))/(3*b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2) - I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_648():
    f = tan(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*a**2 - 8*b**2)*sqrt(tan(c + d*x))/(3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_649():
    f = sqrt(tan(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -2*b*sqrt(tan(c + d*x))/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2*b*(5*a**2 - b**2)*sqrt(tan(c + d*x))/(3*a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_650():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(tan(c + d*x)))
    F = -atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + 2*b**2*sqrt(tan(c + d*x))/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 4*b**2*(4*a**2 + b**2)*sqrt(tan(c + d*x))/(3*a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_651():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(3)/2))
    F = -I*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + I*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))) - 2*b*(3*a**2 + 4*b**2)*sqrt(tan(c + d*x))/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) - 2*b*(3*a**4 + 17*a**2*b**2 + 8*b**4)*sqrt(tan(c + d*x))/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_652():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2))
    F = atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 4*b/(a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(tan(c + d*x))) + 2*b**2*(7*a**2 + 8*b**2)*sqrt(tan(c + d*x))/(3*a**3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)) + 4*b**2*(4*a**4 + 15*a**2*b**2 + 8*b**4)*sqrt(tan(c + d*x))/(3*a**4*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_653():
    f = 1/(sqrt(3*tan(c + d*x) + 2)*sqrt(tan(c + d*x)))
    F = atanh(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) + 2))/(d*sqrt(3 - 2*I)) + atanh(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) + 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_654():
    f = 1/(sqrt(3*tan(c + d*x) - 2)*sqrt(tan(c + d*x)))
    F = atanh(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) - 2))/(d*sqrt(3 - 2*I)) + atanh(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) - 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_655():
    f = 1/(sqrt(2 - 3*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = atan(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(2 - 3*tan(c + d*x)))/(d*sqrt(3 - 2*I)) + atan(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(2 - 3*tan(c + d*x)))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_656():
    f = 1/(sqrt(-3*tan(c + d*x) - 2)*sqrt(tan(c + d*x)))
    F = atan(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(-3*tan(c + d*x) - 2))/(d*sqrt(3 - 2*I)) + atan(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(-3*tan(c + d*x) - 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_657():
    f = 1/(sqrt(2*tan(c + d*x) + 3)*sqrt(tan(c + d*x)))
    F = atanh(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) + 3))/(d*sqrt(2 - 3*I)) + atanh(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) + 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_658():
    f = 1/(sqrt(3 - 2*tan(c + d*x))*sqrt(tan(c + d*x)))
    F = atan(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(3 - 2*tan(c + d*x)))/(d*sqrt(2 - 3*I)) + atan(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(3 - 2*tan(c + d*x)))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_659():
    f = 1/(sqrt(2*tan(c + d*x) - 3)*sqrt(tan(c + d*x)))
    F = atanh(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) - 3))/(d*sqrt(2 - 3*I)) + atanh(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) - 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_660():
    f = 1/(sqrt(-2*tan(c + d*x) - 3)*sqrt(tan(c + d*x)))
    F = atan(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(-2*tan(c + d*x) - 3))/(d*sqrt(2 - 3*I)) + atan(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(-2*tan(c + d*x) - 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_661():
    f = sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) + 2)
    F = I*atanh(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) + 2))/(d*sqrt(3 - 2*I)) - I*atanh(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) + 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_662():
    f = sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) - 2)
    F = -I*atanh(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) - 2))/(d*sqrt(3 - 2*I)) + I*atanh(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(3*tan(c + d*x) - 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_663():
    f = sqrt(tan(c + d*x))/sqrt(2 - 3*tan(c + d*x))
    F = -I*atan(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(2 - 3*tan(c + d*x)))/(d*sqrt(3 - 2*I)) + I*atan(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(2 - 3*tan(c + d*x)))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_664():
    f = sqrt(tan(c + d*x))/sqrt(-3*tan(c + d*x) - 2)
    F = I*atan(sqrt(3 - 2*I)*sqrt(tan(c + d*x))/sqrt(-3*tan(c + d*x) - 2))/(d*sqrt(3 - 2*I)) - I*atan(sqrt(3 + 2*I)*sqrt(tan(c + d*x))/sqrt(-3*tan(c + d*x) - 2))/(d*sqrt(3 + 2*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_665():
    f = sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) + 3)
    F = I*atanh(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) + 3))/(d*sqrt(2 - 3*I)) - I*atanh(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) + 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_666():
    f = sqrt(tan(c + d*x))/sqrt(3 - 2*tan(c + d*x))
    F = -I*atan(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(3 - 2*tan(c + d*x)))/(d*sqrt(2 - 3*I)) + I*atan(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(3 - 2*tan(c + d*x)))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_667():
    f = sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) - 3)
    F = -I*atanh(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) - 3))/(d*sqrt(2 - 3*I)) + I*atanh(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(2*tan(c + d*x) - 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_668():
    f = sqrt(tan(c + d*x))/sqrt(-2*tan(c + d*x) - 3)
    F = I*atan(sqrt(2 - 3*I)*sqrt(tan(c + d*x))/sqrt(-2*tan(c + d*x) - 3))/(d*sqrt(2 - 3*I)) - I*atan(sqrt(2 + 3*I)*sqrt(tan(c + d*x))/sqrt(-2*tan(c + d*x) - 3))/(d*sqrt(2 + 3*I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_669():
    f = tan(c + d*x)**(sympy.S(5)/3)/(a + b*tan(c + d*x))
    F = -3*a**(sympy.S(5)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(2*b**(sympy.S(2)/3)*d*(a**2 + b**2)) + a**(sympy.S(5)/3)*log(a + b*tan(c + d*x))/(2*b**(sympy.S(2)/3)*d*(a**2 + b**2)) - sqrt(3)*a**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(b**(sympy.S(2)/3)*d*(a**2 + b**2)) + a*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(4*a**2 + 4*b**2)) - a*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(2*a**2 + 2*b**2)) + sqrt(3)*a*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(d*(2*a**2 + 2*b**2)) + sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) - sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) + b*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(d*(2*a**2 + 2*b**2)) + b*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(d*(2*a**2 + 2*b**2)) + b*atan(tan(c + d*x)**(sympy.S(1)/3))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_670():
    f = tan(c + d*x)**(sympy.S(1)/3)/(a + b*tan(c + d*x))
    F = -3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(d*(2*a**2 + 2*b**2)) + a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*log(a + b*tan(c + d*x))/(d*(2*a**2 + 2*b**2)) + sqrt(3)*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(d*(a**2 + b**2)) + a*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(4*a**2 + 4*b**2)) - a*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(2*a**2 + 2*b**2)) - sqrt(3)*a*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(d*(2*a**2 + 2*b**2)) - sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) + sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) + b*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(d*(2*a**2 + 2*b**2)) + b*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(d*(2*a**2 + 2*b**2)) + b*atan(tan(c + d*x)**(sympy.S(1)/3))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_671():
    f = 1/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(1)/3))
    F = -a*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(4*a**2 + 4*b**2)) + a*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(2*a**2 + 2*b**2)) - sqrt(3)*a*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(d*(2*a**2 + 2*b**2)) - sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) + sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) - b*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(d*(2*a**2 + 2*b**2)) - b*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(d*(2*a**2 + 2*b**2)) - b*atan(tan(c + d*x)**(sympy.S(1)/3))/(d*(a**2 + b**2)) - 3*b**(sympy.S(4)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(2*a**(sympy.S(1)/3)*d*(a**2 + b**2)) + b**(sympy.S(4)/3)*log(a + b*tan(c + d*x))/(2*a**(sympy.S(1)/3)*d*(a**2 + b**2)) - sqrt(3)*b**(sympy.S(4)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(1)/3)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_672():
    f = 1/((a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(5)/3))
    F = -a*log(tan(c + d*x)**(sympy.S(4)/3) - tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(4*a**2 + 4*b**2)) + a*log(tan(c + d*x)**(sympy.S(2)/3) + 1)/(d*(2*a**2 + 2*b**2)) + sqrt(3)*a*atan(sqrt(3)*(1 - 2*tan(c + d*x)**(sympy.S(2)/3))/3)/(d*(2*a**2 + 2*b**2)) - 3*a/(d*(2*a**2 + 2*b**2)*tan(c + d*x)**(sympy.S(2)/3)) + sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) - sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) - sqrt(3)*b*log(tan(c + d*x)**(sympy.S(2)/3) + sqrt(3)*tan(c + d*x)**(sympy.S(1)/3) + 1)/(d*(4*a**2 + 4*b**2)) - b*atan(2*tan(c + d*x)**(sympy.S(1)/3) - sqrt(3))/(d*(2*a**2 + 2*b**2)) - b*atan(2*tan(c + d*x)**(sympy.S(1)/3) + sqrt(3))/(d*(2*a**2 + 2*b**2)) - b*atan(tan(c + d*x)**(sympy.S(1)/3))/(d*(a**2 + b**2)) - 3*b**2/(2*a*d*(a**2 + b**2)*tan(c + d*x)**(sympy.S(2)/3)) - 3*b**(sympy.S(8)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(2*a**(sympy.S(5)/3)*d*(a**2 + b**2)) + b**(sympy.S(8)/3)*log(a + b*tan(c + d*x))/(2*a**(sympy.S(5)/3)*d*(a**2 + b**2)) + sqrt(3)*b**(sympy.S(8)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tan(c + d*x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(a**(sympy.S(5)/3)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_673():
    f = tan(c + d*x)**(sympy.S(4)/3)/sqrt(a + b*tan(c + d*x))
    F = 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, sympy.S.Half, 1, sympy.S(10)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(14*d*sqrt(a + b*tan(c + d*x))) + 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(7)/3)*appellf1(sympy.S(7)/3, sympy.S.Half, 1, sympy.S(10)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(14*d*sqrt(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_674():
    f = tan(c + d*x)**(sympy.S(2)/3)/sqrt(a + b*tan(c + d*x))
    F = 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, sympy.S.Half, 1, sympy.S(8)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(10*d*sqrt(a + b*tan(c + d*x))) + 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(5)/3)*appellf1(sympy.S(5)/3, sympy.S.Half, 1, sympy.S(8)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(10*d*sqrt(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_675():
    f = tan(c + d*x)**(sympy.S(1)/3)/sqrt(a + b*tan(c + d*x))
    F = 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(8*d*sqrt(a + b*tan(c + d*x))) + 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(4)/3)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(8*d*sqrt(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_676():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(1)/3))
    F = 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(4*d*sqrt(a + b*tan(c + d*x))) + 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(4*d*sqrt(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_677():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(2)/3))
    F = 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))) + 3*sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_678():
    f = 1/(sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(4)/3))
    F = -3*sqrt(1 + b*tan(c + d*x)/a)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(1)/3)) - 3*sqrt(1 + b*tan(c + d*x)/a)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_679():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**4
    F = -9*c*(c + d*tan(e + f*x))**(sympy.S(4)/3)*tan(e + f*x)/(35*d**2*f) - x*(c - sqrt(-d**2))**(sympy.S(1)/3)/4 - x*(c + sqrt(-d**2))**(sympy.S(1)/3)/4 + 3*sqrt(-d**2)*(c - sqrt(-d**2))**(sympy.S(1)/3)*log((c - sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*d*f) + sqrt(-d**2)*(c - sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*d*f) - sqrt(3)*sqrt(-d**2)*(c - sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*d*f) - 3*sqrt(-d**2)*(c + sqrt(-d**2))**(sympy.S(1)/3)*log((c + sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*d*f) - sqrt(-d**2)*(c + sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*d*f) + sqrt(3)*sqrt(-d**2)*(c + sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*d*f) + 3*(c + d*tan(e + f*x))**(sympy.S(4)/3)*tan(e + f*x)**2/(10*d*f) + (c + d*tan(e + f*x))**(sympy.S(4)/3)*(27*c**2 - 105*d**2)/(140*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_680():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**3
    F = -9*c*(c + d*tan(e + f*x))**(sympy.S(4)/3)/(28*d**2*f) - I*x*(c - I*d)**(sympy.S(1)/3)/4 + I*x*(c + I*d)**(sympy.S(1)/3)/4 - 3*(c - I*d)**(sympy.S(1)/3)*log((c - I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c - I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c - I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - I*d)**(sympy.S(1)/3))/3)/(2*f) - 3*(c + I*d)**(sympy.S(1)/3)*log((c + I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c + I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c + I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + I*d)**(sympy.S(1)/3))/3)/(2*f) - 3*(c + d*tan(e + f*x))**(sympy.S(1)/3)/f + 3*(c + d*tan(e + f*x))**(sympy.S(4)/3)*tan(e + f*x)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_681():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**2
    F = 3*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log((c - sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) + d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) - sqrt(3)*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) - 3*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log((c + sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) - d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) + sqrt(3)*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) + x*(c - sqrt(-d**2))**(sympy.S(1)/3)/4 + x*(c + sqrt(-d**2))**(sympy.S(1)/3)/4 + 3*(c + d*tan(e + f*x))**(sympy.S(4)/3)/(4*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_682():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)
    F = I*x*(c - I*d)**(sympy.S(1)/3)/4 - I*x*(c + I*d)**(sympy.S(1)/3)/4 + 3*(c - I*d)**(sympy.S(1)/3)*log((c - I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) + (c - I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) - sqrt(3)*(c - I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - I*d)**(sympy.S(1)/3))/3)/(2*f) + 3*(c + I*d)**(sympy.S(1)/3)*log((c + I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) + (c + I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) - sqrt(3)*(c + I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + I*d)**(sympy.S(1)/3))/3)/(2*f) + 3*(c + d*tan(e + f*x))**(sympy.S(1)/3)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_683():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)
    F = -3*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log((c - sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) - d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) + sqrt(3)*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) + 3*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log((c + sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) + d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) - sqrt(3)*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) - x*(c - sqrt(-d**2))**(sympy.S(1)/3)/4 - x*(c + sqrt(-d**2))**(sympy.S(1)/3)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_684():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*cot(e + f*x)
    F = 3*c**(sympy.S(1)/3)*log(c**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(2*f) - c**(sympy.S(1)/3)*log(tan(e + f*x))/(2*f) - sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3))/(3*c**(sympy.S(1)/3)))/f - I*x*(c - I*d)**(sympy.S(1)/3)/4 + I*x*(c + I*d)**(sympy.S(1)/3)/4 - 3*(c - I*d)**(sympy.S(1)/3)*log((c - I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c - I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c - I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - I*d)**(sympy.S(1)/3))/3)/(2*f) - 3*(c + I*d)**(sympy.S(1)/3)*log((c + I*d)**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f) - (c + I*d)**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f) + sqrt(3)*(c + I*d)**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + I*d)**(sympy.S(1)/3))/3)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_685():
    f = (c + d*tan(e + f*x))**(sympy.S(1)/3)*cot(e + f*x)**2
    F = 3*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log((c - sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) + d*(c - sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) - sqrt(3)*d*(c - sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c - sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) - 3*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log((c + sqrt(-d**2))**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(4*f*sqrt(-d**2)) - d*(c + sqrt(-d**2))**(sympy.S(1)/3)*log(cos(e + f*x))/(4*f*sqrt(-d**2)) + sqrt(3)*d*(c + sqrt(-d**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3)/(c + sqrt(-d**2))**(sympy.S(1)/3))/3)/(2*f*sqrt(-d**2)) + x*(c - sqrt(-d**2))**(sympy.S(1)/3)/4 + x*(c + sqrt(-d**2))**(sympy.S(1)/3)/4 - (c + d*tan(e + f*x))**(sympy.S(1)/3)*cot(e + f*x)/f + d*log(c**(sympy.S(1)/3) - (c + d*tan(e + f*x))**(sympy.S(1)/3))/(2*c**(sympy.S(2)/3)*f) - d*log(tan(e + f*x))/(6*c**(sympy.S(2)/3)*f) - sqrt(3)*d*atan(sqrt(3)*(c**(sympy.S(1)/3) + 2*(c + d*tan(e + f*x))**(sympy.S(1)/3))/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_686():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/3)
    F = 3*b*(a + b*tan(c + d*x))**(sympy.S(2)/3)/(2*d) - x*(a - I*b)**(sympy.S(5)/3)/4 - x*(a + I*b)**(sympy.S(5)/3)/4 + 3*I*(a - I*b)**(sympy.S(5)/3)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) + I*(a - I*b)**(sympy.S(5)/3)*log(cos(c + d*x))/(4*d) + sqrt(3)*I*(a - I*b)**(sympy.S(5)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d) - 3*I*(a + I*b)**(sympy.S(5)/3)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) - I*(a + I*b)**(sympy.S(5)/3)*log(cos(c + d*x))/(4*d) - sqrt(3)*I*(a + I*b)**(sympy.S(5)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_687():
    f = (a + b*tan(c + d*x))**(sympy.S(4)/3)
    F = 3*b*(a + b*tan(c + d*x))**(sympy.S(1)/3)/d - x*(a - I*b)**(sympy.S(4)/3)/4 - x*(a + I*b)**(sympy.S(4)/3)/4 + 3*I*(a - I*b)**(sympy.S(4)/3)*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) + I*(a - I*b)**(sympy.S(4)/3)*log(cos(c + d*x))/(4*d) - sqrt(3)*I*(a - I*b)**(sympy.S(4)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d) - 3*I*(a + I*b)**(sympy.S(4)/3)*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d) - I*(a + I*b)**(sympy.S(4)/3)*log(cos(c + d*x))/(4*d) + sqrt(3)*I*(a + I*b)**(sympy.S(4)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_688():
    f = (a + b*tan(c + d*x))**(sympy.S(2)/3)
    F = -3*b*(a - sqrt(-b**2))**(sympy.S(2)/3)*log((a - sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)) - b*(a - sqrt(-b**2))**(sympy.S(2)/3)*log(cos(c + d*x))/(4*d*sqrt(-b**2)) - sqrt(3)*b*(a - sqrt(-b**2))**(sympy.S(2)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)) + 3*b*(a + sqrt(-b**2))**(sympy.S(2)/3)*log((a + sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)) + b*(a + sqrt(-b**2))**(sympy.S(2)/3)*log(cos(c + d*x))/(4*d*sqrt(-b**2)) + sqrt(3)*b*(a + sqrt(-b**2))**(sympy.S(2)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)) - x*(a - sqrt(-b**2))**(sympy.S(2)/3)/4 - x*(a + sqrt(-b**2))**(sympy.S(2)/3)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_689():
    f = (a + b*tan(c + d*x))**(sympy.S(1)/3)
    F = -3*b*(a - sqrt(-b**2))**(sympy.S(1)/3)*log((a - sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)) - b*(a - sqrt(-b**2))**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d*sqrt(-b**2)) + sqrt(3)*b*(a - sqrt(-b**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)) + 3*b*(a + sqrt(-b**2))**(sympy.S(1)/3)*log((a + sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)) + b*(a + sqrt(-b**2))**(sympy.S(1)/3)*log(cos(c + d*x))/(4*d*sqrt(-b**2)) - sqrt(3)*b*(a + sqrt(-b**2))**(sympy.S(1)/3)*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)) - x*(a - sqrt(-b**2))**(sympy.S(1)/3)/4 - x*(a + sqrt(-b**2))**(sympy.S(1)/3)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_690():
    f = (a + b*tan(c + d*x))**(sympy.S(-1)/3)
    F = 3*b*log((a + sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(1)/3)) + b*log(cos(c + d*x))/(4*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(1)/3)) + sqrt(3)*b*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(1)/3)) - 3*b*log((a - sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(1)/3)) - b*log(cos(c + d*x))/(4*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(1)/3)) - sqrt(3)*b*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(1)/3)) - x/(4*(a + sqrt(-b**2))**(sympy.S(1)/3)) - x/(4*(a - sqrt(-b**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_691():
    f = (a + b*tan(c + d*x))**(sympy.S(-2)/3)
    F = 3*b*log((a + sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(2)/3)) + b*log(cos(c + d*x))/(4*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(2)/3)) - sqrt(3)*b*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))**(sympy.S(2)/3)) - 3*b*log((a - sqrt(-b**2))**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(2)/3)) - b*log(cos(c + d*x))/(4*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(2)/3)) + sqrt(3)*b*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - sqrt(-b**2))**(sympy.S(1)/3))/3)/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))**(sympy.S(2)/3)) - x/(4*(a + sqrt(-b**2))**(sympy.S(2)/3)) - x/(4*(a - sqrt(-b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_692():
    f = (a + b*tan(c + d*x))**(sympy.S(-4)/3)
    F = -3*b/(d*(a + b*tan(c + d*x))**(sympy.S(1)/3)*(a**2 + b**2)) - x/(4*(a + I*b)**(sympy.S(4)/3)) - x/(4*(a - I*b)**(sympy.S(4)/3)) - 3*I*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a + I*b)**(sympy.S(4)/3)) - I*log(cos(c + d*x))/(4*d*(a + I*b)**(sympy.S(4)/3)) - sqrt(3)*I*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d*(a + I*b)**(sympy.S(4)/3)) + 3*I*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a - I*b)**(sympy.S(4)/3)) + I*log(cos(c + d*x))/(4*d*(a - I*b)**(sympy.S(4)/3)) + sqrt(3)*I*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d*(a - I*b)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_693():
    f = (a + b*tan(c + d*x))**(sympy.S(-5)/3)
    F = -3*b/(d*(a + b*tan(c + d*x))**(sympy.S(2)/3)*(2*a**2 + 2*b**2)) - x/(4*(a + I*b)**(sympy.S(5)/3)) - x/(4*(a - I*b)**(sympy.S(5)/3)) - 3*I*log((a + I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a + I*b)**(sympy.S(5)/3)) - I*log(cos(c + d*x))/(4*d*(a + I*b)**(sympy.S(5)/3)) + sqrt(3)*I*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a + I*b)**(sympy.S(1)/3))/3)/(2*d*(a + I*b)**(sympy.S(5)/3)) + 3*I*log((a - I*b)**(sympy.S(1)/3) - (a + b*tan(c + d*x))**(sympy.S(1)/3))/(4*d*(a - I*b)**(sympy.S(5)/3)) + I*log(cos(c + d*x))/(4*d*(a - I*b)**(sympy.S(5)/3)) - sqrt(3)*I*atan(sqrt(3)*(1 + 2*(a + b*tan(c + d*x))**(sympy.S(1)/3)/(a - I*b)**(sympy.S(1)/3))/3)/(2*d*(a - I*b)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_694():
    f = (d*tan(e + f*x))**n*(a + b*tan(e + f*x))**4
    F = 2*a*b**3*(d*tan(e + f*x))**(n + 1)*(n + 4)*tan(e + f*x)/(d*f*(n + 2)*(n + 3)) + 4*a*b*(d*tan(e + f*x))**(n + 2)*(a**2 - b**2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(n + 2)) + b**2*(d*tan(e + f*x))**(n + 1)*(a + b*tan(e + f*x))**2/(d*f*(n + 3)) - b**2*(d*tan(e + f*x))**(n + 1)*(-a**2*(5*n + 17) + b**2*(n + 3))/(d*f*(n + 1)*(n + 3)) + (d*tan(e + f*x))**(n + 1)*(a**4 - 6*a**2*b**2 + b**4)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_695():
    f = (d*tan(e + f*x))**n*(a + b*tan(e + f*x))**3
    F = a*b**2*(d*tan(e + f*x))**(n + 1)*(2*n + 5)/(d*f*(n + 1)*(n + 2)) + a*(d*tan(e + f*x))**(n + 1)*(a**2 - 3*b**2)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1)) + b**2*(d*tan(e + f*x))**(n + 1)*(a + b*tan(e + f*x))/(d*f*(n + 2)) + b*(d*tan(e + f*x))**(n + 2)*(3*a**2 - b**2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_696():
    f = (d*tan(e + f*x))**n*(a + b*tan(e + f*x))**2
    F = 2*a*b*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(n + 2)) + b**2*(d*tan(e + f*x))**(n + 1)/(d*f*(n + 1)) + (d*tan(e + f*x))**(n + 1)*(a**2 - b**2)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_697():
    f = (d*tan(e + f*x))**n*(a + b*tan(e + f*x))
    F = a*(d*tan(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1)) + b*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_698():
    f = (d*tan(e + f*x))**n/(a + b*tan(e + f*x))
    F = a*(d*tan(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(a**2 + b**2)*(n + 1)) - b*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(a**2 + b**2)*(n + 2)) + b**2*(d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), -b*tan(e + f*x)/a)/(a*d*f*(a**2 + b**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_699():
    f = (d*tan(e + f*x))**n/(a + b*tan(e + f*x))**2
    F = -2*a*b*(d*tan(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -tan(e + f*x)**2)/(d**2*f*(a**2 + b**2)**2*(n + 2)) + (d*tan(e + f*x))**(n + 1)*(a**2 - b**2)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(a**2 + b**2)**2*(n + 1)) + b**2*(d*tan(e + f*x))**(n + 1)/(a*d*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + b**2*(d*tan(e + f*x))**(n + 1)*(a**2*(2 - n) - b**2*n)*hyper((1, n + 1), (n + 2,), -b*tan(e + f*x)/a)/(a**2*d*f*(a**2 + b**2)**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_700():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**m
    F = a*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1)) + a*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_701():
    f = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**m
    F = sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1)) + sqrt(a + b*tan(c + d*x))*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(1 + b*tan(c + d*x)/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_702():
    f = tan(c + d*x)**m/sqrt(a + b*tan(c + d*x))
    F = sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -b*tan(c + d*x)/a, -I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*(m + 1)) + sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -b*tan(c + d*x)/a, I*tan(c + d*x))/(2*d*sqrt(a + b*tan(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_703():
    f = tan(c + d*x)**m/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a*d*sqrt(a + b*tan(c + d*x))*(m + 1)) + sqrt(1 + b*tan(c + d*x)/a)*tan(c + d*x)**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(2*a*d*sqrt(a + b*tan(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_704():
    f = (d*tan(e + f*x))**n*(a + b*tan(e + f*x))**m
    F = (d*tan(e + f*x))**(n + 1)*(a + b*tan(e + f*x))**m*appellf1(n + 1, 1, -m, n + 2, -I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*d*f*(1 + b*tan(e + f*x)/a)**m*(n + 1)) + (d*tan(e + f*x))**(n + 1)*(a + b*tan(e + f*x))**m*appellf1(n + 1, 1, -m, n + 2, I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*d*f*(1 + b*tan(e + f*x)/a)**m*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_705():
    f = (a + b*tan(c + d*x))**n*tan(c + d*x)**4
    F = -2*a*(a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)/(b**2*d*(n + 2)*(n + 3)) + sqrt(-b**2)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(2*b*d*(a + sqrt(-b**2))*(n + 1)) - sqrt(-b**2)*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(2*b*d*(a - sqrt(-b**2))*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)**2/(b*d*(n + 3)) + (a + b*tan(c + d*x))**(n + 1)*(2*a**2 - b**2*(n + 2)*(n + 3))/(b**3*d*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_706():
    f = (a + b*tan(c + d*x))**n*tan(c + d*x)**3
    F = -a*(a + b*tan(c + d*x))**(n + 1)/(b**2*d*(n + 1)*(n + 2)) + (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*tan(c + d*x)/(b*d*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_707():
    f = (a + b*tan(c + d*x))**n*tan(c + d*x)**2
    F = b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))*(n + 1)) - b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_708():
    f = (a + b*tan(c + d*x))**n*tan(c + d*x)
    F = -(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_709():
    f = (a + b*tan(c + d*x))**n
    F = -b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))*(n + 1)) + b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_710():
    f = (a + b*tan(c + d*x))**n*cot(c + d*x)
    F = (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_711():
    f = (a + b*tan(c + d*x))**n*cot(c + d*x)**2
    F = b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a + sqrt(-b**2))*(n + 1)) - b*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(2*d*sqrt(-b**2)*(a - sqrt(-b**2))*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)/(a*d) - b*n*(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(a**2*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_712():
    f = (a + b*tan(c + d*x))**n*cot(c + d*x)**3
    F = -(a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + I*b))/(d*(2*a + 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - I*b))/(d*(2*a - 2*I*b)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)**2/(2*a*d) + b*(1 - n)*(a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)/(2*a**2*d) + (a + b*tan(c + d*x))**(n + 1)*(2*a**2 + b**2*n*(1 - n))*hyper((1, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(2*a**3*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_713():
    f = (a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)
    F = (a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n) + (a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_714():
    f = (a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))
    F = (a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n) + (a + b*tan(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_715():
    f = (a + b*tan(c + d*x))**n/sqrt(tan(c + d*x))
    F = (a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n) + (a + b*tan(c + d*x))**n*sqrt(tan(c + d*x))*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_716():
    f = (a + b*tan(c + d*x))**n/tan(c + d*x)**(sympy.S(3)/2)
    F = -(a + b*tan(c + d*x))**n*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(tan(c + d*x))) - (a + b*tan(c + d*x))**n*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_717():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*I*a*sqrt(cot(c + d*x))/d + 2*(-1)**(sympy.S(1)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_718():
    f = (I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(cot(c + d*x))/d - 2*(-1)**(sympy.S(3)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_719():
    f = (I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))
    F = -2*(-1)**(sympy.S(1)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_720():
    f = (I*a*tan(c + d*x) + a)/sqrt(cot(c + d*x))
    F = 2*(-1)**(sympy.S(3)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 2*I*a/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_721():
    f = (I*a*tan(c + d*x) + a)/cot(c + d*x)**(sympy.S(3)/2)
    F = 2*(-1)**(sympy.S(1)/4)*a*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 2*a/(d*sqrt(cot(c + d*x))) + 2*I*a/(3*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_722():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*a**2*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 4*I*a**2*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a**2*sqrt(cot(c + d*x))/d + 4*(-1)**(sympy.S(3)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_723():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a**2*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 4*I*a**2*sqrt(cot(c + d*x))/d + 4*(-1)**(sympy.S(1)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_724():
    f = (I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2*sqrt(cot(c + d*x))/d - 4*(-1)**(sympy.S(3)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_725():
    f = (I*a*tan(c + d*x) + a)**2*sqrt(cot(c + d*x))
    F = -4*(-1)**(sympy.S(1)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*a**2/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_726():
    f = (I*a*tan(c + d*x) + a)**2/sqrt(cot(c + d*x))
    F = 4*(-1)**(sympy.S(3)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 4*I*a**2/(d*sqrt(cot(c + d*x))) - 2*a**2/(3*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_727():
    f = (I*a*tan(c + d*x) + a)**2/cot(c + d*x)**(sympy.S(3)/2)
    F = 4*(-1)**(sympy.S(1)/4)*a**2*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 4*a**2/(d*sqrt(cot(c + d*x))) + 4*I*a**2/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - 2*a**2/(5*d*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_728():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(7)/2)
    F = -8*I*a**3*cot(c + d*x)**(sympy.S(3)/2)/(5*d) + 8*a**3*sqrt(cot(c + d*x))/d + 8*(-1)**(sympy.S(3)/4)*a**3*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*(a**3*cot(c + d*x) + I*a**3)*cot(c + d*x)**(sympy.S(3)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_729():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(5)/2)
    F = -16*I*a**3*sqrt(cot(c + d*x))/(3*d) + 8*(-1)**(sympy.S(1)/4)*a**3*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 2*(a**3*cot(c + d*x) + I*a**3)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_730():
    f = (I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(3)/2)
    F = -8*(-1)**(sympy.S(3)/4)*a**3*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - (2*a**3*cot(c + d*x) + 2*I*a**3)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_731():
    f = (I*a*tan(c + d*x) + a)**3*sqrt(cot(c + d*x))
    F = -8*(-1)**(sympy.S(1)/4)*a**3*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d - 16*a**3/(3*d*sqrt(cot(c + d*x))) - (2*a**3*cot(c + d*x) + 2*I*a**3)/(3*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_732():
    f = (I*a*tan(c + d*x) + a)**3/sqrt(cot(c + d*x))
    F = 8*(-1)**(sympy.S(3)/4)*a**3*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/d + 8*I*a**3/(d*sqrt(cot(c + d*x))) - 8*a**3/(5*d*cot(c + d*x)**(sympy.S(3)/2)) - (2*a**3*cot(c + d*x) + 2*I*a**3)/(5*d*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_733():
    f = cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cot(c + d*x) + I*a)) - sqrt(2)*(sympy.S(5)/8 - 3*I/8)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(sympy.S(5)/8 - 3*I/8)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - 5*sqrt(cot(c + d*x))/(2*a*d) + sqrt(2)*(sympy.S(5)/4 + 3*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) + sqrt(2)*(sympy.S(5)/4 + 3*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_734():
    f = sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)
    F = sqrt(cot(c + d*x))/(2*d*(a*cot(c + d*x) + I*a)) - sqrt(2)*(sympy.S(3)/8 + I/8)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(sympy.S(3)/8 + I/8)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(3)/4 - I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(3)/4 - I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_735():
    f = 1/((I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x)))
    F = I*sqrt(cot(c + d*x))/(2*d*(a*cot(c + d*x) + I*a)) + (-1)**(sympy.S(3)/4)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_736():
    f = 1/((I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(cot(c + d*x))/(2*d*(a*cot(c + d*x) + I*a)) - sqrt(2)*(sympy.S(1)/8 + 3*I/8)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) + sqrt(2)*(sympy.S(1)/8 + 3*I/8)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/4 - 3*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(1)/4 - 3*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_737():
    f = 1/((I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2))
    F = -1/(2*d*(a*cot(c + d*x) + I*a)*sqrt(cot(c + d*x))) + sqrt(2)*(sympy.S(3)/8 - 5*I/8)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(3)/8 - 5*I/8)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a*d) - sqrt(2)*(sympy.S(3)/4 + 5*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(sympy.S(3)/4 + 5*I/4)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a*d) - 5*I/(2*a*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_738():
    f = cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**2
    F = cot(c + d*x)**(sympy.S(5)/2)/(4*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*(sympy.S(25)/32 - 21*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(25)/32 - 21*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - 25*sqrt(cot(c + d*x))/(8*a**2*d) + sqrt(2)*(sympy.S(25)/16 + 21*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) + sqrt(2)*(sympy.S(25)/16 + 21*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) + 7*cot(c + d*x)**(sympy.S(3)/2)/(8*a**2*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_739():
    f = sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**2
    F = cot(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*(sympy.S(9)/32 + 5*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(9)/32 + 5*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(9)/16 - 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(9)/16 - 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) + 5*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_740():
    f = 1/((I*a*tan(c + d*x) + a)**2*sqrt(cot(c + d*x)))
    F = sqrt(cot(c + d*x))/(4*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*(sympy.S(1)/32 + 3*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/32 + 3*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/16 - 3*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/16 - 3*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) + 3*I*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_741():
    f = 1/((I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(3)/2))
    F = I*sqrt(cot(c + d*x))/(4*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*(sympy.S(1)/32 - 3*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(1)/32 - 3*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/16 + 3*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) + sqrt(2)*(sympy.S(1)/16 + 3*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) + sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_742():
    f = 1/((I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(cot(c + d*x))/(4*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*(sympy.S(9)/32 - 5*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(9)/32 - 5*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(9)/16 + 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) + sqrt(2)*(sympy.S(9)/16 + 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) + 5*I*sqrt(cot(c + d*x))/(8*a**2*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_743():
    f = 1/((I*a*tan(c + d*x) + a)**2*cot(c + d*x)**(sympy.S(7)/2))
    F = -1/(4*d*(a*cot(c + d*x) + I*a)**2*sqrt(cot(c + d*x))) - sqrt(2)*(sympy.S(25)/32 + 21*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) + sqrt(2)*(sympy.S(25)/32 + 21*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**2*d) - sqrt(2)*(sympy.S(25)/16 - 21*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**2*d) - sqrt(2)*(sympy.S(25)/16 - 21*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**2*d) - 25/(8*a**2*d*sqrt(cot(c + d*x))) + 7*I/(8*a**2*d*(cot(c + d*x) + I)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_744():
    f = sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**3
    F = 5*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + cot(c + d*x)**(sympy.S(5)/2)/(6*d*(a*cot(c + d*x) + I*a)**3) + cot(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*(sympy.S(7)/32 + 5*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) + sqrt(2)*(sympy.S(7)/32 + 5*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(7)/16 - 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**3*d) - sqrt(2)*(sympy.S(7)/16 - 5*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_745():
    f = 1/((I*a*tan(c + d*x) + a)**3*sqrt(cot(c + d*x)))
    F = I*sqrt(cot(c + d*x))/(4*d*(a**3*cot(c + d*x) + I*a**3)) + cot(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cot(c + d*x) + I*a)**3) + sqrt(cot(c + d*x))/(4*a*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*I*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(32*a**3*d) - sqrt(2)*I*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(32*a**3*d) + sqrt(2)*I*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(16*a**3*d) + sqrt(2)*I*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(16*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_746():
    f = 1/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(3)/2))
    F = sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) + sqrt(cot(c + d*x))/(6*d*(a*cot(c + d*x) + I*a)**3) + I*sqrt(cot(c + d*x))/(6*a*d*(a*cot(c + d*x) + I*a)**2) + (-1)**(sympy.S(1)/4)*atan((-1)**(sympy.S(1)/4)*sqrt(cot(c + d*x)))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_747():
    f = 1/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(5)/2))
    F = I*sqrt(cot(c + d*x))/(6*d*(a*cot(c + d*x) + I*a)**3) + sqrt(cot(c + d*x))/(12*a*d*(a*cot(c + d*x) + I*a)**2) - sqrt(2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(32*a**3*d) + sqrt(2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(32*a**3*d) + sqrt(2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(16*a**3*d) + sqrt(2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(16*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_748():
    f = 1/((I*a*tan(c + d*x) + a)**3*cot(c + d*x)**(sympy.S(7)/2))
    F = 5*sqrt(cot(c + d*x))/(8*d*(a**3*cot(c + d*x) + I*a**3)) - sqrt(cot(c + d*x))/(6*d*(a*cot(c + d*x) + I*a)**3) + I*sqrt(cot(c + d*x))/(3*a*d*(a*cot(c + d*x) + I*a)**2) + sqrt(2)*(sympy.S(5)/32 + 7*I/32)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) - sqrt(2)*(sympy.S(5)/32 + 7*I/32)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(2*a**3*d) + sqrt(2)*(sympy.S(5)/16 - 7*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*a**3*d) + sqrt(2)*(sympy.S(5)/16 - 7*I/16)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_749():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)
    F = sqrt(a)*(-1 - I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 2*I*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + 26*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_750():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)
    F = sqrt(a)*(-1 + I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_751():
    f = sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a)*(1 + I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_752():
    f = sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))
    F = sqrt(a)*(1 - I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_753():
    f = sqrt(I*a*tan(c + d*x) + a)/sqrt(cot(c + d*x))
    F = -2*(-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - sqrt(a)*(1 + I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_754():
    f = sqrt(I*a*tan(c + d*x) + a)/cot(c + d*x)**(sympy.S(3)/2)
    F = -(-1)**(sympy.S(1)/4)*sqrt(a)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - sqrt(a)*(1 - I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_755():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = a**(sympy.S(3)/2)*(-2 - 2*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*cot(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(I*a*tan(c + d*x) + a)) - 2*I*a**2*cot(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(I*a*tan(c + d*x) + a)) - 4*I*a*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(5*d) + 12*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_756():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = a**(sympy.S(3)/2)*(-2 + 2*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*I*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d - 2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_757():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*(2 + 2*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_758():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))
    F = 2*(-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(3)/2)*(2 - 2*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_759():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cot(c + d*x))
    F = -3*(-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**(sympy.S(3)/2)*(2 + 2*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + I*a**2/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - a**2/(d*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_760():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = a**(sympy.S(5)/2)*(4 - 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - 6*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2)/(7*d) + 32*a**2*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(21*d) + 104*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_761():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = a**(sympy.S(5)/2)*(-4 - 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 4*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d - 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_762():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = a**(sympy.S(5)/2)*(-4 + 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 4*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d - 2*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_763():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = 2*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 + 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - 2*a**2*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_764():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))
    F = 5*(-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + a**(sympy.S(5)/2)*(4 - 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d - a**2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_765():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cot(c + d*x))
    F = -23*(-1)**(sympy.S(3)/4)*a**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(4*d) - a**(sympy.S(5)/2)*(4 + 4*I)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/d + 9*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(4*d*sqrt(cot(c + d*x))) - a**2*sqrt(I*a*tan(c + d*x) + a)/(2*d*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_766():
    f = cot(c + d*x)**(sympy.S(5)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = cot(c + d*x)**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) + 7*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(3*a*d) - (sympy.S.Half - I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_767():
    f = cot(c + d*x)**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = sqrt(cot(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(a*d) + (sympy.S.Half + I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_768():
    f = sqrt(cot(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = 1/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S.Half - I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_769():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x)))
    F = I/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(-1)/2 - I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_770():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2))
    F = -1/(d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - 2*(-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) - (sympy.S.Half - I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_771():
    f = 1/(sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(5)/2))
    F = -1/(d*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(a*d*sqrt(cot(c + d*x))) - (-1)**(sympy.S(3)/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d) + (sympy.S.Half + I/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_772():
    f = cot(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = cot(c + d*x)**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*cot(c + d*x)**(sympy.S(3)/2)/(2*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(2*a**2*d) + 13*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(2*a**2*d) + (sympy.S(-1)/4 + I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_773():
    f = cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = sqrt(cot(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 11*sqrt(cot(c + d*x))/(6*a*d*sqrt(I*a*tan(c + d*x) + a)) - 25*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(6*a**2*d) + (sympy.S(1)/4 + I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_774():
    f = sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + 7/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/4 - I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_775():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x)))
    F = 1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + I/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) - (sympy.S(1)/4 + I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_776():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = I/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + 1/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(-1)/4 + I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_777():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = -1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + 3*I/(2*a*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + 2*(-1)**(sympy.S(3)/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (sympy.S(1)/4 + I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_778():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2))
    F = -1/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)) + 13*I/(6*a*d*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)) - 7*sqrt(I*a*tan(c + d*x) + a)/(2*a**2*d*sqrt(cot(c + d*x))) - 3*(-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (sympy.S(1)/4 - I/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_779():
    f = cot(c + d*x)**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = cot(c + d*x)**(sympy.S(3)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 7*cot(c + d*x)**(sympy.S(3)/2)/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 89*cot(c + d*x)**(sympy.S(3)/2)/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 361*sqrt(I*a*tan(c + d*x) + a)*cot(c + d*x)**(sympy.S(3)/2)/(60*a**3*d) + 707*I*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(60*a**3*d) + (sympy.S(-1)/8 + I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_780():
    f = cot(c + d*x)**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = sqrt(cot(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 17*sqrt(cot(c + d*x))/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 151*sqrt(cot(c + d*x))/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 317*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))/(60*a**3*d) + (sympy.S(1)/8 + I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_781():
    f = sqrt(cot(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))) + 13/(30*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) + 67/(60*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/8 - I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_782():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x)))
    F = I/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cot(c + d*x))) + I/(10*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cot(c + d*x))) - I/(20*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(-1)/8 - I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_783():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)) + I/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + 1/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(-1)/8 + I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_784():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = I/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)) + 1/(6*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) - I/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + (sympy.S(1)/8 + I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_785():
    f = 1/((I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(7)/2))
    F = -1/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)) + I/(2*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)) + 7/(4*a**2*d*sqrt(I*a*tan(c + d*x) + a)*sqrt(cot(c + d*x))) + 2*(-1)**(sympy.S(1)/4)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan((-1)**(sympy.S(3)/4)*sqrt(a)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (sympy.S(1)/8 - I/8)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(a)*(1 + I)*sqrt(tan(c + d*x))/sqrt(I*a*tan(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_786():
    f = (d*cot(e + f*x))**n*(I*a*tan(e + f*x) + a)**3
    F = I*a**3*d**2*(d*cot(e + f*x))**(n - 2)*(1 - 2*n)/(f*(1 - n)*(2 - n)) - 4*I*a**3*d**2*(d*cot(e + f*x))**(n - 2)*hyper((1, n - 2), (n - 1,), -I*cot(e + f*x))/(f*(2 - n)) + d**2*(d*cot(e + f*x))**(n - 2)*(a**3*cot(e + f*x) + I*a**3)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_787():
    f = (d*cot(e + f*x))**n*(I*a*tan(e + f*x) + a)**2
    F = -2*a**2*d*(d*cot(e + f*x))**(n - 1)*hyper((1, n - 1), (n,), -I*cot(e + f*x))/(f*(1 - n)) + a**2*d*(d*cot(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_788():
    f = (d*cot(e + f*x))**n*(I*a*tan(e + f*x) + a)
    F = -I*a*(d*cot(e + f*x))**n*hyper((1, n), (n + 1,), -I*cot(e + f*x))/(f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_789():
    f = (d*cot(e + f*x))**n/(I*a*tan(e + f*x) + a)
    F = -(d*cot(e + f*x))**(n + 2)/(2*d**2*f*(a*cot(e + f*x) + I*a)) - I*n*(d*cot(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -cot(e + f*x)**2)/(2*a*d**2*f*(n + 2)) + (d*cot(e + f*x))**(n + 3)*(n + 1)*hyper((1, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -cot(e + f*x)**2)/(2*a*d**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_790():
    f = (d*cot(e + f*x))**n/(I*a*tan(e + f*x) + a)**2
    F = -(d*cot(e + f*x))**(n + 3)/(4*d**3*f*(a*cot(e + f*x) + I*a)**2) - I*n*(d*cot(e + f*x))**(n + 3)/(4*a**2*d**3*f*(cot(e + f*x) + I)) + (d*cot(e + f*x))**(n + 3)*(n + 1)**2*hyper((1, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -cot(e + f*x)**2)/(4*a**2*d**3*f*(n + 3)) + I*n*(d*cot(e + f*x))**(n + 4)*(n + 2)*hyper((1, n/2 + 2), (n/2 + 3,), -cot(e + f*x)**2)/(4*a**2*d**4*f*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_791():
    f = (d*cot(e + f*x))**n*(I*a*tan(e + f*x) + a)**m
    F = (d*cot(e + f*x))**n*(I*a*tan(e + f*x) + a)**m*tan(e + f*x)*appellf1(1 - n, 1, 1 - m, 2 - n, I*tan(e + f*x), -I*tan(e + f*x))/(f*(1 - n)*(I*tan(e + f*x) + 1)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_792():
    f = (I*a*tan(c + d*x) + a)**n*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*(I*a*tan(c + d*x) + a)**n*sqrt(cot(c + d*x))*appellf1(sympy.S(-1)/2, 1, 1 - n, sympy.S.Half, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_793():
    f = (I*a*tan(c + d*x) + a)**n*sqrt(cot(c + d*x))
    F = 2*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, I*tan(c + d*x), -I*tan(c + d*x))/(d*(I*tan(c + d*x) + 1)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_794():
    f = (I*a*tan(c + d*x) + a)**n/sqrt(cot(c + d*x))
    F = 2*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S(3)/2, 1, 1 - n, sympy.S(5)/2, I*tan(c + d*x), -I*tan(c + d*x))/(3*d*(I*tan(c + d*x) + 1)**n*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_795():
    f = (I*a*tan(c + d*x) + a)**n/cot(c + d*x)**(sympy.S(3)/2)
    F = 2*(I*a*tan(c + d*x) + a)**n*appellf1(sympy.S(5)/2, 1, 1 - n, sympy.S(7)/2, I*tan(c + d*x), -I*tan(c + d*x))/(5*d*(I*tan(c + d*x) + 1)**n*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_796():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*a*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*sqrt(cot(c + d*x))/d - 2*b*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_797():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*b*sqrt(cot(c + d*x))/d + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_798():
    f = (a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(cot(c + d*x))/d + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_799():
    f = (a + b*tan(c + d*x))*sqrt(cot(c + d*x))
    F = -sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_800():
    f = (a + b*tan(c + d*x))/sqrt(cot(c + d*x))
    F = 2*b/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_801():
    f = (a + b*tan(c + d*x))/cot(c + d*x)**(sympy.S(3)/2)
    F = 2*a/(d*sqrt(cot(c + d*x))) + 2*b/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_802():
    f = (a + b*tan(c + d*x))/cot(c + d*x)**(sympy.S(5)/2)
    F = 2*a/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - 2*b/(d*sqrt(cot(c + d*x))) + 2*b/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_803():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*a**2*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - 4*a*b*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + 4*a*b*sqrt(cot(c + d*x))/d + (2*a**2 - 2*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_804():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*a**2*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 4*a*b*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*a**2 - 2*b**2)*sqrt(cot(c + d*x))/d - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_805():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a**2*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 4*a*b*sqrt(cot(c + d*x))/d + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_806():
    f = (a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2*sqrt(cot(c + d*x))/d + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_807():
    f = (a + b*tan(c + d*x))**2*sqrt(cot(c + d*x))
    F = 2*b**2/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_808():
    f = (a + b*tan(c + d*x))**2/sqrt(cot(c + d*x))
    F = 4*a*b/(d*sqrt(cot(c + d*x))) + 2*b**2/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_809():
    f = (a + b*tan(c + d*x))**2/cot(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*b**2/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + (2*a**2 - 2*b**2)/(d*sqrt(cot(c + d*x))) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_810():
    f = (a + b*tan(c + d*x))**2/cot(c + d*x)**(sympy.S(5)/2)
    F = -4*a*b/(d*sqrt(cot(c + d*x))) + 4*a*b/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*b**2/(7*d*cot(c + d*x)**(sympy.S(7)/2)) + (2*a**2 - 2*b**2)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_811():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(9)/2)
    F = -32*a**2*b*cot(c + d*x)**(sympy.S(5)/2)/(35*d) - 2*a**2*(a*cot(c + d*x) + b)*cot(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*(a**2 - 3*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b*(3*a**2 - b**2)*sqrt(cot(c + d*x))/d - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_812():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(7)/2)
    F = -8*a**2*b*cot(c + d*x)**(sympy.S(3)/2)/(5*d) - 2*a**2*(a*cot(c + d*x) + b)*cot(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*(a**2 - 3*b**2)*sqrt(cot(c + d*x))/d + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_813():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(5)/2)
    F = -16*a**2*b*sqrt(cot(c + d*x))/(3*d) - 2*a**2*(a*cot(c + d*x) + b)*sqrt(cot(c + d*x))/(3*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_814():
    f = (a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a*(a**2 + b**2)*sqrt(cot(c + d*x))/d + 2*b**2*(a*cot(c + d*x) + b)/(d*sqrt(cot(c + d*x))) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_815():
    f = (a + b*tan(c + d*x))**3*sqrt(cot(c + d*x))
    F = 16*a*b**2/(3*d*sqrt(cot(c + d*x))) + 2*b**2*(a*cot(c + d*x) + b)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_816():
    f = (a + b*tan(c + d*x))**3/sqrt(cot(c + d*x))
    F = 8*a*b**2/(5*d*cot(c + d*x)**(sympy.S(3)/2)) + 2*b**2*(a*cot(c + d*x) + b)/(5*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*b*(3*a**2 - b**2)/(d*sqrt(cot(c + d*x))) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_817():
    f = (a + b*tan(c + d*x))**3/cot(c + d*x)**(sympy.S(3)/2)
    F = 32*a*b**2/(35*d*cot(c + d*x)**(sympy.S(5)/2)) + 2*a*(a**2 - 3*b**2)/(d*sqrt(cot(c + d*x))) + 2*b**2*(a*cot(c + d*x) + b)/(7*d*cot(c + d*x)**(sympy.S(7)/2)) + 2*b*(3*a**2 - b**2)/(3*d*cot(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_818():
    f = cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))
    F = sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - 2*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) + 2*b*sqrt(cot(c + d*x))/(a**2*d) - 2*b**(sympy.S(7)/2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(5)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_819():
    f = cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))
    F = -sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - 2*sqrt(cot(c + d*x))/(a*d) + 2*b**(sympy.S(5)/2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(3)/2)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_820():
    f = sqrt(cot(c + d*x))/(a + b*tan(c + d*x))
    F = -sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - 2*b**(sympy.S(3)/2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(a)*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_821():
    f = 1/((a + b*tan(c + d*x))*sqrt(cot(c + d*x)))
    F = 2*sqrt(a)*sqrt(b)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) - sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_822():
    f = 1/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2))
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(b)*d*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a - b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) - sqrt(2)*(a + b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_823():
    f = 1/((a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2))
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)) - sqrt(2)*(a - b)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a - b)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)) + sqrt(2)*(a + b)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)) + 2/(b*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_824():
    f = cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**2
    F = sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b**2*cot(c + d*x)**(sympy.S(5)/2)/(a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - (2*a**2 + 5*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(a**2 + b**2)) + b*(4*a**2 + 5*b**2)*sqrt(cot(c + d*x))/(a**3*d*(a**2 + b**2)) - b**(sympy.S(7)/2)*(9*a**2 + 5*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(7)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_825():
    f = cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**2
    F = -sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + b**2*cot(c + d*x)**(sympy.S(3)/2)/(a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - (2*a**2 + 3*b**2)*sqrt(cot(c + d*x))/(a**2*d*(a**2 + b**2)) + b**(sympy.S(5)/2)*(7*a**2 + 3*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(5)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_826():
    f = sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**2
    F = -sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + b**2*sqrt(cot(c + d*x))/(a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - b**(sympy.S(3)/2)*(5*a**2 + b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(a**(sympy.S(3)/2)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_827():
    f = 1/((a + b*tan(c + d*x))**2*sqrt(cot(c + d*x)))
    F = -b*sqrt(cot(c + d*x))/(d*(a**2 + b**2)*(a*cot(c + d*x) + b)) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(b)*(3*a**2 - b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(a)*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_828():
    f = 1/((a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(a)*(a**2 - 3*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(sqrt(b)*d*(a**2 + b**2)**2) + a*sqrt(cot(c + d*x))/(d*(a**2 + b**2)*(a*cot(c + d*x) + b)) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_829():
    f = 1/((a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(5)/2))
    F = -a**(sympy.S(3)/2)*(a**2 + 5*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d*(a**2 + b**2)**2) - a**2*sqrt(cot(c + d*x))/(b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)) - sqrt(2)*(a**2 - 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 - 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_830():
    f = 1/((a + b*tan(c + d*x))**2*cot(c + d*x)**(sympy.S(7)/2))
    F = a**(sympy.S(5)/2)*(3*a**2 + 7*b**2)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d*(a**2 + b**2)**2) - a**2/(b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)*sqrt(cot(c + d*x))) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 - 2*a*b - b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**2) - sqrt(2)*(a**2 + 2*a*b - b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + sqrt(2)*(a**2 + 2*a*b - b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**2) + (3*a**2 + 2*b**2)/(b**2*d*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_831():
    f = cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**3
    F = sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b**2*cot(c + d*x)**(sympy.S(7)/2)/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) + b**2*(15*a**2 + 7*b**2)*cot(c + d*x)**(sympy.S(5)/2)/(4*a**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - (8*a**4 + 67*a**2*b**2 + 35*b**4)*cot(c + d*x)**(sympy.S(3)/2)/(12*a**3*d*(a**2 + b**2)**2) + b*(24*a**4 + 67*a**2*b**2 + 35*b**4)*sqrt(cot(c + d*x))/(4*a**4*d*(a**2 + b**2)**2) - b**(sympy.S(7)/2)*(99*a**4 + 102*a**2*b**2 + 35*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(9)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_832():
    f = cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**3
    F = sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b**2*cot(c + d*x)**(sympy.S(5)/2)/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) + b**2*(13*a**2 + 5*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(4*a**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - (8*a**4 + 31*a**2*b**2 + 15*b**4)*sqrt(cot(c + d*x))/(4*a**3*d*(a**2 + b**2)**2) + b**(sympy.S(5)/2)*(63*a**4 + 46*a**2*b**2 + 15*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(7)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_833():
    f = sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**3
    F = -sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + b**2*cot(c + d*x)**(sympy.S(3)/2)/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) + b**2*(11*a**2 + 3*b**2)*sqrt(cot(c + d*x))/(4*a**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - b**(sympy.S(3)/2)*(35*a**4 + 6*a**2*b**2 + 3*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(5)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_834():
    f = 1/((a + b*tan(c + d*x))**3*sqrt(cot(c + d*x)))
    F = -sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + b**2*sqrt(cot(c + d*x))/(2*a*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) - b*(9*a**2 + b**2)*sqrt(cot(c + d*x))/(4*a*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) + sqrt(b)*(15*a**4 - 18*a**2*b**2 - b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*a**(sympy.S(3)/2)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_835():
    f = 1/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(3)/2))
    F = -b*sqrt(cot(c + d*x))/(d*(2*a**2 + 2*b**2)*(a*cot(c + d*x) + b)**2) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) + (5*a**2 - 3*b**2)*sqrt(cot(c + d*x))/(4*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - (3*a**4 - 26*a**2*b**2 + 3*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*sqrt(a)*sqrt(b)*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_836():
    f = 1/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(a)*(a**4 + 18*a**2*b**2 - 15*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*b**(sympy.S(3)/2)*d*(a**2 + b**2)**3) + a*sqrt(cot(c + d*x))/(d*(2*a**2 + 2*b**2)*(a*cot(c + d*x) + b)**2) - a*(a**2 - 7*b**2)*sqrt(cot(c + d*x))/(4*b*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_837():
    f = 1/((a + b*tan(c + d*x))**3*cot(c + d*x)**(sympy.S(7)/2))
    F = -a**(sympy.S(3)/2)*(3*a**4 + 6*a**2*b**2 + 35*b**4)*atan(sqrt(a)*sqrt(cot(c + d*x))/sqrt(b))/(4*b**(sympy.S(5)/2)*d*(a**2 + b**2)**3) - a**2*sqrt(cot(c + d*x))/(2*b*d*(a**2 + b**2)*(a*cot(c + d*x) + b)**2) - a**2*(3*a**2 + 11*b**2)*sqrt(cot(c + d*x))/(4*b**2*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b)) - sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(-sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) + sqrt(2)*(a - b)*(a**2 + 4*a*b + b**2)*log(sqrt(2)*sqrt(cot(c + d*x)) + cot(c + d*x) + 1)/(4*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) - 1)/(2*d*(a**2 + b**2)**3) - sqrt(2)*(a + b)*(a**2 - 4*a*b + b**2)*atan(sqrt(2)*sqrt(cot(c + d*x)) + 1)/(2*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_838():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*d) + sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(15*a*d) + sqrt(a + b*tan(c + d*x))*(30*a**2 + 4*b**2)*sqrt(cot(c + d*x))/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_839():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*d) + I*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - 2*b*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_840():
    f = sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/d - sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_841():
    f = sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))
    F = -I*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_842():
    f = sqrt(a + b*tan(c + d*x))/sqrt(cot(c + d*x))
    F = 2*sqrt(b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_843():
    f = sqrt(a + b*tan(c + d*x))/cot(c + d*x)**(sympy.S(3)/2)
    F = a*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d) + sqrt(a + b*tan(c + d*x))/(d*sqrt(cot(c + d*x))) + I*sqrt(I*a - b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*sqrt(I*a + b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_844():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - 16*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(35*d) - (I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(70*a**2 - 6*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(105*a*d) + 4*b*sqrt(a + b*tan(c + d*x))*(70*a**2 + 3*b**2)*sqrt(cot(c + d*x))/(105*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_845():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 4*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(5*d) - I*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + sqrt(a + b*tan(c + d*x))*(10*a**2 - 2*b**2)*sqrt(cot(c + d*x))/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_846():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 8*b*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(3*d) + (I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_847():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/d + I*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_848():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x))
    F = 2*b**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_849():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/sqrt(cot(c + d*x))
    F = 3*a*sqrt(b)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b*sqrt(a + b*tan(c + d*x))/(d*sqrt(cot(c + d*x))) - I*(I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_850():
    f = (a + b*tan(c + d*x))**(sympy.S(3)/2)/cot(c + d*x)**(sympy.S(3)/2)
    F = 3*a*sqrt(a + b*tan(c + d*x))/(4*d*sqrt(cot(c + d*x))) + (a + b*tan(c + d*x))**(sympy.S(3)/2)/(2*d*sqrt(cot(c + d*x))) + (I*a - b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (3*a**2 - 8*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_851():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(11)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(9)/2)/(9*d) - 38*a*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/(63*d) + sqrt(a + b*tan(c + d*x))*(42*a**2 - 50*b**2)*cot(c + d*x)**(sympy.S(5)/2)/(105*d) + (I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + 2*b*sqrt(a + b*tan(c + d*x))*(231*a**2 - 5*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(315*a*d) - sqrt(a + b*tan(c + d*x))*(630*a**4 - 966*a**2*b**2 - 20*b**4)*sqrt(cot(c + d*x))/(315*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_852():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(9)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(7)/2)/(7*d) - 6*a*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(7*d) + sqrt(a + b*tan(c + d*x))*(14*a**2 - 18*b**2)*cot(c + d*x)**(sympy.S(3)/2)/(21*d) + I*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + 2*b*sqrt(a + b*tan(c + d*x))*(49*a**2 - 3*b**2)*sqrt(cot(c + d*x))/(21*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_853():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(7)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2)/(5*d) - 22*a*b*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(15*d) + sqrt(a + b*tan(c + d*x))*(30*a**2 - 46*b**2)*sqrt(cot(c + d*x))/(15*d) - (I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_854():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*d) - 14*a*b*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(3*d) - I*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_855():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2)
    F = -2*a**2*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/d + 2*b**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - (I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_856():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(cot(c + d*x))
    F = 5*a*b**(sympy.S(3)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + b**2*sqrt(a + b*tan(c + d*x))/(d*sqrt(cot(c + d*x))) + I*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + I*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_857():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/sqrt(cot(c + d*x))
    F = 9*a*b*sqrt(a + b*tan(c + d*x))/(4*d*sqrt(cot(c + d*x))) + sqrt(b)*(15*a**2 - 8*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(4*d) + b**2*sqrt(a + b*tan(c + d*x))/(2*d*cot(c + d*x)**(sympy.S(3)/2)) - (I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d + (I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_858():
    f = (a + b*tan(c + d*x))**(sympy.S(5)/2)/cot(c + d*x)**(sympy.S(3)/2)
    F = 13*a*b*sqrt(a + b*tan(c + d*x))/(12*d*cot(c + d*x)**(sympy.S(3)/2)) + 5*a*(a**2 - 8*b**2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(8*sqrt(b)*d) + b**2*sqrt(a + b*tan(c + d*x))/(3*d*cot(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*tan(c + d*x))*(11*a**2 - 8*b**2)/(8*d*sqrt(cot(c + d*x))) - I*(I*a - b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d - I*(I*a + b)**(sympy.S(5)/2)*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_859():
    f = cot(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*tan(c + d*x))
    F = -sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d) + 4*b*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_860():
    f = cot(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*tan(c + d*x))
    F = I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) - 2*sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_861():
    f = sqrt(cot(c + d*x))/sqrt(a + b*tan(c + d*x))
    F = sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_862():
    f = 1/(sqrt(a + b*tan(c + d*x))*sqrt(cot(c + d*x)))
    F = -I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_863():
    f = 1/(sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + 2*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_864():
    f = 1/(sqrt(a + b*tan(c + d*x))*cot(c + d*x)**(sympy.S(5)/2))
    F = -a*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a + b)) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*sqrt(I*a - b)) + sqrt(a + b*tan(c + d*x))/(b*d*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_865():
    f = cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d*sqrt(a + b*tan(c + d*x))) + 8*b*sqrt(cot(c + d*x))/(3*a**2*d*sqrt(a + b*tan(c + d*x))) + 2*b**2*(5*a**2 + 8*b**2)/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_866():
    f = cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = -sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) - 2*sqrt(cot(c + d*x))/(a*d*sqrt(a + b*tan(c + d*x))) - 2*b*(a**2 + 2*b**2)/(a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_867():
    f = sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(3)/2)
    F = I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + 2*b**2/(a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_868():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*sqrt(cot(c + d*x)))
    F = -2*b/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x))) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_869():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*a/(d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x))) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_870():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = -2*a**2/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*sqrt(cot(c + d*x))) - sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + 2*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_871():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(3)/2)*cot(c + d*x)**(sympy.S(7)/2))
    F = -2*a**2/(b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)*cot(c + d*x)**(sympy.S(3)/2)) - 3*a*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(b**(sympy.S(5)/2)*d) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(3)/2)) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(3)/2)) + sqrt(a + b*tan(c + d*x))*(3*a**2 + b**2)/(b**2*d*(a**2 + b**2)*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_872():
    f = cot(c + d*x)**(sympy.S(5)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2*cot(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + 4*b*sqrt(cot(c + d*x))/(a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) + 2*b**2*(7*a**2 + 8*b**2)/(3*a**3*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + 4*b**2*(4*a**4 + 15*a**2*b**2 + 8*b**4)/(3*a**4*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_873():
    f = cot(c + d*x)**(sympy.S(3)/2)/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2*sqrt(cot(c + d*x))/(a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)) - 2*b*(3*a**2 + 4*b**2)/(3*a**2*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) - 2*b*(3*a**4 + 17*a**2*b**2 + 8*b**4)/(3*a**3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_874():
    f = sqrt(cot(c + d*x))/(a + b*tan(c + d*x))**(sympy.S(5)/2)
    F = -sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + 2*b**2/(3*a*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + 4*b**2*(4*a**2 + b**2)/(3*a**2*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_875():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*sqrt(cot(c + d*x)))
    F = -2*b/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*sqrt(cot(c + d*x))) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) - 2*b*(5*a**2 - b**2)/(3*a*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_876():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(3)/2))
    F = 2*a/(d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*sqrt(cot(c + d*x))) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2)) + (4*a**2 - 8*b**2)/(3*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_877():
    f = 1/((a + b*tan(c + d*x))**(sympy.S(5)/2)*cot(c + d*x)**(sympy.S(5)/2))
    F = -2*a**2/(3*b*d*(a + b*tan(c + d*x))**(sympy.S(3)/2)*(a**2 + b**2)*sqrt(cot(c + d*x))) + 2*a*(a**2 + 7*b**2)/(3*b*d*sqrt(a + b*tan(c + d*x))*(a**2 + b**2)**2*sqrt(cot(c + d*x))) - I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atanh(sqrt(I*a + b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a + b)**(sympy.S(5)/2)) + I*sqrt(tan(c + d*x))*sqrt(cot(c + d*x))*atan(sqrt(I*a - b)*sqrt(tan(c + d*x))/sqrt(a + b*tan(c + d*x)))/(d*(I*a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_878():
    f = (d*cot(e + f*x))**n*(a + b*tan(e + f*x))**3
    F = a**2*b*d**2*(d*cot(e + f*x))**(n - 2)*(1 - 2*n)/(f*(1 - n)*(2 - n)) + a**2*d**2*(d*cot(e + f*x))**(n - 2)*(a*cot(e + f*x) + b)/(f*(1 - n)) - a*d*(d*cot(e + f*x))**(n - 1)*(a**2 - 3*b**2)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), -cot(e + f*x)**2)/(f*(1 - n)) - b*d**2*(d*cot(e + f*x))**(n - 2)*(3*a**2 - b**2)*hyper((1, n/2 - 1), (n/2,), -cot(e + f*x)**2)/(f*(2 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_879():
    f = (d*cot(e + f*x))**n*(a + b*tan(e + f*x))**2
    F = a**2*d*(d*cot(e + f*x))**(n - 1)/(f*(1 - n)) - 2*a*b*(d*cot(e + f*x))**n*hyper((1, n/2), (n/2 + 1,), -cot(e + f*x)**2)/(f*n) - d*(d*cot(e + f*x))**(n - 1)*(a**2 - b**2)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), -cot(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_880():
    f = (d*cot(e + f*x))**n*(a + b*tan(e + f*x))
    F = -a*(d*cot(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -cot(e + f*x)**2)/(d*f*(n + 1)) - b*(d*cot(e + f*x))**n*hyper((1, n/2), (n/2 + 1,), -cot(e + f*x)**2)/(f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_881():
    f = (d*cot(e + f*x))**n/(a + b*tan(e + f*x))
    F = -a**2*(d*cot(e + f*x))**(n + 2)*hyper((1, n + 2), (n + 3,), -a*cot(e + f*x)/b)/(b*d**2*f*(a**2 + b**2)*(n + 2)) + a*(d*cot(e + f*x))**(n + 3)*hyper((1, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -cot(e + f*x)**2)/(d**3*f*(a**2 + b**2)*(n + 3)) - b*(d*cot(e + f*x))**(n + 2)*hyper((1, n/2 + 1), (n/2 + 2,), -cot(e + f*x)**2)/(d**2*f*(a**2 + b**2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_882():
    f = (d*cot(e + f*x))**n/(a + b*tan(e + f*x))**2
    F = -a**2*(d*cot(e + f*x))**(n + 3)/(b*d**3*f*(a**2 + b**2)*(a*cot(e + f*x) + b)) + a**2*(d*cot(e + f*x))**(n + 3)*(a**2*(n + 2) + b**2*n)*hyper((1, n + 3), (n + 4,), -a*cot(e + f*x)/b)/(b**2*d**3*f*(a**2 + b**2)**2*(n + 3)) + 2*a*b*(d*cot(e + f*x))**(n + 4)*hyper((1, n/2 + 2), (n/2 + 3,), -cot(e + f*x)**2)/(d**4*f*(a**2 + b**2)**2*(n + 4)) + (d*cot(e + f*x))**(n + 3)*(a**2 - b**2)*hyper((1, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -cot(e + f*x)**2)/(d**3*f*(a**2 + b**2)**2*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_883():
    f = (d*cot(e + f*x))**n*(a + b*tan(e + f*x))**m
    F = (d*cot(e + f*x))**n*(a + b*tan(e + f*x))**m*tan(e + f*x)*appellf1(1 - n, 1, -m, 2 - n, -I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*f*(1 - n)*(1 + b*tan(e + f*x)/a)**m) + (d*cot(e + f*x))**n*(a + b*tan(e + f*x))**m*tan(e + f*x)*appellf1(1 - n, 1, -m, 2 - n, I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*f*(1 - n)*(1 + b*tan(e + f*x)/a)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_884():
    f = (a + b*tan(c + d*x))**n*cot(c + d*x)**(sympy.S(3)/2)
    F = -(a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n) - (a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))*appellf1(sympy.S(-1)/2, 1, -n, sympy.S.Half, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_885():
    f = (a + b*tan(c + d*x))**n*sqrt(cot(c + d*x))
    F = (a + b*tan(c + d*x))**n*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(cot(c + d*x))) + (a + b*tan(c + d*x))**n*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(d*(1 + b*tan(c + d*x)/a)**n*sqrt(cot(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_886():
    f = (a + b*tan(c + d*x))**n/sqrt(cot(c + d*x))
    F = (a + b*tan(c + d*x))**n*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(3)/2)) + (a + b*tan(c + d*x))**n*appellf1(sympy.S(3)/2, 1, -n, sympy.S(5)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(3*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_887():
    f = (a + b*tan(c + d*x))**n/cot(c + d*x)**(sympy.S(3)/2)
    F = (a + b*tan(c + d*x))**n*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, -I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(5)/2)) + (a + b*tan(c + d*x))**n*appellf1(sympy.S(5)/2, 1, -n, sympy.S(7)/2, I*tan(c + d*x), -b*tan(c + d*x)/a)/(5*d*(1 + b*tan(c + d*x)/a)**n*cot(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_888():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)
    F = -I*c*(I*a*tan(e + f*x) + a)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_889():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)
    F = a*c*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_890():
    f = (-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)
    F = I*c/(f*(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_891():
    f = (-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**2
    F = I*c/(2*f*(I*a*tan(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_892():
    f = (-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**3
    F = I*c/(3*f*(I*a*tan(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_893():
    f = (I*a*tan(e + f*x) + a)**4*(-I*c*tan(e + f*x) + c)**2
    F = -I*c**2*(I*a*tan(e + f*x) + a)**4/(2*f) + I*c**2*(I*a*tan(e + f*x) + a)**5/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_894():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**2
    F = a**3*c**2*tan(e + f*x)**3/(3*f) + a**3*c**2*tan(e + f*x)/f + I*a**3*c**2*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_895():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**2
    F = a**2*c**2*tan(e + f*x)**3/(3*f) + a**2*c**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_896():
    f = (-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)
    F = 2*I*c**2/(f*(I*a*tan(e + f*x) + a)) - c**2*x/a - I*c**2*log(cos(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_897():
    f = (-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)**2
    F = c**2*tan(e + f*x)/(f*(I*a*tan(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_898():
    f = (-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)**3
    F = 2*I*c**2/(3*f*(I*a*tan(e + f*x) + a)**3) - I*c**2/(2*a*f*(I*a*tan(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_899():
    f = (-I*c*tan(e + f*x) + c)**2/(I*a*tan(e + f*x) + a)**4
    F = -I*a**2*c**2/(3*f*(I*a**2*tan(e + f*x) + a**2)**3) + I*c**2/(2*f*(I*a*tan(e + f*x) + a)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_900():
    f = (I*a*tan(e + f*x) + a)**5*(-I*c*tan(e + f*x) + c)**3
    F = -4*I*c**3*(I*a*tan(e + f*x) + a)**5/(5*f) + 2*I*c**3*(I*a*tan(e + f*x) + a)**6/(3*a*f) - I*c**3*(I*a*tan(e + f*x) + a)**7/(7*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_901():
    f = (I*a*tan(e + f*x) + a)**4*(-I*c*tan(e + f*x) + c)**3
    F = a**4*c**3*tan(e + f*x)**5/(5*f) + 2*a**4*c**3*tan(e + f*x)**3/(3*f) + a**4*c**3*tan(e + f*x)/f + I*a**4*c**3*sec(e + f*x)**6/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_902():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**3
    F = a**3*c**3*tan(e + f*x)**5/(5*f) + 2*a**3*c**3*tan(e + f*x)**3/(3*f) + a**3*c**3*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_903():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**3
    F = a**2*c**3*tan(e + f*x)**3/(3*f) + a**2*c**3*tan(e + f*x)/f - I*a**2*c**3*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_904():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**3
    F = I*a*(-I*c*tan(e + f*x) + c)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_905():
    f = (-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)
    F = 4*I*c**3/(f*(I*a*tan(e + f*x) + a)) - 4*c**3*x/a - 4*I*c**3*log(cos(e + f*x))/(a*f) + c**3*tan(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_906():
    f = (-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**2
    F = -4*I*c**3/(f*(I*a**2*tan(e + f*x) + a**2)) + 2*I*c**3/(f*(I*a*tan(e + f*x) + a)**2) + c**3*x/a**2 + I*c**3*log(cos(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_907():
    f = (-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**3
    F = I*c**3*(-I*a**2*tan(e + f*x) + a**2)**3/(6*f*(I*a**3*tan(e + f*x) + a**3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_908():
    f = (-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**4
    F = I*c**3/(2*f*(I*a**2*tan(e + f*x) + a**2)**2) + I*c**3/(f*(I*a*tan(e + f*x) + a)**4) - 4*I*c**3/(3*a*f*(I*a*tan(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_909():
    f = (-I*c*tan(e + f*x) + c)**3/(I*a*tan(e + f*x) + a)**5
    F = -I*a**3*c**3/(f*(I*a**2*tan(e + f*x) + a**2)**4) + 4*I*c**3/(5*f*(I*a*tan(e + f*x) + a)**5) + I*c**3/(3*a**2*f*(I*a*tan(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_910():
    f = (I*a*tan(e + f*x) + a)**5*(-I*c*tan(e + f*x) + c)**4
    F = a**5*c**4*tan(e + f*x)**7/(7*f) + 3*a**5*c**4*tan(e + f*x)**5/(5*f) + a**5*c**4*tan(e + f*x)**3/f + a**5*c**4*tan(e + f*x)/f + I*a**5*c**4*sec(e + f*x)**8/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_911():
    f = (I*a*tan(e + f*x) + a)**4*(-I*c*tan(e + f*x) + c)**4
    F = a**4*c**4*tan(e + f*x)**7/(7*f) + 3*a**4*c**4*tan(e + f*x)**5/(5*f) + a**4*c**4*tan(e + f*x)**3/f + a**4*c**4*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_912():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**4
    F = a**3*c**4*tan(e + f*x)**5/(5*f) + 2*a**3*c**4*tan(e + f*x)**3/(3*f) + a**3*c**4*tan(e + f*x)/f - I*a**3*c**4*sec(e + f*x)**6/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_913():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**4
    F = I*a**2*(-I*c*tan(e + f*x) + c)**4/(2*f) - I*a**2*(-I*c*tan(e + f*x) + c)**5/(5*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_914():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**4
    F = I*a*(-I*c*tan(e + f*x) + c)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_915():
    f = (-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)
    F = 8*I*c**4/(f*(I*a*tan(e + f*x) + a)) - 12*c**4*x/a - 12*I*c**4*log(cos(e + f*x))/(a*f) - I*c**4*tan(e + f*x)**2/(2*a*f) + 5*c**4*tan(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_916():
    f = (-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**2
    F = -12*I*c**4/(f*(I*a**2*tan(e + f*x) + a**2)) + 4*I*c**4/(f*(I*a*tan(e + f*x) + a)**2) + 6*c**4*x/a**2 + 6*I*c**4*log(cos(e + f*x))/(a**2*f) - c**4*tan(e + f*x)/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_917():
    f = (-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**3
    F = 6*I*c**4/(f*(I*a**3*tan(e + f*x) + a**3)) + 8*I*c**4/(3*f*(I*a*tan(e + f*x) + a)**3) - 6*I*c**4/(a*f*(I*a*tan(e + f*x) + a)**2) - c**4*x/a**3 - I*c**4*log(cos(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_918():
    f = (-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**4
    F = I*c**4*(-I*a**2*tan(e + f*x) + a**2)**4/(8*f*(I*a**3*tan(e + f*x) + a**3)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_919():
    f = (-I*c*tan(e + f*x) + c)**4/(I*a*tan(e + f*x) + a)**5
    F = I*c**4*(-I*tan(e + f*x) + 1)**4/(10*f*(I*a*tan(e + f*x) + a)**5) + I*c**4*(-I*a*tan(e + f*x) + a)**4/(80*a**5*f*(I*a*tan(e + f*x) + a)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_920():
    f = (I*a*tan(e + f*x) + a)**4/(-I*c*tan(e + f*x) + c)
    F = -8*I*a**4/(f*(-I*c*tan(e + f*x) + c)) - 12*a**4*x/c + 12*I*a**4*log(cos(e + f*x))/(c*f) + I*a**4*tan(e + f*x)**2/(2*c*f) + 5*a**4*tan(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_921():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)
    F = -4*I*a**3/(f*(-I*c*tan(e + f*x) + c)) - 4*a**3*x/c + 4*I*a**3*log(cos(e + f*x))/(c*f) + a**3*tan(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_922():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)
    F = -2*I*a**2/(f*(-I*c*tan(e + f*x) + c)) - a**2*x/c + I*a**2*log(cos(e + f*x))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_923():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)
    F = -I*a/(f*(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_924():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c))
    F = x/(2*a*c) + sin(e + f*x)*cos(e + f*x)/(2*a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_925():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c))
    F = 3*x/(8*a**2*c) + sin(e + f*x)*cos(e + f*x)**3/(4*a**2*c*f) + 3*sin(e + f*x)*cos(e + f*x)/(8*a**2*c*f) + I*cos(e + f*x)**4/(4*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_926():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c))
    F = I*c**2/(12*a**3*f*(I*c*tan(e + f*x) + c)**3) + I*c/(8*a**3*f*(I*c*tan(e + f*x) + c)**2) + 3*I/(16*a**3*f*(I*c*tan(e + f*x) + c)) - I/(16*a**3*f*(-I*c*tan(e + f*x) + c)) + x/(4*a**3*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_927():
    f = (I*a*tan(e + f*x) + a)**4/(-I*c*tan(e + f*x) + c)**2
    F = 12*I*a**4/(f*(-I*c**2*tan(e + f*x) + c**2)) - 4*I*a**4/(f*(-I*c*tan(e + f*x) + c)**2) + 6*a**4*x/c**2 - 6*I*a**4*log(cos(e + f*x))/(c**2*f) - a**4*tan(e + f*x)/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_928():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**2
    F = 4*I*a**3/(f*(-I*c**2*tan(e + f*x) + c**2)) - 2*I*a**3/(f*(-I*c*tan(e + f*x) + c)**2) + a**3*x/c**2 - I*a**3*log(cos(e + f*x))/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_929():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**2
    F = a**2*tan(e + f*x)/(f*(-I*c*tan(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_930():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**2
    F = -I*a/(2*f*(-I*c*tan(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_931():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**2)
    F = I/(8*a*f*(I*c**2*tan(e + f*x) + c**2)) - I/(4*a*f*(-I*c**2*tan(e + f*x) + c**2)) - I/(8*a*f*(-I*c*tan(e + f*x) + c)**2) + 3*x/(8*a*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_932():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**2)
    F = 3*x/(8*a**2*c**2) + sin(e + f*x)*cos(e + f*x)**3/(4*a**2*c**2*f) + 3*sin(e + f*x)*cos(e + f*x)/(8*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_933():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**2)
    F = 5*x/(16*a**3*c**2) + sin(e + f*x)*cos(e + f*x)**5/(6*a**3*c**2*f) + 5*sin(e + f*x)*cos(e + f*x)**3/(24*a**3*c**2*f) + 5*sin(e + f*x)*cos(e + f*x)/(16*a**3*c**2*f) + I*cos(e + f*x)**6/(6*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_934():
    f = (I*a*tan(e + f*x) + a)**6/(-I*c*tan(e + f*x) + c)**3
    F = -80*I*a**6/(f*(-I*c**3*tan(e + f*x) + c**3)) - 32*I*a**6/(3*f*(-I*c*tan(e + f*x) + c)**3) + 40*I*a**6/(c*f*(-I*c*tan(e + f*x) + c)**2) - 40*a**6*x/c**3 + 40*I*a**6*log(cos(e + f*x))/(c**3*f) + I*a**6*tan(e + f*x)**2/(2*c**3*f) + 9*a**6*tan(e + f*x)/(c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_935():
    f = (I*a*tan(e + f*x) + a)**5/(-I*c*tan(e + f*x) + c)**3
    F = 16*I*a**5*c**5/(f*(-I*c**4*tan(e + f*x) + c**4)**2) - 24*I*a**5/(f*(-I*c**3*tan(e + f*x) + c**3)) - 16*I*a**5/(3*f*(-I*c*tan(e + f*x) + c)**3) - 8*a**5*x/c**3 + 8*I*a**5*log(cos(e + f*x))/(c**3*f) + a**5*tan(e + f*x)/(c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_936():
    f = (I*a*tan(e + f*x) + a)**4/(-I*c*tan(e + f*x) + c)**3
    F = -6*I*a**4/(f*(-I*c**3*tan(e + f*x) + c**3)) - 8*I*a**4/(3*f*(-I*c*tan(e + f*x) + c)**3) + 6*I*a**4/(c*f*(-I*c*tan(e + f*x) + c)**2) - a**4*x/c**3 + I*a**4*log(cos(e + f*x))/(c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_937():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**3
    F = -I*a**3*(I*c**2*tan(e + f*x) + c**2)**3/(6*f*(-I*c**3*tan(e + f*x) + c**3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_938():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**3
    F = -2*I*a**2/(3*f*(-I*c*tan(e + f*x) + c)**3) + I*a**2/(2*c*f*(-I*c*tan(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_939():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**3
    F = -I*a/(3*f*(-I*c*tan(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_940():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**3)
    F = I/(16*a*f*(I*c**3*tan(e + f*x) + c**3)) - 3*I/(16*a*f*(-I*c**3*tan(e + f*x) + c**3)) - I/(12*a*f*(-I*c*tan(e + f*x) + c)**3) - I/(8*a*c*f*(-I*c*tan(e + f*x) + c)**2) + x/(4*a*c**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_941():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**3)
    F = I/(8*a**2*f*(I*c**3*tan(e + f*x) + c**3)) - 3*I/(16*a**2*f*(-I*c**3*tan(e + f*x) + c**3)) - I/(24*a**2*f*(-I*c*tan(e + f*x) + c)**3) + I/(32*a**2*c*f*(I*c*tan(e + f*x) + c)**2) - 3*I/(32*a**2*c*f*(-I*c*tan(e + f*x) + c)**2) + 5*x/(16*a**2*c**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_942():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**3)
    F = 5*x/(16*a**3*c**3) + sin(e + f*x)*cos(e + f*x)**5/(6*a**3*c**3*f) + 5*sin(e + f*x)*cos(e + f*x)**3/(24*a**3*c**3*f) + 5*sin(e + f*x)*cos(e + f*x)/(16*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_943():
    f = (I*a*tan(e + f*x) + a)**6/(-I*c*tan(e + f*x) + c)**4
    F = 40*I*a**6/(f*(-I*c**4*tan(e + f*x) + c**4)) - 40*I*a**6/(f*(-I*c**2*tan(e + f*x) + c**2)**2) - 8*I*a**6/(f*(-I*c*tan(e + f*x) + c)**4) + 80*I*a**6/(3*c*f*(-I*c*tan(e + f*x) + c)**3) + 10*a**6*x/c**4 - 10*I*a**6*log(cos(e + f*x))/(c**4*f) - a**6*tan(e + f*x)/(c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_944():
    f = (I*a*tan(e + f*x) + a)**5/(-I*c*tan(e + f*x) + c)**4
    F = 32*I*a**5*c**5/(3*f*(-I*c**3*tan(e + f*x) + c**3)**3) + 8*I*a**5/(f*(-I*c**4*tan(e + f*x) + c**4)) - 12*I*a**5/(f*(-I*c**2*tan(e + f*x) + c**2)**2) - 4*I*a**5/(f*(-I*c*tan(e + f*x) + c)**4) + a**5*x/c**4 - I*a**5*log(cos(e + f*x))/(c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_945():
    f = (I*a*tan(e + f*x) + a)**4/(-I*c*tan(e + f*x) + c)**4
    F = -I*a**4*(I*c**2*tan(e + f*x) + c**2)**4/(8*f*(-I*c**3*tan(e + f*x) + c**3)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_946():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**4
    F = -I*a**3/(2*f*(-I*c**2*tan(e + f*x) + c**2)**2) - I*a**3/(f*(-I*c*tan(e + f*x) + c)**4) + 4*I*a**3/(3*c*f*(-I*c*tan(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_947():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**4
    F = I*a**2*c**2/(3*f*(-I*c**2*tan(e + f*x) + c**2)**3) - I*a**2/(2*f*(-I*c*tan(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_948():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**4
    F = -I*a/(4*f*(-I*c*tan(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_949():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**4)
    F = I/(32*a*f*(I*c**4*tan(e + f*x) + c**4)) - I/(8*a*f*(-I*c**4*tan(e + f*x) + c**4)) - 3*I/(32*a*f*(-I*c**2*tan(e + f*x) + c**2)**2) - I/(16*a*f*(-I*c*tan(e + f*x) + c)**4) - I/(12*a*c*f*(-I*c*tan(e + f*x) + c)**3) + 5*x/(32*a*c**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_950():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**4)
    F = 5*I/(64*a**2*f*(I*c**4*tan(e + f*x) + c**4)) - 5*I/(32*a**2*f*(-I*c**4*tan(e + f*x) + c**4)) + I/(64*a**2*f*(I*c**2*tan(e + f*x) + c**2)**2) - 3*I/(32*a**2*f*(-I*c**2*tan(e + f*x) + c**2)**2) - I/(32*a**2*f*(-I*c*tan(e + f*x) + c)**4) - I/(16*a**2*c*f*(-I*c*tan(e + f*x) + c)**3) + 15*x/(64*a**2*c**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_951():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**4)
    F = 15*I/(128*a**3*f*(I*c**4*tan(e + f*x) + c**4)) - 5*I/(32*a**3*f*(-I*c**4*tan(e + f*x) + c**4)) + 5*I/(128*a**3*f*(I*c**2*tan(e + f*x) + c**2)**2) - 5*I/(64*a**3*f*(-I*c**2*tan(e + f*x) + c**2)**2) - I/(64*a**3*f*(-I*c*tan(e + f*x) + c)**4) + I/(96*a**3*c*f*(I*c*tan(e + f*x) + c)**3) - I/(24*a**3*c*f*(-I*c*tan(e + f*x) + c)**3) + 35*x/(128*a**3*c**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_952():
    f = (I*a*tan(e + f*x) + a)**3*sqrt(-I*c*tan(e + f*x) + c)
    F = 8*I*a**3*sqrt(-I*c*tan(e + f*x) + c)/f - 8*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c*f) + 2*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_953():
    f = (I*a*tan(e + f*x) + a)**2*sqrt(-I*c*tan(e + f*x) + c)
    F = 4*I*a**2*sqrt(-I*c*tan(e + f*x) + c)/f - 2*I*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_954():
    f = (I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)
    F = 2*I*a*sqrt(-I*c*tan(e + f*x) + c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_955():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)
    F = sqrt(2)*I*sqrt(c)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(4*a*f) + I*sqrt(-I*c*tan(e + f*x) + c)/(2*a*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_956():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**2
    F = 3*sqrt(2)*I*sqrt(c)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a**2*f) + 3*I*sqrt(-I*c*tan(e + f*x) + c)/(16*a**2*f*(I*tan(e + f*x) + 1)) + I*sqrt(-I*c*tan(e + f*x) + c)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_957():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**3
    F = 5*sqrt(2)*I*sqrt(c)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(128*a**3*f) + 5*I*sqrt(-I*c*tan(e + f*x) + c)/(64*a**3*f*(I*tan(e + f*x) + 1)) + 5*I*sqrt(-I*c*tan(e + f*x) + c)/(48*a**3*f*(I*tan(e + f*x) + 1)**2) + I*sqrt(-I*c*tan(e + f*x) + c)/(6*a**3*f*(I*tan(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_958():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 8*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f) - 8*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c*f) + 2*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_959():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 4*I*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f) - 2*I*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_960():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = 2*I*a*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_961():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)
    F = -sqrt(2)*I*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(2*a*f) + I*c**2*sqrt(-I*c*tan(e + f*x) + c)/(a*f*(I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_962():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**2
    F = -sqrt(2)*I*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(16*a**2*f) + I*c**3*sqrt(-I*c*tan(e + f*x) + c)/(2*a**2*f*(I*c*tan(e + f*x) + c)**2) - I*c**2*sqrt(-I*c*tan(e + f*x) + c)/(8*a**2*f*(I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_963():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**3
    F = -sqrt(2)*I*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(64*a**3*f) + I*c**4*sqrt(-I*c*tan(e + f*x) + c)/(3*a**3*f*(I*c*tan(e + f*x) + c)**3) - I*c**3*sqrt(-I*c*tan(e + f*x) + c)/(24*a**3*f*(I*c*tan(e + f*x) + c)**2) - I*c**2*sqrt(-I*c*tan(e + f*x) + c)/(32*a**3*f*(I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_964():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 8*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f) - 8*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c*f) + 2*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(9)/2)/(9*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_965():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 4*I*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f) - 2*I*a**2*(-I*c*tan(e + f*x) + c)**(sympy.S(7)/2)/(7*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_966():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*I*a*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_967():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)
    F = -3*sqrt(2)*I*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(a*f) + I*c**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(a*f*(I*c*tan(e + f*x) + c)) + 3*I*c**2*sqrt(-I*c*tan(e + f*x) + c)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_968():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**2
    F = 3*sqrt(2)*I*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(8*a**2*f) + I*c**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(2*a**2*f*(I*c*tan(e + f*x) + c)**2) - 3*I*c**3*sqrt(-I*c*tan(e + f*x) + c)/(4*a**2*f*(I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_969():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**3
    F = sqrt(2)*I*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a**3*f) + I*c**4*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*a**3*f*(I*c*tan(e + f*x) + c)**3) - I*c**4*sqrt(-I*c*tan(e + f*x) + c)/(4*a**3*f*(I*c*tan(e + f*x) + c)**2) + I*c**3*sqrt(-I*c*tan(e + f*x) + c)/(16*a**3*f*(I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_970():
    f = (I*a*tan(e + f*x) + a)**3/sqrt(-I*c*tan(e + f*x) + c)
    F = -8*I*a**3/(f*sqrt(-I*c*tan(e + f*x) + c)) - 8*I*a**3*sqrt(-I*c*tan(e + f*x) + c)/(c*f) + 2*I*a**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_971():
    f = (I*a*tan(e + f*x) + a)**2/sqrt(-I*c*tan(e + f*x) + c)
    F = -4*I*a**2/(f*sqrt(-I*c*tan(e + f*x) + c)) - 2*I*a**2*sqrt(-I*c*tan(e + f*x) + c)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_972():
    f = (I*a*tan(e + f*x) + a)/sqrt(-I*c*tan(e + f*x) + c)
    F = -2*I*a/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_973():
    f = 1/((I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    F = -3*I/(4*a*f*sqrt(-I*c*tan(e + f*x) + c)) + I/(2*a*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) + 3*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(8*a*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_974():
    f = 1/((I*a*tan(e + f*x) + a)**2*sqrt(-I*c*tan(e + f*x) + c))
    F = -15*I/(32*a**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 5*I/(16*a**2*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) + I/(4*a**2*f*(I*tan(e + f*x) + 1)**2*sqrt(-I*c*tan(e + f*x) + c)) + 15*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(64*a**2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_975():
    f = 1/((I*a*tan(e + f*x) + a)**3*sqrt(-I*c*tan(e + f*x) + c))
    F = -35*I/(128*a**3*f*sqrt(-I*c*tan(e + f*x) + c)) + 35*I/(192*a**3*f*(I*tan(e + f*x) + 1)*sqrt(-I*c*tan(e + f*x) + c)) + 7*I/(48*a**3*f*(I*tan(e + f*x) + 1)**2*sqrt(-I*c*tan(e + f*x) + c)) + I/(6*a**3*f*(I*tan(e + f*x) + 1)**3*sqrt(-I*c*tan(e + f*x) + c)) + 35*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(256*a**3*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_976():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -8*I*a**3/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*I*a**3/(c*f*sqrt(-I*c*tan(e + f*x) + c)) + 2*I*a**3*sqrt(-I*c*tan(e + f*x) + c)/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_977():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -4*I*a**2/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 2*I*a**2/(c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_978():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*I*a/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_979():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = -5*I/(12*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + I/(2*a*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 5*I/(8*a*c*f*sqrt(-I*c*tan(e + f*x) + c)) + 5*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(16*a*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_980():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = -35*I/(96*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 7*I/(16*a**2*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + I/(4*a**2*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 35*I/(64*a**2*c*f*sqrt(-I*c*tan(e + f*x) + c)) + 35*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(128*a**2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_981():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = -35*I/(128*a**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 21*I/(64*a**3*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 3*I/(16*a**3*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + I/(6*a**3*f*(I*tan(e + f*x) + 1)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 105*I/(256*a**3*c*f*sqrt(-I*c*tan(e + f*x) + c)) + 105*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(512*a**3*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_982():
    f = (I*a*tan(e + f*x) + a)**3/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -8*I*a**3/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 8*I*a**3/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*a**3/(c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_983():
    f = (I*a*tan(e + f*x) + a)**2/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -4*I*a**2/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 2*I*a**2/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_984():
    f = (I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*I*a/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_985():
    f = 1/((I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = -7*I/(20*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + I/(2*a*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 7*I/(24*a*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 7*I/(16*a*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 7*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(32*a*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_986():
    f = 1/((I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = -63*I/(160*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 9*I/(16*a**2*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + I/(4*a**2*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 21*I/(64*a**2*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 63*I/(128*a**2*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 63*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(256*a**2*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_987():
    f = 1/((I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = -231*I/(640*a**3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 33*I/(64*a**3*f*(I*tan(e + f*x) + 1)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 11*I/(48*a**3*f*(I*tan(e + f*x) + 1)**2*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + I/(6*a**3*f*(I*tan(e + f*x) + 1)**3*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 77*I/(256*a**3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 231*I/(512*a**3*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 231*sqrt(2)*I*atanh(sqrt(2)*sqrt(-I*c*tan(e + f*x) + c)/(2*sqrt(c)))/(1024*a**3*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_988():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)
    F = -3*I*a**(sympy.S(5)/2)*sqrt(c)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + 3*I*a**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f) + I*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_989():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)
    F = -2*I*a**(sympy.S(3)/2)*sqrt(c)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + I*a*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_990():
    f = sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)
    F = -2*I*sqrt(a)*sqrt(c)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_991():
    f = sqrt(-I*c*tan(e + f*x) + c)/sqrt(I*a*tan(e + f*x) + a)
    F = I*sqrt(-I*c*tan(e + f*x) + c)/(f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_992():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = I*sqrt(-I*c*tan(e + f*x) + c)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + I*sqrt(-I*c*tan(e + f*x) + c)/(3*a*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_993():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = I*sqrt(-I*c*tan(e + f*x) + c)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + 2*I*sqrt(-I*c*tan(e + f*x) + c)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + 2*I*sqrt(-I*c*tan(e + f*x) + c)/(15*a**2*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_994():
    f = sqrt(-I*c*tan(e + f*x) + c)/(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)
    F = I*sqrt(-I*c*tan(e + f*x) + c)/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + 3*I*sqrt(-I*c*tan(e + f*x) + c)/(35*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + 2*I*sqrt(-I*c*tan(e + f*x) + c)/(35*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + 2*I*sqrt(-I*c*tan(e + f*x) + c)/(35*a**3*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_995():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -I*a**(sympy.S(5)/2)*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + a**2*c*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(2*f) + I*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_996():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -I*a**(sympy.S(3)/2)*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + a*c*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_997():
    f = sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*I*sqrt(a)*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f - I*c*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_998():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = 2*I*c*sqrt(-I*c*tan(e + f*x) + c)/(f*sqrt(I*a*tan(e + f*x) + a)) + 2*I*c**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_999():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1000():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1001():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(35*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(105*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1002():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(9*f*(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)) + I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(21*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(105*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(315*a**3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1003():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -3*I*a**(sympy.S(5)/2)*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(4*f) + 3*a**2*c**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(8*f) + a*c*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1004():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -I*a**(sympy.S(3)/2)*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f + a*c**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)*tan(e + f*x)/(2*f) - I*c*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1005():
    f = sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -3*I*sqrt(a)*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/f - 3*I*c**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*f) - I*c*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1006():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = 2*I*c*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(f*sqrt(I*a*tan(e + f*x) + a)) + 3*I*c**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(a*f) + 6*I*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1007():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*I*c*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) - 2*I*c**2*sqrt(-I*c*tan(e + f*x) + c)/(a*f*sqrt(I*a*tan(e + f*x) + a)) - 2*I*c**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1008():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1009():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(35*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1010():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(9*f*(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(63*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(315*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1011():
    f = (-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(11)/2)
    F = I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(11*f*(I*a*tan(e + f*x) + a)**(sympy.S(11)/2)) + I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(33*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(231*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)) + 2*I*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)/(1155*a**3*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1012():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 15*I*a**(sympy.S(7)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - 15*I*a**3*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c*f) - 5*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(2*c*f) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1013():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 6*I*a**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - 3*I*a**2*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c*f) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1014():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/sqrt(-I*c*tan(e + f*x) + c)
    F = 2*I*a**(sympy.S(3)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(sqrt(c)*f) - 2*I*a*sqrt(I*a*tan(e + f*x) + a)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1015():
    f = sqrt(I*a*tan(e + f*x) + a)/sqrt(-I*c*tan(e + f*x) + c)
    F = -I*sqrt(I*a*tan(e + f*x) + a)/(f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1016():
    f = 1/(sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    F = tan(e + f*x)/(f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1017():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c))
    F = I/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 2*tan(e + f*x)/(3*a*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1018():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c))
    F = I/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)) + I/(5*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 2*tan(e + f*x)/(5*a**2*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1019():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-I*c*tan(e + f*x) + c))
    F = I/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 4*I/(35*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 4*I/(35*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)) + 8*tan(e + f*x)/(35*a**3*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1020():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(9)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -35*I*a**(sympy.S(9)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + 35*I*a**4*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c**2*f) + 35*I*a**3*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(6*c**2*f) + 14*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1021():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -10*I*a**(sympy.S(7)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + 5*I*a**3*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c**2*f) + 10*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1022():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*I*a**(sympy.S(5)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) + 2*I*a**2*sqrt(I*a*tan(e + f*x) + a)/(c*f*sqrt(-I*c*tan(e + f*x) + c)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1023():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -I*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1024():
    f = sqrt(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)
    F = -I*sqrt(I*a*tan(e + f*x) + a)/(3*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - I*sqrt(I*a*tan(e + f*x) + a)/(3*c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1025():
    f = 1/(sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = I/(f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(3*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(3*a*c*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1026():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = tan(e + f*x)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 2*tan(e + f*x)/(3*a*c*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1027():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = I/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 4*tan(e + f*x)/(15*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*tan(e + f*x)/(15*a**2*c*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1028():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    F = I/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + I/(7*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 4*tan(e + f*x)/(21*a**2*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*tan(e + f*x)/(21*a**3*c*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1029():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(11)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 63*I*a**(sympy.S(11)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(5)/2)*f) - 63*I*a**5*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(2*c**3*f) - 21*I*a**4*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-I*c*tan(e + f*x) + c)/(2*c**3*f) - 42*I*a**3*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(5*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 6*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(5*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(9)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1030():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(9)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 14*I*a**(sympy.S(9)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(5)/2)*f) - 7*I*a**4*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c)/(c**3*f) - 14*I*a**3*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 14*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1031():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(7)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = 2*I*a**(sympy.S(7)/2)*atan(sqrt(c)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(-I*c*tan(e + f*x) + c)))/(c**(sympy.S(5)/2)*f) - 2*I*a**3*sqrt(I*a*tan(e + f*x) + a)/(c**2*f*sqrt(-I*c*tan(e + f*x) + c)) + 2*I*a**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(3*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*a*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1032():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -I*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1033():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -I*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - I*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1034():
    f = sqrt(I*a*tan(e + f*x) + a)/(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)
    F = -I*sqrt(I*a*tan(e + f*x) + a)/(5*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(15*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(15*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1035():
    f = 1/(sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = I/(f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 3*I*sqrt(I*a*tan(e + f*x) + a)/(5*a*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(5*a*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 2*I*sqrt(I*a*tan(e + f*x) + a)/(5*a*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1036():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = I/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 4*I/(3*a*f*sqrt(I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 4*I*sqrt(I*a*tan(e + f*x) + a)/(5*a**2*f*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) - 8*I*sqrt(I*a*tan(e + f*x) + a)/(15*a**2*c*f*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) - 8*I*sqrt(I*a*tan(e + f*x) + a)/(15*a**2*c**2*f*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1037():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = tan(e + f*x)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 4*tan(e + f*x)/(15*a*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 8*tan(e + f*x)/(15*a**2*c**2*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1038():
    f = 1/((I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2))
    F = I/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 6*tan(e + f*x)/(35*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(5)/2)) + 8*tan(e + f*x)/(35*a**2*c*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)*(-I*c*tan(e + f*x) + c)**(sympy.S(3)/2)) + 16*tan(e + f*x)/(35*a**3*c**2*f*sqrt(I*a*tan(e + f*x) + a)*sqrt(-I*c*tan(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1039():
    f = (I*a*tan(e + f*x) + a)**4*(-I*c*tan(e + f*x) + c)**n
    F = 8*I*a**4*(-I*c*tan(e + f*x) + c)**n/(f*n) - 12*I*a**4*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1)) + 6*I*a**4*(-I*c*tan(e + f*x) + c)**(n + 2)/(c**2*f*(n + 2)) - I*a**4*(-I*c*tan(e + f*x) + c)**(n + 3)/(c**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1040():
    f = (I*a*tan(e + f*x) + a)**3*(-I*c*tan(e + f*x) + c)**n
    F = 4*I*a**3*(-I*c*tan(e + f*x) + c)**n/(f*n) - 4*I*a**3*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1)) + I*a**3*(-I*c*tan(e + f*x) + c)**(n + 2)/(c**2*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1041():
    f = (I*a*tan(e + f*x) + a)**2*(-I*c*tan(e + f*x) + c)**n
    F = 2*I*a**2*(-I*c*tan(e + f*x) + c)**n/(f*n) - I*a**2*(-I*c*tan(e + f*x) + c)**(n + 1)/(c*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1042():
    f = (I*a*tan(e + f*x) + a)*(-I*c*tan(e + f*x) + c)**n
    F = I*a*(-I*c*tan(e + f*x) + c)**n/(f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1043():
    f = (-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)
    F = I*(-I*c*tan(e + f*x) + c)**n*hyper((2, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(4*a*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1044():
    f = (-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)**2
    F = I*(-I*c*tan(e + f*x) + c)**n*hyper((3, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(8*a**2*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1045():
    f = (-I*c*tan(e + f*x) + c)**n/(I*a*tan(e + f*x) + a)**3
    F = I*(-I*c*tan(e + f*x) + c)**n*hyper((4, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(16*a**3*f*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1046():
    f = (I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**4
    F = -8*I*c**4*(I*a*tan(e + f*x) + a)**m/(f*m) + 12*I*c**4*(I*a*tan(e + f*x) + a)**(m + 1)/(a*f*(m + 1)) - 6*I*c**4*(I*a*tan(e + f*x) + a)**(m + 2)/(a**2*f*(m + 2)) + I*c**4*(I*a*tan(e + f*x) + a)**(m + 3)/(a**3*f*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1047():
    f = (I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**3
    F = -4*I*c**3*(I*a*tan(e + f*x) + a)**m/(f*m) + 4*I*c**3*(I*a*tan(e + f*x) + a)**(m + 1)/(a*f*(m + 1)) - I*c**3*(I*a*tan(e + f*x) + a)**(m + 2)/(a**2*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1048():
    f = (I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)**2
    F = -2*I*c**2*(I*a*tan(e + f*x) + a)**m/(f*m) + I*c**2*(I*a*tan(e + f*x) + a)**(m + 1)/(a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1049():
    f = (I*a*tan(e + f*x) + a)**m*(-I*c*tan(e + f*x) + c)
    F = -I*c*(I*a*tan(e + f*x) + a)**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1050():
    f = (I*a*tan(e + f*x) + a)**m/(-I*c*tan(e + f*x) + c)
    F = -I*(I*a*tan(e + f*x) + a)**m*hyper((2, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(4*c*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1051():
    f = (I*a*tan(e + f*x) + a)**m/(-I*c*tan(e + f*x) + c)**2
    F = -I*(I*a*tan(e + f*x) + a)**m*hyper((3, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(8*c**2*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1052():
    f = (I*a*tan(e + f*x) + a)**m/(-I*c*tan(e + f*x) + c)**3
    F = -I*(I*a*tan(e + f*x) + a)**m*hyper((4, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(16*c**3*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1053():
    f = (I*a*tan(e + f*x) + a)**m/(-I*c*tan(e + f*x) + c)**4
    F = -I*(I*a*tan(e + f*x) + a)**m*hyper((5, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(32*c**4*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1054():
    f = (c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3
    F = 4*a**3*x*(c - I*d) - 2*a**3*(c - I*d)*tan(e + f*x)/f - 4*a**3*(I*c + d)*log(cos(e + f*x))/f + a*(I*c + d)*(I*a*tan(e + f*x) + a)**2/(2*f) + d*(I*a*tan(e + f*x) + a)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1055():
    f = (c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2
    F = 2*a**2*x*(c - I*d) - a**2*(c - I*d)*tan(e + f*x)/f - 2*a**2*(I*c + d)*log(cos(e + f*x))/f + d*(I*a*tan(e + f*x) + a)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1056():
    f = (c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)
    F = I*a*d*tan(e + f*x)/f + a*x*(c - I*d) - a*(I*c + d)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1057():
    f = (c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)
    F = (I*c - d)/(2*f*(I*a*tan(e + f*x) + a)) + x*(c - I*d)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1058():
    f = (c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**2
    F = (I*c - d)/(4*f*(I*a*tan(e + f*x) + a)**2) + (I*c + d)/(4*f*(I*a**2*tan(e + f*x) + a**2)) + x*(c - I*d)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1059():
    f = (c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**3
    F = (I*c - d)/(6*f*(I*a*tan(e + f*x) + a)**3) + (I*c + d)/(8*f*(I*a**3*tan(e + f*x) + a**3)) + (I*c + d)/(8*a*f*(I*a*tan(e + f*x) + a)**2) + x*(c - I*d)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1060():
    f = (c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**3
    F = 4*a**3*x*(c - I*d)**2 - 4*I*a**3*(c - I*d)**2*log(cos(e + f*x))/f - 2*a**3*(c - I*d)**2*tan(e + f*x)/f + I*a*(c - I*d)**2*(I*a*tan(e + f*x) + a)**2/(2*f) + 2*c*d*(I*a*tan(e + f*x) + a)**3/(3*f) - I*d**2*(I*a*tan(e + f*x) + a)**4/(4*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1061():
    f = (c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**2
    F = 2*a**2*x*(c - I*d)**2 - 2*I*a**2*(c - I*d)**2*log(cos(e + f*x))/f - a**2*(c - I*d)**2*tan(e + f*x)/f + c*d*(I*a*tan(e + f*x) + a)**2/f - I*d**2*(I*a*tan(e + f*x) + a)**3/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1062():
    f = (c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)
    F = a*d*(I*c + d)*tan(e + f*x)/f + a*x*(c - I*d)**2 - I*a*(c - I*d)**2*log(cos(e + f*x))/f + I*a*(c + d*tan(e + f*x))**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1063():
    f = (c + d*tan(e + f*x))**2/(I*a*tan(e + f*x) + a)
    F = I*(c + I*d)**2/(2*f*(I*a*tan(e + f*x) + a)) + I*d**2*log(cos(e + f*x))/(a*f) + x*(c**2 - 2*I*c*d + d**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1064():
    f = (c + d*tan(e + f*x))**2/(I*a*tan(e + f*x) + a)**2
    F = I*(c + I*d)**2/(4*f*(I*a*tan(e + f*x) + a)**2) + x*(c - I*d)**2/(4*a**2) + (c + I*d)*(I*c + 3*d)/(4*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1065():
    f = (c + d*tan(e + f*x))**2/(I*a*tan(e + f*x) + a)**3
    F = I*(c - I*d)**2/(8*f*(I*a**3*tan(e + f*x) + a**3)) + I*(c + I*d)**2/(6*f*(I*a*tan(e + f*x) + a)**3) + (c + I*d)*(I*c + 3*d)/(8*a*f*(I*a*tan(e + f*x) + a)**2) + x*(c - I*d)**2/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1066():
    f = (c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)**3
    F = 4*I*a**3*d*(c - I*d)**2*tan(e + f*x)/f + 4*a**3*x*(c - I*d)**3 + 4*I*a**3*(c + d*tan(e + f*x))**3/(3*f) + 2*a**3*(c + d*tan(e + f*x))**2*(I*c + d)/f + 4*a**3*(I*c + d)**3*log(cos(e + f*x))/f + a**3*(c + d*tan(e + f*x))**4*(I*c - 11*d)/(20*d**2*f) - (c + d*tan(e + f*x))**4*(I*a**3*tan(e + f*x) + a**3)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1067():
    f = (c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)**2
    F = 2*I*a**2*d*(c - I*d)**2*tan(e + f*x)/f + 2*a**2*x*(c - I*d)**3 + 2*I*a**2*(c + d*tan(e + f*x))**3/(3*f) + a**2*(c + d*tan(e + f*x))**2*(I*c + d)/f + 2*a**2*(I*c + d)**3*log(cos(e + f*x))/f - a**2*(c + d*tan(e + f*x))**4/(4*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1068():
    f = (c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)
    F = I*a*d*(c - I*d)**2*tan(e + f*x)/f + a*x*(c - I*d)**3 + I*a*(c + d*tan(e + f*x))**3/(3*f) + a*(c + d*tan(e + f*x))**2*(I*c + d)/(2*f) + a*(I*c + d)**3*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1069():
    f = (c + d*tan(e + f*x))**3/(I*a*tan(e + f*x) + a)
    F = (c + d*tan(e + f*x))**2*(I*c - d)/(2*f*(I*a*tan(e + f*x) + a)) - d**2*(c + 3*I*d)*tan(e + f*x)/(2*a*f) + d**2*(3*I*c - d)*log(cos(e + f*x))/(a*f) + x*(c**3 - 3*I*c**2*d + 3*c*d**2 + 3*I*d**3)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1070():
    f = (c + d*tan(e + f*x))**3/(I*a*tan(e + f*x) + a)**2
    F = (c + d*tan(e + f*x))**2*(I*c - d)/(4*f*(I*a*tan(e + f*x) + a)**2) + d**3*log(cos(e + f*x))/(a**2*f) + x*(c**3 - 3*I*c**2*d - 3*c*d**2 - 3*I*d**3)/(4*a**2) + (c + I*d)**2*(I*c + 3*d)/(4*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1071():
    f = (c + d*tan(e + f*x))**3/(I*a*tan(e + f*x) + a)**3
    F = I*(c + d*tan(e + f*x))**3/(6*f*(I*a*tan(e + f*x) + a)**3) + (c + I*d)**2*(I*c + d)/(8*a*f*(I*a*tan(e + f*x) + a)**2) + x*(c - I*d)**3/(8*a**3) + (c - 3*I*d)*(c + I*d)*(I*c + d)/(8*a**3*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1072():
    f = (I*a*tan(e + f*x) + a)**3/(c + d*tan(e + f*x))
    F = 4*a**3*x/(c - I*d) - a**3*(c + I*d)**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(d**2*f*(I*c + d)) - a**3*(I*c - 3*d)*log(cos(e + f*x))/(d**2*f) - (I*a**3*tan(e + f*x) + a**3)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1073():
    f = (I*a*tan(e + f*x) + a)**2/(c + d*tan(e + f*x))
    F = -a**2*c*x*(c + I*d)/(d**2*(c - I*d)) - a**2*(I*c - d)*log(c*cos(e + f*x) + d*sin(e + f*x))/(d*f*(I*c + d)) + a**2*log(cos(e + f*x))/(d*f) + a**2*x*(c + 2*I*d)/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1074():
    f = (I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))
    F = a*x/(c - I*d) + a*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1075():
    f = 1/((c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a))
    F = -1/(f*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) - c*d*x/(a*(c**2 + d**2)*(I*c - d)) - d**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(a*f*(c**2 + d**2)*(I*c - d)) + x/(2*a*(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1076():
    f = 1/((c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2)
    F = -1/(f*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) - d**3*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**2*f*(c - I*d)*(c + I*d)**3) + x*(c**3 + 3*I*c**2*d - 3*c*d**2 + 3*I*d**3)/(4*a**2*(c - I*d)*(c + I*d)**3) + (I*c - 3*d)/(4*a**2*f*(c + I*d)**2*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1077():
    f = 1/((c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3)
    F = -1/(f*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (c**2 + 4*I*c*d - 7*d**2)/(8*f*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + (I*c - 3*d)/(8*a*f*(c + I*d)**2*(I*a*tan(e + f*x) + a)**2) + d**4*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**3*f*(c + I*d)**4*(I*c + d)) + x*(c**4 + 4*I*c**3*d - 6*c**2*d**2 - 4*I*c*d**3 - 7*d**4)/(8*a**3*(c - I*d)*(c + I*d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1078():
    f = (I*a*tan(e + f*x) + a)**3/(c + d*tan(e + f*x))**2
    F = 4*a**3*x/(c - I*d)**2 - a**3*(c - 3*I*d)*(I*c - d)*log(c*cos(e + f*x) + d*sin(e + f*x))/(d**2*f*(c - I*d)**2) + I*a**3*log(cos(e + f*x))/(d**2*f) + (c + I*d)*(I*a**3*tan(e + f*x) + a**3)/(d*f*(c - I*d)*(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1079():
    f = (I*a*tan(e + f*x) + a)**2/(c + d*tan(e + f*x))**2
    F = 2*a**2*x/(c - I*d)**2 - 2*I*a**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c - I*d)**2) + a**2*(I*c - d)/(d*f*(c + d*tan(e + f*x))*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1080():
    f = (I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**2
    F = a*x/(c - I*d)**2 - a/(f*(c + d*tan(e + f*x))*(I*c + d)) - I*a*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c - I*d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1081():
    f = 1/((c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a))
    F = -1/(f*(c + d*tan(e + f*x))*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + d**2*(3*c - I*d)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a*f*(c - I*d)**2*(I*c - d)**3) + d*(c - 3*I*d)/(2*a*f*(c - I*d)*(c + I*d)**2*(c + d*tan(e + f*x))) + x*(c**3 + 3*I*c**2*d + 3*c*d**2 - 3*I*d**3)/(2*a*(c - I*d)**2*(c + I*d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1082():
    f = 1/((c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**2)
    F = -1/(f*(c + d*tan(e + f*x))*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) - d**3*(4*c - 2*I*d)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**2*f*(c - I*d)**2*(c + I*d)**4) + d*(c**2 + 4*I*c*d + 9*d**2)/(4*a**2*f*(c - I*d)*(c + I*d)**3*(c + d*tan(e + f*x))) + x*(c**4 + 4*I*c**3*d - 6*c**2*d**2 + 12*I*c*d**3 + 9*d**4)/(4*a**2*(c - I*d)**2*(c + I*d)**4) + (I*c - 4*d)/(4*a**2*f*(c + I*d)**2*(c + d*tan(e + f*x))*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1083():
    f = 1/((c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**3)
    F = -1/(f*(c + d*tan(e + f*x))*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (c**2 + 5*I*c*d - 12*d**2)/(8*f*(c + d*tan(e + f*x))*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + (3*I*c - 11*d)/(24*a*f*(c + I*d)**2*(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2) + d**4*(5*c - 3*I*d)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**3*f*(c - I*d)**2*(I*c - d)**5) + d*(c**3 + 5*I*c**2*d - 11*c*d**2 + 25*I*d**3)/(8*a**3*f*(c - I*d)*(c + I*d)**4*(c + d*tan(e + f*x))) + x*(c**5 + 5*I*c**4*d - 10*c**3*d**2 - 10*I*c**2*d**3 - 35*c*d**4 + 25*I*d**5)/(8*a**3*(c - I*d)**2*(c + I*d)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1084():
    f = (I*a*tan(e + f*x) + a)**3/(c + d*tan(e + f*x))**3
    F = 4*a**3*x/(c - I*d)**3 - 4*a**3*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(I*c + d)**3) + 2*a**3*(c + I*d)/(d*f*(c - I*d)**2*(c + d*tan(e + f*x))) - a*(I*a*tan(e + f*x) + a)**2/(f*(c + d*tan(e + f*x))**2*(2*I*c + 2*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1085():
    f = (I*a*tan(e + f*x) + a)**2/(c + d*tan(e + f*x))**3
    F = 2*a**2*x/(c - I*d)**3 - 2*a**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(I*c + d)**3) + 2*I*a**2/(f*(c - I*d)**2*(c + d*tan(e + f*x))) + a**2*(I*c - d)/(2*d*f*(c + d*tan(e + f*x))**2*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1086():
    f = (I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**3
    F = a*x/(c - I*d)**3 - a*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(I*c + d)**3) - a/(f*(c + d*tan(e + f*x))**2*(2*I*c + 2*d)) + I*a/(f*(c - I*d)**2*(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1087():
    f = 1/((c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a))
    F = -1/(f*(c + d*tan(e + f*x))**2*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + 2*d**2*(3*c**2 - 2*I*c*d - d**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a*f*(c + I*d)**4*(I*c + d)**3) + d*(c - 2*I*d)/(2*a*f*(c - I*d)*(c + I*d)**2*(c + d*tan(e + f*x))**2) + d*(c**2 - 8*I*c*d - 3*d**2)/(2*a*f*(c - I*d)**2*(c + I*d)**3*(c + d*tan(e + f*x))) + x*(c**4 + 4*I*c**3*d + 6*c**2*d**2 - 12*I*c*d**3 - 3*d**4)/(2*a*(c - I*d)**3*(c + I*d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1088():
    f = 1/((c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)**2)
    F = -1/(f*(c + d*tan(e + f*x))**2*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) - 2*d**3*(5*c**2 - 5*I*c*d - 2*d**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**2*f*(I*c - d)**5*(I*c + d)**3) + d*(c - 3*I*d)*(c**2 + 8*I*c*d + 5*d**2)/(4*a**2*f*(c - I*d)**2*(c + I*d)**4*(c + d*tan(e + f*x))) + d*(c**2 + 5*I*c*d + 8*d**2)/(4*a**2*f*(c - I*d)*(c + I*d)**3*(c + d*tan(e + f*x))**2) + x*(c**5 + 5*I*c**4*d - 10*c**3*d**2 + 30*I*c**2*d**3 + 45*c*d**4 - 15*I*d**5)/(4*a**2*(c - I*d)**3*(c + I*d)**5) + (I*c - 5*d)/(4*a**2*f*(c + I*d)**2*(c + d*tan(e + f*x))**2*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1089():
    f = 1/((c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)**3)
    F = -1/(f*(c + d*tan(e + f*x))**2*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (3*c**2 + 18*I*c*d - 55*d**2)/(24*f*(c + d*tan(e + f*x))**2*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + (3*I*c - 13*d)/(24*a*f*(c + I*d)**2*(c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**2) - d**4*(15*c**2 - 18*I*c*d - 7*d**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(a**3*f*(c + I*d)**6*(I*c + d)**3) + d*(c**3 + 6*I*c**2*d - 17*c*d**2 + 28*I*d**3)/(8*a**3*f*(c - I*d)*(c + I*d)**4*(c + d*tan(e + f*x))**2) + d*(c**4 + 6*I*c**3*d - 16*c**2*d**2 + 94*I*c*d**3 + 55*d**4)/(8*a**3*f*(c - I*d)**2*(c + I*d)**5*(c + d*tan(e + f*x))) + x*(c**6 + 6*I*c**5*d - 15*c**4*d**2 - 20*I*c**3*d**3 - 105*c**2*d**4 + 150*I*c*d**5 + 55*d**6)/(8*a**3*(c - I*d)**3*(c + I*d)**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1090():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3
    F = -8*I*a**3*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 8*I*a**3*sqrt(c + d*tan(e + f*x))/f + 4*a**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - 6*d)/(15*d**2*f) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*I*a**3*tan(e + f*x) + 2*a**3)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1091():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2
    F = -4*I*a**2*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 4*I*a**2*sqrt(c + d*tan(e + f*x))/f - 2*a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1092():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)
    F = -2*I*a*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 2*I*a*sqrt(c + d*tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1093():
    f = sqrt(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)
    F = I*sqrt(c + d*tan(e + f*x))/(2*f*(I*a*tan(e + f*x) + a)) + I*c*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f*sqrt(c + I*d)) - I*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1094():
    f = sqrt(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**2
    F = I*sqrt(c + d*tan(e + f*x))/(4*f*(I*a*tan(e + f*x) + a)**2) - I*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f) + sqrt(c + d*tan(e + f*x))*(2*I*c - d)/(8*a**2*f*(c + I*d)*(I*tan(e + f*x) + 1)) - (2*c*d - I*(2*c**2 + d**2))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f*(c + I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1095():
    f = sqrt(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**3
    F = c*sqrt(c + d*tan(e + f*x))*(2*I*c - 3*d)/(16*f*(c + I*d)**2*(I*a**3*tan(e + f*x) + a**3)) + I*sqrt(c + d*tan(e + f*x))/(6*f*(I*a*tan(e + f*x) + a)**3) + sqrt(c + d*tan(e + f*x))*(3*I*c - 2*d)/(24*a*f*(c + I*d)*(I*a*tan(e + f*x) + a)**2) - I*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f) + (2*I*c**3 - 4*c**2*d - I*c*d**2 - 2*d**3)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*(c + I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1096():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**3
    F = -8*I*a**3*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 8*I*a**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 8*a**3*sqrt(c + d*tan(e + f*x))*(I*c + d)/f + 4*a**3*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*c - 8*d)/(35*d**2*f) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*I*a**3*tan(e + f*x) + 2*a**3)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1097():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2
    F = -4*I*a**2*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 4*I*a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 4*a**2*sqrt(c + d*tan(e + f*x))*(I*c + d)/f - 2*a**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1098():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)
    F = -2*I*a*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 2*I*a*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*a*sqrt(c + d*tan(e + f*x))*(I*c + d)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1099():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)
    F = sqrt(c + d*tan(e + f*x))*(I*c - d)/(2*f*(I*a*tan(e + f*x) + a)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f) + sqrt(c + I*d)*(I*c + 2*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1100():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**2
    F = sqrt(c + d*tan(e + f*x))*(I*c - d)/(4*f*(I*a*tan(e + f*x) + a)**2) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f) + sqrt(c + d*tan(e + f*x))*(2*I*c + 3*d)/(8*a**2*f*(I*tan(e + f*x) + 1)) + (2*c*d + I*(2*c**2 + d**2))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1101():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**3
    F = sqrt(c + d*tan(e + f*x))*(I*c - d)/(6*f*(I*a*tan(e + f*x) + a)**3) - sqrt(c + d*tan(e + f*x))*(2*c**2 - I*c*d + 2*d**2)/(f*(16*I*c - 16*d)*(I*a**3*tan(e + f*x) + a**3)) + sqrt(c + d*tan(e + f*x))*(3*I*c + 4*d)/(24*a*f*(I*a*tan(e + f*x) + a)**2) + I*c*(2*c**2 + 3*d**2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*(c + I*d)**(sympy.S(3)/2)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1102():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**3
    F = -8*I*a**3*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 8*I*a**3*(c - I*d)**2*sqrt(c + d*tan(e + f*x))/f + 8*I*a**3*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 8*a**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c + d)/(3*f) + 4*a**3*(c + d*tan(e + f*x))**(sympy.S(7)/2)*(I*c - 10*d)/(63*d**2*f) - (c + d*tan(e + f*x))**(sympy.S(7)/2)*(2*I*a**3*tan(e + f*x) + 2*a**3)/(9*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1103():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**2
    F = -4*I*a**2*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 4*I*a**2*(c - I*d)**2*sqrt(c + d*tan(e + f*x))/f + 4*I*a**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 4*a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c + d)/(3*f) - 2*a**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1104():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)
    F = -2*I*a*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + 2*I*a*(c - I*d)**2*sqrt(c + d*tan(e + f*x))/f + 2*I*a*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 2*a*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c + d)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1105():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)
    F = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)/(2*f*(I*a*tan(e + f*x) + a)) - d*(c + 5*I*d)*sqrt(c + d*tan(e + f*x))/(2*a*f) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f) + (c + I*d)**(sympy.S(3)/2)*(I*c + 4*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1106():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**2
    F = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)/(4*f*(I*a*tan(e + f*x) + a)**2) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f) + sqrt(c + I*d)*(2*I*c**2 + 6*c*d - 7*I*d**2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f) + (c + I*d)*sqrt(c + d*tan(e + f*x))*(2*I*c + 5*d)/(8*a**2*f*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1107():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**3
    F = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)/(6*f*(I*a*tan(e + f*x) + a)**3) + sqrt(c + d*tan(e + f*x))*(2*I*c**2 + 5*c*d - 4*I*d**2)/(16*f*(I*a**3*tan(e + f*x) + a**3)) + (c + I*d)*sqrt(c + d*tan(e + f*x))*(I*c + 2*d)/(8*a*f*(I*a*tan(e + f*x) + a)**2) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f) + (2*I*c**3 + 4*c**2*d - I*c*d**2 + 2*d**3)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1108():
    f = (I*a*tan(e + f*x) + a)**3/sqrt(c + d*tan(e + f*x))
    F = -8*I*a**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) + 4*a**3*sqrt(c + d*tan(e + f*x))*(I*c - 4*d)/(3*d**2*f) - sqrt(c + d*tan(e + f*x))*(2*I*a**3*tan(e + f*x) + 2*a**3)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1109():
    f = (I*a*tan(e + f*x) + a)**2/sqrt(c + d*tan(e + f*x))
    F = -4*I*a**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) - 2*a**2*sqrt(c + d*tan(e + f*x))/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1110():
    f = (I*a*tan(e + f*x) + a)/sqrt(c + d*tan(e + f*x))
    F = -2*I*a*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1111():
    f = 1/(sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a))
    F = -sqrt(c + d*tan(e + f*x))/(f*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + (I*c - 2*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f*(c + I*d)**(sympy.S(3)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1112():
    f = 1/(sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2)
    F = -sqrt(c + d*tan(e + f*x))/(f*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) + sqrt(c + d*tan(e + f*x))*(2*I*c - 5*d)/(8*a**2*f*(c + I*d)**2*(I*tan(e + f*x) + 1)) + (2*I*c**2 - 6*c*d - 7*I*d**2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1113():
    f = 1/(sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**3)
    F = -sqrt(c + d*tan(e + f*x))/(f*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + sqrt(c + d*tan(e + f*x))*(2*c**2 + 7*I*c*d - 10*d**2)/(16*f*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + sqrt(c + d*tan(e + f*x))*(3*I*c - 8*d)/(24*a*f*(c + I*d)**2*(I*a*tan(e + f*x) + a)**2) + (2*I*c**3 - 8*c**2*d - 13*I*c*d**2 + 12*d**3)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*(c + I*d)**(sympy.S(7)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1114():
    f = (I*a*tan(e + f*x) + a)**3/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 4*a**3*c*sqrt(c + d*tan(e + f*x))/(d**2*f*(I*c + d)) - 8*I*a**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) + (2*c + 2*I*d)*(I*a**3*tan(e + f*x) + a**3)/(d*f*(c - I*d)*sqrt(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1115():
    f = (I*a*tan(e + f*x) + a)**2/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -4*I*a**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) + 2*a**2*(I*c - d)/(d*f*sqrt(c + d*tan(e + f*x))*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1116():
    f = (I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*a/(f*sqrt(c + d*tan(e + f*x))*(I*c + d)) - 2*I*a*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1117():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a))
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + d*(c - 5*I*d)/(2*a*f*(c - I*d)*(c + I*d)**2*sqrt(c + d*tan(e + f*x))) + (I*c - 4*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1118():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2)
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) + d*(2*c**2 + 7*I*c*d + 25*d**2)/(8*a**2*f*(c - I*d)*(c + I*d)**3*sqrt(c + d*tan(e + f*x))) + (2*I*c - 7*d)/(8*a**2*f*(c + I*d)**2*sqrt(c + d*tan(e + f*x))*(I*tan(e + f*x) + 1)) + (2*I*c**2 - 10*c*d - 23*I*d**2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f*(c + I*d)**(sympy.S(7)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1119():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**3)
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (6*c**2 + 27*I*c*d - 56*d**2)/(48*f*sqrt(c + d*tan(e + f*x))*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + (3*I*c - 10*d)/(24*a*f*(c + I*d)**2*sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**2) + d*(2*c**3 + 9*I*c**2*d - 17*c*d**2 + 60*I*d**3)/(16*a**3*f*(c - I*d)*(c + I*d)**4*sqrt(c + d*tan(e + f*x))) + (2*I*c**3 - 12*c**2*d - 33*I*c*d**2 + 58*d**3)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*(c + I*d)**(sympy.S(9)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1120():
    f = (I*a*tan(e + f*x) + a)**3/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -8*I*a**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) + 4*a**3*(c - 4*I*d)*(I*c - d)/(3*d**2*f*(c - I*d)**2*sqrt(c + d*tan(e + f*x))) + (2*c + 2*I*d)*(I*a**3*tan(e + f*x) + a**3)/(d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c - 3*I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1121():
    f = (I*a*tan(e + f*x) + a)**2/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 4*I*a**2/(f*(c - I*d)**2*sqrt(c + d*tan(e + f*x))) - 4*I*a**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) + 2*a**2*(I*c - d)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1122():
    f = (I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*a/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*I*c + 3*d)) + 2*I*a/(f*(c - I*d)**2*sqrt(c + d*tan(e + f*x))) - 2*I*a*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1123():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a))
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + d*(3*I*c + 7*d)/(6*a*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(I*c - d)) + d*(c**2 - 14*I*c*d - 5*d**2)/(2*a*f*(c - I*d)**2*(c + I*d)**3*sqrt(c + d*tan(e + f*x))) + (I*c - 6*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(2*a*f*(c + I*d)**(sympy.S(7)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(2*a*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1124():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**2)
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) + d*(6*c**2 + 27*I*c*d + 49*d**2)/(24*a**2*f*(c - I*d)*(c + I*d)**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)) + d*(2*c**3 + 9*I*c**2*d + 88*c*d**2 - 45*I*d**3)/(8*a**2*f*(c - I*d)**2*(c + I*d)**4*sqrt(c + d*tan(e + f*x))) + (2*I*c - 9*d)/(8*a**2*f*(c + I*d)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*tan(e + f*x) + 1)) + (2*I*c**2 - 14*c*d - 47*I*d**2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(8*a**2*f*(c + I*d)**(sympy.S(9)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(4*a**2*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1125():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**3)
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (2*c**2 + 11*I*c*d - 30*d**2)/(16*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)**3*(I*a**3*tan(e + f*x) + a**3)) + (I*c - 4*d)/(8*a*f*(c + I*d)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**2) + d*(6*c**3 + 33*I*c**2*d - 83*c*d**2 + 154*I*d**3)/(48*a**3*f*(c - I*d)*(c + I*d)**4*(c + d*tan(e + f*x))**(sympy.S(3)/2)) + d*(2*c**4 + 11*I*c**3*d - 26*c**2*d**2 + 253*I*c*d**3 + 150*d**4)/(16*a**3*f*(c - I*d)**2*(c + I*d)**5*sqrt(c + d*tan(e + f*x))) + (2*I*c**3 - 16*c**2*d - 61*I*c*d**2 + 152*d**3)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(16*a**3*f*(c + I*d)**(sympy.S(11)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(8*a**3*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1126():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*(c**2 + 10*I*c*d + 23*d**2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(4*d**(sympy.S(3)/2)*f) + a**2*(c + 9*I*d)*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(4*d*f) - a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1127():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*(I*c + 3*d)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(sqrt(d)*f) + a**2*(c + I*d)*sqrt(c + d*tan(e + f*x))/(d*f*sqrt(I*a*tan(e + f*x) + a)) - a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(d*f*sqrt(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1128():
    f = sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)
    F = -2*(-1)**(sympy.S(1)/4)*sqrt(a)*sqrt(d)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/f - sqrt(2)*I*sqrt(a)*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1129():
    f = sqrt(c + d*tan(e + f*x))/sqrt(I*a*tan(e + f*x) + a)
    F = I*sqrt(c + d*tan(e + f*x))/(f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1130():
    f = sqrt(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c + d*tan(e + f*x))**(sympy.S(3)/2)/(f*(3*I*c - 3*d)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + I*sqrt(c + d*tan(e + f*x))/(2*a*f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1131():
    f = sqrt(c + d*tan(e + f*x))/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = I*sqrt(c + d*tan(e + f*x))/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + sqrt(c + d*tan(e + f*x))*(5*I*c - 3*d)/(30*a*f*(c + I*d)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(c + d*tan(e + f*x))*(20*c*d - I*(15*c**2 + 3*d**2))/(60*a**2*f*(c + I*d)**2*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*sqrt(c - I*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1132():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*(c - 3*I*d)*(c**2 + 18*I*c*d + 15*d**2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(8*d**(sympy.S(3)/2)*f) + a**2*(c + 13*I*d)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)/(12*d*f) - a**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)*sqrt(I*a*tan(e + f*x) + a)/(3*d*f) + a**2*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*(c**2 + 14*I*c*d + 19*d**2)/(8*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1133():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*(3*I*c**2 + 18*c*d - 11*I*d**2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(4*sqrt(d)*f) + a**2*(c + I*d)*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(2*d*f*sqrt(I*a*tan(e + f*x) + a)) - a**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(2*d*f*sqrt(I*a*tan(e + f*x) + a)) + a*sqrt(c + d*tan(e + f*x))*(3*I*c + 5*d)*sqrt(I*a*tan(e + f*x) + a)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1134():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)
    F = -(-1)**(sympy.S(1)/4)*sqrt(a)*sqrt(d)*(3*c - I*d)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/f - sqrt(2)*I*sqrt(a)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f + d*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1135():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = sqrt(c + d*tan(e + f*x))*(I*c - d)/(f*sqrt(I*a*tan(e + f*x) + a)) + 2*(-1)**(sympy.S(3)/4)*d**(sympy.S(3)/2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(sqrt(a)*f) - sqrt(2)*I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1136():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = I*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + sqrt(c + d*tan(e + f*x))*(I*c + d)/(2*a*f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1137():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c + d*tan(e + f*x))**(sympy.S(5)/2)/(f*(5*I*c - 5*d)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + I*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(6*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + sqrt(c + d*tan(e + f*x))*(I*c + d)/(4*a**2*f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1138():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*(5*c**4 + 100*I*c**3*d + 690*c**2*d**2 - 900*I*c*d**3 - 363*d**4)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(64*d**(sympy.S(3)/2)*f) + a**2*(c + 17*I*d)*(c + d*tan(e + f*x))**(sympy.S(5)/2)*sqrt(I*a*tan(e + f*x) + a)/(24*d*f) - a**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)*sqrt(I*a*tan(e + f*x) + a)/(4*d*f) + a**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)*(5*c**2 + 90*I*c*d + 107*d**2)/(96*d*f) + a**2*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)*(5*c**3 + 95*I*c**2*d + 273*c*d**2 - 149*I*d**3)/(64*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1139():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f - (-1)**(sympy.S(1)/4)*a**(sympy.S(3)/2)*(5*I*c**3 + 45*c**2*d - 55*I*c*d**2 - 23*d**3)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(8*sqrt(d)*f) + a**2*(c + I*d)*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(3*d*f*sqrt(I*a*tan(e + f*x) + a)) - a**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(3*d*f*sqrt(I*a*tan(e + f*x) + a)) + a*(c - 3*I*d)*sqrt(c + d*tan(e + f*x))*(5*I*c + 3*d)*sqrt(I*a*tan(e + f*x) + a)/(8*f) + a*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(5*I*c + 7*d)*sqrt(I*a*tan(e + f*x) + a)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1140():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*sqrt(I*a*tan(e + f*x) + a)
    F = -(-1)**(sympy.S(1)/4)*sqrt(a)*sqrt(d)*(15*c**2 - 10*I*c*d - 7*d**2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(4*f) - sqrt(2)*I*sqrt(a)*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/f + d*(c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)/(2*f) + d*sqrt(c + d*tan(e + f*x))*(7*c - I*d)*sqrt(I*a*tan(e + f*x) + a)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1141():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/sqrt(I*a*tan(e + f*x) + a)
    F = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)/(f*sqrt(I*a*tan(e + f*x) + a)) - d*(c + 2*I*d)*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(a*f) + (-1)**(sympy.S(1)/4)*d**(sympy.S(3)/2)*(5*I*c - d)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(sqrt(a)*f) - sqrt(2)*I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1142():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)
    F = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)/(3*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (c + I*d)*sqrt(c + d*tan(e + f*x))*(I*c + 3*d)/(2*a*f*sqrt(I*a*tan(e + f*x) + a)) + 2*(-1)**(sympy.S(1)/4)*d**(sympy.S(5)/2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(a**(sympy.S(3)/2)*f) - sqrt(2)*I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1143():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)
    F = I*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c + d)/(6*a*f*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + I*(c - I*d)**2*sqrt(c + d*tan(e + f*x))/(4*a**2*f*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1144():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/sqrt(c + d*tan(e + f*x))
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*sqrt(c - I*d)) - (-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*(c + 5*I*d)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) - a**2*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1145():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/sqrt(c + d*tan(e + f*x))
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*sqrt(c - I*d)) - 2*(-1)**(sympy.S(3)/4)*a**(sympy.S(3)/2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1146():
    f = sqrt(I*a*tan(e + f*x) + a)/sqrt(c + d*tan(e + f*x))
    F = -sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1147():
    f = 1/(sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a))
    F = 2*d*sqrt(c + d*tan(e + f*x))/(f*(c**2 + d**2)*sqrt(I*a*tan(e + f*x) + a)) - sqrt(c + d*tan(e + f*x))/(f*(I*c + d)*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1148():
    f = 1/(sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    F = -sqrt(c + d*tan(e + f*x))/(f*(3*I*c - 3*d)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + sqrt(c + d*tan(e + f*x))*(3*I*c - 7*d)/(6*a*f*(c + I*d)**2*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1149():
    f = 1/(sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    F = -sqrt(c + d*tan(e + f*x))/(f*(5*I*c - 5*d)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + sqrt(c + d*tan(e + f*x))*(5*I*c - 13*d)/(30*a*f*(c + I*d)**2*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + sqrt(c + d*tan(e + f*x))*(15*c**2 + 50*I*c*d - 67*d**2)/(60*a**2*f*(I*c - d)**3*sqrt(I*a*tan(e + f*x) + a)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1150():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(3)/2)) + 2*(-1)**(sympy.S(1)/4)*a**(sympy.S(5)/2)*atanh((-1)**(sympy.S(3)/4)*sqrt(d)*sqrt(I*a*tan(e + f*x) + a)/(sqrt(a)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) + 2*a**2*(c + I*d)*sqrt(I*a*tan(e + f*x) + a)/(d*f*(c - I*d)*sqrt(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1151():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(3)/2)) - 2*a*sqrt(I*a*tan(e + f*x) + a)/(f*sqrt(c + d*tan(e + f*x))*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1152():
    f = sqrt(I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(3)/2)) - 2*d*sqrt(I*a*tan(e + f*x) + a)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1153():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a))
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(I*c - d)*sqrt(I*a*tan(e + f*x) + a)) + d*(c - 3*I*d)*sqrt(I*a*tan(e + f*x) + a)/(a*f*(c - I*d)*(c + I*d)**2*sqrt(c + d*tan(e + f*x))) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1154():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(3*I*c - 3*d)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (3*I*c - 11*d)/(6*a*f*(c + I*d)**2*sqrt(c + d*tan(e + f*x))*sqrt(I*a*tan(e + f*x) + a)) + d*(c + 5*I*d)*(3*c - 5*I*d)*sqrt(I*a*tan(e + f*x) + a)/(6*a**2*f*(c - I*d)*(c + I*d)**3*sqrt(c + d*tan(e + f*x))) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1155():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    F = -1/(f*sqrt(c + d*tan(e + f*x))*(5*I*c - 5*d)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (5*I*c - 17*d)/(30*a*f*(c + I*d)**2*sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (15*c**2 + 70*I*c*d - 151*d**2)/(60*a**2*f*sqrt(c + d*tan(e + f*x))*(I*c - d)**3*sqrt(I*a*tan(e + f*x) + a)) + d*sqrt(I*a*tan(e + f*x) + a)*(15*c**3 + 65*I*c**2*d - 117*c*d**2 + 317*I*d**3)/(60*a**3*f*(c - I*d)*(c + I*d)**4*sqrt(c + d*tan(e + f*x))) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1156():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(5)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -4*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(5)/2)) + 4*I*a**2*sqrt(I*a*tan(e + f*x) + a)/(f*(c - I*d)**2*sqrt(c + d*tan(e + f*x))) - 2*a*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*I*c + 3*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1157():
    f = (I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(5)/2)) + 2*I*a*sqrt(I*a*tan(e + f*x) + a)/(f*(c - I*d)**2*sqrt(c + d*tan(e + f*x))) - 2*d*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c**2 + 3*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1158():
    f = sqrt(I*a*tan(e + f*x) + a)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(f*(c - I*d)**(sympy.S(5)/2)) - d*(10*c + 2*I*d)*sqrt(I*a*tan(e + f*x) + a)/(3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) - 2*d*sqrt(I*a*tan(e + f*x) + a)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c**2 + 3*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1159():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*sqrt(I*a*tan(e + f*x) + a))
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)*sqrt(I*a*tan(e + f*x) + a)) + d*(c - 7*I*d)*(3*c - I*d)*sqrt(I*a*tan(e + f*x) + a)/(3*a*f*(c - I*d)**2*(c + I*d)**3*sqrt(c + d*tan(e + f*x))) + d*(3*I*c + 5*d)*sqrt(I*a*tan(e + f*x) + a)/(3*a*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(I*c - d)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(2*sqrt(a)*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1160():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2))
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*I*c - 3*d)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (I*c - 5*d)/(2*a*f*(c + I*d)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*sqrt(I*a*tan(e + f*x) + a)) + d*(c - 3*I*d)*sqrt(I*a*tan(e + f*x) + a)*(3*c**2 + 22*I*c*d + 13*d**2)/(6*a**2*f*(c - I*d)**2*(c + I*d)**4*sqrt(c + d*tan(e + f*x))) + d*sqrt(I*a*tan(e + f*x) + a)*(3*c**2 + 14*I*c*d + 21*d**2)/(6*a**2*f*(c - I*d)*(c + I*d)**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1161():
    f = 1/((c + d*tan(e + f*x))**(sympy.S(5)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2))
    F = -1/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(5*I*c - 5*d)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/2)) + (5*I*c - 21*d)/(30*a*f*(c + I*d)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**(sympy.S(3)/2)) + (5*c**2 + 30*I*c*d - 89*d**2)/(20*a**2*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*c - d)**3*sqrt(I*a*tan(e + f*x) + a)) + d*sqrt(I*a*tan(e + f*x) + a)*(15*c**3 + 85*I*c**2*d - 221*c*d**2 + 361*I*d**3)/(60*a**3*f*(c - I*d)*(c + I*d)**4*(c + d*tan(e + f*x))**(sympy.S(3)/2)) + d*sqrt(I*a*tan(e + f*x) + a)*(15*c**4 + 80*I*c**3*d - 182*c**2*d**2 + 1224*I*c*d**3 + 707*d**4)/(60*a**3*f*(c - I*d)**2*(c + I*d)**5*sqrt(c + d*tan(e + f*x))) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sqrt(c + d*tan(e + f*x))/(sqrt(c - I*d)*sqrt(I*a*tan(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1162():
    f = (c + d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**m
    F = -I*(c + d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**m*appellf1(m, 1, -n, m + 1, I*tan(e + f*x)/2 + sympy.S.Half, -d*(I*tan(e + f*x) + 1)/(I*c - d))/(2*f*m*((c + d*tan(e + f*x))/(c + I*d))**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1163():
    f = (c + d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**3
    F = 4*a**3*(c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(f*(n + 1)*(I*c + d)) + a**3*(c + d*tan(e + f*x))**(n + 1)*(I*c - d*(2*n + 5))/(d**2*f*(n + 1)*(n + 2)) - (c + d*tan(e + f*x))**(n + 1)*(I*a**3*tan(e + f*x) + a**3)/(d*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1164():
    f = (c + d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)**2
    F = 2*a**2*(c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(f*(n + 1)*(I*c + d)) - a**2*(c + d*tan(e + f*x))**(n + 1)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1165():
    f = (c + d*tan(e + f*x))**n*(I*a*tan(e + f*x) + a)
    F = a*(c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(f*(n + 1)*(I*c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1166():
    f = (c + d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)
    F = -(c + d*tan(e + f*x))**(n + 1)/(f*(2*I*c - 2*d)*(I*a*tan(e + f*x) + a)) + (c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(4*a*f*(n + 1)*(I*c + d)) + (c + d*tan(e + f*x))**(n + 1)*(I*c + 2*d*n - d)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c + I*d))/(4*a*f*(c + I*d)**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1167():
    f = (c + d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**2
    F = -(c + d*tan(e + f*x))**(n + 1)/(f*(4*I*c - 4*d)*(I*a*tan(e + f*x) + a)**2) + (c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(8*a**2*f*(n + 1)*(I*c + d)) + (c + d*tan(e + f*x))**(n + 1)*(c**2 + 2*I*c*d*(1 - n) - d**2*(2*n**2 - 4*n + 1))*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c + I*d))/(8*a**2*f*(n + 1)*(I*c - d)**3) + (c + d*tan(e + f*x))**(n + 1)*(I*c - d*(2 - n))/(4*a**2*f*(c + I*d)**2*(I*tan(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1168():
    f = (c + d*tan(e + f*x))**n/(I*a*tan(e + f*x) + a)**3
    F = -(c + d*tan(e + f*x))**(n + 1)/(f*(6*I*c - 6*d)*(I*a*tan(e + f*x) + a)**3) + (c + d*tan(e + f*x))**(n + 1)*(3*I*c**2 - 3*c*d*(3 - n) - I*d**2*(2*n**2 - 9*n + 10))/(24*f*(c + I*d)**3*(I*a**3*tan(e + f*x) + a**3)) + (c + d*tan(e + f*x))**(n + 1)*(3*I*c - d*(7 - 2*n))/(24*a*f*(c + I*d)**2*(I*a*tan(e + f*x) + a)**2) + (c + d*tan(e + f*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c - I*d))/(16*a**3*f*(n + 1)*(I*c + d)) + (c + d*tan(e + f*x))**(n + 1)*(3*I*c**3 - c**2*d*(9 - 6*n) - 3*I*c*d**2*(2*n**2 - 6*n + 3) + d**3*(-4*n**3 + 18*n**2 - 20*n + 3))*hyper((1, n + 1), (n + 2,), (c + d*tan(e + f*x))/(c + I*d))/(48*a**3*f*(c + I*d)**4*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1169():
    f = (c + d*tan(e + f*x))**3*(I*a*tan(e + f*x) + a)**m
    F = d*(c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**m/(f*(m + 2)) - 2*d*(I*a*tan(e + f*x) + a)**m*(-c**2*(m + 3) + I*c*d*m + d**2)/(f*m*(m + 2)) + (I*c + d)**3*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m) - d**2*(I*a*tan(e + f*x) + a)**(m + 1)*(I*c*(m + 4) + d*m)/(a*f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1170():
    f = (c + d*tan(e + f*x))**2*(I*a*tan(e + f*x) + a)**m
    F = 2*c*d*(I*a*tan(e + f*x) + a)**m/(f*m) - I*(c - I*d)**2*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m) - I*d**2*(I*a*tan(e + f*x) + a)**(m + 1)/(a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1171():
    f = (c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**m
    F = d*(I*a*tan(e + f*x) + a)**m/(f*m) - (I*c + d)*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1172():
    f = (I*a*tan(e + f*x) + a)**m/(c + d*tan(e + f*x))
    F = -d*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), -d*(I*tan(e + f*x) + 1)/(I*c - d))/(f*m*(c**2 + d**2)) + (I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(f*m*(2*I*c + 2*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1173():
    f = (I*a*tan(e + f*x) + a)**m/(c + d*tan(e + f*x))**2
    F = -d*(I*a*tan(e + f*x) + a)**m/(f*(c + d*tan(e + f*x))*(c**2 + d**2)) - d*(c*(2 - m) + I*d*m)*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), -d*(I*tan(e + f*x) + 1)/(I*c - d))/(f*m*(c**2 + d**2)**2) - I*(I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m*(c - I*d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1174():
    f = (I*a*tan(e + f*x) + a)**m/(c + d*tan(e + f*x))**3
    F = -d*(c*(4 - m) + I*d*m)*(I*a*tan(e + f*x) + a)**m/(2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) - d*(I*a*tan(e + f*x) + a)**m/(f*(c + d*tan(e + f*x))**2*(2*c**2 + 2*d**2)) - d*(I*a*tan(e + f*x) + a)**m*(c**2*(m**2 - 5*m + 6) + 2*I*c*d*m*(3 - m) - d**2*(m**2 - m + 2))*hyper((1, m), (m + 1,), -d*(I*tan(e + f*x) + 1)/(I*c - d))/(2*f*m*(c**2 + d**2)**3) - (I*a*tan(e + f*x) + a)**m*hyper((1, m), (m + 1,), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m*(I*c + d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1175():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(I*a*tan(e + f*x) + a)**m
    F = -sqrt(c + d*tan(e + f*x))*(I*c - d)*(I*a*tan(e + f*x) + a)**m*appellf1(m, sympy.S(-3)/2, 1, m + 1, -d*(I*tan(e + f*x) + 1)/(I*c - d), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m*sqrt((c + d*tan(e + f*x))/(c + I*d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1176():
    f = sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**m
    F = -I*sqrt(c + d*tan(e + f*x))*(I*a*tan(e + f*x) + a)**m*appellf1(m, sympy.S(-1)/2, 1, m + 1, -d*(I*tan(e + f*x) + 1)/(I*c - d), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m*sqrt((c + d*tan(e + f*x))/(c + I*d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1177():
    f = (I*a*tan(e + f*x) + a)**m/sqrt(c + d*tan(e + f*x))
    F = -I*sqrt((c + d*tan(e + f*x))/(c + I*d))*(I*a*tan(e + f*x) + a)**m*appellf1(m, sympy.S.Half, 1, m + 1, -d*(I*tan(e + f*x) + 1)/(I*c - d), I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*m*sqrt(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1178():
    f = (I*a*tan(e + f*x) + a)**m/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = sqrt((c + d*tan(e + f*x))/(c + I*d))*(I*a*tan(e + f*x) + a)**m*appellf1(m, 1, sympy.S(3)/2, m + 1, I*tan(e + f*x)/2 + sympy.S.Half, -d*(I*tan(e + f*x) + 1)/(I*c - d))/(f*m*sqrt(c + d*tan(e + f*x))*(2*I*c - 2*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1179():
    f = (I*a*tan(e + f*x) + a)**m/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*sqrt((c + d*tan(e + f*x))/(c + I*d))*(I*a*tan(e + f*x) + a)**m*appellf1(m, 1, sympy.S(5)/2, m + 1, I*tan(e + f*x)/2 + sympy.S.Half, -d*(I*tan(e + f*x) + 1)/(I*c - d))/(2*f*m*(c + I*d)**2*sqrt(c + d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1180():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))
    F = b*(a**2*d + 2*a*b*c - b**2*d)*tan(e + f*x)/f + d*(a + b*tan(e + f*x))**3/(3*f) + x*(a**3*c - 3*a**2*b*d - 3*a*b**2*c + b**3*d) + (a + b*tan(e + f*x))**2*(a*d + b*c)/(2*f) - (a**3*d + 3*a**2*b*c - 3*a*b**2*d - b**3*c)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1181():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))
    F = b*(a*d + b*c)*tan(e + f*x)/f + d*(a + b*tan(e + f*x))**2/(2*f) + x*(a**2*c - 2*a*b*d - b**2*c) - (a**2*d + 2*a*b*c - b**2*d)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1182():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))
    F = b*d*tan(e + f*x)/f + x*(a*c - b*d) - (a*d + b*c)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1183():
    f = (c + d*tan(e + f*x))/(a + b*tan(e + f*x))
    F = x*(a*c + b*d)/(a**2 + b**2) + (-a*d + b*c)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1184():
    f = (c + d*tan(e + f*x))/(a + b*tan(e + f*x))**2
    F = x*(a**2*c + 2*a*b*d - b**2*c)/(a**2 + b**2)**2 + (-a**2*d + 2*a*b*c + b**2*d)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2) - (-a*d + b*c)/(f*(a + b*tan(e + f*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1185():
    f = (c + d*tan(e + f*x))/(a + b*tan(e + f*x))**3
    F = x*(a**3*c + 3*a**2*b*d - 3*a*b**2*c - b**3*d)/(a**2 + b**2)**3 + (-a**3*d + 3*a**2*b*c + 3*a*b**2*d - b**3*c)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3) - (-a**2*d + 2*a*b*c + b**2*d)/(f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - (-a*d + b*c)/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1186():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**2
    F = 2*b*(a*c - b*d)*(a*d + b*c)*tan(e + f*x)/f + 2*c*d*(a + b*tan(e + f*x))**3/(3*f) - x*(-a**3*(c**2 - d**2) + 6*a**2*b*c*d + 3*a*b**2*(c**2 - d**2) - 2*b**3*c*d) + (a + b*tan(e + f*x))**2*(2*a*c*d + b*(c**2 - d**2))/(2*f) - (2*a**3*c*d + 3*a**2*b*(c**2 - d**2) - 6*a*b**2*c*d - b**3*(c**2 - d**2))*log(cos(e + f*x))/f + d**2*(a + b*tan(e + f*x))**4/(4*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1187():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**2
    F = b*(2*a*c*d + b*(c**2 - d**2))*tan(e + f*x)/f + c*d*(a + b*tan(e + f*x))**2/f + x*(a*c - a*d - b*c - b*d)*(a*c + a*d + b*c - b*d) - (a*c - b*d)*(2*a*d + 2*b*c)*log(cos(e + f*x))/f + d**2*(a + b*tan(e + f*x))**3/(3*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1188():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**2
    F = b*(c + d*tan(e + f*x))**2/(2*f) + d*(a*d + b*c)*tan(e + f*x)/f - x*(-a*(c**2 - d**2) + 2*b*c*d) - (2*a*c*d + b*(c**2 - d**2))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1189():
    f = (c + d*tan(e + f*x))**2/(a + b*tan(e + f*x))
    F = a*x*(-a*d + b*c)**2/(b**2*(a**2 + b**2)) - d**2*log(cos(e + f*x))/(b*f) + (-a*d + b*c)**2*log(a*cos(e + f*x) + b*sin(e + f*x))/(b*f*(a**2 + b**2)) + d*x*(-a*d + 2*b*c)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1190():
    f = (c + d*tan(e + f*x))**2/(a + b*tan(e + f*x))**2
    F = -x*(a*(c - d) + b*(c + d))*(-a*(c + d) + b*(c - d))/(a**2 + b**2)**2 + (a*c + b*d)*(-2*a*d + 2*b*c)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2) - (-a*d + b*c)**2/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1191():
    f = (c + d*tan(e + f*x))**2/(a + b*tan(e + f*x))**3
    F = x*(a**3*(c**2 - d**2) + 6*a**2*b*c*d - 3*a*b**2*(c**2 - d**2) - 2*b**3*c*d)/(a**2 + b**2)**3 - (2*a**3*c*d - 3*a**2*b*(c**2 - d**2) - 6*a*b**2*c*d + b**3*(c**2 - d**2))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3) - (a*c + b*d)*(-2*a*d + 2*b*c)/(f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - (-a*d + b*c)**2/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1192():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**3
    F = b**2*(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**4/(5*d*f) - b**2*(c + d*tan(e + f*x))**4*(-11*a*d + b*c)/(20*d**2*f) + b*(3*a**2 - b**2)*(c + d*tan(e + f*x))**3/(3*f) + d*(2*a**3*c*d + 3*a**2*b*(c**2 - d**2) - 6*a*b**2*c*d - b**3*(c**2 - d**2))*tan(e + f*x)/f - x*(a*c - b*d)*(-a**2*(c**2 - 3*d**2) + 8*a*b*c*d + b**2*(3*c**2 - d**2)) + (c + d*tan(e + f*x))**2*(a**3*d + 3*a**2*b*c - 3*a*b**2*d - b**3*c)/(2*f) + (a*d + b*c)*(-a**2*(3*c**2 - d**2) + 8*a*b*c*d + b**2*(c**2 - 3*d**2))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1193():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3
    F = 2*a*b*(c + d*tan(e + f*x))**3/(3*f) + b**2*(c + d*tan(e + f*x))**4/(4*d*f) + 2*d*(a*c - b*d)*(a*d + b*c)*tan(e + f*x)/f - x*(-a**2*(c**3 - 3*c*d**2) + 2*a*b*d*(3*c**2 - d**2) + b**2*c*(c**2 - 3*d**2)) + (c + d*tan(e + f*x))**2*(a**2*d + 2*a*b*c - b**2*d)/(2*f) - (a**2*(3*c**2*d - d**3) + 2*a*b*c*(c**2 - 3*d**2) - b**2*d*(3*c**2 - d**2))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1194():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**3
    F = b*(c + d*tan(e + f*x))**3/(3*f) + d*(2*a*c*d + b*(c**2 - d**2))*tan(e + f*x)/f - x*(-a*(c**3 - 3*c*d**2) + b*d*(3*c**2 - d**2)) + (c + d*tan(e + f*x))**2*(a*d + b*c)/(2*f) - (3*a*c**2*d - a*d**3 + b*c**3 - 3*b*c*d**2)*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1195():
    f = (c + d*tan(e + f*x))**3/(a + b*tan(e + f*x))
    F = x*(a*c**3 - 3*a*c*d**2 + 3*b*c**2*d - b*d**3)/(a**2 + b**2) + (-3*a*c**2*d + a*d**3 + b*c**3 - 3*b*c*d**2)*log(cos(e + f*x))/(f*(a**2 + b**2)) + d**2*(c + d*tan(e + f*x))/(b*f) + (-a*d + b*c)**3*log(a + b*tan(e + f*x))/(b**2*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1196():
    f = (c + d*tan(e + f*x))**3/(a + b*tan(e + f*x))**2
    F = -x*(-a**2*(c**3 - 3*c*d**2) - 2*a*b*d*(3*c**2 - d**2) + b**2*c*(c**2 - 3*d**2))/(a**2 + b**2)**2 + (-a**2*(3*c**2*d - d**3) + 2*a*b*c*(c**2 - 3*d**2) + b**2*d*(3*c**2 - d**2))*log(cos(e + f*x))/(f*(a**2 + b**2)**2) - (c + d*tan(e + f*x))*(-a*d + b*c)**2/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + (-a*d + b*c)**2*(a**2*d + 2*a*b*c + 3*b**2*d)*log(a + b*tan(e + f*x))/(b**2*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1197():
    f = (c + d*tan(e + f*x))**3/(a + b*tan(e + f*x))**3
    F = x*(a*c + b*d)*(a**2*(c**2 - 3*d**2) + 8*a*b*c*d - b**2*(3*c**2 - d**2))/(a**2 + b**2)**3 + (-a*d + b*c)*(a**2*(3*c**2 - d**2) + 8*a*b*c*d - b**2*(c**2 - 3*d**2))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3) - (c + d*tan(e + f*x))*(-a*d + b*c)**2/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) - (-a*d + b*c)**2*(a**2*d + 4*a*b*c + 5*b**2*d)/(2*b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1198():
    f = (a + b*tan(e + f*x))**4/(c + d*tan(e + f*x))
    F = -b**3*(-3*a*d + b*c)*tan(e + f*x)/(d**2*f) + b**2*(a + b*tan(e + f*x))**2/(2*d*f) + x*(a**4*c + 4*a**3*b*d - 6*a**2*b**2*c - 4*a*b**3*d + b**4*c)/(c**2 + d**2) - (-a**4*d + 4*a**3*b*c + 6*a**2*b**2*d - 4*a*b**3*c - b**4*d)*log(cos(e + f*x))/(f*(c**2 + d**2)) + (-a*d + b*c)**4*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1199():
    f = (a + b*tan(e + f*x))**3/(c + d*tan(e + f*x))
    F = b**2*(a + b*tan(e + f*x))/(d*f) + x*(a**3*c + 3*a**2*b*d - 3*a*b**2*c - b**3*d)/(c**2 + d**2) - (-a**3*d + 3*a**2*b*c + 3*a*b**2*d - b**3*c)*log(cos(e + f*x))/(f*(c**2 + d**2)) - (-a*d + b*c)**3*log(c + d*tan(e + f*x))/(d**2*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1200():
    f = (a + b*tan(e + f*x))**2/(c + d*tan(e + f*x))
    F = -b**2*log(cos(e + f*x))/(d*f) - b*x*(-2*a*d + b*c)/d**2 + c*x*(-a*d + b*c)**2/(d**2*(c**2 + d**2)) + (-a*d + b*c)**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(d*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1201():
    f = (a + b*tan(e + f*x))/(c + d*tan(e + f*x))
    F = x*(a*c + b*d)/(c**2 + d**2) - (-a*d + b*c)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1202():
    f = 1/((a + b*tan(e + f*x))*(c + d*tan(e + f*x)))
    F = b**2*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c)) - d**2*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)) + x*(a*c - b*d)/((a**2 + b**2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1203():
    f = 1/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x)))
    F = b**2*(-3*a**2*d + 2*a*b*c - b**2*d)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**2) - b**2/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c)) + d**3*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)**2) + x*(a**2*c - 2*a*b*d - b**2*c)/((a**2 + b**2)**2*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1204():
    f = 1/((a + b*tan(e + f*x))**3*(c + d*tan(e + f*x)))
    F = -b**2*(-6*a**4*d**2 + 8*a**3*b*c*d - 3*a**2*b**2*(c**2 + d**2) + b**4*(c**2 - d**2))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3*(-a*d + b*c)**3) - b**2*(-3*a**2*d + 2*a*b*c - b**2*d)/(f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)**2) - b**2/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)*(-a*d + b*c)) - d**4*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)**3) + x*(a**3*c - 3*a**2*b*d - 3*a*b**2*c + b**3*d)/((a**2 + b**2)**3*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1205():
    f = (a + b*tan(e + f*x))**4/(c + d*tan(e + f*x))**2
    F = -b**2*(a*d*(-a*d + 2*b*c) - b**2*(2*c**2 + d**2))*tan(e + f*x)/(d**2*f*(c**2 + d**2)) + x*(a**4*(c**2 - d**2) + 8*a**3*b*c*d - 6*a**2*b**2*(c**2 - d**2) - 8*a*b**3*c*d + b**4*(c**2 - d**2))/(c**2 + d**2)**2 - (2*a**2*c + 4*a*b*d - 2*b**2*c)*(-a**2*d + 2*a*b*c + b**2*d)*log(cos(e + f*x))/(f*(c**2 + d**2)**2) - (a + b*tan(e + f*x))**2*(-a*d + b*c)**2/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2)) - 2*(-a*d + b*c)**3*(a*c*d + b*(c**2 + 2*d**2))*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1206():
    f = (a + b*tan(e + f*x))**3/(c + d*tan(e + f*x))**2
    F = x*(a**3*(c**2 - d**2) + 6*a**2*b*c*d - 3*a*b**2*(c**2 - d**2) - 2*b**3*c*d)/(c**2 + d**2)**2 + (2*a**3*c*d - 3*a**2*b*(c**2 - d**2) - 6*a*b**2*c*d + b**3*(c**2 - d**2))*log(cos(e + f*x))/(f*(c**2 + d**2)**2) - (a + b*tan(e + f*x))*(-a*d + b*c)**2/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2)) + (-a*d + b*c)**2*(2*a*c*d + b*(c**2 + 3*d**2))*log(c + d*tan(e + f*x))/(d**2*f*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1207():
    f = (a + b*tan(e + f*x))**2/(c + d*tan(e + f*x))**2
    F = -x*(a*(c - d) + b*(c + d))*(-a*(c + d) + b*(c - d))/(c**2 + d**2)**2 - (a*c + b*d)*(-2*a*d + 2*b*c)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2) - (-a*d + b*c)**2/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1208():
    f = (a + b*tan(e + f*x))/(c + d*tan(e + f*x))**2
    F = x*(a*(c**2 - d**2) + 2*b*c*d)/(c**2 + d**2)**2 + (2*a*c*d - b*(c**2 - d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2) + (-a*d + b*c)/(f*(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1209():
    f = 1/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**2)
    F = b**3*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c)**2) + d**2*(2*a*c*d - b*(3*c**2 + d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**2) + d**2/(f*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) - x*(-a*(c**2 - d**2) + 2*b*c*d)/((a**2 + b**2)*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1210():
    f = 1/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**2)
    F = 2*b**3*(-2*a**2*d + a*b*c - b**2*d)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**3) - b**2/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))*(-a*d + b*c)) - 2*d**3*(a*c*d - b*(2*c**2 + d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**3) - d*(a**2*d**2 + b**2*(c**2 + 2*d**2))/(f*(a**2 + b**2)*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) + x*(a*(c - d) - b*(c + d))*(a*(c + d) + b*(c - d))/((a**2 + b**2)**2*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1211():
    f = 1/((a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**2)
    F = -b**3*(-10*a**4*d**2 + 10*a**3*b*c*d - 3*a**2*b**2*(c**2 + 3*d**2) + 2*a*b**3*c*d + b**4*(c**2 - 3*d**2))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3*(-a*d + b*c)**4) - b**2*(-7*a**2*d + 4*a*b*c - 3*b**2*d)/(2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(c + d*tan(e + f*x))*(-a*d + b*c)**2) - b**2/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)*(c + d*tan(e + f*x))*(-a*d + b*c)) - d**4*(-2*a*c*d + 5*b*c**2 + 3*b*d**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**4) + d*(a**4*d**3 + 2*a**2*b**2*d*(2*c**2 + 3*d**2) - 2*a*b**3*c*(c**2 + d**2) + b**4*d*(2*c**2 + 3*d**2))/(f*(a**2 + b**2)**2*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**3) - x*(-a**3*(c**2 - d**2) + 6*a**2*b*c*d + 3*a*b**2*(c**2 - d**2) - 2*b**3*c*d)/((a**2 + b**2)**3*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1212():
    f = (a + b*tan(e + f*x))**4/(c + d*tan(e + f*x))**3
    F = -x*(-a**4*(c**3 - 3*c*d**2) - 4*a**3*b*d*(3*c**2 - d**2) + 6*a**2*b**2*c*(c**2 - 3*d**2) + 4*a*b**3*d*(3*c**2 - d**2) - b**4*c*(c**2 - 3*d**2))/(c**2 + d**2)**3 - (-a**4*(3*c**2*d - d**3) + 4*a**3*b*c*(c**2 - 3*d**2) + 6*a**2*b**2*d*(3*c**2 - d**2) - 4*a*b**3*c*(c**2 - 3*d**2) - b**4*d*(3*c**2 - d**2))*log(cos(e + f*x))/(f*(c**2 + d**2)**3) - (a + b*tan(e + f*x))**2*(-a*d + b*c)**2/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2)) + (-a*d + b*c)**2*(a**2*d**2*(3*c**2 - d**2) + 2*a*b*c*d*(c**2 + 5*d**2) + b**2*(c**4 + 3*c**2*d**2 + 6*d**4))*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2)**3) + (-a*d + b*c)**3*(2*a*c*d + b*(c**2 + 3*d**2))/(d**3*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1213():
    f = (a + b*tan(e + f*x))**3/(c + d*tan(e + f*x))**3
    F = x*(a*c + b*d)*(a**2*(c**2 - 3*d**2) + 8*a*b*c*d - b**2*(3*c**2 - d**2))/(c**2 + d**2)**3 - (-a*d + b*c)*(a**2*(3*c**2 - d**2) + 8*a*b*c*d - b**2*(c**2 - 3*d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3) - (a + b*tan(e + f*x))*(-a*d + b*c)**2/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2)) - (-a*d + b*c)**2*(4*a*c*d + b*(c**2 + 5*d**2))/(2*d**2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1214():
    f = (a + b*tan(e + f*x))**2/(c + d*tan(e + f*x))**3
    F = -x*(-a**2*(c**3 - 3*c*d**2) - 2*a*b*d*(3*c**2 - d**2) + b**2*c*(c**2 - 3*d**2))/(c**2 + d**2)**3 - (-a**2*(3*c**2*d - d**3) + 2*a*b*c*(c**2 - 3*d**2) + b**2*d*(3*c**2 - d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3) + (a*c + b*d)*(-2*a*d + 2*b*c)/(f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) - (-a*d + b*c)**2/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1215():
    f = (a + b*tan(e + f*x))/(c + d*tan(e + f*x))**3
    F = x*(a*c**3 - 3*a*c*d**2 + 3*b*c**2*d - b*d**3)/(c**2 + d**2)**3 + (a*d*(3*c**2 - d**2) - b*(c**3 - 3*c*d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3) - (2*a*c*d - b*(c**2 - d**2))/(f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) + (-a*d + b*c)/(f*(c + d*tan(e + f*x))**2*(2*c**2 + 2*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1216():
    f = 1/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**3)
    F = b**4*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c)**3) + d**2*(-a**2*d**2*(3*c**2 - d**2) + 8*a*b*c**3*d - b**2*(6*c**4 + 3*c**2*d**2 + d**4))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3*(-a*d + b*c)**3) - d**2*(2*a*c*d - b*(3*c**2 + d**2))/(f*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + d**2/(f*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-2*a*d + 2*b*c)) - x*(-a*(c**3 - 3*c*d**2) + b*d*(3*c**2 - d**2))/((a**2 + b**2)*(c**2 + d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1217():
    f = 1/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3)
    F = b**4*(-5*a**2*d + 2*a*b*c - 3*b**2*d)*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**4) - b**2/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**2*(-a*d + b*c)) + d**3*(a**2*d**2*(3*c**2 - d**2) - 2*a*b*c*d*(5*c**2 + d**2) + b**2*(10*c**4 + 9*c**2*d**2 + 3*d**4))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3*(-a*d + b*c)**4) - d*(a**2*d**2 + b**2*(2*c**2 + 3*d**2))/(f*(2*a**2 + 2*b**2)*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-a*d + b*c)**2) + d*(2*a**3*c*d**3 - 2*a**2*b*d**2*(2*c**2 + d**2) + 2*a*b**2*c*d**3 - b**3*(c**4 + 6*c**2*d**2 + 3*d**4))/(f*(a**2 + b**2)*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) - x*(-a**2*(c**3 - 3*c*d**2) + a*b*(6*c**2*d - 2*d**3) + b**2*c*(c**2 - 3*d**2))/((a**2 + b**2)**2*(c**2 + d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1218():
    f = (a + b*tan(e + f*x))**3*sqrt(c + d*tan(e + f*x))
    F = 2*b**2*(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(5*d*f) - 4*b**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-6*a*d + b*c)/(15*d**2*f) + 2*b*(3*a**2 - b**2)*sqrt(c + d*tan(e + f*x))/f + sqrt(c - I*d)*(I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f - sqrt(c + I*d)*(I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1219():
    f = (a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))
    F = 4*a*b*sqrt(c + d*tan(e + f*x))/f + 2*b**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f) - I*(a - I*b)**2*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + I*(a + I*b)**2*sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1220():
    f = (a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))
    F = 2*b*sqrt(c + d*tan(e + f*x))/f - sqrt(c - I*d)*(I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + sqrt(c + I*d)*(I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1221():
    f = sqrt(c + d*tan(e + f*x))/(a + b*tan(e + f*x))
    F = -2*sqrt(b)*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)) + sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)) - sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1222():
    f = sqrt(c + d*tan(e + f*x))/(a + b*tan(e + f*x))**2
    F = -sqrt(b)*(-3*a**2*d + 4*a*b*c + b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*sqrt(-a*d + b*c)) - b*sqrt(c + d*tan(e + f*x))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)) + I*sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - I*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1223():
    f = sqrt(c + d*tan(e + f*x))/(a + b*tan(e + f*x))**3
    F = sqrt(b)*(-15*a**4*d**2 + 40*a**3*b*c*d - 6*a**2*b**2*(4*c**2 - 3*d**2) - 24*a*b**3*c*d + b**4*(8*c**2 + d**2))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*f*(a**2 + b**2)**3*(-a*d + b*c)**(sympy.S(3)/2)) - b*sqrt(c + d*tan(e + f*x))*(-7*a**2*d + 8*a*b*c + b**2*d)/(4*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)) - b*sqrt(c + d*tan(e + f*x))/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)) - sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1224():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*b**2*(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(7*d*f) - 4*b**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(-8*a*d + b*c)/(35*d**2*f) + 2*b*(3*a**2 - b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + (c - I*d)**(sympy.S(3)/2)*(I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f - (c + I*d)**(sympy.S(3)/2)*(I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + sqrt(c + d*tan(e + f*x))*(2*a**3*d + 6*a**2*b*c - 6*a*b**2*d - 2*b**3*c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1225():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 4*a*b*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*b**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*d*f) - I*(a - I*b)**2*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + I*(a + I*b)**2*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + sqrt(c + d*tan(e + f*x))*(2*a**2*d + 4*a*b*c - 2*b**2*d)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1226():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*b*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) - (c - I*d)**(sympy.S(3)/2)*(I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + I*d)**(sympy.S(3)/2)*(I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + sqrt(c + d*tan(e + f*x))*(2*a*d + 2*b*c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1227():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))
    F = (c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)) - (c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)) - 2*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(sqrt(b)*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1228():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**2
    F = -sqrt(c + d*tan(e + f*x))*(-a*d + b*c)/(f*(a + b*tan(e + f*x))*(a**2 + b**2)) + I*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2) - sqrt(-a*d + b*c)*(-a**2*d + 4*a*b*c + 3*b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(sqrt(b)*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1229():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**3
    F = -(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + (c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3) - sqrt(c + d*tan(e + f*x))*(-3*a**2*d + 8*a*b*c + 5*b**2*d)/(4*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(-a*d + b*c)/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)) + (-3*a**4*d**2 + 24*a**3*b*c*d - 2*a**2*b**2*(12*c**2 - 13*d**2) - 40*a*b**3*c*d + b**4*(8*c**2 - 3*d**2))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*sqrt(b)*f*(a**2 + b**2)**3*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1230():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*b**2*(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(9*d*f) - 4*b**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)*(-10*a*d + b*c)/(63*d**2*f) + 2*b*(3*a**2 - b**2)*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + (c - I*d)**(sympy.S(5)/2)*(I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f - (c + I*d)**(sympy.S(5)/2)*(I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*a**3*d + 6*a**2*b*c - 6*a*b**2*d - 2*b**3*c)/(3*f) + sqrt(c + d*tan(e + f*x))*(4*a**3*c*d + 6*a**2*b*(c**2 - d**2) - 12*a*b**2*c*d - 2*b**3*(c**2 - d**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1231():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 4*a*b*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 2*b**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f) - I*(a - I*b)**2*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + I*(a + I*b)**2*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*a**2*d + 4*a*b*c - 2*b**2*d)/(3*f) + sqrt(c + d*tan(e + f*x))*(a*c - b*d)*(4*a*d + 4*b*c)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1232():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*b*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) - (c - I*d)**(sympy.S(5)/2)*(I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + I*d)**(sympy.S(5)/2)*(I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*a*d + 2*b*c)/(3*f) + sqrt(c + d*tan(e + f*x))*(4*a*c*d + 2*b*(c**2 - d**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1233():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))
    F = (c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)) - (c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)) + 2*d**2*sqrt(c + d*tan(e + f*x))/(b*f) - 2*(-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1234():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**2
    F = I*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2) - sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) - (-a*d + b*c)**(sympy.S(3)/2)*(a**2*d + 4*a*b*c + 5*b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1235():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**3
    F = -(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + (c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3) - sqrt(c + d*tan(e + f*x))*(-a*d + b*c)*(a**2*d + 8*a*b*c + 9*b**2*d)/(4*b*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) + sqrt(-a*d + b*c)*(a**4*d**2 + 8*a**3*b*c*d - 6*a**2*b**2*(4*c**2 - 3*d**2) - 56*a*b**3*c*d + b**4*(8*c**2 - 15*d**2))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*f*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1236():
    f = (a + b*tan(e + f*x))**4/sqrt(c + d*tan(e + f*x))
    F = -4*b**3*sqrt(c + d*tan(e + f*x))*(-7*a*d + 2*b*c)*tan(e + f*x)/(15*d**2*f) + 2*b**2*(a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))/(5*d*f) - 2*b**2*sqrt(c + d*tan(e + f*x))*(-87*a**2*d**2 + 40*a*b*c*d - b**2*(8*c**2 - 15*d**2))/(15*d**3*f) - I*(a - I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) + I*(a + I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1237():
    f = (a + b*tan(e + f*x))**3/sqrt(c + d*tan(e + f*x))
    F = 2*b**2*(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/(3*d*f) - 4*b**2*sqrt(c + d*tan(e + f*x))*(-4*a*d + b*c)/(3*d**2*f) - (I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) + (I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1238():
    f = (a + b*tan(e + f*x))**2/sqrt(c + d*tan(e + f*x))
    F = 2*b**2*sqrt(c + d*tan(e + f*x))/(d*f) - I*(a - I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) + I*(a + I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1239():
    f = (a + b*tan(e + f*x))/sqrt(c + d*tan(e + f*x))
    F = (I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) - (I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1240():
    f = 1/((a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x)))
    F = -2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)*sqrt(-a*d + b*c)) - atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)*(I*a - b)) + atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)*(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1241():
    f = 1/((a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x)))
    F = -b**(sympy.S(3)/2)*(-5*a**2*d + 4*a*b*c - b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(3)/2)) - b**2*sqrt(c + d*tan(e + f*x))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c)) + I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*sqrt(c + I*d)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1242():
    f = (a + b*tan(e + f*x))**4/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*sqrt(c + d*tan(e + f*x))*(3*a*d*(-a*d + 2*b*c) - b**2*(4*c**2 + d**2))*tan(e + f*x)/(3*d**2*f*(c**2 + d**2)) - 2*b*sqrt(c + d*tan(e + f*x))*(-6*a**3*d**3 + 15*a**2*b*c*d**2 - 12*a*b**2*d*(2*c**2 + d**2) + b**3*(8*c**3 + 5*c*d**2))/(3*d**3*f*(c**2 + d**2)) - I*(a - I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) + I*(a + I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - 2*(a + b*tan(e + f*x))**2*(-a*d + b*c)**2/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1243():
    f = (a + b*tan(e + f*x))**3/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*b*sqrt(c + d*tan(e + f*x))*(a*d*(-a*d + 2*b*c) - b**2*(2*c**2 + d**2))/(d**2*f*(c**2 + d**2)) - (I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) + (I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) - 2*(a + b*tan(e + f*x))*(-a*d + b*c)**2/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1244():
    f = (a + b*tan(e + f*x))**2/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -I*(a - I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) + I*(a + I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - 2*(-a*d + b*c)**2/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1245():
    f = (a + b*tan(e + f*x))/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = (-2*a*d + 2*b*c)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)) + (I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - (I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1246():
    f = 1/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)*(-a*d + b*c)**(sympy.S(3)/2)) + 2*d**2/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) - atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)*(I*a - b)) + atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)*(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1247():
    f = 1/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -b**(sympy.S(5)/2)*(-7*a**2*d + 4*a*b*c - 3*b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(5)/2)) - b**2/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) - d*(2*a**2*d**2 + b**2*(c**2 + 3*d**2))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) + I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*(c + I*d)**(sympy.S(3)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1248():
    f = (a + b*tan(e + f*x))**4/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*b**2*sqrt(c + d*tan(e + f*x))*(a*d*(-a*d + 2*b*c) - b**2*(4*c**2 + 3*d**2))/(3*d**3*f*(c**2 + d**2)) - I*(a - I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**4*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - 2*(a + b*tan(e + f*x))**2*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) + 4*(-a*d + b*c)**3*(3*a*c*d + 2*b*c**2 + 5*b*d**2)/(3*d**3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1249():
    f = (a + b*tan(e + f*x))**3/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -(I*a - b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) + (I*a + b)**3*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) - 2*(a + b*tan(e + f*x))*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - 4*(-a*d + b*c)**2*(3*a*c*d + b*(c**2 + 4*d**2))/(3*d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1250():
    f = (a + b*tan(e + f*x))**2/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*(a - I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) + (a*c + b*d)*(-4*a*d + 4*b*c)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) - 2*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1251():
    f = (a + b*tan(e + f*x))/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -(4*a*c*d - 2*b*(c**2 - d**2))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) + (-2*a*d + 2*b*c)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c**2 + 3*d**2)) + (I*a - b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - (I*a + b)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1252():
    f = 1/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -2*b**(sympy.S(7)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)*(-a*d + b*c)**(sympy.S(5)/2)) - 2*d**2*(2*a*c*d - b*(3*c**2 + d**2))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + 2*d**2/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-3*a*d + 3*b*c)) - atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)*(I*a - b)) + atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)*(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1253():
    f = 1/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -b**(sympy.S(7)/2)*(-9*a**2*d + 4*a*b*c - 5*b**2*d)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(7)/2)) - b**2/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - d*(2*a**2*d**2 + b**2*(3*c**2 + 5*d**2))/(f*(3*a**2 + 3*b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-a*d + b*c)**2) + d*(4*a**3*c*d**3 - 4*a**2*b*d**2*(2*c**2 + d**2) + 4*a*b**2*c*d**3 - b**3*(c**4 + 10*c**2*d**2 + 5*d**4))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) + I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1254():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x))
    F = sqrt(b)*(15*a**2*d**2 + 10*a*b*c*d - b**2*(c**2 + 8*d**2))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*d**(sympy.S(3)/2)*f) + b**2*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(2*d*f) - b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-9*a*d + b*c)/(4*d*f) - I*(a - I*b)**(sympy.S(5)/2)*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*(a + I*b)**(sympy.S(5)/2)*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1255():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))
    F = sqrt(b)*(3*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(d)*f) + b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/f - I*(a - I*b)**(sympy.S(3)/2)*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*(a + I*b)**(sympy.S(3)/2)*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1256():
    f = sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))
    F = 2*sqrt(b)*sqrt(d)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/f - I*sqrt(a - I*b)*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*sqrt(a + I*b)*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1257():
    f = sqrt(c + d*tan(e + f*x))/sqrt(a + b*tan(e + f*x))
    F = I*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - I*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1258():
    f = sqrt(c + d*tan(e + f*x))/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*b*sqrt(c + d*tan(e + f*x))/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)) + I*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - I*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1259():
    f = sqrt(c + d*tan(e + f*x))/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*b*sqrt(c + d*tan(e + f*x))*(-5*a**2*d + 6*a*b*c + b**2*d)/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)) - 2*b*sqrt(c + d*tan(e + f*x))/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + I*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - I*sqrt(c - I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1260():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = b*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(2*f) - I*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(5*a*d + 3*b*c)/(4*f) + (3*a**2*d**2 + 18*a*b*c*d + b**2*(3*c**2 - 8*d**2))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*sqrt(b)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1261():
    f = sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = d*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/f - I*sqrt(a - I*b)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*sqrt(a + I*b)*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f + sqrt(d)*(a*d + 3*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1262():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/sqrt(a + b*tan(e + f*x))
    F = I*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)) + 2*d**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1263():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(c + d*tan(e + f*x))*(-2*a*d + 2*b*c)/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)) + I*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1264():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(c + d*tan(e + f*x))*(-4*a**2*d + 12*a*b*c + 8*b**2*d)/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(-2*a*d + 2*b*c)/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + I*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1265():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**(sympy.S(7)/2)
    F = sqrt(c + d*tan(e + f*x))*(-16*a**4*d**2 + 100*a**3*b*c*d - 2*a**2*b**2*(45*c**2 - 49*d**2) - 140*a*b**3*c*d + 6*b**4*(5*c**2 - d**2))/(15*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**3*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(-8*a**2*d + 20*a*b*c + 12*b**2*d)/(15*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(-2*a*d + 2*b*c)/(f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(5*a**2 + 5*b**2)) + I*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(7)/2)) - I*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1266():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f + d**2*(a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x))/(3*b*f) + d*(a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))*(-a*d + 13*b*c)/(12*b*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-a**2*d**2 + 14*a*b*c*d + b**2*(11*c**2 - 8*d**2))/(8*b*f) + (-a**3*d**3 + 15*a**2*b*c*d**2 + 3*a*b**2*d*(15*c**2 - 8*d**2) + 5*b**3*(c**3 - 8*c*d**2))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(8*b**(sympy.S(3)/2)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1267():
    f = sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*sqrt(a - I*b)*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + I*sqrt(a + I*b)*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f + d**2*(a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))/(2*b*f) + d*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-a*d + 9*b*c)/(4*b*f) + sqrt(d)*(-a**2*d**2 + 10*a*b*c*d + b**2*(15*c**2 - 8*d**2))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1268():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/sqrt(a + b*tan(e + f*x))
    F = I*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)) + d**2*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/(b*f) + d**(sympy.S(3)/2)*(-a*d + 5*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1269():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = I*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)) - 2*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2/(b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)) + 2*d**(sympy.S(5)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1270():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = I*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)) - sqrt(c + d*tan(e + f*x))*(-2*a*d + 2*b*c)*(a**2*d + 6*a*b*c + 7*b**2*d)/(3*b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2) - 2*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2/(3*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1271():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**(sympy.S(7)/2)
    F = I*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(7)/2)) - I*(c - I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(7)/2)) + sqrt(c + d*tan(e + f*x))*(4*a**4*d**2 + 40*a**3*b*c*d - 6*a**2*b**2*(15*c**2 - 13*d**2) - 200*a*b**3*c*d + 2*b**4*(15*c**2 - 23*d**2))/(15*b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**3) - sqrt(c + d*tan(e + f*x))*(-2*a*d + 2*b*c)*(a**2*d + 10*a*b*c + 11*b**2*d)/(15*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**2) - 2*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2/(5*b*f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1272():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)/sqrt(c + d*tan(e + f*x))
    F = -b**(sympy.S(3)/2)*(-5*a*d + b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) + b**2*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/(d*f) - I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) + I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1273():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)/sqrt(c + d*tan(e + f*x))
    F = 2*b**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(d)*f) - I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) + I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1274():
    f = sqrt(a + b*tan(e + f*x))/sqrt(c + d*tan(e + f*x))
    F = -I*sqrt(a - I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) + I*sqrt(a + I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1275():
    f = 1/(sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x)))
    F = I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*sqrt(c + I*d)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1276():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x)))
    F = -2*b**2*sqrt(c + d*tan(e + f*x))/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c)) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*sqrt(c + I*d)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1277():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x)))
    F = -4*b**2*sqrt(c + d*tan(e + f*x))*(-4*a**2*d + 3*a*b*c - b**2*d)/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)**2) - 2*b**2*sqrt(c + d*tan(e + f*x))/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*(-a*d + b*c)) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)*sqrt(c + I*d)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1278():
    f = (a + b*tan(e + f*x))**(sympy.S(7)/2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -b**(sympy.S(5)/2)*(-7*a*d + 3*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f) - b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(2*a*d*(-a*d + 2*b*c) - b**2*(3*c**2 + d**2))/(d**2*f*(c**2 + d**2)) - I*(a - I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) + I*(a + I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) - 2*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)**2/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1279():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*b**(sympy.S(5)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) - I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) + I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) - 2*sqrt(a + b*tan(e + f*x))*(-a*d + b*c)**2/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1280():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) + I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) + sqrt(a + b*tan(e + f*x))*(-2*a*d + 2*b*c)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1281():
    f = sqrt(a + b*tan(e + f*x))/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*d*sqrt(a + b*tan(e + f*x))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)) - I*sqrt(a - I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) + I*sqrt(a + I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1282():
    f = 1/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = 2*d**2*sqrt(a + b*tan(e + f*x))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*(c + I*d)**(sympy.S(3)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1283():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -2*b**2/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) - 2*d*sqrt(a + b*tan(e + f*x))*(a**2*d**2 + b**2*(c**2 + 2*d**2))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(3)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1284():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(5)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -4*b**2*(-5*a**2*d + 3*a*b*c - 2*b**2*d)/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2) - 2*b**2/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) + 2*d*sqrt(a + b*tan(e + f*x))*(3*a**4*d**3 + a**2*b**2*d*(11*c**2 + 17*d**2) - 6*a*b**3*c*(c**2 + d**2) + b**4*d*(5*c**2 + 8*d**2))/(3*f*(a**2 + b**2)**2*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**3) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)*(c + I*d)**(sympy.S(3)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1285():
    f = (a + b*tan(e + f*x))**(sympy.S(9)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -b**(sympy.S(7)/2)*(-9*a*d + 5*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(7)/2)*f) + b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(4*a**3*c*d**3 - 4*a**2*b*d**2*(c**2 - 2*d**2) - 4*a*b**2*c*d*(c**2 + 4*d**2) + b**3*(5*c**4 + 10*c**2*d**2 + d**4))/(d**3*f*(c**2 + d**2)**2) - I*(a - I*b)**(sympy.S(9)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**(sympy.S(9)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) - 2*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - 2*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)**2*(6*a*c*d + 5*b*c**2 + 11*b*d**2)/(3*d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1286():
    f = (a + b*tan(e + f*x))**(sympy.S(7)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*b**(sympy.S(7)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f) - I*(a - I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**(sympy.S(7)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) - 2*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - 2*sqrt(a + b*tan(e + f*x))*(-a*d + b*c)**2*(2*a*c*d + b*(c**2 + 3*d**2))/(d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1287():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*(a - I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) + sqrt(a + b*tan(e + f*x))*(-2*a*d + 2*b*c)*(6*a*c*d + b*(c**2 + 7*d**2))/(3*d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) - 2*sqrt(a + b*tan(e + f*x))*(-a*d + b*c)**2/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1288():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -I*(a - I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + I*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) + sqrt(a + b*tan(e + f*x))*(-12*a*c*d + 4*b*c**2 - 8*b*d**2)/(3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) + sqrt(a + b*tan(e + f*x))*(-2*a*d + 2*b*c)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c**2 + 3*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1289():
    f = sqrt(a + b*tan(e + f*x))/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*d*sqrt(a + b*tan(e + f*x))*(6*a*c*d - b*(5*c**2 - d**2))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-3*a*d + 3*b*c)) - 2*d*sqrt(a + b*tan(e + f*x))/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*c**2 + 3*d**2)) - I*sqrt(a - I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + I*sqrt(a + I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1290():
    f = 1/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -4*d**2*sqrt(a + b*tan(e + f*x))*(3*a*c*d - b*(4*c**2 + d**2))/(3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + 2*d**2*sqrt(a + b*tan(e + f*x))/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-3*a*d + 3*b*c)) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1291():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -2*b**2/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - 2*d*sqrt(a + b*tan(e + f*x))*(a**2*d**2 + b**2*(3*c**2 + 4*d**2))/(f*(3*a**2 + 3*b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-a*d + b*c)**2) + sqrt(a + b*tan(e + f*x))*(12*a**3*c*d**4 - 2*a**2*b*d**3*(11*c**2 + 5*d**2) + 12*a*b**2*c*d**4 - 2*b**3*(3*c**4*d + 17*c**2*d**3 + 8*d**5))/(f*(3*a**2 + 3*b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1292():
    f = 1/((a + b*tan(e + f*x))**(sympy.S(5)/2)*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -4*b**2*(-2*a**2*d + a*b*c - b**2*d)/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)**2) - 2*b**2/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - 4*d*sqrt(a + b*tan(e + f*x))*(3*a**5*c*d**4 - a**4*b*d**3*(7*c**2 + 4*d**2) + 6*a**3*b**2*c*d**4 - a**2*b**3*d*(7*c**4 + 28*c**2*d**2 + 15*d**4) + 3*a*b**4*c*(c**4 + 2*c**2*d**2 + 2*d**4) - b**5*d*(4*c**4 + 15*c**2*d**2 + 8*d**4))/(3*f*(a**2 + b**2)**2*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**4) + 2*d*sqrt(a + b*tan(e + f*x))*(a**4*d**3 + a**2*b**2*d*(13*c**2 + 15*d**2) - 6*a*b**3*c*(c**2 + d**2) + b**4*d*(7*c**2 + 8*d**2))/(3*f*(a**2 + b**2)**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-a*d + b*c)**3) + I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)*(c + I*d)**(sympy.S(5)/2)) - I*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1293():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n
    F = (a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**n*appellf1(m + 1, 1, -n, m + 2, (a + b*tan(e + f*x))/(a - I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(b*(c + d*tan(e + f*x))/(-a*d + b*c))**n*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**n*appellf1(m + 1, 1, -n, m + 2, (a + b*tan(e + f*x))/(a + I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(b*(c + d*tan(e + f*x))/(-a*d + b*c))**n*(m + 1)*(2*I*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1294():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**2
    F = (a + b*tan(e + f*x))**(m + 1)*(c - I*d)**2*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(c + I*d)**2*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(m + 1)*(2*I*a - 2*b)) + d**2*(a + b*tan(e + f*x))**(m + 1)/(b*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1295():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))
    F = (a + b*tan(e + f*x))**(m + 1)*(c - I*d)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*(I*c - d)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1296():
    f = (a + b*tan(e + f*x))**m
    F = -b*(a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + sqrt(-b**2)))/(2*f*sqrt(-b**2)*(a + sqrt(-b**2))*(m + 1)) + b*(a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - sqrt(-b**2)))/(2*f*sqrt(-b**2)*(a - sqrt(-b**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1297():
    f = (a + b*tan(e + f*x))**m/(c + d*tan(e + f*x))
    F = d**2*(a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(c**2 + d**2)*(m + 1)*(-a*d + b*c)) + (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(c - I*d)*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1)*(I*c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1298():
    f = (a + b*tan(e + f*x))**m/(c + d*tan(e + f*x))**2
    F = -d**2*(a + b*tan(e + f*x))**(m + 1)*(2*a*c*d - b*(c**2*(2 - m) - d**2*m))*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(c**2 + d**2)**2*(m + 1)*(-a*d + b*c)**2) + d**2*(a + b*tan(e + f*x))**(m + 1)/(f*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) - (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(c + I*d)**2*(m + 1)*(2*I*a - 2*b)) + (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(c - I*d)**2*(m + 1)*(2*I*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1299():
    f = (a + b*tan(e + f*x))**m/(c + d*tan(e + f*x))**3
    F = d**2*(a + b*tan(e + f*x))**(m + 1)*(2*a**2*d**2*(3*c**2 - d**2) - 4*a*b*c*d*(c**2*(3 - m) - d**2*(m + 1)) - b**2*(-c**4*(m**2 - 5*m + 6) + 2*c**2*d**2*(-m**2 + 3*m + 1) + d**4*m*(1 - m)))*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(2*f*(c**2 + d**2)**3*(m + 1)*(-a*d + b*c)**3) - d**2*(a + b*tan(e + f*x))**(m + 1)*(4*a*c*d - b*(c**2*(5 - m) + d**2*(1 - m)))/(2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + d**2*(a + b*tan(e + f*x))**(m + 1)/(f*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-2*a*d + 2*b*c)) + (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(c - I*d)**3*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1)*(I*c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1300():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -(a + b*tan(e + f*x))**(m + 1)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a + I*b))/(b*f*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(m + 1)*(2*I*a - 2*b)) + (a + b*tan(e + f*x))**(m + 1)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)*appellf1(m + 1, sympy.S(-3)/2, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a - I*b))/(2*b*f*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(m + 1)*(I*a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1301():
    f = (a + b*tan(e + f*x))**m*sqrt(c + d*tan(e + f*x))
    F = (a + b*tan(e + f*x))**(m + 1)*sqrt(c + d*tan(e + f*x))*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a - I*b))/(f*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*sqrt(c + d*tan(e + f*x))*appellf1(m + 1, sympy.S(-1)/2, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a + I*b))/(f*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(m + 1)*(2*I*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1302():
    f = (a + b*tan(e + f*x))**m/sqrt(c + d*tan(e + f*x))
    F = sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a - I*b))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a + 2*b)) - sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, sympy.S.Half, 1, m + 2, -d*(a + b*tan(e + f*x))/(-a*d + b*c), (a + b*tan(e + f*x))/(a + I*b))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1303():
    f = (a + b*tan(e + f*x))**m/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = b*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, (a + b*tan(e + f*x))/(a - I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a + 2*b)*(-a*d + b*c)) - b*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, 1, sympy.S(3)/2, m + 2, (a + b*tan(e + f*x))/(a + I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a - 2*b)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1304():
    f = (a + b*tan(e + f*x))**m/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = b**2*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, 1, sympy.S(5)/2, m + 2, (a + b*tan(e + f*x))/(a - I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a + 2*b)*(-a*d + b*c)**2) - b**2*sqrt(b*(c + d*tan(e + f*x))/(-a*d + b*c))*(a + b*tan(e + f*x))**(m + 1)*appellf1(m + 1, 1, sympy.S(5)/2, m + 2, (a + b*tan(e + f*x))/(a + I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*sqrt(c + d*tan(e + f*x))*(m + 1)*(2*I*a - 2*b)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1305():
    f = (c*(d*tan(e + f*x))**p)**n*(I*a*tan(e + f*x) + a)**m
    F = (c*(d*tan(e + f*x))**p)**n*(I*a*tan(e + f*x) + a)**m*tan(e + f*x)*appellf1(n*p + 1, 1, 1 - m, n*p + 2, I*tan(e + f*x), -I*tan(e + f*x))/(f*(I*tan(e + f*x) + 1)**m*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1306():
    f = (c*(d*tan(e + f*x))**p)**n*(I*a*tan(e + f*x) + a)**3
    F = -I*a**3*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2/(f*(n*p + 2)) + 4*a**3*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), I*tan(e + f*x))/(f*(n*p + 1)) - 3*a**3*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1307():
    f = (c*(d*tan(e + f*x))**p)**n*(I*a*tan(e + f*x) + a)**2
    F = 2*a**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), I*tan(e + f*x))/(f*(n*p + 1)) - a**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1308():
    f = (c*(d*tan(e + f*x))**p)**n*(I*a*tan(e + f*x) + a)
    F = a*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), I*tan(e + f*x))/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1309():
    f = (c*(d*tan(e + f*x))**p)**n/(I*a*tan(e + f*x) + a)
    F = -I*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2*hyper((2, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(a*f*(n*p + 2)) + (c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((2, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(a*f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1310():
    f = (c*(d*tan(e + f*x))**p)**n/(I*a*tan(e + f*x) + a)**2
    F = (c*(d*tan(e + f*x))**p)**n*(2*n**2*p**2 - 4*n*p + 1)*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), -I*tan(e + f*x))/(8*a**2*f*(n*p + 1)) + (c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), I*tan(e + f*x))/(8*a**2*f*(n*p + 1)) + (c*(d*tan(e + f*x))**p)**n*(-n*p + 2)*tan(e + f*x)/(4*a**2*f*(I*tan(e + f*x) + 1)) + (c*(d*tan(e + f*x))**p)**n*tan(e + f*x)/(4*a**2*f*(I*tan(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1311():
    f = (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))**m
    F = (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))**m*tan(e + f*x)*appellf1(n*p + 1, 1, -m, n*p + 2, -I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*f*(1 + b*tan(e + f*x)/a)**m*(n*p + 1)) + (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))**m*tan(e + f*x)*appellf1(n*p + 1, 1, -m, n*p + 2, I*tan(e + f*x), -b*tan(e + f*x)/a)/(2*f*(1 + b*tan(e + f*x)/a)**m*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1312():
    f = (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))**3
    F = 3*a*b**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1)) + a*(c*(d*tan(e + f*x))**p)**n*(a**2 - 3*b**2)*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1)) + b**3*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2/(f*(n*p + 2)) + b*(c*(d*tan(e + f*x))**p)**n*(3*a**2 - b**2)*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1313():
    f = (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))**2
    F = 2*a*b*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(n*p + 2)) + b**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1)) + (c*(d*tan(e + f*x))**p)**n*(a**2 - b**2)*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1314():
    f = (c*(d*tan(e + f*x))**p)**n*(a + b*tan(e + f*x))
    F = a*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(n*p + 1)) + b*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1315():
    f = (c*(d*tan(e + f*x))**p)**n/(a + b*tan(e + f*x))
    F = a*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(a**2 + b**2)*(n*p + 1)) - b*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(a**2 + b**2)*(n*p + 2)) + b**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), -b*tan(e + f*x)/a)/(a*f*(a**2 + b**2)*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_2_1_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_1316():
    f = (c*(d*tan(e + f*x))**p)**n/(a + b*tan(e + f*x))**2
    F = -2*a*b*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)**2*hyper((1, n*p/2 + 1), (n*p/2 + 2,), -tan(e + f*x)**2)/(f*(a**2 + b**2)**2*(n*p + 2)) + 2*b**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((1, n*p + 1), (n*p + 2,), -b*tan(e + f*x)/a)/(f*(a**2 + b**2)**2*(n*p + 1)) + (c*(d*tan(e + f*x))**p)**n*(a**2 - b**2)*tan(e + f*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(a**2 + b**2)**2*(n*p + 1)) + b**2*(c*(d*tan(e + f*x))**p)**n*tan(e + f*x)*hyper((2, n*p + 1), (n*p + 2,), -b*tan(e + f*x)/a)/(a**2*f*(a**2 + b**2)*(n*p + 1))
    assert integrate(f, x) == F

