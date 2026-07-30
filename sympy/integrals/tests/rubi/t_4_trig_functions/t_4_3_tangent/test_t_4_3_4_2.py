"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.4.2 (a+b tan)^m (c+d tan)^n (A+B tan+C tan^2).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, m, n = symbols('A B C a b c d e f m n')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_1():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)
    F = C*b*tan(c + d*x)**3/(3*d) - x*(B*a - C*b) + (B*a - C*b)*tan(c + d*x)/d + (B*b + C*a)*log(cos(c + d*x))/d + (B*b + C*a)*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_2():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)
    F = B*b*tan(c + d*x)/d + C*(a + b*tan(c + d*x))**2/(2*b*d) - x*(B*b + C*a) - (B*a - C*b)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_3():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)
    F = C*b*tan(c + d*x)/d + x*(B*a - C*b) - (B*b + C*a)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_4():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2
    F = B*a*log(sin(c + d*x))/d - C*b*log(cos(c + d*x))/d + x*(B*b + C*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_5():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3
    F = -B*a*cot(c + d*x)/d - x*(B*a - C*b) + (B*b + C*a)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_6():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**4
    F = -B*a*cot(c + d*x)**2/(2*d) - x*(B*b + C*a) - (B*a - C*b)*log(sin(c + d*x))/d - (B*b + C*a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_7():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**5
    F = -B*a*cot(c + d*x)**3/(3*d) + x*(B*a - C*b) + (B*a - C*b)*cot(c + d*x)/d - (B*b + C*a)*log(sin(c + d*x))/d - (B*b + C*a)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_8():
    f = (a + b*tan(c + d*x))*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**6
    F = -B*a*cot(c + d*x)**4/(4*d) + x*(B*b + C*a) + (B*a - C*b)*log(sin(c + d*x))/d + (B*a - C*b)*cot(c + d*x)**2/(2*d) - (B*b + C*a)*cot(c + d*x)**3/(3*d) + (B*b + C*a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_9():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)
    F = -C*(a + b*tan(c + d*x))**2/(2*d) + C*(a + b*tan(c + d*x))**3*tan(c + d*x)/(4*b*d) - b*(B*b + C*a)*tan(c + d*x)/d - x*(B*a**2 - B*b**2 - 2*C*a*b) + (2*B*a*b + C*a**2 - C*b**2)*log(cos(c + d*x))/d + (a + b*tan(c + d*x))**3*(4*B*b - C*a)/(12*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_10():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)
    F = B*(a + b*tan(c + d*x))**2/(2*d) + C*(a + b*tan(c + d*x))**3/(3*b*d) + b*(B*a - C*b)*tan(c + d*x)/d - x*(2*B*a*b + C*a**2 - C*b**2) - (B*a**2 - B*b**2 - 2*C*a*b)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_11():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)
    F = C*(a + b*tan(c + d*x))**2/(2*d) + b*(B*b + C*a)*tan(c + d*x)/d + x*(B*a**2 - B*b**2 - 2*C*a*b) - (2*B*a*b + C*a**2 - C*b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_12():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2
    F = B*a**2*log(sin(c + d*x))/d + C*b**2*tan(c + d*x)/d - b*(B*b + 2*C*a)*log(cos(c + d*x))/d + x*(2*B*a*b + C*a**2 - C*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_13():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3
    F = -B*a**2*cot(c + d*x)/d - C*b**2*log(cos(c + d*x))/d + a*(2*B*b + C*a)*log(sin(c + d*x))/d - x*(B*a**2 - B*b**2 - 2*C*a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_14():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**4
    F = -B*a**2*cot(c + d*x)**2/(2*d) - a*(2*B*b + C*a)*cot(c + d*x)/d + x*(C*b**2 - a*(2*B*b + C*a)) - (B*a**2 - B*b**2 - 2*C*a*b)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_15():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**5
    F = -B*a**2*cot(c + d*x)**3/(3*d) - a*(2*B*b + C*a)*cot(c + d*x)**2/(2*d) + x*(B*a**2 - B*b**2 - 2*C*a*b) + (C*b**2 - a*(2*B*b + C*a))*log(sin(c + d*x))/d + (B*a**2 - B*b**2 - 2*C*a*b)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_16():
    f = (a + b*tan(c + d*x))**2*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**6
    F = -B*a**2*cot(c + d*x)**4/(4*d) - a*(2*B*b + C*a)*cot(c + d*x)**3/(3*d) + x*(2*B*a*b + C*a**2 - C*b**2) - (C*b**2 - a*(2*B*b + C*a))*cot(c + d*x)/d + (B*a**2 - B*b**2 - 2*C*a*b)*log(sin(c + d*x))/d + (B*a**2 - B*b**2 - 2*C*a*b)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_17():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)
    F = B*(a + b*tan(c + d*x))**3/(3*d) + C*(a + b*tan(c + d*x))**4/(4*b*d) + b*(B*a**2 - B*b**2 - 2*C*a*b)*tan(c + d*x)/d - x*(3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2) + (a + b*tan(c + d*x))**2*(B*a - C*b)/(2*d) - (B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_18():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)
    F = C*(a + b*tan(c + d*x))**3/(3*d) + b*(2*B*a*b + C*a**2 - C*b**2)*tan(c + d*x)/d + x*(B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3) + (a + b*tan(c + d*x))**2*(B*b + C*a)/(2*d) - (3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_19():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2
    F = B*a**3*log(sin(c + d*x))/d + C*b*(a + b*tan(c + d*x))**2/(2*d) + b**2*(B*b + 2*C*a)*tan(c + d*x)/d - b*(3*B*a*b + 3*C*a**2 - C*b**2)*log(cos(c + d*x))/d + x*(3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_20():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3
    F = -B*a*(a + b*tan(c + d*x))**2*cot(c + d*x)/d + a**2*(3*B*b + C*a)*log(sin(c + d*x))/d + b**2*(B*a + C*b)*tan(c + d*x)/d - b**2*(B*b + 3*C*a)*log(cos(c + d*x))/d - x*(B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_21():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**4
    F = -B*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**2/(2*d) - C*b**3*log(cos(c + d*x))/d - a**2*(2*B*b + C*a)*cot(c + d*x)/d - a*(B*a**2 - 3*B*b**2 - 3*C*a*b)*log(sin(c + d*x))/d - x*(3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_22():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**5
    F = -B*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**3/(3*d) - a**2*(5*B*b + 3*C*a)*cot(c + d*x)**2/(6*d) + a*(3*B*a**2 - 8*B*b**2 - 9*C*a*b)*cot(c + d*x)/(3*d) + x*(B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3) - (3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_23():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**6
    F = -B*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**4/(4*d) - a**2*(3*B*b + 2*C*a)*cot(c + d*x)**3/(6*d) + a*(2*B*a**2 - 5*B*b**2 - 6*C*a*b)*cot(c + d*x)**2/(4*d) + x*(3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2) + (B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3)*log(sin(c + d*x))/d + (3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_24():
    f = (a + b*tan(c + d*x))**3*(B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**7
    F = -B*a*(a + b*tan(c + d*x))**2*cot(c + d*x)**5/(5*d) - a**2*(7*B*b + 5*C*a)*cot(c + d*x)**4/(20*d) + a*(5*B*a**2 - 12*B*b**2 - 15*C*a*b)*cot(c + d*x)**3/(15*d) - x*(B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3) - (B*a**3 - 3*B*a*b**2 - 3*C*a**2*b + C*b**3)*cot(c + d*x)/d + (3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)*log(sin(c + d*x))/d + (3*B*a**2*b - B*b**3 + C*a**3 - 3*C*a*b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_25():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**2/(a + b*tan(c + d*x))
    F = C*tan(c + d*x)**2/(2*b*d) - a**3*(B*b - C*a)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)) - x*(B*b - C*a)/(a**2 + b**2) + (B*a + C*b)*log(cos(c + d*x))/(d*(a**2 + b**2)) + (B*b - C*a)*tan(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_26():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)/(a + b*tan(c + d*x))
    F = C*tan(c + d*x)/(b*d) + a**2*(B*b - C*a)*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)) - x*(B*a + C*b)/(a**2 + b**2) - (B*b - C*a)*log(cos(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_27():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)/(a + b*tan(c + d*x))
    F = -a*(B*b - C*a)*log(a + b*tan(c + d*x))/(b*d*(a**2 + b**2)) + x*(B*b - C*a)/(a**2 + b**2) - (B*a + C*b)*log(cos(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_28():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)/(a + b*tan(c + d*x))
    F = x*(B*a + C*b)/(a**2 + b**2) + (B*b - C*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_29():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2/(a + b*tan(c + d*x))
    F = B*log(sin(c + d*x))/(a*d) - x*(B*b - C*a)/(a**2 + b**2) - b*(B*b - C*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_30():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3/(a + b*tan(c + d*x))
    F = -B*cot(c + d*x)/(a*d) - x*(B*a + C*b)/(a**2 + b**2) + b**2*(B*b - C*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)) - (B*b - C*a)*log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_31():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**4/(a + b*tan(c + d*x))
    F = -B*cot(c + d*x)**2/(2*a*d) + x*(B*b - C*a)/(a**2 + b**2) + (B*b - C*a)*cot(c + d*x)/(a**2*d) - b**3*(B*b - C*a)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)) - (B*a**2 - B*b**2 + C*a*b)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_32():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = a**2*(B*a**2*b + 3*B*b**3 - 2*C*a**3 - 4*C*a*b**2)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**2) + a*(B*b - C*a)*tan(c + d*x)**2/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - x*(2*B*a*b - C*a**2 + C*b**2)/(a**2 + b**2)**2 + (B*a**2 - B*b**2 + 2*C*a*b)*log(cos(c + d*x))/(d*(a**2 + b**2)**2) - (B*a*b - 2*C*a**2 - C*b**2)*tan(c + d*x)/(b**2*d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_33():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)/(a + b*tan(c + d*x))**2
    F = -a**2*(B*b - C*a)/(b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - a*(2*B*b**3 - C*a**3 - 3*C*a*b**2)*log(a + b*tan(c + d*x))/(b**2*d*(a**2 + b**2)**2) - x*(B*a**2 - B*b**2 + 2*C*a*b)/(a**2 + b**2)**2 - (2*B*a*b - C*a**2 + C*b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_34():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)/(a + b*tan(c + d*x))**2
    F = a*(B*b - C*a)/(b*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + x*(2*B*a*b - C*a**2 + C*b**2)/(a**2 + b**2)**2 - (B*a**2 - B*b**2 + 2*C*a*b)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_35():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)/(a + b*tan(c + d*x))**2
    F = x*(B*a**2 - B*b**2 + 2*C*a*b)/(a**2 + b**2)**2 + (2*B*a*b - C*a**2 + C*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) - (B*b - C*a)/(d*(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_36():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = B*log(sin(c + d*x))/(a**2*d) - x*(2*B*a*b - C*a**2 + C*b**2)/(a**2 + b**2)**2 + b*(B*b - C*a)/(a*d*(a + b*tan(c + d*x))*(a**2 + b**2)) - b*(3*B*a**2*b + B*b**3 - 2*C*a**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_37():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = -B*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))) - x*(B*a**2 - B*b**2 + 2*C*a*b)/(a**2 + b**2)**2 - b*(B*a**2 + 2*B*b**2 - C*a*b)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)) + b**2*(4*B*a**2*b + 2*B*b**3 - 3*C*a**3 - C*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**2) - (2*B*b - C*a)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_38():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = a**2*(B*a**4*b + 3*B*a**2*b**3 + 6*B*b**5 - 3*C*a**5 - 9*C*a**3*b**2 - 10*C*a*b**4)*log(a + b*tan(c + d*x))/(b**4*d*(a**2 + b**2)**3) + a*(B*b - C*a)*tan(c + d*x)**3/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(B*a**2*b + 5*B*b**3 - 3*C*a**3 - 7*C*a*b**2)*tan(c + d*x)**2/(2*b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + x*(B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)/(a**2 + b**2)**3 + (3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)*log(cos(c + d*x))/(d*(a**2 + b**2)**3) - (B*a**3*b + 3*B*a*b**3 - 3*C*a**4 - 6*C*a**2*b**2 - C*b**4)*tan(c + d*x)/(b**3*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_39():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -a**2*(2*B*b**3 - C*a**3 - 3*C*a*b**2)/(b**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + a*(B*b - C*a)*tan(c + d*x)**2/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(B*a**2*b**3 - 3*B*b**5 + C*a**5 + 3*C*a**3*b**2 + 6*C*a*b**4)*log(a + b*tan(c + d*x))/(b**3*d*(a**2 + b**2)**3) - x*(3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)/(a**2 + b**2)**3 + (B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)*log(cos(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_40():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)/(a + b*tan(c + d*x))**3
    F = -a**2*(B*b - C*a)/(2*b**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + a*(2*B*b**3 - C*a**3 - 3*C*a*b**2)/(b**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - x*(B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)/(a**2 + b**2)**3 - (3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_41():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)/(a + b*tan(c + d*x))**3
    F = a*(B*b - C*a)/(2*b*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + x*(3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)/(a**2 + b**2)**3 - (B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + (B*a**2 - B*b**2 + 2*C*a*b)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_42():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)/(a + b*tan(c + d*x))**3
    F = x*(B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)/(a**2 + b**2)**3 + (3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) - (2*B*a*b - C*a**2 + C*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - (B*b - C*a)/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_43():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = B*log(sin(c + d*x))/(a**3*d) - x*(3*B*a**2*b - B*b**3 - C*a**3 + 3*C*a*b**2)/(a**2 + b**2)**3 + b*(B*b - C*a)/(2*a*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) + b*(3*B*a**2*b + B*b**3 - 2*C*a**3)/(a**2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b*(6*B*a**4*b + 3*B*a**2*b**3 + B*b**5 - 3*C*a**5 + C*a**3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**3*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_44():
    f = (B*tan(c + d*x) + C*tan(c + d*x)**2)*cot(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = -B*cot(c + d*x)/(a*d*(a + b*tan(c + d*x))**2) - x*(B*a**3 - 3*B*a*b**2 + 3*C*a**2*b - C*b**3)/(a**2 + b**2)**3 - b*(2*B*a**2 + 3*B*b**2 - C*a*b)/(2*a**2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)) - b*(B*a**4 + 6*B*a**2*b**2 + 3*B*b**4 - 3*C*a**3*b - C*a*b**3)/(a**3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + b**2*(10*B*a**4*b + 9*B*a**2*b**3 + 3*B*b**5 - 6*C*a**5 - 3*C*a**3*b**2 - C*a*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(a**4*d*(a**2 + b**2)**3) - (3*B*b - C*a)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_45():
    f = (b*tan(c + d*x))**n*(A + B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**2
    F = B*(b*tan(c + d*x))**(n + 4)*hyper((1, n/2 + 2), (n/2 + 3,), -tan(c + d*x)**2)/(b**4*d*(n + 4)) + C*(b*tan(c + d*x))**(n + 3)/(b**3*d*(n + 3)) + (b*tan(c + d*x))**(n + 3)*(A - C)*hyper((1, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -tan(c + d*x)**2)/(b**3*d*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_46():
    f = (b*tan(c + d*x))**n*(A + B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**m
    F = B*(b*tan(c + d*x))**n*tan(c + d*x)**(m + 2)*hyper((1, m/2 + n/2 + 1), (m/2 + n/2 + 2,), -tan(c + d*x)**2)/(d*(m + n + 2)) + C*(b*tan(c + d*x))**n*tan(c + d*x)**(m + 1)/(d*(m + n + 1)) + (b*tan(c + d*x))**n*(A - C)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_47():
    f = sqrt(b*tan(c + d*x))*(A + B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**m
    F = 2*B*sqrt(b*tan(c + d*x))*tan(c + d*x)**(m + 2)*hyper((1, m/2 + sympy.S(5)/4), (m/2 + sympy.S(9)/4,), -tan(c + d*x)**2)/(d*(2*m + 5)) + 2*C*sqrt(b*tan(c + d*x))*tan(c + d*x)**(m + 1)/(d*(2*m + 3)) + sqrt(b*tan(c + d*x))*(2*A - 2*C)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S(3)/4), (m/2 + sympy.S(7)/4,), -tan(c + d*x)**2)/(d*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_48():
    f = (A + B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**m/sqrt(b*tan(c + d*x))
    F = 2*B*tan(c + d*x)**(m + 2)*hyper((1, m/2 + sympy.S(3)/4), (m/2 + sympy.S(7)/4,), -tan(c + d*x)**2)/(d*sqrt(b*tan(c + d*x))*(2*m + 3)) + 2*C*tan(c + d*x)**(m + 1)/(d*sqrt(b*tan(c + d*x))*(2*m + 1)) + (2*A - 2*C)*tan(c + d*x)**(m + 1)*hyper((1, m/2 + sympy.S(1)/4), (m/2 + sympy.S(5)/4,), -tan(c + d*x)**2)/(d*sqrt(b*tan(c + d*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_49():
    f = (A + B*tan(c + d*x) + C*tan(c + d*x)**2)*tan(c + d*x)**m/sqrt(a + b*tan(c + d*x))
    F = 2*C*sqrt(a + b*tan(c + d*x))*tan(c + d*x)**m*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), 1 + b*tan(c + d*x)/a)/(b*d*(-b*tan(c + d*x)/a)**m) - sqrt(a + b*tan(c + d*x))*(B*b - sqrt(-b**2)*(A - C))*tan(c + d*x)**m*appellf1(sympy.S.Half, 1, -m, sympy.S(3)/2, (a + b*tan(c + d*x))/(a + sqrt(-b**2)), 1 + b*tan(c + d*x)/a)/(b*d*(-b*tan(c + d*x)/a)**m*(a + sqrt(-b**2))) - sqrt(a + b*tan(c + d*x))*(B*b + sqrt(-b**2)*(A - C))*tan(c + d*x)**m*appellf1(sympy.S.Half, 1, -m, sympy.S(3)/2, (a + b*tan(c + d*x))/(a - sqrt(-b**2)), 1 + b*tan(c + d*x)/a)/(b*d*(-b*tan(c + d*x)/a)**m*(a - sqrt(-b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_50():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*d*(a + b*tan(e + f*x))**4*tan(e + f*x)/(5*b*f) + b*(a**2*(B*c + d*(A - C)) + 2*a*b*(A*c - B*d - C*c) - b**2*(B*c + d*(A - C)))*tan(e + f*x)/f + x*(a**3*(A*c - B*d - C*c) - 3*a**2*b*(B*c + d*(A - C)) - 3*a*b**2*(A*c - B*d - C*c) + b**3*(B*c + d*(A - C))) + (a + b*tan(e + f*x))**3*(B*c + d*(A - C))/(3*f) + (a + b*tan(e + f*x))**2*(A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c)/(2*f) - (a**3*(B*c + d*(A - C)) + 3*a**2*b*(A*c - B*d - C*c) - 3*a*b**2*(B*c + d*(A - C)) - b**3*(A*c - B*d - C*c))*log(cos(e + f*x))/f - (a + b*tan(e + f*x))**4*(C*a*d - 5*b*(B*d + C*c))/(20*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_51():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*d*(a + b*tan(e + f*x))**3*tan(e + f*x)/(4*b*f) + b*(A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c)*tan(e + f*x)/f + x*(a**2*(A*c - B*d - C*c) - 2*a*b*(B*c + d*(A - C)) - b**2*(A*c - B*d - C*c)) + (a + b*tan(e + f*x))**2*(B*c + d*(A - C))/(2*f) - (a**2*(B*c + d*(A - C)) + 2*a*b*(A*c - B*d - C*c) - b**2*(B*c + d*(A - C)))*log(cos(e + f*x))/f - (a + b*tan(e + f*x))**3*(C*a*d - 4*b*(B*d + C*c))/(12*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_52():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*b*(c + d*tan(e + f*x))**2*tan(e + f*x)/(3*d*f) + d*(A*b + B*a - C*b)*tan(e + f*x)/f + x*(a*(A*c - B*d - C*c) - b*(B*c + d*(A - C))) - (A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c)*log(cos(e + f*x))/f - (c + d*tan(e + f*x))**2*(-3*B*b*d - 3*C*a*d + C*b*c)/(6*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_53():
    f = (c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = B*d*tan(e + f*x)/f + C*(c + d*tan(e + f*x))**2/(2*d*f) + x*(A*c - B*d - C*c) - (B*c + d*(A - C))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_54():
    f = (c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = C*d*tan(e + f*x)/(b*f) + x*(a*(A*c - B*d - C*c) + b*(B*c + d*(A - C)))/(a**2 + b**2) + (-A*a*d + A*b*c - B*a*c - B*b*d + C*a*d - C*b*c)*log(cos(e + f*x))/(f*(a**2 + b**2)) + (A*b**2 - a*(B*b - C*a))*(-a*d + b*c)*log(a + b*tan(e + f*x))/(b**2*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_55():
    f = (c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = x*(a**2*(A*c - B*d - C*c) + 2*a*b*(B*c + d*(A - C)) - b**2*(A*c - B*d - C*c))/(a**2 + b**2)**2 + (-a**2*(B*c + d*(A - C)) + 2*a*b*(A*c - B*d - C*c) + b**2*(B*c + d*(A - C)))*log(cos(e + f*x))/(f*(a**2 + b**2)**2) + (C*a**4*d - a**2*b**2*(B*c + d*(A - 3*C)) + 2*a*b**3*(A*c - B*d - C*c) + b**4*(A*d + B*c))*log(a + b*tan(e + f*x))/(b**2*f*(a**2 + b**2)**2) - (A*b**2 - a*(B*b - C*a))*(-a*d + b*c)/(b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_56():
    f = (c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = x*(a**3*(A*c - B*d - C*c) + 3*a**2*b*(B*c + d*(A - C)) - 3*a*b**2*(A*c - B*d - C*c) - b**3*(B*c + d*(A - C)))/(a**2 + b**2)**3 + (-a**3*(B*c + d*(A - C)) + 3*a**2*b*(A*c - B*d - C*c) + 3*a*b**2*(B*c + d*(A - C)) - b**3*(A*c - B*d - C*c))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3) - (C*a**4*d - a**2*b**2*(B*c + d*(A - 3*C)) + 2*a*b**3*(A*c - B*d - C*c) + b**4*(A*d + B*c))/(b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - (A*b**2 - a*(B*b - C*a))*(-a*d + b*c)/(2*b**2*f*(a + b*tan(e + f*x))**2*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_57():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**3/(6*d*f) + b*(c + d*tan(e + f*x))**3*(5*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-2*B*b*d - C*a*d + C*b*c))*tan(e + f*x)/(20*d**3*f) + d*(a**3*(B*c + d*(A - C)) + 3*a**2*b*(A*c - B*d - C*c) - 3*a*b**2*(B*c + d*(A - C)) - b**3*(A*c - B*d - C*c))*tan(e + f*x)/f - x*(a**3*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + 3*a**2*b*(B*(c**2 - d**2) + 2*c*d*(A - C)) - 3*a*b**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - b**3*(B*(c**2 - d**2) + 2*c*d*(A - C))) + (c + d*tan(e + f*x))**2*(B*a**3 - 3*B*a*b**2 + 3*a**2*b*(A - C) - b**3*(A - C))/(2*f) + (-a**3*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 3*a**2*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + 3*a*b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) - b**3*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/f - (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3*(-2*B*b*d - C*a*d + C*b*c)/(10*d**2*f) + (c + d*tan(e + f*x))**3*(4*C*a**3*d**3 - 3*a**2*b*d**2*(-16*B*d + 3*C*c) + 3*a*b**2*d*(-5*B*c*d + 2*C*c**2 + d**2*(20*A - 20*C)) - b**3*(-2*B*c**2*d + 20*B*d**3 + C*c**3 + 5*c*d**2*(A - C)))/(60*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_58():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3/(5*d*f) - b*(c + d*tan(e + f*x))**3*(-5*B*b*d - 2*C*a*d + 2*C*b*c)*tan(e + f*x)/(20*d**2*f) + d*(a**2*(B*c + d*(A - C)) + 2*a*b*(A*c - B*d - C*c) - b**2*(B*c + d*(A - C)))*tan(e + f*x)/f - x*(a**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + 2*a*b*(B*(c**2 - d**2) + 2*c*d*(A - C)) - b**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2)) + (c + d*tan(e + f*x))**2*(B*a**2 - B*b**2 + 2*a*b*(A - C))/(2*f) + (-a**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*a*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)))*log(cos(e + f*x))/f + (c + d*tan(e + f*x))**3*(8*C*a**2*d**2 - 10*a*b*d*(-4*B*d + C*c) + b**2*(-5*B*c*d + 2*C*c**2 + d**2*(20*A - 20*C)))/(60*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_59():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*b*(c + d*tan(e + f*x))**3*tan(e + f*x)/(4*d*f) + d*(A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c)*tan(e + f*x)/f - x*(a*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + b*(B*(c**2 - d**2) + 2*c*d*(A - C))) + (c + d*tan(e + f*x))**2*(A*b + B*a - C*b)/(2*f) - (A*(2*a*c*d + b*(c**2 - d**2)) + a*(B*c**2 - B*d**2 - 2*C*c*d) - b*(2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/f - (c + d*tan(e + f*x))**3*(-4*B*b*d - 4*C*a*d + C*b*c)/(12*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_60():
    f = (c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = B*(c + d*tan(e + f*x))**2/(2*f) + C*(c + d*tan(e + f*x))**3/(3*d*f) + d*(B*c + d*(A - C))*tan(e + f*x)/f - x*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - (B*(c**2 - d**2) + 2*c*d*(A - C))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_61():
    f = (c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = C*(c + d*tan(e + f*x))**2/(2*b*f) - x*(a*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - b*(B*(c**2 - d**2) + 2*c*d*(A - C)))/(a**2 + b**2) - (A*(2*a*c*d - b*(c**2 - d**2)) + a*(B*c**2 - B*d**2 - 2*C*c*d) + b*(2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/(f*(a**2 + b**2)) + d*(B*b*d - C*a*d + C*b*c)*tan(e + f*x)/(b**2*f) + (A*b**2 - a*(B*b - C*a))*(-a*d + b*c)**2*log(a + b*tan(e + f*x))/(b**3*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_62():
    f = (c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = -x*(a**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - 2*a*b*(B*(c**2 - d**2) + 2*c*d*(A - C)) - b**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2))/(a**2 + b**2)**2 - (a**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*a*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)))*log(cos(e + f*x))/(f*(a**2 + b**2)**2) - (c + d*tan(e + f*x))**2*(A*b**2 - a*(B*b - C*a))/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + d**2*(A*b**2 - B*a*b + 2*C*a**2 + C*b**2)*tan(e + f*x)/(b**2*f*(a**2 + b**2)) - (-a*d + b*c)*(B*a**3*b*d - 2*C*a**4*d + a**2*b**2*(B*c - 4*C*d) - a*b**3*(2*A*c - 3*B*d - 2*C*c) - b**4*(2*A*d + B*c))*log(a + b*tan(e + f*x))/(b**3*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_63():
    f = (c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = -x*(a**3*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - 3*a**2*b*(B*(c**2 - d**2) + 2*c*d*(A - C)) - 3*a*b**2*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + b**3*(B*(c**2 - d**2) + 2*c*d*(A - C)))/(a**2 + b**2)**3 - (a**3*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 3*a**2*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) - 3*a*b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) - b**3*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/(f*(a**2 + b**2)**3) - (c + d*tan(e + f*x))**2*(A*b**2 - a*(B*b - C*a))/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) + (C*a**6*d**2 + 3*C*a**4*b**2*d**2 - a**3*b**3*(B*(c**2 - d**2) + 2*c*d*(A - C)) - 3*a**2*b**4*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - 2*C*d**2) + 3*a*b**5*(B*(c**2 - d**2) + 2*c*d*(A - C)) + b**6*(-A*(c**2 - d**2) + c*(2*B*d + C*c)))*log(a + b*tan(e + f*x))/(b**3*f*(a**2 + b**2)**3) - (-a*d + b*c)*(C*a**4*d - a**2*b**2*(B*c + d*(A - 3*C)) + 2*a*b**3*(A*c - B*d - C*c) + b**4*(A*d + B*c))/(b**3*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_64():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**4/(6*d*f) - b*(c + d*tan(e + f*x))**4*(-3*B*b*d - C*a*d + C*b*c)*tan(e + f*x)/(15*d**2*f) - d*(-a**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*a*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)))*tan(e + f*x)/f + x*(a**2*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2) - 2*a*b*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + b**2*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2)) + (c + d*tan(e + f*x))**3*(B*a**2 - B*b**2 + 2*a*b*(A - C))/(3*f) + (c + d*tan(e + f*x))**2*(a**2*(B*c + d*(A - C)) + 2*a*b*(A*c - B*d - C*c) - b**2*(B*c + d*(A - C)))/(2*f) + (-a**2*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 2*a*b*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) + b**2*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))*log(cos(e + f*x))/f + (c + d*tan(e + f*x))**4*(5*C*a**2*d**2 - 6*a*b*d*(-5*B*d + C*c) + b**2*(-3*B*c*d + C*c**2 + d**2*(15*A - 15*C)))/(60*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_65():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*b*(c + d*tan(e + f*x))**4*tan(e + f*x)/(5*d*f) + d*(A*(2*a*c*d + b*(c**2 - d**2)) + a*(B*c**2 - B*d**2 - 2*C*c*d) - b*(2*B*c*d + C*c**2 - C*d**2))*tan(e + f*x)/f + x*(a*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2) - b*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2))) + (c + d*tan(e + f*x))**3*(A*b + B*a - C*b)/(3*f) + (c + d*tan(e + f*x))**2*(A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c)/(2*f) - (A*(3*a*c**2*d - a*d**3 + b*c**3 - 3*b*c*d**2) + a*(B*c**3 - 3*B*c*d**2 - 3*C*c**2*d + C*d**3) - b*(3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2))*log(cos(e + f*x))/f - (c + d*tan(e + f*x))**4*(-5*B*b*d - 5*C*a*d + C*b*c)/(20*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_66():
    f = (c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = B*(c + d*tan(e + f*x))**3/(3*f) + C*(c + d*tan(e + f*x))**4/(4*d*f) + d*(B*(c**2 - d**2) + 2*c*d*(A - C))*tan(e + f*x)/f - x*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) + (c + d*tan(e + f*x))**2*(B*c + d*(A - C))/(2*f) - (B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2))*log(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_67():
    f = (c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = C*(c + d*tan(e + f*x))**3/(3*b*f) - x*(a*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) - b*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))/(a**2 + b**2) - (A*(a*d*(3*c**2 - d**2) - b*(c**3 - 3*c*d**2)) + a*(B*c**3 - 3*B*c*d**2 - 3*C*c**2*d + C*d**3) + b*(3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2))*log(cos(e + f*x))/(f*(a**2 + b**2)) + (c + d*tan(e + f*x))**2*(B*b*d - C*a*d + C*b*c)/(2*b**2*f) + d*(b**2*d*(B*c + d*(A - C)) + (-a*d + b*c)*(B*b*d - C*a*d + C*b*c))*tan(e + f*x)/(b**3*f) + (A*b**2 - a*(B*b - C*a))*(-a*d + b*c)**3*log(a + b*tan(e + f*x))/(b**4*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_68():
    f = (c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = -x*(a**2*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) - 2*a*b*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + b**2*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2))/(a**2 + b**2)**2 + (-a**2*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 2*a*b*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2) + b**2*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))*log(cos(e + f*x))/(f*(a**2 + b**2)**2) - (c + d*tan(e + f*x))**3*(A*b**2 - a*(B*b - C*a))/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + d*(c + d*tan(e + f*x))**2*(2*A*b**2 - 2*B*a*b + 3*C*a**2 + C*b**2)/(2*b**2*f*(a**2 + b**2)) - d**2*(-A*b**2*(-a*d + b*c) + 3*C*a**3*d - a**2*b*(2*B*d + 3*C*c) + a*b**2*(B*c + 2*C*d) - b**3*(B*d + 2*C*c))*tan(e + f*x)/(b**3*f*(a**2 + b**2)) - (-a*d + b*c)**2*(2*B*a**3*b*d - 3*C*a**4*d + a**2*b**2*(B*c - d*(A + 5*C)) - 2*a*b**3*(A*c - 2*B*d - C*c) - b**4*(3*A*d + B*c))*log(a + b*tan(e + f*x))/(b**4*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_69():
    f = (c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = -x*(a**3*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) - 3*a**2*b*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 3*a*b**2*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2) + b**3*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))/(a**2 + b**2)**3 - (a**3*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 3*a**2*b*(-A*(c**3 - 3*c*d**2) + 3*B*c**2*d - B*d**3 + C*c**3 - 3*C*c*d**2) - 3*a*b**2*(B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + b**3*(A*c**3 - 3*A*c*d**2 - 3*B*c**2*d + B*d**3 - C*c**3 + 3*C*c*d**2))*log(cos(e + f*x))/(f*(a**2 + b**2)**3) - (c + d*tan(e + f*x))**3*(A*b**2 - a*(B*b - C*a))/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) + (c + d*tan(e + f*x))**2*(B*a**3*b*d - 3*C*a**4*d + a**2*b**2*(2*B*c + d*(A - 7*C)) - a*b**3*(4*A*c - 5*B*d - 4*C*c) - b**4*(3*A*d + 2*B*c))/(2*b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - d**2*(B*a**3*b*d - 3*C*a**4*d + a**2*b**2*(B*c - 6*C*d) - a*b**3*(2*A*c - 3*B*d - 2*C*c) - b**4*(B*c + d*(2*A + C)))*tan(e + f*x)/(b**3*f*(a**2 + b**2)**2) - (-a*d + b*c)*(B*a**5*b*d**2 + B*a**3*b**3*(c**2 + 3*d**2) - 3*C*a**6*d**2 + a**4*b**2*d*(B*c - 9*C*d) + a**2*b**4*(-A*(3*c**2 - d**2) + 6*B*c*d + 3*C*c**2 - 10*C*d**2) - a*b**5*(3*B*(c**2 - 2*d**2) + 8*c*d*(A - C)) - b**6*(-A*(c**2 - 3*d**2) + c*(3*B*d + C*c)))*log(a + b*tan(e + f*x))/(b**4*f*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_70():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))
    F = C*(a + b*tan(e + f*x))**3/(3*d*f) + b*(b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-B*b*d - C*a*d + C*b*c))*tan(e + f*x)/(d**3*f) + x*(a**3*(A*c + B*d - C*c) - 3*a**2*b*(B*c - d*(A - C)) - 3*a*b**2*(A*c + B*d - C*c) + b**3*(B*c - d*(A - C)))/(c**2 + d**2) - (a**3*(B*c - d*(A - C)) + 3*a**2*b*(A*c + B*d - C*c) - 3*a*b**2*(B*c - d*(A - C)) - b**3*(A*c + B*d - C*c))*log(cos(e + f*x))/(f*(c**2 + d**2)) - (a + b*tan(e + f*x))**2*(-B*b*d - C*a*d + C*b*c)/(2*d**2*f) - (-a*d + b*c)**3*(A*d**2 - B*c*d + C*c**2)*log(c + d*tan(e + f*x))/(d**4*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_71():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))
    F = C*(a + b*tan(e + f*x))**2/(2*d*f) - b*(-B*b*d - C*a*d + C*b*c)*tan(e + f*x)/(d**2*f) + x*(a**2*(A*c + B*d - C*c) - 2*a*b*(B*c - d*(A - C)) - b**2*(A*c + B*d - C*c))/(c**2 + d**2) - (a**2*(B*c - d*(A - C)) + 2*a*b*(A*c + B*d - C*c) - b**2*(B*c - d*(A - C)))*log(cos(e + f*x))/(f*(c**2 + d**2)) + (-a*d + b*c)**2*(A*d**2 - B*c*d + C*c**2)*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_72():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))
    F = C*b*tan(e + f*x)/(d*f) + x*(a*(A*c + B*d - C*c) - b*(B*c - d*(A - C)))/(c**2 + d**2) - (-A*a*d + A*b*c + B*a*c + B*b*d + C*a*d - C*b*c)*log(cos(e + f*x))/(f*(c**2 + d**2)) - (-a*d + b*c)*(A*d**2 - B*c*d + C*c**2)*log(c + d*tan(e + f*x))/(d**2*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_73():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))
    F = x*(A*c + B*d - C*c)/(c**2 + d**2) - (B*c - d*(A - C))*log(cos(e + f*x))/(f*(c**2 + d**2)) + (A*d**2 - B*c*d + C*c**2)*log(c + d*tan(e + f*x))/(d*f*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_74():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*(c + d*tan(e + f*x)))
    F = x*(a*(A*c + B*d - C*c) + b*(B*c - d*(A - C)))/((a**2 + b**2)*(c**2 + d**2)) - (A*d**2 - B*c*d + C*c**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)) + (A*b**2 - a*(B*b - C*a))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_75():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x)))
    F = d*(A*d**2 - B*c*d + C*c**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)**2) + x*(a**2*(A*c + B*d - C*c) + 2*a*b*(B*c - d*(A - C)) - b**2*(A*c + B*d - C*c))/((a**2 + b**2)**2*(c**2 + d**2)) + (2*B*a**3*b*d - C*a**4*d - a**2*b**2*(3*A*d + B*c - C*d) + 2*a*b**3*c*(A - C) + b**4*(-A*d + B*c))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**2) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_76():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**3*(c + d*tan(e + f*x)))
    F = -d**2*(A*d**2 - B*c*d + C*c**2)*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)*(-a*d + b*c)**3) + x*(a**3*(A*c + B*d - C*c) + 3*a**2*b*(B*c - d*(A - C)) - 3*a*b**2*(A*c + B*d - C*c) - b**3*(B*c - d*(A - C)))/((a**2 + b**2)**3*(c**2 + d**2)) + (-3*B*a**5*b*d**2 + 3*B*a*b**5*c**2 + C*a**6*d**2 + 3*a**4*b**2*d*(2*A*d + B*c - C*d) - a**3*b**3*(B*(c**2 - d**2) + 8*c*d*(A - C)) - 3*a**2*b**4*(-A*(c**2 + d**2) + c*(2*B*d + C*c)) + b**6*(-A*(c**2 - d**2) + c*(-B*d + C*c)))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3*(-a*d + b*c)**3) - (2*B*a**3*b*d - C*a**4*d - a**2*b**2*(3*A*d + B*c - C*d) + 2*a*b**3*c*(A - C) + b**4*(-A*d + B*c))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)**2) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_77():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**2
    F = b**2*(a*d*(-B*c*d + 3*C*c**2 + d**2*(A + 2*C)) - b*(-2*B*c**2*d - B*d**3 + 3*C*c**3 + c*d**2*(A + 2*C)))*tan(e + f*x)/(d**3*f*(c**2 + d**2)) + b*(a + b*tan(e + f*x))**2*(-2*B*c*d + 3*C*c**2 + d**2*(2*A + C))/(2*d**2*f*(c**2 + d**2)) - x*(a**3*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - 3*a**2*b*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - 3*a*b**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) + b**3*(-B*(c**2 - d**2) + 2*c*d*(A - C)))/(c**2 + d**2)**2 + (a**3*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 3*a**2*b*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - 3*a*b**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - b**3*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/(f*(c**2 + d**2)**2) - (a + b*tan(e + f*x))**3*(A*d**2 - B*c*d + C*c**2)/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2)) + (-a*d + b*c)**2*(a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(3*A*d**4 - 2*B*c**3*d - 4*B*c*d**3 + 3*C*c**4 + c**2*d**2*(A + 5*C)))*log(c + d*tan(e + f*x))/(d**4*f*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_78():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**2
    F = b**2*(-B*c*d + 2*C*c**2 + d**2*(A + C))*tan(e + f*x)/(d**2*f*(c**2 + d**2)) - x*(a**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - 2*a*b*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - b**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2))/(c**2 + d**2)**2 + (a**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*a*b*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - b**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)))*log(cos(e + f*x))/(f*(c**2 + d**2)**2) - (a + b*tan(e + f*x))**2*(A*d**2 - B*c*d + C*c**2)/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2)) - (-a*d + b*c)*(a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(2*A*d**4 - B*c**3*d - 3*B*c*d**3 + 2*C*c**4 + 4*C*c**2*d**2))*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_79():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**2
    F = -x*(a*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - b*(-B*(c**2 - d**2) + 2*c*d*(A - C)))/(c**2 + d**2)**2 - (-A*(2*a*c*d - b*(c**2 - d**2)) + a*(B*c**2 - B*d**2 + 2*C*c*d) - b*(-2*B*c*d + C*c**2 - C*d**2))*log(cos(e + f*x))/(f*(c**2 + d**2)**2) + (a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(A*d**4 - 2*B*c*d**3 + C*c**4 - c**2*d**2*(A - 3*C)))*log(c + d*tan(e + f*x))/(d**2*f*(c**2 + d**2)**2) + (-a*d + b*c)*(A*d**2 - B*c*d + C*c**2)/(d**2*f*(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_80():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**2
    F = -x*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2)/(c**2 + d**2)**2 + (-B*(c**2 - d**2) + 2*c*d*(A - C))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2) - (A*d**2 - B*c*d + C*c**2)/(d*f*(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_81():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**2)
    F = b*(A*b**2 - a*(B*b - C*a))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c)**2) - x*(a*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) + b*(-B*(c**2 - d**2) + 2*c*d*(A - C)))/((a**2 + b**2)*(c**2 + d**2)**2) - (-a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(A*d**4 - 2*B*c**3*d + C*c**4 + c**2*d**2*(3*A - C)))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**2) + (A*d**2 - B*c*d + C*c**2)/(f*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_82():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**2)
    F = b*(3*B*a**3*b*d - 2*C*a**4*d - a**2*b**2*(4*A*d + B*c) + a*b**3*(2*A*c + B*d - 2*C*c) + b**4*(-2*A*d + B*c))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**3) + d*(-a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(4*A*c**2*d**2 + 2*A*d**4 - 3*B*c**3*d - B*c*d**3 + 2*C*c**4))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**3) - d*(A*(a**2*d**2 + b**2*(c**2 + 2*d**2)) - B*a*b*(c**2 + d**2) + a**2*(-B*c*d + 2*C*c**2 + C*d**2) + b**2*c*(-B*d + C*c))/(f*(a**2 + b**2)*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) - x*(a**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) + 2*a*b*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - b**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2))/((a**2 + b**2)**2*(c**2 + d**2)**2) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_83():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**2)
    F = -b*(6*B*a**5*b*d**2 - 3*C*a**6*d**2 - a**4*b**2*d*(4*B*c + d*(10*A - C)) + a**3*b**3*(B*(c**2 + 3*d**2) + 10*c*d*(A - C)) + 3*a**2*b**4*(-A*(c**2 + 3*d**2) + c*(2*B*d + C*c)) + a*b**5*(-B*(3*c**2 - d**2) + 2*c*d*(A - C)) - b**6*(-A*(c**2 - 3*d**2) + c*(-2*B*d + C*c)))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**3*(-a*d + b*c)**4) - d**2*(-a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(3*A*d**4 - 4*B*c**3*d - 2*B*c*d**3 + 3*C*c**4 + c**2*d**2*(5*A + C)))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**2*(-a*d + b*c)**4) - d*(3*B*a**3*b*d*(c**2 + d**2) - a**4*d*(-B*c*d + 3*C*c**2 + d**2*(A + 2*C)) - a**2*b**2*(4*A*c**2*d + 6*A*d**3 + B*c**3 - B*c*d**2 + 2*C*c**2*d) + a*b**3*(c**2 + d**2)*(2*A*c + B*d - 2*C*c) - b**4*(-B*(c**3 + 2*c*d**2) + d*(2*A*c**2 + 3*A*d**2 + C*c**2)))/(f*(a**2 + b**2)**2*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**3) - x*(a**3*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) + 3*a**2*b*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - 3*a*b**2*(-A*(c**2 - d**2) - 2*B*c*d + C*c**2 - C*d**2) - b**3*(-B*(c**2 - d**2) + 2*c*d*(A - C)))/((a**2 + b**2)**3*(c**2 + d**2)**2) - (5*B*a**3*b*d - 3*C*a**4*d - a**2*b**2*(2*B*c + d*(7*A - C)) + a*b**3*(4*A*c + B*d - 4*C*c) + b**4*(-3*A*d + 2*B*c))/(2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(c + d*tan(e + f*x))*(-a*d + b*c)**2) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))**2*(2*a**2 + 2*b**2)*(c + d*tan(e + f*x))*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_84():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**3
    F = b**2*(a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(-B*c**3*d - 3*B*c*d**3 + 3*C*c**4 + 6*C*c**2*d**2 + d**4*(2*A + C)))*tan(e + f*x)/(d**3*f*(c**2 + d**2)**2) - x*(a**3*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2) - 3*a**2*b*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 3*a*b**2*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2) + b**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))/(c**2 + d**2)**3 - (-a**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 3*a**2*b*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2) + 3*a*b**2*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) - b**3*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2))*log(cos(e + f*x))/(f*(c**2 + d**2)**3) - (a + b*tan(e + f*x))**3*(A*d**2 - B*c*d + C*c**2)/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2)) - (a + b*tan(e + f*x))**2*(2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(3*A*d**4 - B*c**3*d - 5*B*c*d**3 + 3*C*c**4 - c**2*d**2*(A - 7*C)))/(2*d**2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) - (-a*d + b*c)*(a**2*d**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + a*b*d**2*(-B*(c**4 + 6*c**2*d**2 - 3*d**4) + 8*c*d**3*(A - C)) + b**2*(3*A*d**6 - B*c**5*d - 3*B*c**3*d**3 - 6*B*c*d**5 + 3*C*c**6 + 9*C*c**4*d**2 - c**2*d**4*(A - 10*C)))*log(c + d*tan(e + f*x))/(d**4*f*(c**2 + d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_85():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**3
    F = -x*(a**2*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2) - 2*a*b*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + b**2*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2))/(c**2 + d**2)**3 - (-a**2*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 2*a*b*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2) + b**2*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))*log(cos(e + f*x))/(f*(c**2 + d**2)**3) - (a + b*tan(e + f*x))**2*(A*d**2 - B*c*d + C*c**2)/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2)) - (-a**2*d**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + 2*a*b*d**3*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2) - b**2*(A*d**6 + B*c**3*d**3 - 3*B*c*d**5 + C*c**6 + 3*C*c**4*d**2 - 3*c**2*d**4*(A - 2*C)))*log(c + d*tan(e + f*x))/(d**3*f*(c**2 + d**2)**3) + (-a*d + b*c)*(a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(A*d**4 - 2*B*c*d**3 + C*c**4 - c**2*d**2*(A - 3*C)))/(d**3*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_86():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**3
    F = -x*(a*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2) - b*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))/(c**2 + d**2)**3 + (A*(a*d*(3*c**2 - d**2) - b*(c**3 - 3*c*d**2)) - a*(B*c**3 - 3*B*c*d**2 + 3*C*c**2*d - C*d**3) + b*(-3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3) - (a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(A*d**4 - 2*B*c*d**3 + C*c**4 - c**2*d**2*(A - 3*C)))/(d**2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) + (-a*d + b*c)*(A*d**2 - B*c*d + C*c**2)/(2*d**2*f*(c + d*tan(e + f*x))**2*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_87():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**3
    F = -x*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2)/(c**2 + d**2)**3 + (-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3) - (-B*(c**2 - d**2) + 2*c*d*(A - C))/(f*(c + d*tan(e + f*x))*(c**2 + d**2)**2) - (A*d**2 - B*c*d + C*c**2)/(2*d*f*(c + d*tan(e + f*x))**2*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_88():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**3)
    F = b**2*(A*b**2 - a*(B*b - C*a))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)*(-a*d + b*c)**3) - x*(a*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2) + b*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)))/((a**2 + b**2)*(c**2 + d**2)**3) - (a**2*d**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) - a*b*d**2*(-B*(3*c**4 - 6*c**2*d**2 - d**4) + 8*c**3*d*(A - C)) + b**2*(3*A*c**2*d**4 + A*d**6 - 3*B*c**5*d + B*c**3*d**3 + C*c**6 + 3*c**4*d**2*(2*A - C)))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3*(-a*d + b*c)**3) + (-a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(A*d**4 - 2*B*c**3*d + C*c**4 + c**2*d**2*(3*A - C)))/(f*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + (A*d**2 - B*c*d + C*c**2)/(f*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-2*a*d + 2*b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_89():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**3)
    F = b**2*(4*B*a**3*b*d - 3*C*a**4*d - a**2*b**2*(B*c + d*(5*A + C)) + 2*a*b**3*(A*c + B*d - C*c) + b**4*(-3*A*d + B*c))*log(a*cos(e + f*x) + b*sin(e + f*x))/(f*(a**2 + b**2)**2*(-a*d + b*c)**4) + d*(a**2*d**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) - 2*a*b*d**2*(-B*(2*c**4 - 3*c**2*d**2 - d**4) + c*d*(A - C)*(5*c**2 + d**2)) + b**2*(9*A*c**2*d**4 + 3*A*d**6 - 6*B*c**5*d - 3*B*c**3*d**3 - B*c*d**5 + 3*C*c**6 + c**4*d**2*(10*A - C)))*log(c*cos(e + f*x) + d*sin(e + f*x))/(f*(c**2 + d**2)**3*(-a*d + b*c)**4) - d*(A*(a**2*d**2 + b**2*(2*c**2 + 3*d**2)) - 2*B*a*b*(c**2 + d**2) + a**2*(-B*c*d + 3*C*c**2 + 2*C*d**2) + b**2*c*(-B*d + C*c))/(f*(2*a**2 + 2*b**2)*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-a*d + b*c)**2) - d*(-A*(2*a**3*c*d**3 - 2*a**2*b*d**2*(2*c**2 + d**2) + 2*a*b**2*c*d**3 - b**3*(c**4 + 6*c**2*d**2 + 3*d**4)) + a**3*d**2*(B*(c**2 - d**2) + 2*C*c*d) + a**2*b*(-3*B*c**3*d - B*c*d**3 + 3*C*c**4 + 2*C*c**2*d**2 + C*d**4) + a*b**2*(-B*(c**4 + c**2*d**2 + 2*d**4) + 2*C*c*d**3) + b**3*c*(-3*B*c**2*d - B*d**3 + 2*C*c**3))/(f*(a**2 + b**2)*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) - x*(a**2*(-A*(c**3 - 3*c*d**2) - 3*B*c**2*d + B*d**3 + C*c**3 - 3*C*c*d**2) + 2*a*b*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) + b**2*(A*c**3 - 3*A*c*d**2 + 3*B*c**2*d - B*d**3 - C*c**3 + 3*C*c*d**2))/((a**2 + b**2)**2*(c**2 + d**2)**3) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_90():
    f = (a + b*tan(e + f*x))**3*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*(a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(9*d*f) + 2*b*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(21*b*d**2*(A*b + B*a - C*b) + (-4*a*d + 4*b*c)*(-3*B*b*d - 2*C*a*d + 2*C*b*c))*tan(e + f*x)/(105*d**3*f) - (a - I*b)**3*sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (a + I*b)**3*sqrt(c + I*d)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + sqrt(c + d*tan(e + f*x))*(2*B*a**3 - 6*B*a*b**2 + 6*a**2*b*(A - C) - 2*b**3*(A - C))/f - (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-6*B*b*d - 4*C*a*d + 4*C*b*c)/(21*d**2*f) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(80*C*a**3*d**3 - 12*a**2*b*d**2*(-45*B*d + 16*C*c) + 18*a*b**2*d*(-14*B*c*d + 8*C*c**2 + d**2*(35*A - 35*C)) - 2*b**3*(-24*B*c**2*d + 105*B*d**3 + 16*C*c**3 + 42*c*d**2*(A - C)))/(315*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_91():
    f = (a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*(a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(7*d*f) - 2*b*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-7*B*b*d - 4*C*a*d + 4*C*b*c)*tan(e + f*x)/(35*d**2*f) - (B - I*(A - C))*(a + I*b)**2*sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f - (B + I*(A - C))*(a - I*b)**2*sqrt(c - I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + sqrt(c + d*tan(e + f*x))*(2*B*a**2 - 2*B*b**2 + 4*a*b*(A - C))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(40*C*a**2*d**2 - 28*a*b*d*(-5*B*d + 2*C*c) + 2*b**2*(-14*B*c*d + 8*C*c**2 + d**2*(35*A - 35*C)))/(105*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_92():
    f = (a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*b*(c + d*tan(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)/(5*d*f) - sqrt(c - I*d)*(I*a + b)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + sqrt(c + I*d)*(I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + sqrt(c + d*tan(e + f*x))*(2*A*b + 2*B*a - 2*C*b)/f - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(-10*B*b*d - 10*C*a*d + 4*C*b*c)/(15*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_93():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*B*sqrt(c + d*tan(e + f*x))/f + 2*C*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f) - (B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_94():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = 2*C*sqrt(c + d*tan(e + f*x))/(b*f) + sqrt(c + I*d)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)) - (2*A*b**2 - 2*a*(B*b - C*a))*sqrt(-a*d + b*c)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_95():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = -(B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2) - sqrt(c + d*tan(e + f*x))*(A*b**2 - a*(B*b - C*a))/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) - (B*a**3*b*d + C*a**4*d - a**2*b**2*(3*A*d + 2*B*c - 5*C*d) + a*b**3*(4*A*c - 3*B*d - 4*C*c) + b**4*(A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*f*(a**2 + b**2)**2*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_96():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = -sqrt(c - I*d)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + sqrt(c + I*d)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3) - sqrt(c + d*tan(e + f*x))*(3*B*a**3*b*d + C*a**4*d - a**2*b**2*(7*A*d + 4*B*c - 9*C*d) + a*b**3*(8*A*c - 5*B*d - 8*C*c) + b**4*(A*d + 4*B*c))/(4*b*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(A*b**2 - a*(B*b - C*a))/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) + (3*B*a**5*b*d**2 + C*a**6*d**2 - 3*a**4*b**2*d*(5*A*d + 4*B*c - 6*C*d) + 2*a**3*b**3*(B*(4*c**2 - 13*d**2) + 20*c*d*(A - C)) - 3*a**2*b**4*(8*A*c**2 - 6*A*d**2 - 16*B*c*d - 8*C*c**2 + 5*C*d**2) - 3*a*b**5*(B*(8*c**2 - d**2) + 8*c*d*(A - C)) - b**6*(-A*(8*c**2 + d**2) + 4*c*(B*d + 2*C*c)))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*f*(a**2 + b**2)**3*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_97():
    f = (a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*(a + b*tan(e + f*x))**3*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(11*d*f) + 2*b*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(99*b*d**2*(A*b + B*a - C*b) + (-4*a*d + 4*b*c)*(-11*B*b*d - 6*C*a*d + 6*C*b*c))*tan(e + f*x)/(693*d**3*f) + (a + I*b)**3*(c + I*d)**(sympy.S(3)/2)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c - I*d)**(sympy.S(3)/2)*(I*a + b)**3*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*B*a**3 - 6*B*a*b**2 + 6*a**2*b*(A - C) - 2*b**3*(A - C))/(3*f) + sqrt(c + d*tan(e + f*x))*(2*a**3*(B*c + d*(A - C)) + 6*a**2*b*(A*c - B*d - C*c) - 6*a*b**2*(B*c + d*(A - C)) - 2*b**3*(A*c - B*d - C*c))/f - (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(-22*B*b*d - 12*C*a*d + 12*C*b*c)/(99*d**2*f) + (c + d*tan(e + f*x))**(sympy.S(5)/2)*(336*C*a**3*d**3 - 4*a**2*b*d**2*(-847*B*d + 192*C*c) + 66*a*b**2*d*(-18*B*c*d + 8*C*c**2 + d**2*(63*A - 63*C)) - 2*b**3*(-88*B*c**2*d + 693*B*d**3 + 48*C*c**3 + 198*c*d**2*(A - C)))/(3465*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_98():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*(a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(9*d*f) - 2*b*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(-9*B*b*d - 4*C*a*d + 4*C*b*c)*tan(e + f*x)/(63*d**2*f) - (B + I*(A - C))*(a - I*b)**2*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (a + I*b)**2*(c + I*d)**(sympy.S(3)/2)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*B*a**2 - 2*B*b**2 + 4*a*b*(A - C))/(3*f) + sqrt(c + d*tan(e + f*x))*(2*a**2*(B*c + d*(A - C)) + 4*a*b*(A*c - B*d - C*c) - 2*b**2*(B*c + d*(A - C)))/f + (c + d*tan(e + f*x))**(sympy.S(5)/2)*(56*C*a**2*d**2 - 36*a*b*d*(-7*B*d + 2*C*c) + 2*b**2*(-18*B*c*d + 8*C*c**2 + d**2*(63*A - 63*C)))/(315*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_99():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*b*(c + d*tan(e + f*x))**(sympy.S(5)/2)*tan(e + f*x)/(7*d*f) - (c - I*d)**(sympy.S(3)/2)*(I*a + b)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + I*d)**(sympy.S(3)/2)*(I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*A*b + 2*B*a - 2*C*b)/(3*f) + sqrt(c + d*tan(e + f*x))*(2*A*a*d + 2*A*b*c + 2*B*a*c - 2*B*b*d - 2*C*a*d - 2*C*b*c)/f - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(-14*B*b*d - 14*C*a*d + 4*C*b*c)/(35*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_100():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*B*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*f) + 2*C*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*d*f) - (B - I*(A - C))*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + sqrt(c + d*tan(e + f*x))*(2*B*c + 2*d*(A - C))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_101():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = 2*C*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*b*f) - (c + I*d)**(sympy.S(3)/2)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)) + sqrt(c + d*tan(e + f*x))*(2*B*b*d - 2*C*a*d + 2*C*b*c)/(b**2*f) - (2*A*b**2 - 2*a*(B*b - C*a))*(-a*d + b*c)**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(5)/2)*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_102():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A*b**2 - a*(B*b - C*a))/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + d*sqrt(c + d*tan(e + f*x))*(A*b**2 - B*a*b + 3*C*a**2 + 2*C*b**2)/(b**2*f*(a**2 + b**2)) + sqrt(-a*d + b*c)*(B*a**3*b*d - 3*C*a**4*d + a**2*b**2*(2*B*c + d*(A - 7*C)) - a*b**3*(4*A*c - 5*B*d - 4*C*c) - b**4*(3*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(5)/2)*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_103():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = -(c - I*d)**(sympy.S(3)/2)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + (c + I*d)**(sympy.S(3)/2)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A*b**2 - a*(B*b - C*a))/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) - sqrt(c + d*tan(e + f*x))*(B*a**3*b*d + 3*C*a**4*d - a**2*b**2*(5*A*d + 4*B*c - 11*C*d) + a*b**3*(8*A*c - 7*B*d - 8*C*c) + b**4*(3*A*d + 4*B*c))/(4*b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - (B*a**5*b*d**2 + 3*C*a**6*d**2 + a**4*b**2*d*(4*B*c + d*(3*A + 6*C)) - 2*a**3*b**3*(B*(4*c**2 - 9*d**2) + 12*c*d*(A - C)) + a**2*b**4*(24*A*c**2 - 26*A*d**2 - 48*B*c*d - 24*C*c**2 + 35*C*d**2) + a*b**5*(3*B*(8*c**2 - 5*d**2) + 40*c*d*(A - C)) - b**6*(8*A*c**2 - 3*A*d**2 - 12*B*c*d - 8*C*c**2))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*b**(sympy.S(5)/2)*f*(a**2 + b**2)**3*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_104():
    f = (a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*(a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(11*d*f) - 2*b*(c + d*tan(e + f*x))**(sympy.S(7)/2)*(-11*B*b*d - 4*C*a*d + 4*C*b*c)*tan(e + f*x)/(99*d**2*f) - (a - I*b)**2*(c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (a + I*b)**2*(c + I*d)**(sympy.S(5)/2)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*B*a**2 - 2*B*b**2 + 4*a*b*(A - C))/(5*f) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*a**2*(B*c + d*(A - C)) + 4*a*b*(A*c - B*d - C*c) - 2*b**2*(B*c + d*(A - C)))/(3*f) - sqrt(c + d*tan(e + f*x))*(-2*a**2*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 4*a*b*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - C*d**2) + 2*b**2*(B*(c**2 - d**2) + 2*c*d*(A - C)))/f + (c + d*tan(e + f*x))**(sympy.S(7)/2)*(72*C*a**2*d**2 - 44*a*b*d*(-9*B*d + 2*C*c) + 2*b**2*(-22*B*c*d + 8*C*c**2 + d**2*(99*A - 99*C)))/(693*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_105():
    f = (a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*C*b*(c + d*tan(e + f*x))**(sympy.S(7)/2)*tan(e + f*x)/(9*d*f) - (c - I*d)**(sympy.S(5)/2)*(I*a + b)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + I*d)**(sympy.S(5)/2)*(I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f + (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*A*b + 2*B*a - 2*C*b)/(5*f) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*A*a*d + 2*A*b*c + 2*B*a*c - 2*B*b*d - 2*C*a*d - 2*C*b*c)/(3*f) + sqrt(c + d*tan(e + f*x))*(2*A*(2*a*c*d + b*(c**2 - d**2)) + 2*a*(B*c**2 - B*d**2 - 2*C*c*d) - 2*b*(2*B*c*d + C*c**2 - C*d**2))/f - (c + d*tan(e + f*x))**(sympy.S(7)/2)*(-18*B*b*d - 18*C*a*d + 4*C*b*c)/(63*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_106():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = 2*B*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*f) + 2*C*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f) - (B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/f - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/f + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*B*c + 2*d*(A - C))/(3*f) + sqrt(c + d*tan(e + f*x))*(2*B*(c**2 - d**2) + 4*c*d*(A - C))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_107():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))
    F = 2*C*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(5*b*f) + (c + I*d)**(sympy.S(5)/2)*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*B*b*d - 2*C*a*d + 2*C*b*c)/(3*b**2*f) + sqrt(c + d*tan(e + f*x))*(2*b**2*d*(B*c + d*(A - C)) + 2*(-a*d + b*c)*(B*b*d - C*a*d + C*b*c))/(b**3*f) - (2*A*b**2 - 2*a*(B*b - C*a))*(-a*d + b*c)**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(7)/2)*f*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_108():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**2
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A*b**2 - a*(B*b - C*a))/(b*f*(a + b*tan(e + f*x))*(a**2 + b**2)) + d*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(3*A*b**2 - 3*B*a*b + 5*C*a**2 + 2*C*b**2)/(3*b**2*f*(a**2 + b**2)) - d*sqrt(c + d*tan(e + f*x))*(-A*b**2*(-a*d + b*c) + 5*C*a**3*d - a**2*b*(3*B*d + 5*C*c) + a*b**2*(B*c + 4*C*d) - 2*b**3*(B*d + 2*C*c))/(b**3*f*(a**2 + b**2)) + (-a*d + b*c)**(sympy.S(3)/2)*(3*B*a**3*b*d - 5*C*a**4*d + a**2*b**2*(2*B*c - d*(A + 9*C)) - a*b**3*(4*A*c - 7*B*d - 4*C*c) - b**4*(5*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(b**(sympy.S(7)/2)*f*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_109():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**3
    F = -(c - I*d)**(sympy.S(5)/2)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(I*a + b)**3) + (c + I*d)**(sympy.S(5)/2)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(I*a - b)**3) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A*b**2 - a*(B*b - C*a))/(2*b*f*(a + b*tan(e + f*x))**2*(a**2 + b**2)) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(B*a**3*b*d - 5*C*a**4*d + a**2*b**2*(3*A*d + 4*B*c - 13*C*d) - a*b**3*(8*A*c - 9*B*d - 8*C*c) - b**4*(5*A*d + 4*B*c))/(4*b**2*f*(a + b*tan(e + f*x))*(a**2 + b**2)**2) - d*sqrt(c + d*tan(e + f*x))*(3*B*a**3*b*d - 15*C*a**4*d + a**2*b**2*(4*B*c + d*(A - 31*C)) - a*b**3*(8*A*c - 11*B*d - 8*C*c) - b**4*(7*A*d + 4*B*c + 8*C*d))/(4*b**3*f*(a**2 + b**2)**2) + sqrt(-a*d + b*c)*(3*B*a**5*b*d**2 - 15*C*a**6*d**2 + a**4*b**2*d*(4*B*c + d*(A - 46*C)) + 2*a**3*b**3*(B*(4*c**2 + 3*d**2) + 4*c*d*(A - C)) - 3*a**2*b**4*(8*A*c**2 - 6*A*d**2 - 16*B*c*d - 8*C*c**2 + 21*C*d**2) - a*b**5*(B*(24*c**2 - 35*d**2) + 56*c*d*(A - C)) - b**6*(-A*(8*c**2 - 15*d**2) + 4*c*(5*B*d + 2*C*c)))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(4*b**(sympy.S(7)/2)*f*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_110():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = 2*C*(a + b*tan(e + f*x))**3*sqrt(c + d*tan(e + f*x))/(7*d*f) + 2*b*sqrt(c + d*tan(e + f*x))*(35*b*d**2*(A*b + B*a - C*b) + (-4*a*d + 4*b*c)*(-7*B*b*d - 6*C*a*d + 6*C*b*c))*tan(e + f*x)/(105*d**3*f) - (I*a - b)**3*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) + (I*a + b)**3*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) - (a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))*(-14*B*b*d - 12*C*a*d + 12*C*b*c)/(35*d**2*f) + sqrt(c + d*tan(e + f*x))*(144*C*a**3*d**3 - 12*a**2*b*d**2*(-49*B*d + 32*C*c) + 42*a*b**2*d*(-10*B*c*d + 8*C*c**2 + d**2*(15*A - 15*C)) - 2*b**3*(-56*B*c**2*d + 105*B*d**3 + 48*C*c**3 + 70*c*d**2*(A - C)))/(105*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_111():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = 2*C*(a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))/(5*d*f) - 2*b*sqrt(c + d*tan(e + f*x))*(-5*B*b*d - 4*C*a*d + 4*C*b*c)*tan(e + f*x)/(15*d**2*f) - (B + I*(A - C))*(a - I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) + (a + I*b)**2*(I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) + sqrt(c + d*tan(e + f*x))*(24*C*a**2*d**2 - 20*a*b*d*(-3*B*d + 2*C*c) + 2*b**2*(-10*B*c*d + 8*C*c**2 + d**2*(15*A - 15*C)))/(15*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_112():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = 2*C*b*sqrt(c + d*tan(e + f*x))*tan(e + f*x)/(3*d*f) + (I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) - (I*a + b)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d)) - sqrt(c + d*tan(e + f*x))*(-6*B*b*d - 6*C*a*d + 4*C*b*c)/(3*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_113():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = 2*C*sqrt(c + d*tan(e + f*x))/(d*f) - (B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_114():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x)))
    F = -(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*sqrt(c + I*d)*(I*a - b)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)*sqrt(c - I*d)) - (2*A*b**2 - 2*a*(B*b - C*a))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(sqrt(b)*f*(a**2 + b**2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_115():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x)))
    F = -(B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*sqrt(c + I*d)) - sqrt(c + d*tan(e + f*x))*(A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*sqrt(c - I*d)) - (3*B*a**3*b*d - C*a**4*d - a**2*b**2*(5*A*d + 2*B*c - 3*C*d) + a*b**3*(4*A*c - B*d - 4*C*c) + b**4*(-A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(sqrt(b)*f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_116():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*sqrt(c + d*tan(e + f*x))*(-5*d**2*(B*(a*c + b*d) + (A - C)*(-a*d + b*c)) + (-4*a*d + 4*b*c)*(-5*B*c*d + 6*C*c**2 + d**2*(5*A + C)))*tan(e + f*x)/(15*d**3*f*(c**2 + d**2)) + 2*b*(a + b*tan(e + f*x))**2*sqrt(c + d*tan(e + f*x))*(-5*B*c*d + 6*C*c**2 + d**2*(5*A + C))/(5*d**2*f*(c**2 + d**2)) + 2*b*sqrt(c + d*tan(e + f*x))*(6*a**2*d**2*(-5*B*c*d + 12*C*c**2 + d**2*(5*A + 7*C)) - 15*a*b*d*(-6*B*c**2*d - 3*B*d**3 + 8*C*c**3 + c*d**2*(3*A + 5*C)) + b**2*(-40*B*c**3*d - 25*B*c*d**3 + 48*C*c**4 + 6*c**2*d**2*(5*A + 3*C) + d**4*(15*A - 15*C)))/(15*d**4*f*(c**2 + d**2)) - (a - I*b)**3*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) - (I*a - b)**3*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - (a + b*tan(e + f*x))**3*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_117():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*b**2*sqrt(c + d*tan(e + f*x))*(-3*B*c*d + 4*C*c**2 + d**2*(3*A + C))*tan(e + f*x)/(3*d**2*f*(c**2 + d**2)) + 2*b*sqrt(c + d*tan(e + f*x))*(6*a*d*(-B*c*d + 2*C*c**2 + d**2*(A + C)) - b*(-6*B*c**2*d - 3*B*d**3 + 8*C*c**3 + c*d**2*(3*A + 5*C)))/(3*d**3*f*(c**2 + d**2)) - (B - I*(A - C))*(a + I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - (a - I*b)**2*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) - (a + b*tan(e + f*x))**2*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_118():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*C*b*sqrt(c + d*tan(e + f*x))/(d**2*f) + (I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - (I*a + b)*(A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) + (-2*a*d + 2*b*c)*(A*d**2 - B*c*d + C*c**2)/(d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_119():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -(B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(3)/2)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)) - (2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_120():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -2*sqrt(b)*(A*b**2 - a*(B*b - C*a))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)*(-a*d + b*c)**(sympy.S(3)/2)) + (2*A*d**2 - 2*B*c*d + 2*C*c**2)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) + (A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(3)/2)*(I*a + b)) + (I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)*(c + I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_121():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -sqrt(b)*(5*B*a**3*b*d - 3*C*a**4*d - a**2*b**2*(2*B*c + d*(7*A - C)) + a*b**3*(4*A*c + B*d - 4*C*c) + b**4*(-3*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(5)/2)) - d*(A*(2*a**2*d**2 + b**2*(c**2 + 3*d**2)) - B*a*b*(c**2 + d**2) + a**2*(-2*B*c*d + 3*C*c**2 + C*d**2) + 2*b**2*c*(-B*d + C*c))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) - (B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*(c + I*d)**(sympy.S(3)/2)) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_122():
    f = (a + b*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*b**2*sqrt(c + d*tan(e + f*x))*(3*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(-4*B*c**3*d - 10*B*c*d**3 + 8*C*c**4 + c**2*d**2*(A + 15*C) + d**4*(7*A + C)))*tan(e + f*x)/(3*d**3*f*(c**2 + d**2)**2) + 2*b*sqrt(c + d*tan(e + f*x))*(6*a**2*d**3*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 3*a*b*d*(-2*B*c**3*d - 8*B*c*d**3 + 8*C*c**4 - c**2*d**2*(A - 17*C) + d**4*(5*A + 3*C)) - b**2*(-8*B*c**4*d - 17*B*c**2*d**3 - 3*B*d**5 + 16*C*c**5 + 2*c**3*d**2*(A + 15*C) + 8*c*d**4*(A + C)))/(3*d**4*f*(c**2 + d**2)**2) - (a - I*b)**3*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) - (I*a - b)**3*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - (a + b*tan(e + f*x))**3*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - (a + b*tan(e + f*x))**2*(2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(2*A*d**4 - B*c**3*d - 3*B*c*d**3 + 2*C*c**4 + 4*C*c**2*d**2))/(d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_123():
    f = (a + b*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*b**2*sqrt(c + d*tan(e + f*x))*(-B*c*d + 4*C*c**2 + d**2*(A + 3*C))/(3*d**3*f*(c**2 + d**2)) - (B - I*(A - C))*(a + I*b)**2*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - (a - I*b)**2*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) - (a + b*tan(e + f*x))**2*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) + (-2*a*d + 2*b*c)*(3*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(4*A*d**4 - B*c**3*d - 7*B*c*d**3 + 4*C*c**4 - 2*c**2*d**2*(A - 5*C)))/(3*d**3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_124():
    f = (a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -(a - I*b)*(I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) + (I*a - b)*(A + I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - (2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(A*d**4 - 2*B*c*d**3 + C*c**4 - c**2*d**2*(A - 3*C)))/(d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) + (-2*a*d + 2*b*c)*(A*d**2 - B*c*d + C*c**2)/(3*d**2*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_125():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -(B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(c + I*d)**(sympy.S(5)/2)) - (-2*B*(c**2 - d**2) + 4*c*d*(A - C))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)) - (2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_126():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -2*b**(sympy.S(3)/2)*(A*b**2 - a*(B*b - C*a))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)*(-a*d + b*c)**(sympy.S(5)/2)) + (-2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(A*d**4 - 2*B*c**3*d + C*c**4 + c**2*d**2*(3*A - C)))/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + (2*A*d**2 - 2*B*c*d + 2*C*c**2)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-3*a*d + 3*b*c)) + (A - I*B - C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(c - I*d)**(sympy.S(5)/2)*(I*a + b)) + (I*A - B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)*(c + I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_127():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**2*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -b**(sympy.S(3)/2)*(7*B*a**3*b*d - 5*C*a**4*d - a**2*b**2*(2*B*c + d*(9*A + C)) + a*b**3*(4*A*c + 3*B*d - 4*C*c) + b**4*(-5*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*tan(e + f*x))/sqrt(-a*d + b*c))/(f*(a**2 + b**2)**2*(-a*d + b*c)**(sympy.S(7)/2)) - d*(A*(2*a**2*d**2 + b**2*(3*c**2 + 5*d**2)) - 3*B*a*b*(c**2 + d**2) + a**2*(-2*B*c*d + 5*C*c**2 + 3*C*d**2) + 2*b**2*c*(-B*d + C*c))/(f*(3*a**2 + 3*b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-a*d + b*c)**2) - d*(-A*(4*a**3*c*d**3 - 4*a**2*b*d**2*(2*c**2 + d**2) + 4*a*b**2*c*d**3 - b**3*(c**4 + 10*c**2*d**2 + 5*d**4)) + 2*a**3*d**2*(B*c**2 - B*d**2 + 2*C*c*d) + a**2*b*(-6*B*c**3*d - 2*B*c*d**3 + 5*C*c**4 + 2*C*c**2*d**2 + C*d**4) - a*b**2*(B*c**4 + 3*B*d**4 - 4*C*c*d**3) + 2*b**3*c*(-3*B*c**2*d - B*d**3 + 2*C*c**3))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) - (B - I*(A - C))*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c + I*d))/(f*(a + I*b)**2*(c + I*d)**(sympy.S(5)/2)) - (A*b**2 - a*(B*b - C*a))/(f*(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(c + d*tan(e + f*x))/sqrt(c - I*d))/(f*(a - I*b)**2*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_128():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(4*d*f) - (B - I*(A - C))*(a + I*b)**(sympy.S(5)/2)*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - (a - I*b)**(sympy.S(5)/2)*sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f - (a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-8*B*b*d - 5*C*a*d + 5*C*b*c)/(24*d**2*f) + sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(16*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-8*B*b*d - 5*C*a*d + 5*C*b*c))/(32*d**3*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(64*b*d**3*(B*a**2 - B*b**2 + 2*a*b*(A - C)) - (-a*d + b*c)*(16*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-8*B*b*d - 5*C*a*d + 5*C*b*c)))/(64*b*d**3*f) - (5*C*a**4*d**4 - 20*a**3*b*d**3*(2*B*d + C*c) + 30*a**2*b**2*d**2*(-4*B*c*d + C*c**2 - d**2*(8*A - 8*C)) - 20*a*b**3*d*(-2*B*c**2*d - 16*B*d**3 + C*c**3 + 8*c*d**2*(A - C)) + b**4*(-8*B*c**3*d + 64*B*c*d**3 + 5*C*c**4 + 16*c**2*d**2*(A - C) + d**4*(128*A - 128*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(64*b**(sympy.S(3)/2)*d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_129():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f) - (a - I*b)**(sympy.S(3)/2)*sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + (a + I*b)**(sympy.S(3)/2)*sqrt(c + I*d)*(I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-2*B*b*d - C*a*d + C*b*c)/(4*d**2*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(8*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-2*B*b*d - C*a*d + C*b*c))/(8*b*d**2*f) - (C*a**3*d**3 - 3*a**2*b*d**2*(2*B*d + C*c) + 3*a*b**2*d*(-4*B*c*d + C*c**2 - d**2*(8*A - 8*C)) - b**3*(-2*B*c**2*d - 16*B*d**3 + C*c**3 + 8*c*d**2*(A - C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(8*b**(sympy.S(3)/2)*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_130():
    f = sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(2*d*f) - (B - I*(A - C))*sqrt(a + I*b)*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a - I*b)*sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-4*B*b*d - C*a*d + C*b*c)/(4*b*d*f) - (C*a**2*d**2 - 2*a*b*d*(2*B*d + C*c) + b**2*(-4*B*c*d + C*c**2 - d**2*(8*A - 8*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*b**(sympy.S(3)/2)*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_131():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(a + b*tan(e + f*x))
    F = C*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/(b*f) - (B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)) + (2*B*b*d - C*a*d + C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(3)/2)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_132():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*C*sqrt(d)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(3)/2)*f) - (B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)) - sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*a*(B*b - C*a))/(b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_133():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = -(B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)) - sqrt(c + d*tan(e + f*x))*(4*B*a**3*b*d + 2*C*a**4*d - 2*a**2*b**2*(5*A*d + 3*B*c - 7*C*d) + 4*a*b**3*(3*A*c - 2*B*d - 3*C*c) + 2*b**4*(A*d + 3*B*c))/(3*b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_134():
    f = sqrt(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(7)/2)
    F = -(B - I*(A - C))*sqrt(c + I*d)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(7)/2)) - sqrt(c - I*d)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(7)/2)) + sqrt(c + d*tan(e + f*x))*(16*B*a**5*b*d**2 + 4*C*a**6*d**2 - 2*a**4*b**2*d*(33*A*d + 25*B*c - 39*C*d) + 2*a**3*b**3*(B*(15*c**2 - 49*d**2) + 80*c*d*(A - C)) - 2*a**2*b**4*(45*A*c**2 - 29*A*d**2 - 90*B*c*d - 45*C*c**2 + 23*C*d**2) - 2*a*b**5*(B*(45*c**2 - 3*d**2) + 40*c*d*(A - C)) - 2*b**6*(-A*(15*c**2 + 2*d**2) + 5*c*(B*d + 3*C*c)))/(15*b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**3*(-a*d + b*c)**2) - sqrt(c + d*tan(e + f*x))*(8*B*a**3*b*d + 2*C*a**4*d - 2*a**2*b**2*(9*A*d + 5*B*c - 11*C*d) + 4*a*b**3*(5*A*c - 3*B*d - 5*C*c) + 2*b**4*(A*d + 5*B*c))/(15*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**2*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_135():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(4*d*f) - (B - I*(A - C))*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - (B + I*(A - C))*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(-8*B*b*d - 3*C*a*d + 3*C*b*c)/(24*d**2*f) + sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(48*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-8*B*b*d - 3*C*a*d + 3*C*b*c))/(96*b*d**2*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(64*b*d**3*(B*a**2 - B*b**2 + 2*a*b*(A - C)) + (-a*d + b*c)*(48*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-8*B*b*d - 3*C*a*d + 3*C*b*c)))/(64*b**2*d**2*f) + (3*C*a**4*d**4 - 4*a**3*b*d**3*(2*B*d + 3*C*c) + 6*a**2*b**2*d**2*(12*B*c*d + 3*C*c**2 + d**2*(8*A - 8*C)) - 12*a*b**3*d*(-6*B*c**2*d + 16*B*d**3 + C*c**3 - 24*c*d**2*(A - C)) + b**4*(-8*B*c**3*d - 192*B*c*d**3 + 3*C*c**4 + 48*c**2*d**2*(A - C) - d**4*(128*A - 128*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(64*b**(sympy.S(5)/2)*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_136():
    f = sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(3*d*f) - (B - I*(A - C))*sqrt(a + I*b)*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a - I*b)*(c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-6*B*b*d - C*a*d + C*b*c)/(12*b*d*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(8*b*d**2*(A*b + B*a - C*b) - (-a*d + b*c)*(-6*B*b*d - C*a*d + C*b*c))/(8*b**2*d*f) + (C*a**3*d**3 - a**2*b*d**2*(2*B*d + 3*C*c) + a*b**2*d*(12*B*c*d + 3*C*c**2 + d**2*(8*A - 8*C)) - b**3*(-6*B*c**2*d + 16*B*d**3 + C*c**3 - 24*c*d**2*(A - C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(8*b**(sympy.S(5)/2)*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_137():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(a + b*tan(e + f*x))
    F = C*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)/(2*b*f) + (c + I*d)**(sympy.S(3)/2)*(I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(4*B*b*d - 3*C*a*d + 3*C*b*c)/(4*b**2*f) + (3*C*a**2*d**2 - 2*a*b*d*(2*B*d + 3*C*c) + b**2*(12*B*c*d + 3*C*c**2 + d**2*(8*A - 8*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*b**(sympy.S(5)/2)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_138():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)) + d*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*B*a*b + 3*C*a**2 + C*b**2)/(b**2*f*(a**2 + b**2)) + sqrt(d)*(2*B*b*d - 3*C*a*d + 3*C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_139():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*C*d**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(5)/2)*f) - (B - I*(A - C))*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)) - sqrt(c + d*tan(e + f*x))*(2*C*a**4*d - 2*a**2*b**2*(B*c + d*(A - 3*C)) + 4*a*b**3*(A*c - B*d - C*c) + 2*b**4*(A*d + B*c))/(b**2*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_140():
    f = (c + d*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(7)/2)
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(7)/2)) - (c - I*d)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(7)/2)) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(a**2 + b**2)) - sqrt(c + d*tan(e + f*x))*(4*B*a**5*b*d**2 + 6*C*a**6*d**2 + 2*a**4*b**2*d*(10*B*c + d*(8*A + C)) - 2*a**3*b**3*(B*(15*c**2 - 39*d**2) + 50*c*d*(A - C)) + 2*a**2*b**4*(45*A*c**2 - 49*A*d**2 - 90*B*c*d - 45*C*c**2 + 58*C*d**2) + 2*a*b**5*(B*(45*c**2 - 23*d**2) + 70*c*d*(A - C)) + 2*b**6*(-3*A*(5*c**2 - d**2) + 5*c*(4*B*d + 3*C*c)))/(15*b**2*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**3*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(4*B*a**3*b*d + 6*C*a**4*d - 2*a**2*b**2*(7*A*d + 5*B*c - 13*C*d) + 4*a*b**3*(5*A*c - 4*B*d - 5*C*c) + 2*b**4*(3*A*d + 5*B*c))/(15*b**2*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_141():
    f = sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(7)/2)/(4*d*f) - sqrt(a - I*b)*(c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/f + sqrt(a + I*b)*(c + I*d)**(sympy.S(5)/2)*(I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/f - sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)*(-8*B*b*d - C*a*d + C*b*c)/(24*b*d*f) + sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(48*b*d**2*(A*b + B*a - C*b) - (-5*a*d + 5*b*c)*(-8*B*b*d - C*a*d + C*b*c))/(96*b**2*d*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(64*b**2*d**2*(A*a*d + A*b*c + B*a*c - B*b*d - C*a*d - C*b*c) + (-a*d + b*c)*(48*b*d**2*(A*b + B*a - C*b) - (-5*a*d + 5*b*c)*(-8*B*b*d - C*a*d + C*b*c)))/(64*b**3*d*f) - (5*C*a**4*d**4 - 4*a**3*b*d**3*(2*B*d + 5*C*c) + 2*a**2*b**2*d**2*(20*B*c*d + 15*C*c**2 + d**2*(8*A - 8*C)) - 4*a*b**3*d*(30*B*c**2*d - 16*B*d**3 + 5*C*c**3 + 40*c*d**2*(A - C)) + b**4*(-40*B*c**3*d + 320*B*c*d**3 + 5*C*c**4 - 240*c**2*d**2*(A - C) + d**4*(128*A - 128*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(64*b**(sympy.S(7)/2)*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_142():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(a + b*tan(e + f*x))
    F = C*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2)/(3*b*f) - (B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)) + sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(6*B*b*d - 5*C*a*d + 5*C*b*c)/(12*b**2*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(8*b**2*d*(B*c + d*(A - C)) + (-a*d + b*c)*(6*B*b*d - 5*C*a*d + 5*C*b*c))/(8*b**3*f) - (5*C*a**3*d**3 - 3*a**2*b*d**2*(2*B*d + 5*C*c) + a*b**2*d*(20*B*c*d + 15*C*c**2 + d**2*(8*A - 8*C)) - b**3*(30*B*c**2*d - 16*B*d**3 + 5*C*c**3 + 40*c*d**2*(A - C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(8*b**(sympy.S(7)/2)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_143():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(3)/2)
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(b*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)) + d*sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(4*A*b**2 - 4*B*a*b + 5*C*a**2 + C*b**2)/(2*b**2*f*(a**2 + b**2)) - d*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-8*A*b**2*(-a*d + b*c) + 15*C*a**3*d - 3*a**2*b*(4*B*d + 5*C*c) + a*b**2*(8*B*c + 7*C*d) - b**3*(4*B*d + 7*C*c))/(4*b**3*f*(a**2 + b**2)) + sqrt(d)*(15*C*a**2*d**2 - 6*a*b*d*(2*B*d + 5*C*c) + b**2*(20*B*c*d + 15*C*c**2 + d**2*(8*A - 8*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*b**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_144():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(5)/2)
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)) + (c + d*tan(e + f*x))**(sympy.S(3)/2)*(4*B*a**3*b*d - 10*C*a**4*d + 2*a**2*b**2*(3*B*c + d*(A - 11*C)) - 4*a*b**3*(3*A*c - 4*B*d - 3*C*c) - 2*b**4*(5*A*d + 3*B*c))/(3*b**2*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2) - d*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(2*B*a**3*b*d - 5*C*a**4*d + 2*a**2*b**2*(B*c - 5*C*d) - 2*a*b**3*(2*A*c - 3*B*d - 2*C*c) - b**4*(2*B*c + d*(4*A + C)))/(b**3*f*(a**2 + b**2)**2) + d**(sympy.S(3)/2)*(2*B*b*d - 5*C*a*d + 5*C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_145():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(7)/2)
    F = 2*C*d**(sympy.S(5)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(b**(sympy.S(7)/2)*f) - (B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(7)/2)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(7)/2)) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(a**2 + b**2)) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(2*C*a**4*d - 2*a**2*b**2*(B*c + d*(A - 3*C)) + 4*a*b**3*(A*c - B*d - C*c) + 2*b**4*(A*d + B*c))/(3*b**2*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(2*C*a**6*d**2 + 6*C*a**4*b**2*d**2 - 2*a**3*b**3*(B*(c**2 - d**2) + 2*c*d*(A - C)) - 6*a**2*b**4*(-A*(c**2 - d**2) + 2*B*c*d + C*c**2 - 2*C*d**2) + 6*a*b**5*(B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b**6*(-A*(c**2 - d**2) + c*(2*B*d + C*c)))/(b**3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_146():
    f = (c + d*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(a + b*tan(e + f*x))**(sympy.S(9)/2)
    F = -(B - I*(A - C))*(c + I*d)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(9)/2)) - (c - I*d)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(9)/2)) - (c + d*tan(e + f*x))**(sympy.S(5)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(7*b*f*(a + b*tan(e + f*x))**(sympy.S(7)/2)*(a**2 + b**2)) - (c + d*tan(e + f*x))**(sympy.S(3)/2)*(4*B*a**3*b*d + 10*C*a**4*d - 2*a**2*b**2*(9*A*d + 7*B*c - 19*C*d) + 4*a*b**3*(7*A*c - 6*B*d - 7*C*c) + 2*b**4*(5*A*d + 7*B*c))/(35*b**2*f*(a + b*tan(e + f*x))**(sympy.S(5)/2)*(a**2 + b**2)**2) - sqrt(c + d*tan(e + f*x))*(12*B*a**7*b*d**3 + 30*C*a**8*d**3 + 4*a**6*b**2*d**2*(4*A*d + 7*B*c + 26*C*d) + 4*a**5*b**3*d*(B*(35*c**2 - 12*d**2) + 56*c*d*(A - C)) - 2*a**4*b**4*(525*A*c**2*d - 311*A*d**3 + 105*B*c**3 - 749*B*c*d**2 - 525*C*c**2*d + 221*C*d**3) - 4*a**3*b**5*(-42*A*(5*c**3 - 19*c*d**2) + 700*B*c**2*d - 317*B*d**3 + 210*C*c**3 - 798*C*c*d**2) + 4*a**2*b**6*(875*A*c**2*d - 261*A*d**3 + 315*B*c**3 - 812*B*c*d**2 - 875*C*c**2*d + 291*C*d**3) - 4*a*b**7*(210*A*c**3 - 406*A*c*d**2 - 525*B*c**2*d + 88*B*d**3 - 210*C*c**3 + 406*C*c*d**2) - 2*b**8*(7*B*(15*c**3 - 23*c*d**2) + 5*d*(49*A*c**2 - 3*A*d**2 - 49*C*c**2)))/(105*b**3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**4*(-a*d + b*c)) - sqrt(c + d*tan(e + f*x))*(12*B*a**5*b*d**2 + 30*C*a**6*d**2 + 2*a**4*b**2*d*(8*A*d + 14*B*c + 37*C*d) - 2*a**3*b**3*(B*(35*c**2 - 75*d**2) + 98*c*d*(A - C)) + 6*a**2*b**4*(35*A*c**2 - 39*A*d**2 - 70*B*c*d - 35*C*c**2 + 54*C*d**2) + 2*a*b**5*(B*(105*c**2 - 71*d**2) + 182*c*d*(A - C)) + 2*b**6*(-5*A*(7*c**2 - 3*d**2) + 7*c*(8*B*d + 5*C*c)))/(105*b**3*f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_147():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = C*(a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x))/(3*d*f) - (B - I*(A - C))*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d)) - (a - I*b)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) - (a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))*(-6*B*b*d - 5*C*a*d + 5*C*b*c)/(12*d**2*f) + sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(8*b*d**2*(A*b + B*a - C*b) + (-a*d + b*c)*(-6*B*b*d - 5*C*a*d + 5*C*b*c))/(8*d**3*f) + (5*C*a**3*d**3 - 15*a**2*b*d**2*(-2*B*d + C*c) + 5*a*b**2*d*(-4*B*c*d + 3*C*c**2 + d**2*(8*A - 8*C)) - b**3*(-6*B*c**2*d + 16*B*d**3 + 5*C*c**3 + 8*c*d**2*(A - C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(8*sqrt(b)*d**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_148():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = C*(a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))/(2*d*f) - (a - I*b)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) + (a + I*b)**(sympy.S(3)/2)*(I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d)) - sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-4*B*b*d - 3*C*a*d + 3*C*b*c)/(4*d**2*f) + (3*C*a**2*d**2 - 6*a*b*d*(-2*B*d + C*c) + b**2*(-4*B*c*d + 3*C*c**2 + d**2*(8*A - 8*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*sqrt(b)*d**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_149():
    f = sqrt(a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/sqrt(c + d*tan(e + f*x))
    F = C*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))/(d*f) - sqrt(a - I*b)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c - I*d)) + sqrt(a + I*b)*(I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(c + I*d)) - (-2*B*b*d - C*a*d + C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(b)*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_150():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x)))
    F = 2*C*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(sqrt(b)*sqrt(d)*f) - (B + I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*sqrt(c - I*d)) + (I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*sqrt(c + I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_151():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x)))
    F = -(B - I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*sqrt(c + I*d)) - sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*a*(B*b - C*a))/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_152():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*tan(e + f*x)))
    F = -(B - I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)*sqrt(c + I*d)) - sqrt(c + d*tan(e + f*x))*(10*B*a**3*b*d - 4*C*a**4*d - 2*a**2*b**2*(8*A*d + 3*B*c - 4*C*d) + 2*a*b**3*(6*A*c - B*d - 6*C*c) + 2*b**4*(-2*A*d + 3*B*c))/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*(-a*d + b*c)**2) - sqrt(c + d*tan(e + f*x))*(2*A*b**2 - 2*a*(B*b - C*a))/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)*sqrt(c - I*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_153():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = sqrt(b)*(15*C*a**2*d**2 - 10*a*b*d*(-2*B*d + 3*C*c) + b**2*(-12*B*c*d + 15*C*c**2 + d**2*(8*A - 8*C)))*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(4*d**(sympy.S(7)/2)*f) + b*(a + b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*tan(e + f*x))*(-4*B*c*d + 5*C*c**2 + d**2*(4*A + C))/(2*d**2*f*(c**2 + d**2)) - b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-4*d**2*(B*(a*c + b*d) + (A - C)*(-a*d + b*c)) + (-3*a*d + 3*b*c)*(-4*B*c*d + 5*C*c**2 + d**2*(4*A + C)))/(4*d**3*f*(c**2 + d**2)) - (B - I*(A - C))*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) - (a - I*b)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) - (a + b*tan(e + f*x))**(sympy.S(5)/2)*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_154():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(b)*(-2*B*b*d - 3*C*a*d + 3*C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f) + b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(-2*B*c*d + 3*C*c**2 + d**2*(2*A + C))/(d**2*f*(c**2 + d**2)) - (B - I*(A - C))*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) - (a - I*b)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) - (a + b*tan(e + f*x))**(sympy.S(3)/2)*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_155():
    f = sqrt(a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*C*sqrt(b)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(3)/2)*f) - (B - I*(A - C))*sqrt(a + I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(3)/2)) - sqrt(a - I*b)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(3)/2)) - sqrt(a + b*tan(e + f*x))*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_156():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -(B + I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*(c - I*d)**(sympy.S(3)/2)) + sqrt(a + b*tan(e + f*x))*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) + (I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*(c + I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_157():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -2*d*sqrt(a + b*tan(e + f*x))*(A*(a**2*d**2 + b**2*(c**2 + 2*d**2)) - B*a*b*(c**2 + d**2) + a**2*(-B*c*d + 2*C*c**2 + C*d**2) + b**2*c*(-B*d + C*c))/(f*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**2) - (B - I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(3)/2)) - (2*A*b**2 - 2*a*(B*b - C*a))/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_158():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**(sympy.S(5)/2)*(c + d*tan(e + f*x))**(sympy.S(3)/2))
    F = -2*d*sqrt(a + b*tan(e + f*x))*(8*B*a**3*b*d*(c**2 + d**2) - a**4*d*(-3*B*c*d + 8*C*c**2 + d**2*(3*A + 5*C)) - a**2*b**2*(11*A*c**2*d + 17*A*d**3 + 3*B*c**3 - 3*B*c*d**2 + 5*C*c**2*d - C*d**3) + 2*a*b**3*(c**2 + d**2)*(3*A*c + B*d - 3*C*c) - b**4*(-3*B*(c**3 + 2*c*d**2) + d*(5*A*c**2 + 8*A*d**2 + 3*C*c**2)))/(3*f*(a**2 + b**2)**2*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)**3) - (B - I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(5)/2)*(c + I*d)**(sympy.S(3)/2)) - (14*B*a**3*b*d - 8*C*a**4*d - 2*a**2*b**2*(3*B*c + d*(10*A - 2*C)) + 2*a*b**3*(6*A*c + B*d - 6*C*c) + 2*b**4*(-4*A*d + 3*B*c))/(3*f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)**2*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)**2) - (2*A*b**2 - 2*a*(B*b - C*a))/(f*(a + b*tan(e + f*x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)*sqrt(c + d*tan(e + f*x))*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(5)/2)*(c - I*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_159():
    f = (a + b*tan(e + f*x))**(sympy.S(5)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -b**(sympy.S(3)/2)*(-2*B*b*d - 5*C*a*d + 5*C*b*c)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(7)/2)*f) + b*sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))*(2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + b*(-2*B*c**3*d - 6*B*c*d**3 + 5*C*c**4 + 10*C*c**2*d**2 + d**4*(4*A + C)))/(d**3*f*(c**2 + d**2)**2) - (B - I*(A - C))*(a + I*b)**(sympy.S(5)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) - (a - I*b)**(sympy.S(5)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) - (a + b*tan(e + f*x))**(sympy.S(5)/2)*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - (a + b*tan(e + f*x))**(sympy.S(3)/2)*(6*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(5*A*d**4 - 2*B*c**3*d - 8*B*c*d**3 + 5*C*c**4 - c**2*d**2*(A - 11*C)))/(3*d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_160():
    f = (a + b*tan(e + f*x))**(sympy.S(3)/2)*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = 2*C*b**(sympy.S(3)/2)*atanh(sqrt(d)*sqrt(a + b*tan(e + f*x))/(sqrt(b)*sqrt(c + d*tan(e + f*x))))/(d**(sympy.S(5)/2)*f) - (B - I*(A - C))*(a + I*b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) - (a - I*b)**(sympy.S(3)/2)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) - (a + b*tan(e + f*x))**(sympy.S(3)/2)*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)) - sqrt(a + b*tan(e + f*x))*(2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(A*d**4 - 2*B*c*d**3 + C*c**4 - c**2*d**2*(A - 3*C)))/(d**2*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_161():
    f = sqrt(a + b*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**(sympy.S(5)/2)
    F = -(B - I*(A - C))*sqrt(a + I*b)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c + I*d)**(sympy.S(5)/2)) - sqrt(a - I*b)*(I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(c - I*d)**(sympy.S(5)/2)) + sqrt(a + b*tan(e + f*x))*(6*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(A*d**4 + 2*B*c**3*d - 4*B*c*d**3 + C*c**4 - c**2*d**2*(5*A - 7*C)))/(3*d*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)) - sqrt(a + b*tan(e + f*x))*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(3*d*f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_162():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -(B + I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a - I*b)*(c - I*d)**(sympy.S(5)/2)) + sqrt(a + b*tan(e + f*x))*(-6*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) + 2*b*(2*A*d**4 - 5*B*c**3*d + B*c*d**3 + 2*C*c**4 + 4*c**2*d**2*(2*A - C)))/(3*f*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + sqrt(a + b*tan(e + f*x))*(2*A*d**2 - 2*B*c*d + 2*C*c**2)/(f*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-3*a*d + 3*b*c)) + (I*A - B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*sqrt(a + I*b)*(c + I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_163():
    f = (A + B*tan(e + f*x) + C*tan(e + f*x)**2)/((a + b*tan(e + f*x))**(sympy.S(3)/2)*(c + d*tan(e + f*x))**(sympy.S(5)/2))
    F = -2*d*sqrt(a + b*tan(e + f*x))*(-A*(6*a**3*c*d**3 - a**2*b*d**2*(11*c**2 + 5*d**2) + 6*a*b**2*c*d**3 - b**3*(3*c**4 + 17*c**2*d**2 + 8*d**4)) + 3*a**3*d**2*(B*(c**2 - d**2) + 2*C*c*d) + a**2*b*(-8*B*c**3*d - 2*B*c*d**3 + 8*C*c**4 + 5*C*c**2*d**2 + 3*C*d**4) + 3*a*b**2*(-B*(c**4 + c**2*d**2 + 2*d**4) + 2*C*c*d**3) + b**3*c*(-8*B*c**2*d - 2*B*d**3 + 5*C*c**3 - C*c*d**2))/(f*(3*a**2 + 3*b**2)*sqrt(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**3) - 2*d*sqrt(a + b*tan(e + f*x))*(A*(a**2*d**2 + b**2*(3*c**2 + 4*d**2)) - 3*B*a*b*(c**2 + d**2) + a**2*(-B*c*d + 4*C*c**2 + 3*C*d**2) + b**2*c*(-B*d + C*c))/(f*(3*a**2 + 3*b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(c**2 + d**2)*(-a*d + b*c)**2) - (B - I*(A - C))*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c + I*d)/(sqrt(a + I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a + I*b)**(sympy.S(3)/2)*(c + I*d)**(sympy.S(5)/2)) - (2*A*b**2 - 2*a*(B*b - C*a))/(f*sqrt(a + b*tan(e + f*x))*(a**2 + b**2)*(c + d*tan(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - (I*A + B - I*C)*atanh(sqrt(a + b*tan(e + f*x))*sqrt(c - I*d)/(sqrt(a - I*b)*sqrt(c + d*tan(e + f*x))))/(f*(a - I*b)**(sympy.S(3)/2)*(c - I*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_164():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**n*hyper((-n, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(b*f*(b*(c + d*tan(e + f*x))/(-a*d + b*c))**n*(m + 1)) - (B + I*(A - C))*(a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**n*appellf1(m + 1, 1, -n, m + 2, (a + b*tan(e + f*x))/(a - I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(b*(c + d*tan(e + f*x))/(-a*d + b*c))**n*(2*a - 2*I*b)*(m + 1)) - (a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**n*(A + I*B - C)*appellf1(m + 1, 1, -n, m + 2, (a + b*tan(e + f*x))/(a + I*b), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(b*(c + d*tan(e + f*x))/(-a*d + b*c))**n*(m + 1)*(2*I*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_165():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**3*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**3/(b*f*(m + 4)) + (a + b*tan(e + f*x))**(m + 1)*(c - I*d)**3*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(c + I*d)**3*(A + I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(m + 1)*(2*I*a - 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**2*(3*C*a*d - b*(B*d*(m + 4) + 3*C*c))/(b**2*f*(m + 3)*(m + 4)) + d*(a + b*tan(e + f*x))**(m + 1)*(b**2*d*(m + 3)*(m + 4)*(B*c + d*(A - C)) - (-2*a*d + 2*b*c)*(3*C*a*d - b*(B*d*(m + 4) + 3*C*c)))*tan(e + f*x)/(b**3*f*(m + 2)*(m + 3)*(m + 4)) + (a + b*tan(e + f*x))**(m + 1)*(b*c*(m + 2)*(b**2*d*(m + 3)*(m + 4)*(B*c + d*(A - C)) - (-2*a*d + 2*b*c)*(3*C*a*d - b*(B*d*(m + 4) + 3*C*c))) + d*(-a*(b**2*d*(m + 3)*(m + 4)*(B*c + d*(A - C)) - (-2*a*d + 2*b*c)*(3*C*a*d - b*(B*d*(m + 4) + 3*C*c))) + b**3*(m + 2)*(m + 3)*(m + 4)*(B*(c**2 - d**2) + 2*c*d*(A - C))))/(b**4*f*(m + 1)*(m + 2)*(m + 3)*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_166():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**2*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(m + 1)*(c + d*tan(e + f*x))**2/(b*f*(m + 3)) + (a + b*tan(e + f*x))**(m + 1)*(c - I*d)**2*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*(c + I*d)**2*(I*A - B - I*C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1)) - d*(a + b*tan(e + f*x))**(m + 1)*(2*C*a*d - b*(B*d*(m + 3) + 2*C*c))*tan(e + f*x)/(b**2*f*(m + 2)*(m + 3)) + (a + b*tan(e + f*x))**(m + 1)*(2*C*a**2*d**2 - a*b*d*(m + 3)*(B*d + 2*C*c) + b**2*(m + 2)*(2*B*c*d*(m + 3) + 2*C*c**2 + d**2*(A - C)*(m + 3)))/(b**3*f*(m + 1)*(m + 2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_167():
    f = (a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*d*(a + b*tan(e + f*x))**(m + 1)*tan(e + f*x)/(b*f*(m + 2)) + (a + b*tan(e + f*x))**(m + 1)*(c - I*d)*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(c + I*d)*(A + I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(m + 1)*(2*I*a - 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(C*a*d - b*(m + 2)*(B*d + C*c))/(b**2*f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_168():
    f = (a + b*tan(e + f*x))**m*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)
    F = C*(a + b*tan(e + f*x))**(m + 1)/(b*f*(m + 1)) + (a + b*tan(e + f*x))**(m + 1)*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*(I*A - B - I*C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_169():
    f = (a + b*tan(e + f*x))**m*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))
    F = (a + b*tan(e + f*x))**(m + 1)*(A*d**2 - B*c*d + C*c**2)*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(c**2 + d**2)*(m + 1)*(-a*d + b*c)) - (a + b*tan(e + f*x))**(m + 1)*(A + I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(c + I*d)*(m + 1)*(2*I*a - 2*b)) - (a + b*tan(e + f*x))**(m + 1)*(I*A + B - I*C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(2*a - 2*I*b)*(c - I*d)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_170():
    f = (a + b*tan(e + f*x))**m*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**2
    F = -(a + b*tan(e + f*x))**(m + 1)*(a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - b*(A*d**2*(c**2*(2 - m) - d**2*m) - B*c*d*(c**2*(1 - m) - d**2*(m + 1)) - C*c**2*(c**2*m + d**2*(m + 2))))*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(f*(c**2 + d**2)**2*(m + 1)*(-a*d + b*c)**2) + (a + b*tan(e + f*x))**(m + 1)*(A*d**2 - B*c*d + C*c**2)/(f*(c + d*tan(e + f*x))*(c**2 + d**2)*(-a*d + b*c)) + (a + b*tan(e + f*x))**(m + 1)*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(c - I*d)**2*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*(I*A - B - I*C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(c + I*d)**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_4_2_a_plus_b_tan_pow_m_c_plus_d_tan_pow_n_A_plus_B_tan_plus_C_tan_pow_2_171():
    f = (a + b*tan(e + f*x))**m*(A + B*tan(e + f*x) + C*tan(e + f*x)**2)/(c + d*tan(e + f*x))**3
    F = (a + b*tan(e + f*x))**(m + 1)*(2*a**2*d**3*(-B*(c**3 - 3*c*d**2) + d*(A - C)*(3*c**2 - d**2)) - 2*a*b*d**2*(B*(-c**4*(2 - m) + 6*c**2*d**2 - d**4*m) + 2*c*d*(A - C)*(c**2*(3 - m) - d**2*(m + 1))) - b**2*(A*d**2*(-c**4*(m**2 - 5*m + 6) + 2*c**2*d**2*(-m**2 + 3*m + 1) + d**4*m*(1 - m)) + B*c*d*(c**4*(m**2 - 3*m + 2) - 2*c**2*d**2*(-m**2 + m + 3) + d**4*m*(m + 1)) + C*c**2*(c**4*m*(1 - m) + 2*c**2*d**2*(-m**2 - m + 3) - d**4*(m**2 + 3*m + 2))))*hyper((1, m + 1), (m + 2,), -d*(a + b*tan(e + f*x))/(-a*d + b*c))/(2*f*(c**2 + d**2)**3*(m + 1)*(-a*d + b*c)**3) - (a + b*tan(e + f*x))**(m + 1)*(2*a*d**2*(-B*(c**2 - d**2) + 2*c*d*(A - C)) - b*(A*d**4*(1 - m) - B*c**3*d*(3 - m) + B*c*d**3*(m + 1) + C*c**4*(1 - m) + c**2*d**2*(A*(5 - m) - C*(m + 3))))/(2*f*(c + d*tan(e + f*x))*(c**2 + d**2)**2*(-a*d + b*c)**2) + (a + b*tan(e + f*x))**(m + 1)*(A*d**2 - B*c*d + C*c**2)/(f*(c + d*tan(e + f*x))**2*(c**2 + d**2)*(-2*a*d + 2*b*c)) + (a + b*tan(e + f*x))**(m + 1)*(A - I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a - I*b))/(f*(c - I*d)**3*(m + 1)*(2*I*a + 2*b)) + (a + b*tan(e + f*x))**(m + 1)*(A + I*B - C)*hyper((1, m + 1), (m + 2,), (a + b*tan(e + f*x))/(a + I*b))/(f*(2*a + 2*I*b)*(m + 1)*(I*c - d)**3)
    assert integrate(f, x) == F

