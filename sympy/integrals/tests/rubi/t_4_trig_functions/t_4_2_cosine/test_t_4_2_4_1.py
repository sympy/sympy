"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.4.1 (a+b cos)^m (A+B cos+C cos^2).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, m, n = symbols('A B C a b c d e f m n')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_1():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**7
    F = C*sin(c + d*x)**9/(9*d) + (A + C)*sin(c + d*x)/d - (A + 4*C)*sin(c + d*x)**7/(7*d) - (3*A + 4*C)*sin(c + d*x)**3/(3*d) + (3*A + 6*C)*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_2():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**5
    F = -C*sin(c + d*x)**7/(7*d) + (A + C)*sin(c + d*x)/d + (A + 3*C)*sin(c + d*x)**5/(5*d) - (2*A + 3*C)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_3():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**3
    F = C*sin(c + d*x)**5/(5*d) + (A + C)*sin(c + d*x)/d - (A + 2*C)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_4():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = -C*sin(c + d*x)**3/(3*d) + (A + C)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_5():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = A*atanh(sin(c + d*x))/d + C*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_6():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = A*tan(c + d*x)*sec(c + d*x)/(2*d) + (A + 2*C)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_7():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = A*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (3*A + 4*C)*tan(c + d*x)*sec(c + d*x)/(8*d) + (3*A + 4*C)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_8():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**7
    F = A*tan(c + d*x)*sec(c + d*x)**5/(6*d) + (5*A + 6*C)*tan(c + d*x)*sec(c + d*x)**3/(24*d) + (5*A + 6*C)*tan(c + d*x)*sec(c + d*x)/(16*d) + (5*A + 6*C)*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_9():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**6
    F = C*sin(c + d*x)*cos(c + d*x)**7/(8*d) + x*(40*A + 35*C)/128 + (8*A + 7*C)*sin(c + d*x)*cos(c + d*x)**5/(48*d) + (40*A + 35*C)*sin(c + d*x)*cos(c + d*x)**3/(192*d) + (40*A + 35*C)*sin(c + d*x)*cos(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_10():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**4
    F = C*sin(c + d*x)*cos(c + d*x)**5/(6*d) + x*(6*A + 5*C)/16 + (6*A + 5*C)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (6*A + 5*C)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_11():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = C*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(4*A + 3*C)/8 + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_12():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = A*tan(c + d*x)/d + C*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_13():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = A*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (2*A + 3*C)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_14():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**6
    F = A*tan(c + d*x)*sec(c + d*x)**4/(5*d) + (4*A + 5*C)*tan(c + d*x)**3/(15*d) + (4*A + 5*C)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_15():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**8
    F = A*tan(c + d*x)*sec(c + d*x)**6/(7*d) + (6*A + 7*C)*tan(c + d*x)**5/(35*d) + (6*A + 7*C)*tan(c + d*x)/(7*d) + (12*A + 14*C)*tan(c + d*x)**3/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_16():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d) + 2*b**2*sqrt(b*cos(c + d*x))*(9*A + 7*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 2*b*(b*cos(c + d*x))**(sympy.S(3)/2)*(9*A + 7*C)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_17():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) + 2*b**2*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 2*b*sqrt(b*cos(c + d*x))*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_18():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_19():
    f = (A + C*cos(c + d*x)**2)/sqrt(b*cos(c + d*x))
    F = 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_20():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_21():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_22():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*A*sin(c + d*x)/(5*b*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*sin(c + d*x)/(5*b**3*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_23():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(9)/2)
    F = 2*A*sin(c + d*x)/(7*b*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + (10*A + 14*C)*sin(c + d*x)/(21*b**3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**4*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_24():
    f = (3 - 5*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = -2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_25():
    f = (1 - 3*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = -2*sin(c + d*x)*sqrt(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_26():
    f = (b*sec(c + d*x))**(sympy.S(9)/2)*(A + C*cos(c + d*x)**2)
    F = 2*A*b**2*(b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*d) + 2*b**4*sqrt(b*sec(c + d*x))*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b**3*(b*sec(c + d*x))**(sympy.S(3)/2)*(5*A + 7*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_27():
    f = (b*sec(c + d*x))**(sympy.S(7)/2)*(A + C*cos(c + d*x)**2)
    F = 2*A*b**2*(b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) - 2*b**4*(3*A + 5*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*b**3*sqrt(b*sec(c + d*x))*(3*A + 5*C)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_28():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)
    F = 2*A*b**2*sqrt(b*sec(c + d*x))*tan(c + d*x)/(3*d) + 2*b**2*sqrt(b*sec(c + d*x))*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_29():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)
    F = 2*A*b**2*tan(c + d*x)/(d*sqrt(b*sec(c + d*x))) - 2*b**2*(A - C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_30():
    f = sqrt(b*sec(c + d*x))*(A + C*cos(c + d*x)**2)
    F = 2*C*b**2*tan(c + d*x)/(3*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + sqrt(b*sec(c + d*x))*(6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_31():
    f = (A + C*cos(c + d*x)**2)/sqrt(b*sec(c + d*x))
    F = 2*C*b**2*tan(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + (10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_32():
    f = (A + C*cos(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*C*b**2*tan(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + (14*A + 10*C)*sin(c + d*x)/(21*b*d*sqrt(b*sec(c + d*x))) + sqrt(b*sec(c + d*x))*(14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_33():
    f = (A + C*cos(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*C*b**2*tan(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(9)/2)) + (18*A + 14*C)*sin(c + d*x)/(45*b*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + (18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_34():
    f = (b*cos(c + d*x))**m*(A + C*cos(c + d*x)**2)
    F = C*(b*cos(c + d*x))**(m + 1)*sin(c + d*x)/(b*d*(m + 2)) - (b*cos(c + d*x))**(m + 1)*(A*(m + 2) + C*(m + 1))*sin(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(m + 1)*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_35():
    f = (b*cos(c + d*x))**m*(-C*(m + 1)/(m + 2) + C*cos(c + d*x)**2)
    F = C*(b*cos(c + d*x))**(m + 1)*sin(c + d*x)/(b*d*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_36():
    f = (b*cos(c + d*x))**m*(A - A*(m + 2)*cos(c + d*x)**2/(m + 1))
    F = -A*(b*cos(c + d*x))**(m + 1)*sin(c + d*x)/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_37():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**3*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_38():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**2*d) + 2*b*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_39():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_40():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_41():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*A*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_42():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_43():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*b*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_44():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*b**2*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_45():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**2*d) + 2*b*sqrt(b*cos(c + d*x))*(9*A + 7*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_46():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) + 2*b**2*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 2*b*sqrt(b*cos(c + d*x))*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_47():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*b*sqrt(b*cos(c + d*x))*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_48():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b**2*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_49():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*b*sqrt(b*cos(c + d*x))*(A - C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_50():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**2*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_51():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*b**2*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 2*b*sqrt(b*cos(c + d*x))*(3*A + 5*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_52():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**6
    F = 2*A*b**5*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*b**3*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**2*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_53():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d) + 2*b**2*sqrt(b*cos(c + d*x))*(9*A + 7*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 2*b*(b*cos(c + d*x))**(sympy.S(3)/2)*(9*A + 7*C)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_54():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + 2*b**3*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 2*b**2*sqrt(b*cos(c + d*x))*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_55():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*C*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*b**2*sqrt(b*cos(c + d*x))*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_56():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b**3*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_57():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*(A - C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_58():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**3*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_59():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**6
    F = 2*A*b**5*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*b**3*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*(3*A + 5*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_60():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**7
    F = 2*A*b**6*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*b**4*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**3*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_61():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**4/sqrt(b*cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**5*d) + (110*A + 90*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(110*A + 90*C)*sin(c + d*x)/(231*b*d) + (b*cos(c + d*x))**(sympy.S(5)/2)*(22*A + 18*C)*sin(c + d*x)/(77*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_62():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**4*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_63():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**3*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_64():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_65():
    f = (A + C*cos(c + d*x)**2)/sqrt(b*cos(c + d*x))
    F = 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_66():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_67():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*A*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_68():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 2*A*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_69():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**4/sqrt(b*cos(c + d*x))
    F = 2*A*b**3*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*b*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_70():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**5/sqrt(b*cos(c + d*x))
    F = 2*A*b**4*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 2*b**2*(7*A + 9*C)*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (14*A + 18*C)*sin(c + d*x)/(15*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(14*A + 18*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_71():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**5*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_72():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**4*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_73():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_74():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_75():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_76():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_77():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*sin(c + d*x)/(5*b*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_78():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*b**2*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + (10*A + 14*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_79():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**5/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**6*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**3*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_80():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**5*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_81():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_82():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_83():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_84():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_85():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*sin(c + d*x)/(5*b**2*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_86():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*b*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + (10*A + 14*C)*sin(c + d*x)/(21*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_87():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*A*sin(c + d*x)/(5*b*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*sin(c + d*x)/(5*b**3*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_88():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(9)/2)
    F = 2*A*sin(c + d*x)/(7*b*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + (10*A + 14*C)*sin(c + d*x)/(21*b**3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**4*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_89():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)
    F = C*sqrt(b*cos(c + d*x))*sin(c + d*x)**5/(5*d*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - sqrt(b*cos(c + d*x))*(A + 2*C)*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_90():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = C*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_91():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = -C*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_92():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = A*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_93():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = A*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_94():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + C*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_95():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_96():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_97():
    f = sqrt(b*cos(c + d*x))*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_98():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**5/(5*d*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - b*sqrt(b*cos(c + d*x))*(A + 2*C)*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_99():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + b*x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_100():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = -C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_101():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = A*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_102():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_103():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + C*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_104():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + b*sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_105():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + b*sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_106():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(13)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + b*sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b*sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_107():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**5/(5*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - b**2*sqrt(b*cos(c + d*x))*(A + 2*C)*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_108():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + b**2*x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_109():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = -C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(A + C)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_110():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_111():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_112():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + C*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_113():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + b**2*sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_114():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(13)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + b**2*sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_115():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(15)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + b**2*sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b**2*sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_116():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/sqrt(b*cos(c + d*x))
    F = C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_117():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/sqrt(b*cos(c + d*x))
    F = -C*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x))) + (A + C)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_118():
    f = (A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    F = A*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x)) + C*x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_119():
    f = (A + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_120():
    f = (A + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_121():
    f = (A + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_122():
    f = (A + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = A*sin(c + d*x)/(3*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + (2*A + 3*C)*sin(c + d*x)/(3*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_123():
    f = (A + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(9)/2))
    F = A*sin(c + d*x)/(4*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + (3*A + 4*C)*sin(c + d*x)/(8*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_124():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*b*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_125():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = -C*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b*d*sqrt(b*cos(c + d*x))) + (A + C)*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_126():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_127():
    f = (A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_128():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_129():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_130():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(3*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + (2*A + 3*C)*sin(c + d*x)/(3*b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_131():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = A*sin(c + d*x)/(4*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + (3*A + 4*C)*sin(c + d*x)/(8*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_132():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(9)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b**2*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*b**2*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_133():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = -C*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b**2*d*sqrt(b*cos(c + d*x))) + (A + C)*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_134():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_135():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_136():
    f = (A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_137():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_138():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(3*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + (2*A + 3*C)*sin(c + d*x)/(3*b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_139():
    f = (A + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(4*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + (3*A + 4*C)*sin(c + d*x)/(8*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_140():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = 3*C*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)/(13*b**3*d) - (b*cos(c + d*x))**(sympy.S(10)/3)*(39*A + 30*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(130*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_141():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)/(10*b**2*d) - (b*cos(c + d*x))**(sympy.S(7)/3)*(30*A + 21*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(70*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_142():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*b*d) - (b*cos(c + d*x))**(sympy.S(4)/3)*(21*A + 12*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(28*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_143():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*d) - (b*cos(c + d*x))**(sympy.S(1)/3)*(12*A + 3*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_144():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 3*A*b*sin(c + d*x)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)) + (b*cos(c + d*x))**(sympy.S(4)/3)*(3*A - 6*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(8*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_145():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)) - (b*cos(c + d*x))**(sympy.S(1)/3)*(6*A + 15*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(5*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_146():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = 3*C*(b*cos(c + d*x))**(sympy.S(11)/3)*sin(c + d*x)/(14*b**3*d) - (b*cos(c + d*x))**(sympy.S(11)/3)*(42*A + 33*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(17)/6,), cos(c + d*x)**2)/(154*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_147():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)/(11*b**2*d) - (b*cos(c + d*x))**(sympy.S(8)/3)*(33*A + 24*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(88*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_148():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_149():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_150():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 3*A*b*sin(c + d*x)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_151():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_152():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = 3*C*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)/(16*b**3*d) - (b*cos(c + d*x))**(sympy.S(13)/3)*(48*A + 39*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(208*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_153():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)/(13*b**2*d) - (b*cos(c + d*x))**(sympy.S(10)/3)*(39*A + 30*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(130*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_154():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)/(10*b*d) - (b*cos(c + d*x))**(sympy.S(7)/3)*(30*A + 21*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(70*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_155():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*d) - (b*cos(c + d*x))**(sympy.S(4)/3)*(21*A + 12*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(28*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_156():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 3*C*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*d) - 3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*(4*A + C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_157():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)) + (b*cos(c + d*x))**(sympy.S(4)/3)*(3*A - 6*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(8*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_158():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)/(11*b**3*d) - (b*cos(c + d*x))**(sympy.S(8)/3)*(33*A + 24*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(88*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_159():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b**2*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_160():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_161():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*sin(c + d*x)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_162():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*b*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_163():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*b**2*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)) + (12*A + 21*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_164():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)/(10*b**3*d) - (b*cos(c + d*x))**(sympy.S(7)/3)*(30*A + 21*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(70*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_165():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*b**2*d) - (b*cos(c + d*x))**(sympy.S(4)/3)*(21*A + 12*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(28*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_166():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*b*d) - (b*cos(c + d*x))**(sympy.S(1)/3)*(12*A + 3*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(4*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_167():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*sin(c + d*x)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)) + (b*cos(c + d*x))**(sympy.S(4)/3)*(3*A - 6*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_168():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)) - (b*cos(c + d*x))**(sympy.S(1)/3)*(6*A + 15*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_169():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*b**2*sin(c + d*x)/(8*d*(b*cos(c + d*x))**(sympy.S(8)/3)) + (15*A + 24*C)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(16*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_170():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b**3*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_171():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b**2*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_172():
    f = (A + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_173():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_174():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*b*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)) + (12*A + 21*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(7*b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_175():
    f = (A + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*b**2*sin(c + d*x)/(10*d*(b*cos(c + d*x))**(sympy.S(10)/3)) + (21*A + 30*C)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(40*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_176():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = 3*C*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)/(d*(3*m + 10)) - 3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*(A*(3*m + 10) + C*(3*m + 7))*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*(3*m + 10)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_177():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(3*m + 8)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A*(3*m + 8) + 3*C*(3*m + 5))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(3*m + 5)*(3*m + 8)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_178():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = 3*C*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(3*m + 7)) - (b*cos(c + d*x))**(sympy.S(1)/3)*(3*A*(3*m + 7) + 3*C*(3*m + 4))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(3*m + 4)*(3*m + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_179():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)) - (3*A*(3*m + 5) + 3*C*(3*m + 2))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*(3*m + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_180():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 4)) - (3*A*(3*m + 4) + 9*C*m + 3*C)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/6), (m/2 + sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 1)*(3*m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_181():
    f = (A + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*C*sin(c + d*x)*cos(c + d*x)**m/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)) - (-3*A*(3*m + 2) + 3*C*(1 - 3*m))*sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2 + sympy.S(-1)/6), (m/2 + sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(1 - 3*m)*(3*m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_182():
    f = (a*cos(c + d*x))**m*(b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)
    F = C*(a*cos(c + d*x))**(m + 1)*(b*cos(c + d*x))**n*sin(c + d*x)/(a*d*(m + n + 2)) - (a*cos(c + d*x))**(m + 1)*(b*cos(c + d*x))**n*(A*(m + n + 2) + C*(m + n + 1))*sin(c + d*x)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(a*d*(m + n + 1)*(m + n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_183():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = C*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)/(b**3*d*(n + 4)) - (b*cos(c + d*x))**(n + 3)*(A*(n + 4) + C*(n + 3))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*(n + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_184():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*cos(c + d*x)
    F = C*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)/(b**2*d*(n + 3)) - (b*cos(c + d*x))**(n + 2)*(A*(n + 3) + C*(n + 2))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*(n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_185():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)
    F = C*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)/(b*d*(n + 2)) - (b*cos(c + d*x))**(n + 1)*(A*(n + 2) + C*(n + 1))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_186():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*sec(c + d*x)
    F = C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(n + 1)) - (b*cos(c + d*x))**n*(A*n + A + C*n)*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_187():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = C*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)/(d*n) - b*(b*cos(c + d*x))**(n - 1)*(-A*n + C*(1 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*n*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_188():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = -C*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)/(d*(1 - n)) + b**2*(b*cos(c + d*x))**(n - 2)*(A*(1 - n) + C*(2 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(1 - n)*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_189():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = -C*b**3*(b*cos(c + d*x))**(n - 3)*sin(c + d*x)/(d*(2 - n)) + b**3*(b*cos(c + d*x))**(n - 3)*(A*(2 - n) + C*(3 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/2), (n/2 + sympy.S(-1)/2,), cos(c + d*x)**2)/(d*(2 - n)*(3 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_190():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(d*(2*n + 9)) - (b*cos(c + d*x))**n*(2*A*(2*n + 9) + 2*C*(2*n + 7))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*(2*n + 9)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_191():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(d*(2*n + 7)) - (b*cos(c + d*x))**n*(2*A*(2*n + 7) + 2*C*(2*n + 5))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*(2*n + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_192():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 5)) - (b*cos(c + d*x))**n*(2*A*(2*n + 5) + 2*C*(2*n + 3))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*(2*n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_193():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(2*n + 3)) - (b*cos(c + d*x))**n*(2*A*(2*n + 3) + 4*C*n + 2*C)*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_194():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(2*n + 1)*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**n*(4*A*n + 2*A - 2*C*(1 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 4*n**2)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_195():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = -2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(1 - 2*n)*cos(c + d*x)**(sympy.S(3)/2)) + (b*cos(c + d*x))**n*(-4*A*n + 2*A + 2*C*(3 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_196():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = -2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(3 - 2*n)*cos(c + d*x)**(sympy.S(5)/2)) + (b*cos(c + d*x))**n*(2*A*(3 - 2*n) + 2*C*(5 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_197():
    f = (b*cos(c + d*x))**n*(A + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = -2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(5 - 2*n)*cos(c + d*x)**(sympy.S(7)/2)) + (b*cos(c + d*x))**n*(2*A*(5 - 2*n) + 2*C*(7 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-7)/4), (n/2 + sympy.S(-3)/4,), cos(c + d*x)**2)/(d*(5 - 2*n)*(7 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_198():
    f = (A + C*cos(e + f*x)**2)*(a*cos(e + f*x) + a)**m
    F = 2**(m + sympy.S.Half)*(A*(m**2 + 3*m + 2) + C*(m**2 + m + 1))*(a*cos(e + f*x) + a)**m*(cos(e + f*x) + 1)**(-m + sympy.S(-1)/2)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - cos(e + f*x)/2)/(f*(m + 1)*(m + 2)) - C*(a*cos(e + f*x) + a)**m*sin(e + f*x)/(f*(m**2 + 3*m + 2)) + C*(a*cos(e + f*x) + a)**(m + 1)*sin(e + f*x)/(a*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_199():
    f = (A + C*cos(c + d*x)**2)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)
    F = -9*C*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)/(40*d) + 3*C*(a*cos(c + d*x) + a)**(sympy.S(5)/3)*sin(c + d*x)/(8*a*d) + 2**(sympy.S(1)/6)*(40*A + 19*C)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(20*d*(cos(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_200():
    f = (A + C*cos(c + d*x)**2)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)
    F = -9*C*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)/(28*d) + 3*C*(a*cos(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)/(7*a*d) + 2**(sympy.S(5)/6)*(28*A + 13*C)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(28*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_201():
    f = (A + C*cos(c + d*x)**2)/(a*cos(c + d*x) + a)**(sympy.S(1)/3)
    F = -9*C*sin(c + d*x)/(10*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)) + 3*C*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)/(5*a*d) + 2**(sympy.S(1)/6)*(10*A + 7*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(10*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*(cos(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_202():
    f = (A + C*cos(c + d*x)**2)/(a*cos(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*C*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)/(4*a*d) + (3*A + 3*C)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(5)/6)*(4*A + 7*C)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(4*a*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_203():
    f = (A + C*cos(c + d*x)**2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*sqrt(2)*C*a*(a + b)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + 3*C*(a + b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b*d) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(3*C*a**2 + b**2*(8*A + 5*C))*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_204():
    f = (A + C*cos(c + d*x)**2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*sqrt(2)*C*a*(a + b)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + 3*C*(a + b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*b*d) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(3*C*a**2 + b**2*(7*A + 4*C))*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_205():
    f = (A + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*sqrt(2)*C*a*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + 3*C*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b*d) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*(3*C*a**2 + b**2*(5*A + 2*C))*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_206():
    f = (A + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*sqrt(2)*C*a*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(4*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + 3*C*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*b*d) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*(3*C*a**2 + b**2*(4*A + C))*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - cos(c + d*x)/2, b*(1 - cos(c + d*x))/(a + b))/(4*b**2*d*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_207():
    f = (a + b*cos(e + f*x))**m*(-A*cos(e + f*x)**2 + A)
    F = -4*sqrt(2)*A*(a + b*cos(e + f*x))**m*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(f*((a + b*cos(e + f*x))/(a + b))**m*sqrt(cos(e + f*x) + 1)) + 4*sqrt(2)*A*(a + b*cos(e + f*x))**m*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(f*((a + b*cos(e + f*x))/(a + b))**m*sqrt(cos(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_208():
    f = (A + C*cos(e + f*x)**2)*(a + b*cos(e + f*x))**m
    F = -sqrt(2)*C*a*(a + b)*(a + b*cos(e + f*x))**m*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1)) + C*(a + b*cos(e + f*x))**(m + 1)*sin(e + f*x)/(b*f*(m + 2)) + sqrt(2)*(a + b*cos(e + f*x))**m*(C*a**2 + b**2*(A*(m + 2) + C*(m + 1)))*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_209():
    f = (a*cos(e + f*x))**m*(B*cos(e + f*x) + C*cos(e + f*x)**2)
    F = -B*(a*cos(e + f*x))**(m + 2)*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(a**2*f*(m + 2)*sqrt(sin(e + f*x)**2)) - C*(a*cos(e + f*x))**(m + 3)*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), cos(e + f*x)**2)/(a**3*f*(m + 3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_210():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*sqrt(sin(c + d*x)**2)) - 3*C*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(5)/3), (m/2 + sympy.S(8)/3,), cos(c + d*x)**2)/(d*(3*m + 10)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_211():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(4)/3), (m/2 + sympy.S(7)/3,), cos(c + d*x)**2)/(d*(3*m + 8)*sqrt(sin(c + d*x)**2)) - 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(11)/6), (m/2 + sympy.S(17)/6,), cos(c + d*x)**2)/(d*(3*m + 11)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_212():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(5)/3), (m/2 + sympy.S(8)/3,), cos(c + d*x)**2)/(d*(3*m + 10)*sqrt(sin(c + d*x)**2)) - 3*C*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 4)*hyper((sympy.S.Half, m/2 + sympy.S(13)/6), (m/2 + sympy.S(19)/6,), cos(c + d*x)**2)/(d*(3*m + 13)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_213():
    f = (B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)*sqrt(sin(c + d*x)**2)) - 3*C*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(4)/3), (m/2 + sympy.S(7)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 8)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_214():
    f = (B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 4)*sqrt(sin(c + d*x)**2)) - 3*C*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_215():
    f = (B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*sqrt(sin(c + d*x)**2)) - 3*C*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_216():
    f = (a*cos(c + d*x))**m*(b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -B*(a*cos(c + d*x))**(m + 2)*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, m/2 + n/2 + 1), (m/2 + n/2 + 2,), cos(c + d*x)**2)/(a**2*d*(m + n + 2)*sqrt(sin(c + d*x)**2)) - C*(a*cos(c + d*x))**(m + 3)*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S(3)/2), (m/2 + n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(a**3*d*(m + n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_217():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = -B*(b*cos(c + d*x))**(n + 4)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 2), (n/2 + 3,), cos(c + d*x)**2)/(b**4*d*(n + 4)*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**(n + 5)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(5)/2), (n/2 + sympy.S(7)/2,), cos(c + d*x)**2)/(b**5*d*(n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_218():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = -B*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**(n + 4)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 2), (n/2 + 3,), cos(c + d*x)**2)/(b**4*d*(n + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_219():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -B*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_220():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = -B*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_221():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = -B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_222():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = B*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2)) - C*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_223():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = B*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2)) + C*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_224():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)*hyper((sympy.S.Half, n/2 + sympy.S(9)/4), (n/2 + sympy.S(13)/4,), cos(c + d*x)**2)/(d*(2*n + 9)*sqrt(sin(c + d*x)**2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(11)/2)*hyper((sympy.S.Half, n/2 + sympy.S(11)/4), (n/2 + sympy.S(15)/4,), cos(c + d*x)**2)/(d*(2*n + 11)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_225():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)*hyper((sympy.S.Half, n/2 + sympy.S(9)/4), (n/2 + sympy.S(13)/4,), cos(c + d*x)**2)/(d*(2*n + 9)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_226():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_227():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_228():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_229():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x))) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_230():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_231():
    f = (b*cos(c + d*x))**n*(B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_232():
    f = (B*cos(e + f*x) + C*cos(e + f*x)**2)*(a*cos(e + f*x) + a)**m
    F = 2**(m + sympy.S.Half)*(a*cos(e + f*x) + a)**m*(B*m*(m + 2) + C*(m**2 + m + 1))*(cos(e + f*x) + 1)**(-m + sympy.S(-1)/2)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - cos(e + f*x)/2)/(f*(m + 1)*(m + 2)) + C*(a*cos(e + f*x) + a)**(m + 1)*sin(e + f*x)/(a*f*(m + 2)) - (-B*(m + 2) + C)*(a*cos(e + f*x) + a)**m*sin(e + f*x)/(f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_233():
    f = (a + b*cos(e + f*x))**m*(B*cos(e + f*x) + C*cos(e + f*x)**2)
    F = C*(a + b*cos(e + f*x))**(m + 1)*sin(e + f*x)/(b*f*(m + 2)) - sqrt(2)*(a + b)*(a + b*cos(e + f*x))**m*(-B*b*(m + 2) + C*a)*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1)) + sqrt(2)*(a + b*cos(e + f*x))**m*(-B*a*b*(m + 2) + C*a**2 + C*b**2*(m + 1))*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_234():
    f = (a + b*cos(c + d*x))**(sympy.S(2)/3)*(B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b*d) + sqrt(2)*(a + b)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(8*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) - sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(8*B*a*b - 3*C*a**2 - 5*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_235():
    f = (a + b*cos(c + d*x))**(sympy.S(1)/3)*(B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*b*d) + sqrt(2)*(a + b)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(7*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) - sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(7*B*a*b - 3*C*a**2 - 4*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_236():
    f = (B*cos(c + d*x) + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b*d) - sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*(5*B*a*b - 3*C*a**2 - 2*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(5*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_237():
    f = (B*cos(c + d*x) + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*b*d) - sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*(4*B*a*b - 3*C*a**2 - C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - cos(c + d*x)/2, b*(1 - cos(c + d*x))/(a + b))/(4*b**2*d*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(4*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(4*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_238():
    f = (a*cos(e + f*x))**m*(A + B*cos(e + f*x) + C*cos(e + f*x)**2)
    F = -B*(a*cos(e + f*x))**(m + 2)*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(a**2*f*(m + 2)*sqrt(sin(e + f*x)**2)) + C*(a*cos(e + f*x))**(m + 1)*sin(e + f*x)/(a*f*(m + 2)) - (a*cos(e + f*x))**(m + 1)*(A*(m + 2) + C*(m + 1))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*(m + 1)*(m + 2)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_239():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = 10*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**2*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**3*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_240():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**2*d) + 2*b*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_241():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_242():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_243():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*A*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) + 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_244():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*b*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_245():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*b*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_246():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*b*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*b**2*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_247():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = 10*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**2*d) + 2*b*sqrt(b*cos(c + d*x))*(9*A + 7*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_248():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 6*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) + 2*b**2*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 2*b*sqrt(b*cos(c + d*x))*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_249():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*b*sqrt(b*cos(c + d*x))*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_250():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b**2*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_251():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) + 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) - 2*b*sqrt(b*cos(c + d*x))*(A - C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_252():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*b**2*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_253():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*b**2*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 2*b*sqrt(b*cos(c + d*x))*(3*A + 5*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_254():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**6
    F = 2*A*b**5*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*b**4*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*b**2*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*b**3*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**2*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_255():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 10*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d) + 2*b**2*sqrt(b*cos(c + d*x))*(9*A + 7*C)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 2*b*(b*cos(c + d*x))**(sympy.S(3)/2)*(9*A + 7*C)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_256():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = 6*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + 2*b**3*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 2*b**2*sqrt(b*cos(c + d*x))*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_257():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*C*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*b**2*sqrt(b*cos(c + d*x))*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_258():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 2*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 2*b**3*(3*A + C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_259():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) + 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*(A - C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_260():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*b**3*(A + 3*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_261():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**6
    F = 2*A*b**5*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*b**3*(3*A + 5*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*(3*A + 5*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_262():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**7
    F = 2*A*b**6*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*b**5*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*b**3*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*b**4*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**3*(5*A + 7*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_263():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**3*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**4*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_264():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**3*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_265():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_266():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(b*cos(c + d*x))
    F = 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x))) + 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_267():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_268():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*A*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x))) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_269():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 2*A*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + (6*A + 10*C)*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_270():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4/sqrt(b*cos(c + d*x))
    F = 2*A*b**3*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*b*(5*A + 7*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_271():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**2*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**4*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**5*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_272():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**4*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_273():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d) + 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_274():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x))) + 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_275():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_276():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x))) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_277():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + (6*A + 10*C)*sin(c + d*x)/(5*b*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_278():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*b**2*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*sin(c + d*x)/(5*b*d*sqrt(b*cos(c + d*x))) - 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + (10*A + 14*C)*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_279():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**5/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**3*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**5*d) + 2*C*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**6*d) + sqrt(b*cos(c + d*x))*(18*A + 14*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**3*d*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**(sympy.S(3)/2)*(18*A + 14*C)*sin(c + d*x)/(45*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_280():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d) + 2*C*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**5*d) + (14*A + 10*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x))) + sqrt(b*cos(c + d*x))*(14*A + 10*C)*sin(c + d*x)/(21*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_281():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d) + 2*C*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d) + sqrt(b*cos(c + d*x))*(10*A + 6*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_282():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x))) + 2*C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d) + (6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_283():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_284():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x))) + (2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_285():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + (6*A + 10*C)*sin(c + d*x)/(5*b**2*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_286():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*b*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 2*B*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*B*sin(c + d*x)/(5*b**2*d*sqrt(b*cos(c + d*x))) - 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + (10*A + 14*C)*sin(c + d*x)/(21*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + (10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_287():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*A*sin(c + d*x)/(5*b*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 2*B*sin(c + d*x)/(3*b**2*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d*sqrt(b*cos(c + d*x))) + (6*A + 10*C)*sin(c + d*x)/(5*b**3*d*sqrt(b*cos(c + d*x))) - sqrt(b*cos(c + d*x))*(6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_288():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 3*B*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d) - sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)**3/(15*d*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_289():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = -B*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_290():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = B*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(b*cos(c + d*x))*(3*A + 2*C)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_291():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = A*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_292():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = A*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_293():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_294():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_295():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_296():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_297():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 3*B*b*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d) - b*sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)**3/(15*d*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_298():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = -B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + b*x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_299():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = B*b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d) + b*sqrt(b*cos(c + d*x))*(3*A + 2*C)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_300():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = A*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_301():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_302():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_303():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + b*sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_304():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_305():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(13)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + b*sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b*sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_306():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = 3*B*b**2*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d) - b**2*sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)**3/(15*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(5*A + 4*C)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_307():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = -B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + b**2*x*sqrt(b*cos(c + d*x))*(4*A + 3*C)/(8*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(4*A + 3*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_308():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = B*b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d) + b**2*sqrt(b*cos(c + d*x))*(3*A + 2*C)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_309():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + C*b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_310():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + C*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_311():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + C*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_312():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + b**2*sqrt(b*cos(c + d*x))*(A + 2*C)*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_313():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(13)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*(2*A + 3*C)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_314():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(15)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + b**2*sqrt(b*cos(c + d*x))*(3*A + 4*C)*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b**2*sqrt(b*cos(c + d*x))*(3*A + 4*C)*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_315():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/sqrt(b*cos(c + d*x))
    F = -B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_316():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/sqrt(b*cos(c + d*x))
    F = B*x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*sqrt(b*cos(c + d*x))) + (3*A + 2*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_317():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    F = A*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x)) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_318():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x)) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_319():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_320():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_321():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = A*sin(c + d*x)/(3*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x))) + (2*A + 3*C)*sin(c + d*x)/(3*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_322():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(9)/2))
    F = A*sin(c + d*x)/(4*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + B*sin(c + d*x)**3/(3*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (3*A + 4*C)*sin(c + d*x)/(8*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_323():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = -B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*b*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_324():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = B*x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*b*d*sqrt(b*cos(c + d*x))) + (3*A + 2*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_325():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_326():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_327():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_328():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_329():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(3*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x))) + (2*A + 3*C)*sin(c + d*x)/(3*b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_330():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = A*sin(c + d*x)/(4*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + B*sin(c + d*x)**3/(3*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (3*A + 4*C)*sin(c + d*x)/(8*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_331():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(9)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = -B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b**2*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b**2*d*sqrt(b*cos(c + d*x))) + x*(4*A + 3*C)*sqrt(cos(c + d*x))/(8*b**2*sqrt(b*cos(c + d*x))) + (4*A + 3*C)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_332():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = B*x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + (3*A + 2*C)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_333():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_334():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x))) + C*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_335():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + C*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_336():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (A + 2*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_337():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(3*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x))) + (2*A + 3*C)*sin(c + d*x)/(3*b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_338():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(4*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + B*sin(c + d*x)**3/(3*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + B*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + (3*A + 4*C)*sin(c + d*x)/(8*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + (3*A + 4*C)*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_339():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(11)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(17)/6,), cos(c + d*x)**2)/(11*b**3*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)/(11*b**2*d) - (b*cos(c + d*x))**(sympy.S(8)/3)*(33*A + 24*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(88*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_340():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_341():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_342():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = 3*A*b*sin(c + d*x)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)) - 3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*sqrt(sin(c + d*x)**2)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_343():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) + 3*B*b*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_344():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 3*A*b**3*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)) + 3*B*b**2*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) + 3*b*(4*A + 7*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_345():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*b**3*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)/(13*b**2*d) - (b*cos(c + d*x))**(sympy.S(10)/3)*(39*A + 30*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(130*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_346():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)/(10*b*d) - (b*cos(c + d*x))**(sympy.S(7)/3)*(30*A + 21*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(70*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_347():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*d) - (b*cos(c + d*x))**(sympy.S(4)/3)*(21*A + 12*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(28*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_348():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = -3*B*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2)) + 3*C*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*d) - 3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*(4*A + C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_349():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)) - 3*B*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2)) + (b*cos(c + d*x))**(sympy.S(4)/3)*(3*A - 6*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(8*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_350():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = 3*A*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)) + 3*B*b**2*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) - 3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*(2*A + 5*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(5*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_351():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(11)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(17)/6,), cos(c + d*x)**2)/(11*b**4*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)/(11*b**3*d) - (b*cos(c + d*x))**(sympy.S(8)/3)*(33*A + 24*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(88*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_352():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**3*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b**2*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_353():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**2*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_354():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*sin(c + d*x)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)) - 3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b*d*sqrt(sin(c + d*x)**2)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_355():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*b*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) + 3*B*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_356():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*A*b**2*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)) + 3*B*b*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) + (12*A + 21*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_357():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(11)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(17)/6,), cos(c + d*x)**2)/(11*b**5*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)/(11*b**4*d) - (b*cos(c + d*x))**(sympy.S(8)/3)*(33*A + 24*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(88*b**4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_358():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**4*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b**3*d) - (b*cos(c + d*x))**(sympy.S(5)/3)*(24*A + 15*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(40*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_359():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*B*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b**2*d) - (b*cos(c + d*x))**(sympy.S(2)/3)*(15*A + 6*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_360():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)) - 3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2)) + (b*cos(c + d*x))**(sympy.S(5)/3)*(6*A - 3*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_361():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)) + 3*B*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A + 12*C)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_362():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*b*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)) + 3*B*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) + (12*A + 21*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(7*b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_363():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(5)/3), (m/2 + sympy.S(8)/3,), cos(c + d*x)**2)/(d*(3*m + 10)*sqrt(sin(c + d*x)**2)) + 3*C*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)/(d*(3*m + 10)) - 3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*(A*(3*m + 10) + C*(3*m + 7))*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*(3*m + 10)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_364():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(4)/3), (m/2 + sympy.S(7)/3,), cos(c + d*x)**2)/(d*(3*m + 8)*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(3*m + 8)) - (b*cos(c + d*x))**(sympy.S(2)/3)*(3*A*(3*m + 8) + 3*C*(3*m + 5))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(3*m + 5)*(3*m + 8)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_365():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m
    F = -3*B*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*sqrt(sin(c + d*x)**2)) + 3*C*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(3*m + 7)) - (b*cos(c + d*x))**(sympy.S(1)/3)*(3*A*(3*m + 7) + 3*C*(3*m + 4))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(3*m + 4)*(3*m + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_366():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)*sqrt(sin(c + d*x)**2)) + 3*C*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)) - (3*A*(3*m + 5) + 3*C*(3*m + 2))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*(3*m + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_367():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 4)*sqrt(sin(c + d*x)**2)) + 3*C*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 4)) - (3*A*(3*m + 4) + 9*C*m + 3*C)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/6), (m/2 + sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 1)*(3*m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_368():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*B*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*sqrt(sin(c + d*x)**2)) + 3*C*sin(c + d*x)*cos(c + d*x)**m/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)) - (-3*A*(3*m + 2) + 3*C*(1 - 3*m))*sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2 + sympy.S(-1)/6), (m/2 + sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(1 - 3*m)*(3*m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_369():
    f = (a*cos(c + d*x))**m*(b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -B*(a*cos(c + d*x))**(m + 2)*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, m/2 + n/2 + 1), (m/2 + n/2 + 2,), cos(c + d*x)**2)/(a**2*d*(m + n + 2)*sqrt(sin(c + d*x)**2)) + C*(a*cos(c + d*x))**(m + 1)*(b*cos(c + d*x))**n*sin(c + d*x)/(a*d*(m + n + 2)) - (a*cos(c + d*x))**(m + 1)*(b*cos(c + d*x))**n*(A*(m + n + 2) + C*(m + n + 1))*sin(c + d*x)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(a*d*(m + n + 1)*(m + n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_370():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**2
    F = -B*(b*cos(c + d*x))**(n + 4)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 2), (n/2 + 3,), cos(c + d*x)**2)/(b**4*d*(n + 4)*sqrt(sin(c + d*x)**2)) + C*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)/(b**3*d*(n + 4)) - (b*cos(c + d*x))**(n + 3)*(A*(n + 4) + C*(n + 3))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*(n + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_371():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)
    F = -B*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2)) + C*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)/(b**2*d*(n + 3)) - (b*cos(c + d*x))**(n + 2)*(A*(n + 3) + C*(n + 2))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*(n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_372():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = -B*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2)) + C*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)/(b*d*(n + 2)) - (b*cos(c + d*x))**(n + 1)*(A*(n + 2) + C*(n + 1))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_373():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)
    F = -B*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2)) + C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(n + 1)) - (b*cos(c + d*x))**n*(A*n + A + C*n)*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_374():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**2
    F = -B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2)) + C*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)/(d*n) - b*(b*cos(c + d*x))**(n - 1)*(-A*n + C*(1 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*n*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_375():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**3
    F = B*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2)) - C*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)/(d*(1 - n)) + b**2*(b*cos(c + d*x))**(n - 2)*(A*(1 - n) + C*(2 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(1 - n)*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_376():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sec(c + d*x)**4
    F = B*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2)) - C*b**3*(b*cos(c + d*x))**(n - 3)*sin(c + d*x)/(d*(2 - n)) + b**3*(b*cos(c + d*x))**(n - 3)*(A*(2 - n) + C*(3 - n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/2), (n/2 + sympy.S(-1)/2,), cos(c + d*x)**2)/(d*(2 - n)*(3 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_377():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(d*(2*n + 7)) - (b*cos(c + d*x))**n*(2*A*(2*n + 7) + 2*C*(2*n + 5))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*(2*n + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_378():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)*sqrt(cos(c + d*x))
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 5)) - (b*cos(c + d*x))**n*(2*A*(2*n + 5) + 2*C*(2*n + 3))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*(2*n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_379():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/sqrt(cos(c + d*x))
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(2*n + 3)) - (b*cos(c + d*x))**n*(2*A*(2*n + 3) + 4*C*n + 2*C)*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_380():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(3)/2)
    F = -2*B*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2)) + 2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(2*n + 1)*sqrt(cos(c + d*x))) + (b*cos(c + d*x))**n*(4*A*n + 2*A - 2*C*(1 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 4*n**2)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_381():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x))) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(1 - 2*n)*cos(c + d*x)**(sympy.S(3)/2)) + (b*cos(c + d*x))**n*(-4*A*n + 2*A + 2*C*(3 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_382():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)) - 2*C*(b*cos(c + d*x))**n*sin(c + d*x)/(d*(3 - 2*n)*cos(c + d*x)**(sympy.S(5)/2)) + (b*cos(c + d*x))**n*(2*A*(3 - 2*n) + 2*C*(5 - 2*n))*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_383():
    f = (a*cos(e + f*x) + a)**m*(A + B*cos(e + f*x) + C*cos(e + f*x)**2)
    F = 2**(m + sympy.S.Half)*(a*cos(e + f*x) + a)**m*(cos(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(A*(m**2 + 3*m + 2) + B*m*(m + 2) + C*(m**2 + m + 1))*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - cos(e + f*x)/2)/(f*(m + 1)*(m + 2)) + C*(a*cos(e + f*x) + a)**(m + 1)*sin(e + f*x)/(a*f*(m + 2)) - (-B*(m + 2) + C)*(a*cos(e + f*x) + a)**m*sin(e + f*x)/(f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_384():
    f = (a*cos(c + d*x) + a)**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a*cos(c + d*x) + a)**(sympy.S(5)/3)*sin(c + d*x)/(8*a*d) + (24*B - 9*C)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)/(40*d) + 2**(sympy.S(1)/6)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*(40*A + 16*B + 19*C)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(20*d*(cos(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_385():
    f = (a*cos(c + d*x) + a)**(sympy.S(1)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a*cos(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)/(7*a*d) + (21*B - 9*C)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)/(28*d) + 2**(sympy.S(5)/6)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*(28*A + 7*B + 13*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(28*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_386():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(a*cos(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*C*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)/(5*a*d) + (15*B - 9*C)*sin(c + d*x)/(10*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/6)*(10*A - 5*B + 7*C)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(10*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*(cos(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_387():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(a*cos(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*C*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)/(4*a*d) + (3*A - 3*B + 3*C)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(5)/6)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*(4*A - 8*B + 7*C)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(4*a*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_388():
    f = (a + b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)/(8*b*d) + sqrt(2)*(a + b)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(8*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(8*A*b**2 - 8*B*a*b + 3*C*a**2 + 5*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(8*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_389():
    f = (a + b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x) + C*cos(c + d*x)**2)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)/(7*b*d) + sqrt(2)*(a + b)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(7*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(7*A*b**2 - 7*B*a*b + 3*C*a**2 + 4*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(7*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_390():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)/(5*b*d) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*(5*A*b**2 - 5*B*a*b + 3*C*a**2 + 2*C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(5*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(5*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_391():
    f = (A + B*cos(c + d*x) + C*cos(c + d*x)**2)/(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*C*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)/(4*b*d) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*(4*A*b**2 - 4*B*a*b + 3*C*a**2 + C*b**2)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - cos(c + d*x)/2, b*(1 - cos(c + d*x))/(a + b))/(4*b**2*d*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(4*B*b - 3*C*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(4*b**2*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_392():
    f = (a + b*cos(e + f*x))**m*(A + C*cos(e + f*x)**2 + (A + C)*cos(e + f*x))
    F = 4*sqrt(2)*C*(a + b*cos(e + f*x))**m*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(f*((a + b*cos(e + f*x))/(a + b))**m*sqrt(cos(e + f*x) + 1)) + 2*sqrt(2)*(A - C)*(a + b*cos(e + f*x))**m*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(f*((a + b*cos(e + f*x))/(a + b))**m*sqrt(cos(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_4_1_a_plus_b_cos_pow_m_A_plus_B_cos_plus_C_cos_pow_2_393():
    f = (a + b*cos(e + f*x))**m*(A + B*cos(e + f*x) + C*cos(e + f*x)**2)
    F = C*(a + b*cos(e + f*x))**(m + 1)*sin(e + f*x)/(b*f*(m + 2)) - sqrt(2)*(a + b)*(a + b*cos(e + f*x))**m*(-B*b*(m + 2) + C*a)*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1)) + sqrt(2)*(a + b*cos(e + f*x))**m*(A*b**2*(m + 2) - B*a*b*(m + 2) + C*a**2 + C*b**2*(m + 1))*sin(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - cos(e + f*x)/2, b*(1 - cos(e + f*x))/(a + b))/(b**2*f*((a + b*cos(e + f*x))/(a + b))**m*(m + 2)*sqrt(cos(e + f*x) + 1))
    assert integrate(f, x) == F

