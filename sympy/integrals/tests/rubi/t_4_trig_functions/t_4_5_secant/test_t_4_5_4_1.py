"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.4.1 (a+b sec)^m (A+B sec+C sec^2).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, b, c, d, e, f, m = symbols('A B C b c d e f m')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_1():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)**6
    F = C*tan(c + d*x)*sec(c + d*x)**6/(7*d) + (7*A + 6*C)*tan(c + d*x)**5/(35*d) + (7*A + 6*C)*tan(c + d*x)/(7*d) + (14*A + 12*C)*tan(c + d*x)**3/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_2():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)**5
    F = C*tan(c + d*x)*sec(c + d*x)**5/(6*d) + (6*A + 5*C)*tan(c + d*x)*sec(c + d*x)**3/(24*d) + (6*A + 5*C)*tan(c + d*x)*sec(c + d*x)/(16*d) + (6*A + 5*C)*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_3():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)**4
    F = C*tan(c + d*x)*sec(c + d*x)**4/(5*d) + (5*A + 4*C)*tan(c + d*x)**3/(15*d) + (5*A + 4*C)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_4():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)**3
    F = C*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (4*A + 3*C)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*A + 3*C)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_5():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)**2
    F = C*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (3*A + 2*C)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_6():
    f = (A + C*sec(c + d*x)**2)*sec(c + d*x)
    F = C*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*A + C)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_7():
    f = A + C*sec(c + d*x)**2
    F = A*x + C*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_8():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)
    F = A*sin(c + d*x)/d + C*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_9():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)**2
    F = A*sin(c + d*x)*cos(c + d*x)/(2*d) + x*(A/2 + C)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_10():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)**3
    F = -A*sin(c + d*x)**3/(3*d) + (A + C)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_11():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)**4
    F = A*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(3*A/8 + C/2) + (3*A + 4*C)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_12():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)**5
    F = A*sin(c + d*x)**5/(5*d) + (A + C)*sin(c + d*x)/d - (2*A + C)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_13():
    f = (A + C*sec(c + d*x)**2)*cos(c + d*x)**6
    F = A*sin(c + d*x)*cos(c + d*x)**5/(6*d) + x*(5*A/16 + 3*C/8) + (5*A + 6*C)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (5*A + 6*C)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_14():
    f = (-C*m/(m + 1) + C*sec(c + d*x)**2)*sec(c + d*x)**m
    F = C*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_15():
    f = (A - A*(m + 1)*sec(c + d*x)**2/m)*sec(c + d*x)**m
    F = -A*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_16():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*(A + C*sec(c + d*x)**2)
    F = 2*C*(b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*d) + 2*b**2*sqrt(b*sec(c + d*x))*(7*A + 5*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b*(b*sec(c + d*x))**(sympy.S(3)/2)*(7*A + 5*C)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_17():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*(A + C*sec(c + d*x)**2)
    F = 2*C*(b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) - 2*b**2*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*b*sqrt(b*sec(c + d*x))*(5*A + 3*C)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_18():
    f = sqrt(b*sec(c + d*x))*(A + C*sec(c + d*x)**2)
    F = 2*C*sqrt(b*sec(c + d*x))*tan(c + d*x)/(3*d) + sqrt(b*sec(c + d*x))*(6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_19():
    f = (A + C*sec(c + d*x)**2)/sqrt(b*sec(c + d*x))
    F = 2*C*tan(c + d*x)/(d*sqrt(b*sec(c + d*x))) + (2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_20():
    f = (A + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*A*tan(c + d*x)/(3*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + sqrt(b*sec(c + d*x))*(2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_21():
    f = (A + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*A*tan(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + (6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_22():
    f = (A + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(7)/2)
    F = 2*A*tan(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + (10*A + 14*C)*sin(c + d*x)/(21*b**3*d*sqrt(b*sec(c + d*x))) + sqrt(b*sec(c + d*x))*(10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_23():
    f = (A + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(9)/2)
    F = 2*A*tan(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(9)/2)) + (14*A + 18*C)*sin(c + d*x)/(45*b**3*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + (14*A + 18*C)*elliptic_e(c/2 + d*x/2, 2)/(15*b**4*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_24():
    f = (3*sec(c + d*x)**2 + 3)/sqrt(sec(c + d*x))
    F = 6*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_25():
    f = (m - (m + 1)*sec(e + f*x)**2)*sec(e + f*x)**m
    F = -sin(e + f*x)*sec(e + f*x)**(m + 1)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_26():
    f = (5 - 6*sec(e + f*x)**2)*sec(e + f*x)**5
    F = -tan(e + f*x)*sec(e + f*x)**5/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_27():
    f = (4 - 5*sec(e + f*x)**2)*sec(e + f*x)**4
    F = -tan(e + f*x)*sec(e + f*x)**4/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_28():
    f = (3 - 4*sec(e + f*x)**2)*sec(e + f*x)**3
    F = -tan(e + f*x)*sec(e + f*x)**3/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_29():
    f = (2 - 3*sec(e + f*x)**2)*sec(e + f*x)**2
    F = -tan(e + f*x)*sec(e + f*x)**2/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_30():
    f = (1 - 2*sec(e + f*x)**2)*sec(e + f*x)
    F = -tan(e + f*x)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_31():
    f = -sec(e + f*x)**2
    F = -tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_32():
    f = -cos(e + f*x)
    F = -sin(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_33():
    f = (sec(e + f*x)**2 - 2)*cos(e + f*x)**2
    F = -sin(e + f*x)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_34():
    f = (2*sec(e + f*x)**2 - 3)*cos(e + f*x)**3
    F = -sin(e + f*x)*cos(e + f*x)**2/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_35():
    f = (3*sec(e + f*x)**2 - 4)*cos(e + f*x)**4
    F = -sin(e + f*x)*cos(e + f*x)**3/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_36():
    f = (4*sec(e + f*x)**2 - 5)*cos(e + f*x)**5
    F = -sin(e + f*x)*cos(e + f*x)**4/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_37():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)**3
    F = B*tan(c + d*x)**3/(3*d) + B*tan(c + d*x)/d + C*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*C*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*C*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_38():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)**2
    F = B*tan(c + d*x)*sec(c + d*x)/(2*d) + B*atanh(sin(c + d*x))/(2*d) + C*tan(c + d*x)**3/(3*d) + C*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_39():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)
    F = B*tan(c + d*x)/d + C*tan(c + d*x)*sec(c + d*x)/(2*d) + C*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_40():
    f = B*sec(c + d*x) + C*sec(c + d*x)**2
    F = B*atanh(sin(c + d*x))/d + C*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_41():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)
    F = B*x + C*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_42():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**2
    F = B*sin(c + d*x)/d + C*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_43():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**3
    F = B*x/2 + B*sin(c + d*x)*cos(c + d*x)/(2*d) + C*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_44():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**4
    F = -B*sin(c + d*x)**3/(3*d) + B*sin(c + d*x)/d + C*x/2 + C*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_45():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**5
    F = 3*B*x/8 + B*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*B*sin(c + d*x)*cos(c + d*x)/(8*d) - C*sin(c + d*x)**3/(3*d) + C*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_46():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**6
    F = B*sin(c + d*x)**5/(5*d) - 2*B*sin(c + d*x)**3/(3*d) + B*sin(c + d*x)/d + 3*C*x/8 + C*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*C*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_47():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*(B*sec(c + d*x) + C*sec(c + d*x)**2)
    F = 2*B*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*B*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d) - 6*C*b**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*C*b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*C*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_48():
    f = sqrt(b*sec(c + d*x))*(B*sec(c + d*x) + C*sec(c + d*x)**2)
    F = -2*B*b*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sin(c + d*x)/d + 2*C*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*C*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_49():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)/sqrt(b*sec(c + d*x))
    F = 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d) - 2*C*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*C*sqrt(b*sec(c + d*x))*sin(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_50():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*B*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*C*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_51():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*B*sin(c + d*x)/(3*b**2*d*sqrt(b*sec(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d) + 2*C*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_52():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(7)/2)
    F = 2*B*sin(c + d*x)/(5*b**2*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*B*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*C*sin(c + d*x)/(3*b**3*d*sqrt(b*sec(c + d*x))) + 2*C*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_53():
    f = (B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(9)/2)
    F = 2*B*sin(c + d*x)/(7*b**2*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*B*sin(c + d*x)/(21*b**4*d*sqrt(b*sec(c + d*x))) + 10*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**5*d) + 2*C*sin(c + d*x)/(5*b**3*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*C*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_54():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)**4
    F = B*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*B*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*B*atanh(sin(c + d*x))/(8*d) + C*tan(c + d*x)*sec(c + d*x)**4/(5*d) + (5*A + 4*C)*tan(c + d*x)**3/(15*d) + (5*A + 4*C)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_55():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)**3
    F = B*tan(c + d*x)**3/(3*d) + B*tan(c + d*x)/d + C*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (4*A + 3*C)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*A + 3*C)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_56():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)**2
    F = B*tan(c + d*x)*sec(c + d*x)/(2*d) + B*atanh(sin(c + d*x))/(2*d) + C*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (3*A + 2*C)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_57():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*sec(c + d*x)
    F = B*tan(c + d*x)/d + C*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*A + C)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_58():
    f = A + B*sec(c + d*x) + C*sec(c + d*x)**2
    F = A*x + B*atanh(sin(c + d*x))/d + C*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_59():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)
    F = A*sin(c + d*x)/d + B*x + C*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_60():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**2
    F = A*sin(c + d*x)*cos(c + d*x)/(2*d) + B*sin(c + d*x)/d + x*(A/2 + C)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_61():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**3
    F = -A*sin(c + d*x)**3/(3*d) + B*x/2 + B*sin(c + d*x)*cos(c + d*x)/(2*d) + (A + C)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_62():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**4
    F = A*sin(c + d*x)*cos(c + d*x)**3/(4*d) - B*sin(c + d*x)**3/(3*d) + B*sin(c + d*x)/d + x*(3*A/8 + C/2) + (3*A + 4*C)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_63():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**5
    F = A*sin(c + d*x)**5/(5*d) + 3*B*x/8 + B*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*B*sin(c + d*x)*cos(c + d*x)/(8*d) + (A + C)*sin(c + d*x)/d - (2*A + C)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_64():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)*cos(c + d*x)**6
    F = A*sin(c + d*x)*cos(c + d*x)**5/(6*d) + B*sin(c + d*x)**5/(5*d) - 2*B*sin(c + d*x)**3/(3*d) + B*sin(c + d*x)/d + x*(5*A/16 + 3*C/8) + (5*A + 6*C)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (5*A + 6*C)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_65():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*(A + B*sec(c + d*x) + C*sec(c + d*x)**2)
    F = 2*B*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*B*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d) + 2*C*(b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) - 2*b**2*(5*A + 3*C)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*b*sqrt(b*sec(c + d*x))*(5*A + 3*C)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_66():
    f = sqrt(b*sec(c + d*x))*(A + B*sec(c + d*x) + C*sec(c + d*x)**2)
    F = -2*B*b*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sin(c + d*x)/d + 2*C*sqrt(b*sec(c + d*x))*tan(c + d*x)/(3*d) + sqrt(b*sec(c + d*x))*(6*A + 2*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_67():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)/sqrt(b*sec(c + d*x))
    F = 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d) + 2*C*tan(c + d*x)/(d*sqrt(b*sec(c + d*x))) + (2*A - 2*C)*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_68():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*A*tan(c + d*x)/(3*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 2*B*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(b*sec(c + d*x))*(2*A + 6*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_69():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*A*tan(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 2*B*sin(c + d*x)/(3*b**2*d*sqrt(b*sec(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d) + (6*A + 10*C)*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_4_1_a_plus_b_sec_pow_m_A_plus_B_sec_plus_C_sec_pow_2_70():
    f = (A + B*sec(c + d*x) + C*sec(c + d*x)**2)/(b*sec(c + d*x))**(sympy.S(7)/2)
    F = 2*A*tan(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 2*B*sin(c + d*x)/(5*b**2*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*B*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + (10*A + 14*C)*sin(c + d*x)/(21*b**3*d*sqrt(b*sec(c + d*x))) + sqrt(b*sec(c + d*x))*(10*A + 14*C)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**4*d)
    assert integrate(f, x) == F

