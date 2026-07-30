"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.3.1 (a+b sec)^m (d sec)^n (A+B sec).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_1():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*(A + B*sec(c + d*x))
    F = 2*A*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*A*b*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d) - 6*B*b**3*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*B*b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*B*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_2():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*(A + B*sec(c + d*x))
    F = -2*A*b**2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*A*b*sqrt(b*sec(c + d*x))*sin(c + d*x)/d + 2*B*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*B*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_3():
    f = sqrt(b*sec(c + d*x))*(A + B*sec(c + d*x))
    F = 2*A*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d - 2*B*b*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_4():
    f = (A + B*sec(c + d*x))/sqrt(b*sec(c + d*x))
    F = 2*A*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_5():
    f = (A + B*sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(3*b*d*sqrt(b*sec(c + d*x))) + 2*A*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d) + 2*B*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_6():
    f = (A + B*sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(5*b*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*A*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*B*sin(c + d*x)/(3*b**2*d*sqrt(b*sec(c + d*x))) + 2*B*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_7():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))*sec(c + d*x)**2
    F = 3*A*(b*sec(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S(-4)/3, sympy.S.Half), (sympy.S(-1)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_8():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))*sec(c + d*x)
    F = 3*A*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_9():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))
    F = -3*A*b*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_10():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))*cos(c + d*x)
    F = -3*A*b**2*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) - 3*B*b*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_11():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))*cos(c + d*x)**2
    F = -3*A*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2)) - 3*B*b**2*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_12():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))*sec(c + d*x)**2
    F = 3*A*(b*sec(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S(-5)/3, sympy.S.Half), (sympy.S(-2)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_13():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))*sec(c + d*x)
    F = 3*A*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_14():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))
    F = 3*A*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_15():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))*cos(c + d*x)
    F = -3*A*b**2*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) + 3*B*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_16():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))*cos(c + d*x)**2
    F = -3*A*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2)) - 3*B*b**2*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_17():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = 3*A*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_18():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_19():
    f = (A + B*sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = -3*A*b*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_20():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_21():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = 3*A*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_22():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_23():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_24():
    f = (A + B*sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*b*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_25():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_26():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_27():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*(A + B*sec(c + d*x))*sec(c + d*x)**m
    F = 3*A*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/6), (sympy.S(5)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(3*m + 1)*sqrt(sin(c + d*x)**2)) + 3*B*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, -m/2 + sympy.S(-2)/3), (sympy.S(1)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m + 1)/(d*(3*m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_28():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*(A + B*sec(c + d*x))*sec(c + d*x)**m
    F = -3*A*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/6 - m/2), (sympy.S(7)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - 3*m)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/3), (sympy.S(2)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(3*m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_29():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*(A + B*sec(c + d*x))*sec(c + d*x)**m
    F = -3*A*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/3 - m/2), (sympy.S(4)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(2 - 3*m)*sqrt(sin(c + d*x)**2)) + 3*B*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/6), (sympy.S(5)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(3*m + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_30():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3 - m/2), (sympy.S(5)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*(4 - 3*m)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/6 - m/2), (sympy.S(7)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*(1 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_31():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6 - m/2), (sympy.S(11)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(b*sec(c + d*x))**(sympy.S(2)/3)*(5 - 3*m)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/3 - m/2), (sympy.S(4)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(b*sec(c + d*x))**(sympy.S(2)/3)*(2 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_32():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*A*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6 - m/2), (sympy.S(13)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*(7 - 3*m)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3 - m/2), (sympy.S(5)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*(4 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_33():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*sec(c + d*x)**m
    F = -A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -m/2 - n/2 + sympy.S.Half), (-m/2 - n/2 + sympy.S(3)/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(-m - n + 1)*sqrt(sin(c + d*x)**2)) + B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -m/2 - n/2), (-m/2 - n/2 + 1,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(m + n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_34():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*sec(c + d*x)**2
    F = A*(b*sec(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2)) + B*(b*sec(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, -n/2 - 1), (-n/2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_35():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*sec(c + d*x)
    F = A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2)) + B*(b*sec(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_36():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))
    F = -A*b*(b*sec(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2)) + B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_37():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*cos(c + d*x)
    F = -A*b**2*(b*sec(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2)) - B*b*(b*sec(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_38():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*cos(c + d*x)**2
    F = -A*b**3*(b*sec(c + d*x))**(n - 3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), cos(c + d*x)**2)/(d*(3 - n)*sqrt(sin(c + d*x)**2)) - B*b**2*(b*sec(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_39():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/4), (sympy.S(3)/4 - n/2,), cos(c + d*x)**2)*sqrt(sec(c + d*x))/(d*(2*n + 1)*sqrt(sin(c + d*x)**2)) + 2*B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-3)/4), (sympy.S(1)/4 - n/2,), cos(c + d*x)**2)*sec(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_40():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))*sqrt(sec(c + d*x))
    F = -2*A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/4 - n/2), (sympy.S(5)/4 - n/2,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(sec(c + d*x))) + 2*B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/4), (sympy.S(3)/4 - n/2,), cos(c + d*x)**2)*sqrt(sec(c + d*x))/(d*(2*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_41():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))/sqrt(sec(c + d*x))
    F = -2*A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/4 - n/2), (sympy.S(7)/4 - n/2,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(3)/2)) - 2*B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/4 - n/2), (sympy.S(5)/4 - n/2,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_42():
    f = (b*sec(c + d*x))**n*(A + B*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = -2*A*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/4 - n/2), (sympy.S(9)/4 - n/2,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(5)/2)) - 2*B*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/4 - n/2), (sympy.S(7)/4 - n/2,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_43():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)**4
    F = B*a*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a*(A + B)*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*(A + B)*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*(A + B)*atanh(sin(c + d*x))/(8*d) + a*(5*A + 4*B)*tan(c + d*x)**3/(15*d) + a*(5*A + 4*B)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_44():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)**3
    F = B*a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a*(A + B)*tan(c + d*x)**3/(3*d) + a*(A + B)*tan(c + d*x)/d + a*(4*A + 3*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a*(4*A + 3*B)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_45():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)**2
    F = B*a*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a*(A + B)*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(A + B)*atanh(sin(c + d*x))/(2*d) + a*(3*A + 2*B)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_46():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)
    F = B*a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(A + B)*tan(c + d*x)/d + a*(2*A + B)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_47():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)
    F = A*a*x + B*a*tan(c + d*x)/d + a*(A + B)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_48():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)
    F = A*a*sin(c + d*x)/d + B*a*atanh(sin(c + d*x))/d + a*x*(A + B)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_49():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**2
    F = A*a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*x*(A + 2*B)/2 + a*(A + B)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_50():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**3
    F = A*a*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a*x*(A + B)/2 + a*(A + B)*sin(c + d*x)*cos(c + d*x)/(2*d) + a*(2*A + 3*B)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_51():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**4
    F = A*a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a*x*(3*A + 4*B)/8 - a*(A + B)*sin(c + d*x)**3/(3*d) + a*(A + B)*sin(c + d*x)/d + a*(3*A + 4*B)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_52():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**5
    F = A*a*sin(c + d*x)*cos(c + d*x)**4/(5*d) + 3*a*x*(A + B)/8 + a*(A + B)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*(A + B)*sin(c + d*x)*cos(c + d*x)/(8*d) - a*(4*A + 5*B)*sin(c + d*x)**3/(15*d) + a*(4*A + 5*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_53():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sec(c + d*x)**3
    F = B*(a**2*sec(c + d*x) + a**2)*tan(c + d*x)*sec(c + d*x)**3/(5*d) + a**2*(5*A + 6*B)*tan(c + d*x)*sec(c + d*x)**3/(20*d) + a**2*(7*A + 6*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a**2*(7*A + 6*B)*atanh(sin(c + d*x))/(8*d) + a**2*(10*A + 9*B)*tan(c + d*x)**3/(15*d) + a**2*(10*A + 9*B)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_54():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sec(c + d*x)**2
    F = B*(a*sec(c + d*x) + a)**3*tan(c + d*x)/(4*a*d) + a**2*(8*A + 7*B)*tan(c + d*x)*sec(c + d*x)/(24*d) + a**2*(8*A + 7*B)*tan(c + d*x)/(6*d) + a**2*(8*A + 7*B)*atanh(sin(c + d*x))/(8*d) + (4*A - B)*(a*sec(c + d*x) + a)**2*tan(c + d*x)/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_55():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sec(c + d*x)
    F = B*(a*sec(c + d*x) + a)**2*tan(c + d*x)/(3*d) + a**2*(3*A + 2*B)*tan(c + d*x)*sec(c + d*x)/(6*d) + 2*a**2*(3*A + 2*B)*tan(c + d*x)/(3*d) + a**2*(3*A + 2*B)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_56():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2
    F = A*a**2*x + B*(a**2*sec(c + d*x) + a**2)*tan(c + d*x)/(2*d) + a**2*(2*A + 3*B)*tan(c + d*x)/(2*d) + a**2*(4*A + 3*B)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_57():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)
    F = B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/d + a**2*x*(2*A + B) + a**2*(A - B)*sin(c + d*x)/d + a**2*(A + 2*B)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_58():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**2
    F = A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)/(2*d) + B*a**2*atanh(sin(c + d*x))/d + a**2*x*(3*A + 4*B)/2 + a**2*(3*A + 2*B)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_59():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**3
    F = A*(a*sec(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**2*x*(2*A + 3*B)/2 + a**2*(2*A + 3*B)*sin(c + d*x)*cos(c + d*x)/(6*d) + 2*a**2*(2*A + 3*B)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_60():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**4
    F = A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a**2*x*(7*A + 8*B)/8 + a**2*(4*A + 5*B)*sin(c + d*x)/(3*d) + a**2*(5*A + 4*B)*sin(c + d*x)*cos(c + d*x)**2/(12*d) + a**2*(7*A + 8*B)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_61():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**5
    F = A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a**2*x*(6*A + 7*B)/8 + a**2*(6*A + 5*B)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + a**2*(6*A + 7*B)*sin(c + d*x)*cos(c + d*x)/(8*d) - a**2*(9*A + 10*B)*sin(c + d*x)**3/(15*d) + a**2*(9*A + 10*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_62():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*sec(c + d*x)**3
    F = B*a*(a*sec(c + d*x) + a)**2*tan(c + d*x)*sec(c + d*x)**3/(6*d) + a**3*(19*A + 17*B)*tan(c + d*x)**3/(15*d) + a**3*(19*A + 17*B)*tan(c + d*x)/(5*d) + a**3*(22*A + 21*B)*tan(c + d*x)*sec(c + d*x)**3/(40*d) + a**3*(26*A + 23*B)*tan(c + d*x)*sec(c + d*x)/(16*d) + a**3*(26*A + 23*B)*atanh(sin(c + d*x))/(16*d) + (3*A + 4*B)*(a**3*sec(c + d*x) + a**3)*tan(c + d*x)*sec(c + d*x)**3/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_63():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*sec(c + d*x)**2
    F = B*(a*sec(c + d*x) + a)**4*tan(c + d*x)/(5*a*d) + a**3*(15*A + 13*B)*tan(c + d*x)**3/(60*d) + 3*a**3*(15*A + 13*B)*tan(c + d*x)*sec(c + d*x)/(40*d) + a**3*(15*A + 13*B)*tan(c + d*x)/(5*d) + a**3*(15*A + 13*B)*atanh(sin(c + d*x))/(8*d) + (5*A - B)*(a*sec(c + d*x) + a)**3*tan(c + d*x)/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_64():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*sec(c + d*x)
    F = B*(a*sec(c + d*x) + a)**3*tan(c + d*x)/(4*d) + a**3*(4*A + 3*B)*tan(c + d*x)**3/(12*d) + 3*a**3*(4*A + 3*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a**3*(4*A + 3*B)*tan(c + d*x)/d + 5*a**3*(4*A + 3*B)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_65():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3
    F = A*a**3*x + B*a*(a*sec(c + d*x) + a)**2*tan(c + d*x)/(3*d) + 5*a**3*(A + B)*tan(c + d*x)/(2*d) + a**3*(7*A + 5*B)*atanh(sin(c + d*x))/(2*d) + (3*A + 5*B)*(a**3*sec(c + d*x) + a**3)*tan(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_66():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)
    F = -5*B*a**3*sin(c + d*x)/(2*d) + B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(2*d) + a**3*x*(3*A + B) + a**3*(6*A + 7*B)*atanh(sin(c + d*x))/(2*d) + (A + 2*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_67():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)**2
    F = 5*A*a**3*sin(c + d*x)/(2*d) + A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)/(2*d) + a**3*x*(7*A + 6*B)/2 + a**3*(A + 3*B)*atanh(sin(c + d*x))/d - (A - 2*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_68():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)**3
    F = A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**2/(3*d) + B*a**3*atanh(sin(c + d*x))/d + a**3*x*(5*A + 7*B)/2 + 5*a**3*(A + B)*sin(c + d*x)/(2*d) + (5*A + 3*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_69():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)**4
    F = A*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 5*a**3*x*(3*A + 4*B)/8 - a**3*(3*A + 4*B)*sin(c + d*x)**3/(12*d) + 3*a**3*(3*A + 4*B)*sin(c + d*x)*cos(c + d*x)/(8*d) + a**3*(3*A + 4*B)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_70():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)**5
    F = A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a**3*x*(13*A + 15*B)/8 + a**3*(13*A + 15*B)*sin(c + d*x)*cos(c + d*x)/(8*d) + a**3*(38*A + 45*B)*sin(c + d*x)/(15*d) + a**3*(43*A + 45*B)*sin(c + d*x)*cos(c + d*x)**2/(60*d) + (7*A + 5*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_71():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*cos(c + d*x)**6
    F = A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + a**3*x*(23*A + 26*B)/16 - a**3*(17*A + 19*B)*sin(c + d*x)**3/(15*d) + a**3*(17*A + 19*B)*sin(c + d*x)/(5*d) + a**3*(21*A + 22*B)*sin(c + d*x)*cos(c + d*x)**3/(40*d) + a**3*(23*A + 26*B)*sin(c + d*x)*cos(c + d*x)/(16*d) + (4*A + 3*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**4/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_72():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*sec(c + d*x)**2
    F = B*(a*sec(c + d*x) + a)**5*tan(c + d*x)/(6*a*d) + 2*a**4*(8*A + 7*B)*tan(c + d*x)**3/(15*d) + a**4*(8*A + 7*B)*tan(c + d*x)*sec(c + d*x)**3/(40*d) + 27*a**4*(8*A + 7*B)*tan(c + d*x)*sec(c + d*x)/(80*d) + 4*a**4*(8*A + 7*B)*tan(c + d*x)/(5*d) + 7*a**4*(8*A + 7*B)*atanh(sin(c + d*x))/(16*d) + (6*A - B)*(a*sec(c + d*x) + a)**4*tan(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_73():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*sec(c + d*x)
    F = B*(a*sec(c + d*x) + a)**4*tan(c + d*x)/(5*d) + 4*a**4*(5*A + 4*B)*tan(c + d*x)**3/(15*d) + a**4*(5*A + 4*B)*tan(c + d*x)*sec(c + d*x)**3/(20*d) + 27*a**4*(5*A + 4*B)*tan(c + d*x)*sec(c + d*x)/(40*d) + 8*a**4*(5*A + 4*B)*tan(c + d*x)/(5*d) + 7*a**4*(5*A + 4*B)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_74():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4
    F = A*a**4*x + B*a*(a*sec(c + d*x) + a)**3*tan(c + d*x)/(4*d) + 5*a**4*(8*A + 7*B)*tan(c + d*x)/(8*d) + a**4*(48*A + 35*B)*atanh(sin(c + d*x))/(8*d) + (4*A + 7*B)*(a**2*sec(c + d*x) + a**2)**2*tan(c + d*x)/(12*d) + (32*A + 35*B)*(a**4*sec(c + d*x) + a**4)*tan(c + d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_75():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)
    F = B*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)/(3*d) + a**4*x*(4*A + B) - 5*a**4*(A + 2*B)*sin(c + d*x)/(2*d) + a**4*(13*A + 12*B)*atanh(sin(c + d*x))/(2*d) + (A + 2*B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)/(2*d) + (9*A + 11*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_76():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**2
    F = A*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)/(2*d) + a**4*x*(13*A + 8*B)/2 + 5*a**4*(A - B)*sin(c + d*x)/(2*d) + a**4*(8*A + 13*B)*atanh(sin(c + d*x))/(2*d) - (A - B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)/(2*d) + (A + 6*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_77():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**3
    F = A*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**4*x*(12*A + 13*B)/2 + a**4*(A + 4*B)*atanh(sin(c + d*x))/d + 5*a**4*(2*A + B)*sin(c + d*x)/(2*d) + (2*A + B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)/(2*d) - (8*A - 3*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_78():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**4
    F = A*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + B*a**4*atanh(sin(c + d*x))/d + a**4*x*(35*A + 48*B)/8 + 5*a**4*(7*A + 8*B)*sin(c + d*x)/(8*d) + (7*A + 4*B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)**2/(12*d) + (35*A + 32*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)*cos(c + d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_79():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**5
    F = A*(a*sec(c + d*x) + a)**4*sin(c + d*x)*cos(c + d*x)**4/(5*d) + 7*a**4*x*(4*A + 5*B)/8 - 4*a**4*(4*A + 5*B)*sin(c + d*x)**3/(15*d) + a**4*(4*A + 5*B)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + 27*a**4*(4*A + 5*B)*sin(c + d*x)*cos(c + d*x)/(40*d) + 8*a**4*(4*A + 5*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_80():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**6
    F = A*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 7*a**4*x*(7*A + 8*B)/16 + 7*a**4*(7*A + 8*B)*sin(c + d*x)*cos(c + d*x)/(16*d) + a**4*(72*A + 83*B)*sin(c + d*x)/(15*d) + a**4*(159*A + 176*B)*sin(c + d*x)*cos(c + d*x)**2/(120*d) + (3*A + 2*B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)**4/(10*d) + (73*A + 72*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)*cos(c + d*x)**3/(120*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_81():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**4*cos(c + d*x)**7
    F = A*a*(a*sec(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**6/(7*d) + a**4*x*(44*A + 49*B)/16 + a**4*(44*A + 49*B)*sin(c + d*x)*cos(c + d*x)/(16*d) - a**4*(227*A + 252*B)*sin(c + d*x)**3/(105*d) + a**4*(227*A + 252*B)*sin(c + d*x)/(35*d) + a**4*(276*A + 301*B)*sin(c + d*x)*cos(c + d*x)**3/(280*d) + (7*A + 7*B)*(a**4*sec(c + d*x) + a**4)*sin(c + d*x)*cos(c + d*x)**4/(15*d) + (10*A + 7*B)*(a**2*sec(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)**5/(42*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_82():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**3/(d*(a*sec(c + d*x) + a)) - (3*A - 4*B)*tan(c + d*x)**3/(3*a*d) - (3*A - 4*B)*tan(c + d*x)/(a*d) + (3*A - 3*B)*tan(c + d*x)*sec(c + d*x)/(2*a*d) + (3*A - 3*B)*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_83():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**2/(d*(a*sec(c + d*x) + a)) - (2*A - 3*B)*tan(c + d*x)*sec(c + d*x)/(2*a*d) - (2*A - 3*B)*atanh(sin(c + d*x))/(2*a*d) + (2*A - 2*B)*tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_84():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)
    F = B*tan(c + d*x)/(a*d) - (A - B)*tan(c + d*x)/(d*(a*sec(c + d*x) + a)) + (A - B)*atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_85():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)
    F = B*atanh(sin(c + d*x))/(a*d) + (A - B)*tan(c + d*x)/(d*(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_86():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)
    F = A*x/a - (A - B)*tan(c + d*x)/(d*(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_87():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)) - x*(A - B)/a + (2*A - B)*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_88():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(d*(a*sec(c + d*x) + a)) + x*(3*A - 2*B)/(2*a) - (2*A - 2*B)*sin(c + d*x)/(a*d) + (3*A - 2*B)*sin(c + d*x)*cos(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_89():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**2/(d*(a*sec(c + d*x) + a)) - x*(3*A - 3*B)/(2*a) - (3*A - 3*B)*sin(c + d*x)*cos(c + d*x)/(2*a*d) - (4*A - 3*B)*sin(c + d*x)**3/(3*a*d) + (4*A - 3*B)*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_90():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**4/(3*d*(a*sec(c + d*x) + a)**2) + (7*A - 10*B)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + (7*A - 10*B)*atanh(sin(c + d*x))/(2*a**2*d) + (7*A - 10*B)*tan(c + d*x)*sec(c + d*x)**3/(3*a**2*d*(sec(c + d*x) + 1)) - (8*A - 12*B)*tan(c + d*x)**3/(3*a**2*d) - (8*A - 12*B)*tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_91():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**2) - (4*A - 7*B)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) - (4*A - 7*B)*atanh(sin(c + d*x))/(2*a**2*d) + (5*A - 8*B)*tan(c + d*x)*sec(c + d*x)**2/(3*a**2*d*(sec(c + d*x) + 1)) + (10*A - 16*B)*tan(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_92():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**2/(3*d*(a*sec(c + d*x) + a)**2) - (A - 4*B)*tan(c + d*x)/(3*a**2*d) + (A - 2*B)*atanh(sin(c + d*x))/(a**2*d) - (A - 2*B)*tan(c + d*x)/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_93():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = B*atanh(sin(c + d*x))/(a**2*d) - (A - B)*tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + (2*A - 5*B)*tan(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_94():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)**2
    F = (A - B)*tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + (A + 2*B)*tan(c + d*x)/(3*d*(a**2*sec(c + d*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_95():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**2
    F = A*x/a**2 - (A - B)*tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) - (4*A - B)*tan(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_96():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) - x*(2*A - B)/a**2 - (2*A - B)*sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)) + (10*A - 4*B)*sin(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_97():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + x*(7*A - 4*B)/(2*a**2) + (7*A - 4*B)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) - (8*A - 5*B)*sin(c + d*x)*cos(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)) - (16*A - 10*B)*sin(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_98():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**2/(3*d*(a*sec(c + d*x) + a)**2) - x*(10*A - 7*B)/(2*a**2) - (10*A - 7*B)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) - (10*A - 7*B)*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(sec(c + d*x) + 1)) - (12*A - 8*B)*sin(c + d*x)**3/(3*a**2*d) + (12*A - 8*B)*sin(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_99():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**4/(5*d*(a*sec(c + d*x) + a)**3) + (36*A - 76*B)*tan(c + d*x)*sec(c + d*x)**2/(15*d*(a**3*sec(c + d*x) + a**3)) + (6*A - 11*B)*tan(c + d*x)*sec(c + d*x)**3/(15*a*d*(a*sec(c + d*x) + a)**2) - (6*A - 13*B)*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) - (6*A - 13*B)*atanh(sin(c + d*x))/(2*a**3*d) + (72*A - 152*B)*tan(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_100():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)**3
    F = -(A - 3*B)*tan(c + d*x)/(d*(a**3*sec(c + d*x) + a**3)) + (A - B)*tan(c + d*x)*sec(c + d*x)**3/(5*d*(a*sec(c + d*x) + a)**3) + (4*A - 9*B)*tan(c + d*x)*sec(c + d*x)**2/(15*a*d*(a*sec(c + d*x) + a)**2) + (A - 3*B)*atanh(sin(c + d*x))/(a**3*d) - (7*A - 27*B)*tan(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_101():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = B*atanh(sin(c + d*x))/(a**3*d) + (A - B)*tan(c + d*x)*sec(c + d*x)**2/(5*d*(a*sec(c + d*x) + a)**3) + (4*A - 29*B)*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - (2*A - 7*B)*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_102():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = -(A - B)*tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) + (3*A + 7*B)*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) + (3*A - 8*B)*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_103():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)**3
    F = (A - B)*tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) + (2*A + 3*B)*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) + (2*A + 3*B)*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_104():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**3
    F = A*x/a**3 - (A - B)*tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - (22*A - 2*B)*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - (7*A - 2*B)*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_105():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - (3*A - B)*sin(c + d*x)/(d*(a**3*sec(c + d*x) + a**3)) - (9*A - 4*B)*sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2) - x*(3*A - B)/a**3 + (72*A - 22*B)*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_106():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - (76*A - 36*B)*sin(c + d*x)*cos(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - (11*A - 6*B)*sin(c + d*x)*cos(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2) + x*(13*A - 6*B)/(2*a**3) + (13*A - 6*B)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) - (152*A - 72*B)*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_107():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**2/(5*d*(a*sec(c + d*x) + a)**3) - (23*A - 13*B)*sin(c + d*x)*cos(c + d*x)**2/(3*d*(a**3*sec(c + d*x) + a**3)) - (13*A - 8*B)*sin(c + d*x)*cos(c + d*x)**2/(15*a*d*(a*sec(c + d*x) + a)**2) - x*(23*A - 13*B)/(2*a**3) - (23*A - 13*B)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) - (136*A - 76*B)*sin(c + d*x)**3/(15*a**3*d) + (136*A - 76*B)*sin(c + d*x)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_108():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**6/(a*sec(c + d*x) + a)**4
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**5/(7*d*(a*sec(c + d*x) + a)**4) + (A - 2*B)*tan(c + d*x)*sec(c + d*x)**4/(5*a*d*(a*sec(c + d*x) + a)**3) - (8*A - 21*B)*tan(c + d*x)*sec(c + d*x)/(2*a**4*d) - (8*A - 21*B)*atanh(sin(c + d*x))/(2*a**4*d) + (52*A - 129*B)*tan(c + d*x)*sec(c + d*x)**3/(105*a**4*d*(sec(c + d*x) + 1)**2) + (332*A - 864*B)*tan(c + d*x)*sec(c + d*x)**2/(105*a**4*d*(sec(c + d*x) + 1)) + (664*A - 1728*B)*tan(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_109():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**5/(a*sec(c + d*x) + a)**4
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**4/(7*d*(a*sec(c + d*x) + a)**4) + (5*A - 12*B)*tan(c + d*x)*sec(c + d*x)**3/(35*a*d*(a*sec(c + d*x) + a)**3) + (A - 4*B)*atanh(sin(c + d*x))/(a**4*d) - (A - 4*B)*tan(c + d*x)/(a**4*d*(sec(c + d*x) + 1)) + (25*A - 88*B)*tan(c + d*x)*sec(c + d*x)**2/(105*a**4*d*(sec(c + d*x) + 1)**2) - (55*A - 244*B)*tan(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_110():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)**4
    F = B*atanh(sin(c + d*x))/(a**4*d) + (A - B)*tan(c + d*x)*sec(c + d*x)**3/(7*d*(a*sec(c + d*x) + a)**4) + (3*A - 10*B)*tan(c + d*x)*sec(c + d*x)**2/(35*a*d*(a*sec(c + d*x) + a)**3) - (6*A - 55*B)*tan(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)**2) + (12*A - 215*B)*tan(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_111():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)**4
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)**3/(7*d*(a*sec(c + d*x) + a)**4) + (4*A + 3*B)*tan(c + d*x)/(15*d*(a**4*sec(c + d*x) + a**4)) - (32*A + 24*B)*tan(c + d*x)/(105*d*(a**2*sec(c + d*x) + a**2)**2) + (4*A + 3*B)*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_112():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)**4
    F = -(A - B)*tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) + (8*A + 13*B)*tan(c + d*x)/(105*d*(a**4*sec(c + d*x) + a**4)) + (8*A + 13*B)*tan(c + d*x)/(105*d*(a**2*sec(c + d*x) + a**2)**2) + (4*A - 11*B)*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_113():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)**4
    F = (A - B)*tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) + (6*A + 8*B)*tan(c + d*x)/(105*d*(a**4*sec(c + d*x) + a**4)) + (6*A + 8*B)*tan(c + d*x)/(105*d*(a**2*sec(c + d*x) + a**2)**2) + (3*A + 4*B)*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_114():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**4
    F = A*x/a**4 - (A - B)*tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - (10*A - 3*B)*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3) - (55*A - 6*B)*tan(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)**2) - (160*A - 6*B)*tan(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_115():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)**4
    F = -(A - B)*sin(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - (12*A - 5*B)*sin(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3) - x*(4*A - B)/a**4 - (4*A - B)*sin(c + d*x)/(a**4*d*(sec(c + d*x) + 1)) - (88*A - 25*B)*sin(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)**2) + (664*A - 160*B)*sin(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_116():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)**4
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - (2*A - B)*sin(c + d*x)*cos(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**3) + x*(21*A - 8*B)/(2*a**4) + (21*A - 8*B)*sin(c + d*x)*cos(c + d*x)/(2*a**4*d) - (129*A - 52*B)*sin(c + d*x)*cos(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)**2) - (864*A - 332*B)*sin(c + d*x)*cos(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)) - (1728*A - 664*B)*sin(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_117():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a*sec(c + d*x) + a)**4
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**2/(7*d*(a*sec(c + d*x) + a)**4) - (16*A - 9*B)*sin(c + d*x)*cos(c + d*x)**2/(35*a*d*(a*sec(c + d*x) + a)**3) - x*(44*A - 21*B)/(2*a**4) - (44*A - 21*B)*sin(c + d*x)*cos(c + d*x)/(2*a**4*d) - (44*A - 21*B)*sin(c + d*x)*cos(c + d*x)**2/(3*a**4*d*(sec(c + d*x) + 1)) - (178*A - 87*B)*sin(c + d*x)*cos(c + d*x)**2/(105*a**4*d*(sec(c + d*x) + 1)**2) - (1816*A - 864*B)*sin(c + d*x)**3/(105*a**4*d) + (1816*A - 864*B)*sin(c + d*x)/(35*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_118():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**4
    F = 2*B*a*tan(c + d*x)*sec(c + d*x)**4/(9*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(9*A + 8*B)*tan(c + d*x)*sec(c + d*x)**3/(63*d*sqrt(a*sec(c + d*x) + a)) + 4*a*(9*A + 8*B)*tan(c + d*x)/(45*d*sqrt(a*sec(c + d*x) + a)) - (72*A + 64*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(315*d) + (36*A + 32*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_119():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**3
    F = 2*B*a*tan(c + d*x)*sec(c + d*x)**3/(7*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(7*A + 6*B)*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) - (28*A + 24*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(105*d) + (14*A + 12*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_120():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**2
    F = 2*B*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*a*d) + 2*a*(5*A + 7*B)*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) + (10*A - 4*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_121():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)
    F = 2*B*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d) + 2*a*(3*A + B)*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_122():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)
    F = 2*A*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*B*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_123():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)
    F = A*a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(A + 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_124():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**2
    F = A*a*sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(3*A + 4*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a*(3*A + 4*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_125():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**3
    F = A*a*sin(c + d*x)*cos(c + d*x)**2/(3*d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(5*A + 6*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a*(5*A + 6*B)*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + a*(5*A + 6*B)*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_126():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**4
    F = A*a*sin(c + d*x)*cos(c + d*x)**3/(4*d*sqrt(a*sec(c + d*x) + a)) + 5*sqrt(a)*(7*A + 8*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a*(7*A + 8*B)*sin(c + d*x)*cos(c + d*x)**2/(24*d*sqrt(a*sec(c + d*x) + a)) + 5*a*(7*A + 8*B)*sin(c + d*x)*cos(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)) + 5*a*(7*A + 8*B)*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_127():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 2*B*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**3/(9*d) + 2*a**2*(9*A + 10*B)*tan(c + d*x)*sec(c + d*x)**3/(63*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(39*A + 34*B)*tan(c + d*x)/(45*d*sqrt(a*sec(c + d*x) + a)) - 4*a*(39*A + 34*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(315*d) + (78*A + 68*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_128():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 2*B*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(7*a*d) + 8*a**2*(21*A + 19*B)*tan(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(21*A + 19*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(105*d) + (14*A - 4*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_129():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*B*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*d) + 8*a**2*(5*A + 3*B)*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(5*A + 3*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_130():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*A*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*B*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d) + 2*a**2*(3*A + 4*B)*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_131():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/d + a**(sympy.S(3)/2)*(3*A + 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**2*(A - 2*B)*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_132():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)/(2*d) + a**(sympy.S(3)/2)*(7*A + 12*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**2*(5*A + 4*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_133():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**(sympy.S(3)/2)*(11*A + 14*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**2*(7*A + 6*B)*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + a**2*(11*A + 14*B)*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_134():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**4
    F = A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a**(sympy.S(3)/2)*(75*A + 88*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a**2*(9*A + 8*B)*sin(c + d*x)*cos(c + d*x)**2/(24*d*sqrt(a*sec(c + d*x) + a)) + a**2*(75*A + 88*B)*sin(c + d*x)*cos(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)) + a**2*(75*A + 88*B)*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_135():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 2*B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)*sec(c + d*x)**3/(11*d) + 2*a**3*(209*A + 194*B)*tan(c + d*x)*sec(c + d*x)**3/(693*d*sqrt(a*sec(c + d*x) + a)) + 2*a**3*(803*A + 710*B)*tan(c + d*x)/(495*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(11*A + 14*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**3/(99*d) - 4*a**2*(803*A + 710*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3465*d) + 2*a*(803*A + 710*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(1155*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_136():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 2*B*(a*sec(c + d*x) + a)**(sympy.S(7)/2)*tan(c + d*x)/(9*a*d) + 64*a**3*(15*A + 13*B)*tan(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*(15*A + 13*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(315*d) + 2*a*(15*A + 13*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(105*d) + (18*A - 4*B)*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_137():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 2*B*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(7*d) + 64*a**3*(7*A + 5*B)*tan(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*(7*A + 5*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(105*d) + 2*a*(7*A + 5*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_138():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*A*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*d) + 2*a**3*(35*A + 32*B)*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(5*A + 8*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_139():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(3*d) + a**(sympy.S(5)/2)*(5*A + 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - a**3*(3*A + 14*B)*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(A + 2*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_140():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(2*d) + a**(sympy.S(5)/2)*(19*A + 20*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**3*(9*A - 4*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)) - a**2*(A - 4*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_141():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3
    F = A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**(sympy.S(5)/2)*(25*A + 38*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**3*(49*A + 54*B)*sin(c + d*x)/(24*d*sqrt(a*sec(c + d*x) + a)) + a**2*(3*A + 2*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_142():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4
    F = A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a**(sympy.S(5)/2)*(163*A + 200*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a**3*(95*A + 104*B)*sin(c + d*x)*cos(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)) + a**3*(163*A + 200*B)*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a)) + a**2*(11*A + 8*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**2/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_143():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5
    F = A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a**(sympy.S(5)/2)*(283*A + 326*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(128*d) + a**3*(157*A + 170*B)*sin(c + d*x)*cos(c + d*x)**2/(240*d*sqrt(a*sec(c + d*x) + a)) + a**3*(283*A + 326*B)*sin(c + d*x)*cos(c + d*x)/(192*d*sqrt(a*sec(c + d*x) + a)) + a**3*(283*A + 326*B)*sin(c + d*x)/(128*d*sqrt(a*sec(c + d*x) + a)) + a**2*(13*A + 10*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**3/(40*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_144():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/sqrt(a*sec(c + d*x) + a)
    F = 2*B*tan(c + d*x)*sec(c + d*x)**3/(7*d*sqrt(a*sec(c + d*x) + a)) + (14*A - 2*B)*tan(c + d*x)*sec(c + d*x)**2/(35*d*sqrt(a*sec(c + d*x) + a)) + (196*A - 148*B)*tan(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)) - (14*A - 62*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(105*a*d) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_145():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/sqrt(a*sec(c + d*x) + a)
    F = 2*B*tan(c + d*x)*sec(c + d*x)**2/(5*d*sqrt(a*sec(c + d*x) + a)) - (20*A - 28*B)*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) + (10*A - 2*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*a*d) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_146():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = 2*B*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*a*d) + (6*A - 4*B)*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_147():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = 2*B*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_148():
    f = (A + B*sec(c + d*x))/sqrt(a*sec(c + d*x) + a)
    F = 2*A*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_149():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = A*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) - (A - 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_150():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = A*sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)) - (A - 4*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d) + (7*A - 4*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_151():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/sqrt(a*sec(c + d*x) + a)
    F = A*sin(c + d*x)*cos(c + d*x)**2/(3*d*sqrt(a*sec(c + d*x) + a)) - (A - 6*B)*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + (7*A - 2*B)*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d) - (9*A - 14*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_152():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**3/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - (5*A - 9*B)*tan(c + d*x)*sec(c + d*x)**2/(10*a*d*sqrt(a*sec(c + d*x) + a)) - (65*A - 93*B)*tan(c + d*x)/(15*a*d*sqrt(a*sec(c + d*x) + a)) + (35*A - 39*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(30*a**2*d) + sqrt(2)*(11*A - 15*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_153():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**2/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (9*A - 13*B)*tan(c + d*x)/(3*a*d*sqrt(a*sec(c + d*x) + a)) - (3*A - 7*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(6*a**2*d) - sqrt(2)*(7*A - 11*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_154():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*tan(c + d*x)/(a*d*sqrt(a*sec(c + d*x) + a)) - (A - B)*tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 7*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_155():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A + 3*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_156():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*A*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - (A - B)*tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(5*A - B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_157():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (3*A - B)*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)) - (3*A - 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + sqrt(2)*(9*A - 5*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_158():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (2*A - B)*sin(c + d*x)*cos(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)) - (7*A - 6*B)*sin(c + d*x)/(4*a*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(13*A - 9*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d) + (19*A - 12*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_159():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**2/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (5*A - 3*B)*sin(c + d*x)*cos(c + d*x)**2/(6*a*d*sqrt(a*sec(c + d*x) + a)) - (13*A - 12*B)*sin(c + d*x)*cos(c + d*x)/(12*a*d*sqrt(a*sec(c + d*x) + a)) + (21*A - 14*B)*sin(c + d*x)/(8*a*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(17*A - 13*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d) - (47*A - 38*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_160():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**3/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (9*A - 17*B)*tan(c + d*x)*sec(c + d*x)**2/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (93*A - 197*B)*tan(c + d*x)/(24*a**2*d*sqrt(a*sec(c + d*x) + a)) - (39*A - 95*B)*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(48*a**3*d) - sqrt(2)*(75*A - 163*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_161():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*tan(c + d*x)*sec(c + d*x)**2/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (5*A - 13*B)*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - (A - 9*B)*tan(c + d*x)/(4*a**2*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(19*A - 75*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_162():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (5*A - 13*B)*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(5*A + 19*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_163():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (3*A + 5*B)*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A + 5*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_164():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*A*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - (A - B)*tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (11*A - 3*B)*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(43*A - 3*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_165():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (15*A - 7*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (35*A - 11*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - (5*A - 2*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + sqrt(2)*(115*A - 43*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_166():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (19*A - 11*B)*sin(c + d*x)*cos(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (31*A - 15*B)*sin(c + d*x)*cos(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - (63*A - 35*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) + (39*A - 20*B)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*a**(sympy.S(5)/2)*d) - sqrt(2)*(219*A - 115*B)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_167():
    f = (A*sec(c + d*x) + A)/sqrt(-a*sec(c + d*x) + a)
    F = 2*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(sqrt(a)*d) - 2*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_168():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)/sqrt(-a*sec(c + d*x) + a)
    F = A*sin(c + d*x)/(d*sqrt(-a*sec(c + d*x) + a)) + 3*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(sqrt(a)*d) - 2*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_169():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**2/sqrt(-a*sec(c + d*x) + a)
    F = A*sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(-a*sec(c + d*x) + a)) + 5*A*sin(c + d*x)/(4*d*sqrt(-a*sec(c + d*x) + a)) + 11*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(4*sqrt(a)*d) - 2*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_170():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**3/sqrt(-a*sec(c + d*x) + a)
    F = A*sin(c + d*x)*cos(c + d*x)**2/(3*d*sqrt(-a*sec(c + d*x) + a)) + 7*A*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(-a*sec(c + d*x) + a)) + 9*A*sin(c + d*x)/(8*d*sqrt(-a*sec(c + d*x) + a)) + 23*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(8*sqrt(a)*d) - 2*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_171():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)/(-a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -A*sin(c + d*x)/(d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*A*sin(c + d*x)/(a*d*sqrt(-a*sec(c + d*x) + a)) + 5*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 7*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_172():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**2/(-a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -A*sin(c + d*x)*cos(c + d*x)/(d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*A*sin(c + d*x)*cos(c + d*x)/(2*a*d*sqrt(-a*sec(c + d*x) + a)) + 13*A*sin(c + d*x)/(4*a*d*sqrt(-a*sec(c + d*x) + a)) + 31*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d) - 11*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_173():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**3/(-a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -A*sin(c + d*x)*cos(c + d*x)**2/(d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 4*A*sin(c + d*x)*cos(c + d*x)**2/(3*a*d*sqrt(-a*sec(c + d*x) + a)) + 25*A*sin(c + d*x)*cos(c + d*x)/(12*a*d*sqrt(-a*sec(c + d*x) + a)) + 35*A*sin(c + d*x)/(8*a*d*sqrt(-a*sec(c + d*x) + a)) + 85*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(8*a**(sympy.S(3)/2)*d) - 15*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_174():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)/(-a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -A*sin(c + d*x)/(2*d*(-a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 11*A*sin(c + d*x)/(8*a*d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 23*A*sin(c + d*x)/(8*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 7*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 79*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(16*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_175():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**2/(-a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -A*sin(c + d*x)*cos(c + d*x)/(2*d*(-a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 15*A*sin(c + d*x)*cos(c + d*x)/(8*a*d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 23*A*sin(c + d*x)*cos(c + d*x)/(8*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 49*A*sin(c + d*x)/(8*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 59*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(4*a**(sympy.S(5)/2)*d) - 167*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(16*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_176():
    f = (A*sec(c + d*x) + A)*cos(c + d*x)**3/(-a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -A*sin(c + d*x)*cos(c + d*x)**2/(2*d*(-a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 19*A*sin(c + d*x)*cos(c + d*x)**2/(8*a*d*(-a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 77*A*sin(c + d*x)*cos(c + d*x)**2/(24*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 119*A*sin(c + d*x)*cos(c + d*x)/(24*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 21*A*sin(c + d*x)/(2*a**2*d*sqrt(-a*sec(c + d*x) + a)) + 203*A*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(8*a**(sympy.S(5)/2)*d) - 287*sqrt(2)*A*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(16*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_177():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 2*a*(A + B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 6*a*(A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*(7*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 2*a*(7*A + 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_178():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*(A + B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a*(5*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 2*a*(5*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_179():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(3*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_180():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*B*a*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*a*(A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_181():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_182():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a*(3*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_183():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(A + B)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*(5*A + 7*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*(5*A + 7*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_184():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d) + 4*a**2*(4*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 4*a**2*(4*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*(7*A + 6*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 4*a**2*(7*A + 6*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*a**2*(7*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_185():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))
    F = 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a**2*(2*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(5*A + 4*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 4*a**2*(5*A + 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a**2*(5*A + 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_186():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sqrt(sec(c + d*x))
    F = -4*B*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(3*A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a**2*(3*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_187():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 4*A*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) - 2*a**2*(A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(2*A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_188():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(4*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a**2*(7*A + 5*B)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_189():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(3*A + 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*(6*A + 7*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a**2*(6*A + 7*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*a**2*(9*A + 7*B)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_190():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 4*a**2*(5*A + 6*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a**2*(5*A + 6*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**2*(8*A + 9*B)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(8*A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 2*a**2*(11*A + 9*B)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_191():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(9*d) + 4*a**3*(13*A + 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 4*a**3*(13*A + 11*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(21*A + 17*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) - 4*a**3*(21*A + 17*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 4*a**3*(24*A + 23*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(105*d) + (18*A + 26*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_192():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 4*a**3*(9*A + 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 4*a**3*(9*A + 7*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(21*A + 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(42*A + 41*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d) + (14*A + 22*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_193():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sqrt(sec(c + d*x))
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(5*A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 4*a**3*(5*A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(20*A + 21*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) + (10*A + 18*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_194():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a**3*(A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**3*(A + 4*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) - (2*A - 2*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_195():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**3*(3*A + 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 4*a**3*(6*A - 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) + 4*a**3*(9*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (18*A + 10*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_196():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a**3*(7*A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(13*A + 21*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(41*A + 42*B)*sin(c + d*x)/(105*d*sqrt(sec(c + d*x))) + (22*A + 14*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_197():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 4*a**3*(11*A + 13*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a**3*(11*A + 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(17*A + 21*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 4*a**3*(23*A + 24*B)*sin(c + d*x)/(105*d*sec(c + d*x)**(sympy.S(3)/2)) + (26*A + 18*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_198():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 4*a**3*(15*A + 17*B)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**3*(15*A + 17*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 20*a**3*(21*A + 22*B)*sin(c + d*x)/(693*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a**3*(105*A + 121*B)*sin(c + d*x)/(231*d*sqrt(sec(c + d*x))) + 4*a**3*(105*A + 121*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(231*d) + (30*A + 22*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(99*d*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_199():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(d*(a*sec(c + d*x) + a)) - (5*A - 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*a*d) + (5*A - 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) + (5*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d) - (15*A - 21*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*a*d) + (15*A - 21*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_200():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(d*(a*sec(c + d*x) + a)) + (-3*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) - (3*A - 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) - (3*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d) + (3*A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_201():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*(a*sec(c + d*x) + a)) - (A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) + (A - 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (A - B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_202():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) - (A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_203():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) - (A - B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (3*A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_204():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + (-3*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (5*A - 3*B)*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) + (5*A - 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_205():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - (5*A - 5*B)*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) - (5*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d) + (7*A - 5*B)*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) + (21*A - 15*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_206():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2))
    F = -(A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + (-21*A + 21*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d) - (7*A - 7*B)*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) + (9*A - 7*B)*sin(c + d*x)/(7*a*d*sec(c + d*x)**(sympy.S(5)/2)) + (45*A - 35*B)*sin(c + d*x)/(21*a*d*sqrt(sec(c + d*x))) + (45*A - 35*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_207():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*d*(a*sec(c + d*x) + a)**2) + (4*A - 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) - (4*A - 7*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + (4*A - 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*a**2*d*(sec(c + d*x) + 1)) - (5*A - 10*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d) - (5*A - 10*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_208():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d*(a*sec(c + d*x) + a)**2) - (A - 4*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) + (A - 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + (2*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + (2*A - 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_209():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - B*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1)) + (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + (A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_210():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**2
    F = -A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + (2*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + (2*A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_211():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + (4*A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - (5*A - 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - (5*A - 2*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_212():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) - (7*A - 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - (7*A - 4*B)*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*sqrt(sec(c + d*x))) + (10*A - 5*B)*sin(c + d*x)/(3*a**2*d*sqrt(sec(c + d*x))) + (10*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_213():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)) - (3*A - 2*B)*sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2)) - (15*A - 10*B)*sin(c + d*x)/(3*a**2*d*sqrt(sec(c + d*x))) - (15*A - 10*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + (56*A - 35*B)*sin(c + d*x)/(15*a**2*d*sec(c + d*x)**(sympy.S(3)/2)) + (56*A - 35*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_214():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)/(a*sec(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(5*d*(a*sec(c + d*x) + a)**3) + (49*A - 119*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(30*d*(a**3*sec(c + d*x) + a**3)) + (A - 2*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*a*d*(a*sec(c + d*x) + a)**2) + (-49*A + 119*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - (13*A - 33*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*a**3*d) - (13*A - 33*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (49*A - 119*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_215():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(5*d*(a*sec(c + d*x) + a)**3) + (3*A - 13*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a**3*sec(c + d*x) + a**3)) + (3*A - 8*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) + (3*A - 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) - (9*A - 49*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d) + (9*A - 49*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_216():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*(a*sec(c + d*x) + a)**3) - (A + 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(10*d*(a**3*sec(c + d*x) + a**3)) + (A - 6*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) + (A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_217():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) + (A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - (A + 4*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) - (A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + (A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_218():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) + (3*A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) + (3*A + 2*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + (3*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) - (9*A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_219():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) - (13*A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - (8*A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) - (13*A - 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (49*A - 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_220():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))) - (119*A - 49*B)*sin(c + d*x)/(30*d*(a**3*sec(c + d*x) + a**3)*sqrt(sec(c + d*x))) - (2*A - B)*sin(c + d*x)/(3*a*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + (-119*A + 49*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + (33*A - 13*B)*sin(c + d*x)/(6*a**3*d*sqrt(sec(c + d*x))) + (33*A - 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_221():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)) - (63*A - 33*B)*sin(c + d*x)/(10*d*(a**3*sec(c + d*x) + a**3)*sec(c + d*x)**(sympy.S(3)/2)) - (12*A - 7*B)*sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)) - (21*A - 11*B)*sin(c + d*x)/(2*a**3*d*sqrt(sec(c + d*x))) - (21*A - 11*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d) + (231*A - 119*B)*sin(c + d*x)/(30*a**3*d*sec(c + d*x)**(sympy.S(3)/2)) + (231*A - 119*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_222():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(6*A + 5*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a*(6*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(12*d*sqrt(a*sec(c + d*x) + a)) + a*(6*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_223():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(4*A + 3*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a*(4*A + 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_224():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))
    F = B*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(a*sec(c + d*x) + a)) + sqrt(a)*(2*A + B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_225():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*A*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a)) + 2*B*sqrt(a)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_226():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_227():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 4*a*(4*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(4*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_228():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 16*a*(6*A + 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 8*a*(6*A + 7*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a*(6*A + 7*B)*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_229():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d) + a**(sympy.S(3)/2)*(88*A + 75*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a**2*(8*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(24*d*sqrt(a*sec(c + d*x) + a)) + a**2*(88*A + 75*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(96*d*sqrt(a*sec(c + d*x) + a)) + a**2*(88*A + 75*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_230():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d) + a**(sympy.S(3)/2)*(14*A + 11*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**2*(6*A + 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(12*d*sqrt(a*sec(c + d*x) + a)) + a**2*(14*A + 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_231():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d) + a**(sympy.S(3)/2)*(12*A + 7*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**2*(4*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_232():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/d + a**(sympy.S(3)/2)*(2*A + 3*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**2*(2*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_233():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*B*a**(sympy.S(3)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**2*(4*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_234():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a**2*(3*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(3*A + 5*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_235():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(8*A + 7*B)*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(52*A + 63*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(52*A + 63*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_236():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**2*(10*A + 9*B)*sin(c + d*x)/(63*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 16*a**2*(34*A + 39*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 8*a**2*(34*A + 39*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a**2*(34*A + 39*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_237():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(5*d) + a**(sympy.S(5)/2)*(326*A + 283*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(128*d) + a**3*(170*A + 157*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(240*d*sqrt(a*sec(c + d*x) + a)) + a**3*(326*A + 283*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(192*d*sqrt(a*sec(c + d*x) + a)) + a**3*(326*A + 283*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(128*d*sqrt(a*sec(c + d*x) + a)) + a**2*(10*A + 13*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(40*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_238():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d) + a**(sympy.S(5)/2)*(200*A + 163*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a**3*(104*A + 95*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(96*d*sqrt(a*sec(c + d*x) + a)) + a**3*(200*A + 163*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*d*sqrt(a*sec(c + d*x) + a)) + a**2*(8*A + 11*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_239():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + a**(sympy.S(5)/2)*(38*A + 25*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**3*(54*A + 49*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(24*d*sqrt(a*sec(c + d*x) + a)) + a**2*(2*A + 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_240():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*d) + a**(sympy.S(5)/2)*(20*A + 19*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**3*(4*A - 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(4*d*sqrt(a*sec(c + d*x) + a)) + a**2*(4*A + 7*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_241():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + a**(sympy.S(5)/2)*(2*A + 5*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**3*(14*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) - a**2*(2*A - 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_242():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*B*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**3*(32*A + 35*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(8*A + 5*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_243():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 64*a**3*(5*A + 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*(5*A + 7*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(105*d*sqrt(sec(c + d*x))) + 2*a*(5*A + 7*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_244():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**3*(124*A + 135*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**3*(292*A + 345*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 2*a**3*(292*A + 345*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a**2*(4*A + 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(21*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_245():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 2*a**3*(194*A + 209*B)*sin(c + d*x)/(693*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 16*a**3*(710*A + 803*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3465*d*sqrt(a*sec(c + d*x) + a)) + 8*a**3*(710*A + 803*B)*sin(c + d*x)/(3465*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a**3*(710*A + 803*B)*sin(c + d*x)/(1155*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(14*A + 11*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(99*d*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_246():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/sqrt(a*sec(c + d*x) + a)
    F = B*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(a*sec(c + d*x) + a)) + (4*A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d) - (4*A - 7*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_247():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/sqrt(a*sec(c + d*x) + a)
    F = B*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d) + (2*A - B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_248():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/sqrt(a*sec(c + d*x) + a)
    F = 2*B*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_249():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = 2*A*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_250():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*A*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) - (2*A - 6*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_251():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*A*sin(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - (2*A - 10*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + (26*A - 10*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_252():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2))
    F = 2*A*sin(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) - (2*A - 14*B)*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + (62*A - 14*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) - (86*A - 182*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_253():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - (A - 2*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*a*d*sqrt(a*sec(c + d*x) + a)) + (6*A - 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*a*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(9*A - 13*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d) - (12*A - 19*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_254():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - (A - 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*a*d*sqrt(a*sec(c + d*x) + a)) + (2*A - 3*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - sqrt(2)*(5*A - 9*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_255():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A - 5*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_256():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A + B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_257():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (5*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*a*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(7*A - 3*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_258():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + (7*A - 3*B)*sin(c + d*x)/(6*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) - (19*A - 15*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*a*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(11*A - 7*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_259():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + (9*A - 5*B)*sin(c + d*x)/(10*a*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - (39*A - 35*B)*sin(c + d*x)/(30*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + (147*A - 95*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(30*a*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(15*A - 11*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_260():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (7*A - 15*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - (11*A - 35*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) + (2*A - 5*B)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - sqrt(2)*(43*A - 115*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_261():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*B*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (3*A - 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 43*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_262():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + (5*A + 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(5*A + 3*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_263():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (13*A - 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (49*A - 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(75*A - 19*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_264():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - (17*A - 9*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + (95*A - 39*B)*sin(c + d*x)/(48*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) - (299*A - 147*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a**2*d*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*(163*A - 75*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_265():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)) - (21*A - 13*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + (157*A - 85*B)*sin(c + d*x)/(80*a**2*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - (787*A - 475*B)*sin(c + d*x)/(240*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + (2671*A - 1495*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(240*a**2*d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*(283*A - 163*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_266():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*sqrt(2)*A*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S(7)/6, sympy.S.Half, 1, sympy.S(13)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(7*d*sqrt(1 - sec(c + d*x))) - 2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 3*B*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_267():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*sqrt(2)*A*tan(c + d*x)*appellf1(sympy.S(1)/6, sympy.S.Half, 1, sympy.S(7)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_268():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(4)/3)
    F = -3*sqrt(2)*A*tan(c + d*x)*appellf1(sympy.S(-5)/6, sympy.S.Half, 1, sympy.S(1)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(5*a*d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(sec(c + d*x) + 1)) - 2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(10*a*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) + 3*B*tan(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_269():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(4)/3)
    F = 3*sqrt(2)*A*a*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(sec(c + d*x) + 1)*tan(c + d*x)*appellf1(sympy.S(11)/6, sympy.S.Half, 1, sympy.S(17)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(11*d*sqrt(1 - sec(c + d*x))) + 15*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*B*a*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) + 5*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*B*a*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(8*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) + 3*B*a*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)/(4*d) - B*a*(15 + 15*sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)/(4*d*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_270():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*sqrt(2)*A*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S(5)/6, sympy.S.Half, 1, sympy.S(11)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(5*d*sqrt(1 - sec(c + d*x))) + 3*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) + 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) - B*(3 + 3*sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)/(d*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_271():
    f = (A + B*sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(2)/3)
    F = -3*sqrt(2)*A*tan(c + d*x)*appellf1(sympy.S(-1)/6, sympy.S.Half, 1, sympy.S(5)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) - 3*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*B*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + 3*B*tan(c + d*x)/(d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + B*(3 + 3*sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)/(d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_272():
    f = (c*sec(e + f*x))**n*(A + B*sec(e + f*x))*(a*sec(e + f*x) + a)**m
    F = -B*(c*sec(e + f*x))**n*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(n, sympy.S.Half, -m + sympy.S(-1)/2, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))) - (c*sec(e + f*x))**n*(A - B)*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(n, sympy.S.Half, sympy.S.Half - m, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_273():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**n*sec(c + d*x)**(-n - 1)
    F = A*(a*sec(c + d*x) + a)**n*sin(c + d*x)/(d*(n + 1)*sec(c + d*x)**n) + ((sec(c + d*x) + 1)/(1 - sec(c + d*x)))**(sympy.S.Half - n)*(a*sec(c + d*x) + a)**n*(A*n + B*n + B)*sin(c + d*x)*hyper((-n, sympy.S.Half - n), (1 - n,), -2*sec(c + d*x)/(1 - sec(c + d*x)))*sec(c + d*x)**(1 - n)/(d*n*(n + 1)*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_274():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sec(c + d*x)**3
    F = B*b*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (4*A*a + 3*B*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*A*a + 3*B*b)*atanh(sin(c + d*x))/(8*d) + (A*b + B*a)*tan(c + d*x)**3/(3*d) + (A*b + B*a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_275():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sec(c + d*x)**2
    F = B*b*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (3*A*a + 2*B*b)*tan(c + d*x)/(3*d) + (A*b + B*a)*tan(c + d*x)*sec(c + d*x)/(2*d) + (A*b + B*a)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_276():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sec(c + d*x)
    F = B*b*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*A*a + B*b)*atanh(sin(c + d*x))/(2*d) + (A*b + B*a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_277():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))
    F = A*a*x + B*b*tan(c + d*x)/d + (A*b + B*a)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_278():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)
    F = A*a*sin(c + d*x)/d + B*b*atanh(sin(c + d*x))/d + x*(A*b + B*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_279():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**2
    F = A*a*sin(c + d*x)*cos(c + d*x)/(2*d) + x*(A*a/2 + B*b) + (A*b + B*a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_280():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**3
    F = A*a*sin(c + d*x)*cos(c + d*x)**2/(3*d) + x*(A*b/2 + B*a/2) + (2*A*a + 3*B*b)*sin(c + d*x)/(3*d) + (A*b + B*a)*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_281():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**4
    F = A*a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(3*A*a/8 + B*b/2) + (3*A*a + 4*B*b)*sin(c + d*x)*cos(c + d*x)/(8*d) - (A*b + B*a)*sin(c + d*x)**3/(3*d) + (A*b + B*a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_282():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sec(c + d*x)**3
    F = B*b*(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)**3/(5*d) + b*(5*A*b + 6*B*a)*tan(c + d*x)*sec(c + d*x)**3/(20*d) + (4*B*b**2 + 5*a*(2*A*b + B*a))*tan(c + d*x)**3/(15*d) + (4*B*b**2 + 5*a*(2*A*b + B*a))*tan(c + d*x)/(5*d) + (4*A*a**2 + 3*A*b**2 + 6*B*a*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*A*a**2 + 3*A*b**2 + 6*B*a*b)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_283():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sec(c + d*x)**2
    F = B*(a + b*sec(c + d*x))**3*tan(c + d*x)/(4*b*d) + (8*A*a*b - 2*B*a**2 + 9*B*b**2)*tan(c + d*x)*sec(c + d*x)/(24*d) + (8*A*a*b + 4*B*a**2 + 3*B*b**2)*atanh(sin(c + d*x))/(8*d) + (a + b*sec(c + d*x))**2*(4*A*b - B*a)*tan(c + d*x)/(12*b*d) + (4*A*a**2*b + 4*A*b**3 - B*a**3 + 8*B*a*b**2)*tan(c + d*x)/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_284():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sec(c + d*x)
    F = B*(a + b*sec(c + d*x))**2*tan(c + d*x)/(3*d) + b*(3*A*b + 2*B*a)*tan(c + d*x)*sec(c + d*x)/(6*d) + (2*A*a**2 + A*b**2 + 2*B*a*b)*atanh(sin(c + d*x))/(2*d) + (6*A*a*b + 2*B*a**2 + 2*B*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_285():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2
    F = A*a**2*x + B*b*(a + b*sec(c + d*x))*tan(c + d*x)/(2*d) + b*(2*A*b + 3*B*a)*tan(c + d*x)/(2*d) + (4*A*a*b + 2*B*a**2 + B*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_286():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)
    F = A*a**2*sin(c + d*x)/d + B*b**2*tan(c + d*x)/d + a*x*(2*A*b + B*a) + b*(A*b + 2*B*a)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_287():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**2
    F = A*a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + B*b**2*atanh(sin(c + d*x))/d + a*(2*A*b + B*a)*sin(c + d*x)/d + x*(A*a**2/2 + A*b**2 + 2*B*a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_288():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**3
    F = A*a**2*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a*(2*A*b + B*a)*sin(c + d*x)*cos(c + d*x)/(2*d) + x*(A*a*b + B*a**2/2 + B*b**2) + (2*A*a**2 + 3*A*b**2 + 6*B*a*b)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_289():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**4
    F = A*a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) - a*(2*A*b + B*a)*sin(c + d*x)**3/(3*d) + x*(3*A*a**2/8 + A*b**2/2 + B*a*b) + (3*A*a**2 + 4*A*b**2 + 8*B*a*b)*sin(c + d*x)*cos(c + d*x)/(8*d) + (2*A*a*b + B*a**2 + B*b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_290():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**5
    F = A*a**2*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a*(2*A*b + B*a)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(3*A*a*b/4 + 3*B*a**2/8 + B*b**2/2) - (4*A*a**2 + 5*A*b**2 + 10*B*a*b)*sin(c + d*x)**3/(15*d) + (4*A*a**2 + 5*A*b**2 + 10*B*a*b)*sin(c + d*x)/(5*d) + (6*A*a*b + 3*B*a**2 + 4*B*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_291():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*sec(c + d*x)**2
    F = B*(a + b*sec(c + d*x))**4*tan(c + d*x)/(5*b*d) + (12*A*a**2*b + 3*A*b**3 + 4*B*a**3 + 9*B*a*b**2)*atanh(sin(c + d*x))/(8*d) + (30*A*a**2*b + 45*A*b**3 - 6*B*a**3 + 71*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(120*d) + (a + b*sec(c + d*x))**3*(5*A*b - B*a)*tan(c + d*x)/(20*b*d) + (a + b*sec(c + d*x))**2*(15*A*a*b - 3*B*a**2 + 16*B*b**2)*tan(c + d*x)/(60*b*d) + (15*A*a**3*b + 60*A*a*b**3 - 3*B*a**4 + 52*B*a**2*b**2 + 16*B*b**4)*tan(c + d*x)/(30*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_292():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*sec(c + d*x)
    F = B*(a + b*sec(c + d*x))**3*tan(c + d*x)/(4*d) + b*(20*A*a*b + 6*B*a**2 + 9*B*b**2)*tan(c + d*x)*sec(c + d*x)/(24*d) + (a + b*sec(c + d*x))**2*(4*A*b + 3*B*a)*tan(c + d*x)/(12*d) + (8*A*a**3 + 12*A*a*b**2 + 12*B*a**2*b + 3*B*b**3)*atanh(sin(c + d*x))/(8*d) + (16*A*a**2*b + 4*A*b**3 + 3*B*a**3 + 12*B*a*b**2)*tan(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_293():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3
    F = A*a**3*x + B*b*(a + b*sec(c + d*x))**2*tan(c + d*x)/(3*d) + b**2*(3*A*b + 5*B*a)*tan(c + d*x)*sec(c + d*x)/(6*d) + b*(9*A*a*b + 8*B*a**2 + 2*B*b**2)*tan(c + d*x)/(3*d) + (6*A*a**2*b + A*b**3 + 2*B*a**3 + 3*B*a*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_294():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*cos(c + d*x)**2
    F = A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)/(2*d) + a**2*(2*A*b + B*a)*sin(c + d*x)/d + a*x*(A*a**2 + 6*A*b**2 + 6*B*a*b)/2 - b**2*(A*a - 2*B*b)*tan(c + d*x)/(2*d) + b**2*(A*b + 3*B*a)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_295():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*cos(c + d*x)**3
    F = A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**2/(3*d) + B*b**3*atanh(sin(c + d*x))/d + a**2*(5*A*b + 3*B*a)*sin(c + d*x)*cos(c + d*x)/(6*d) + a*(2*A*a**2 + 8*A*b**2 + 9*B*a*b)*sin(c + d*x)/(3*d) + x*(3*A*a**2*b/2 + A*b**3 + B*a**3/2 + 3*B*a*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_296():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*cos(c + d*x)**4
    F = A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a**2*(3*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)**2/(6*d) + a*(3*A*a**2 + 10*A*b**2 + 12*B*a*b)*sin(c + d*x)*cos(c + d*x)/(8*d) + x*(3*A*a**3/8 + 3*A*a*b**2/2 + 3*B*a**2*b/2 + B*b**3) + (6*A*a**2*b + 3*A*b**3 + 2*B*a**3 + 9*B*a*b**2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_297():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*cos(c + d*x)**5
    F = A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a**2*(7*A*b + 5*B*a)*sin(c + d*x)*cos(c + d*x)**3/(20*d) - a*(4*A*a**2 + 12*A*b**2 + 15*B*a*b)*sin(c + d*x)**3/(15*d) + x*(9*A*a**2*b/8 + A*b**3/2 + 3*B*a**3/8 + 3*B*a*b**2/2) + (4*A*a**3 + 14*A*a*b**2 + 15*B*a**2*b + 5*B*b**3)*sin(c + d*x)/(5*d) + (9*A*a**2*b + 4*A*b**3 + 3*B*a**3 + 12*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_298():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*sec(c + d*x)**2
    F = B*(a + b*sec(c + d*x))**5*tan(c + d*x)/(6*b*d) + (32*A*a**3*b + 24*A*a*b**3 + 8*B*a**4 + 36*B*a**2*b**2 + 5*B*b**4)*atanh(sin(c + d*x))/(16*d) + (48*A*a**3*b + 232*A*a*b**3 - 8*B*a**4 + 178*B*a**2*b**2 + 75*B*b**4)*tan(c + d*x)*sec(c + d*x)/(240*d) + (a + b*sec(c + d*x))**4*(6*A*b - B*a)*tan(c + d*x)/(30*b*d) + (a + b*sec(c + d*x))**3*(24*A*a*b - 4*B*a**2 + 25*B*b**2)*tan(c + d*x)/(120*b*d) + (a + b*sec(c + d*x))**2*(24*A*a**2*b + 32*A*b**3 - 4*B*a**3 + 53*B*a*b**2)*tan(c + d*x)/(120*b*d) + (24*A*a**4*b + 224*A*a**2*b**3 + 32*A*b**5 - 4*B*a**5 + 121*B*a**3*b**2 + 128*B*a*b**4)*tan(c + d*x)/(60*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_299():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*sec(c + d*x)
    F = B*(a + b*sec(c + d*x))**4*tan(c + d*x)/(5*d) + b*(130*A*a**2*b + 45*A*b**3 + 24*B*a**3 + 116*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(120*d) + (a + b*sec(c + d*x))**3*(5*A*b + 4*B*a)*tan(c + d*x)/(20*d) + (a + b*sec(c + d*x))**2*(35*A*a*b + 12*B*a**2 + 16*B*b**2)*tan(c + d*x)/(60*d) + (8*A*a**4 + 24*A*a**2*b**2 + 3*A*b**4 + 16*B*a**3*b + 12*B*a*b**3)*atanh(sin(c + d*x))/(8*d) + (95*A*a**3*b + 80*A*a*b**3 + 12*B*a**4 + 112*B*a**2*b**2 + 16*B*b**4)*tan(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_300():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4
    F = A*a**4*x + B*b*(a + b*sec(c + d*x))**3*tan(c + d*x)/(4*d) + b**2*(32*A*a*b + 26*B*a**2 + 9*B*b**2)*tan(c + d*x)*sec(c + d*x)/(24*d) + b*(a + b*sec(c + d*x))**2*(4*A*b + 7*B*a)*tan(c + d*x)/(12*d) + b*(34*A*a**2*b + 4*A*b**3 + 19*B*a**3 + 16*B*a*b**2)*tan(c + d*x)/(6*d) + (32*A*a**3*b + 16*A*a*b**3 + 8*B*a**4 + 24*B*a**2*b**2 + 3*B*b**4)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_301():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)/d + a**3*x*(4*A*b + B*a) - b**2*(6*A*a**2 - 3*A*b**2 - 8*B*a*b)*tan(c + d*x)*sec(c + d*x)/(6*d) - b*(a + b*sec(c + d*x))**2*(3*A*a - B*b)*tan(c + d*x)/(3*d) - b*(6*A*a**3 - 12*A*a*b**2 - 17*B*a**2*b - 2*B*b**3)*tan(c + d*x)/(3*d) + b*(12*A*a**2*b + A*b**3 + 8*B*a**3 + 4*B*a*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_302():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)**2
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)*cos(c + d*x)/(2*d) + a**2*x*(A*a**2 + 12*A*b**2 + 8*B*a*b)/2 + a*(a + b*sec(c + d*x))**2*(5*A*b + 2*B*a)*sin(c + d*x)/(2*d) - b**2*(6*A*a*b + 2*B*a**2 - B*b**2)*tan(c + d*x)*sec(c + d*x)/(2*d) + b**2*(8*A*a*b + 12*B*a**2 + B*b**2)*atanh(sin(c + d*x))/(2*d) - b*(13*A*a**2*b - 2*A*b**3 + 4*B*a**3 - 8*B*a*b**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_303():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)**3
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**2*(2*A*a**2 + 9*A*b**2 + 9*B*a*b)*sin(c + d*x)/(3*d) + a*x*(4*A*a**2*b + 8*A*b**3 + B*a**3 + 12*B*a*b**2)/2 + a*(a + b*sec(c + d*x))**2*(2*A*b + B*a)*sin(c + d*x)*cos(c + d*x)/(2*d) + b**3*(A*b + 4*B*a)*atanh(sin(c + d*x))/d - b**2*(8*A*a*b + 3*B*a**2 - 6*B*b**2)*tan(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_304():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)**4
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + B*b**4*atanh(sin(c + d*x))/d + a**2*(9*A*a**2 + 26*A*b**2 + 32*B*a*b)*sin(c + d*x)*cos(c + d*x)/(24*d) + a*(a + b*sec(c + d*x))**2*(7*A*b + 4*B*a)*sin(c + d*x)*cos(c + d*x)**2/(12*d) + a*(16*A*a**2*b + 19*A*b**3 + 4*B*a**3 + 34*B*a*b**2)*sin(c + d*x)/(6*d) + x*(3*A*a**4/8 + 3*A*a**2*b**2 + A*b**4 + 2*B*a**3*b + 4*B*a*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_305():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)**5
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**4/(5*d) + a**2*(8*A*a**2 + 18*A*b**2 + 25*B*a*b)*sin(c + d*x)*cos(c + d*x)**2/(30*d) + a*(a + b*sec(c + d*x))**2*(8*A*b + 5*B*a)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + a*(60*A*a**2*b + 56*A*b**3 + 15*B*a**3 + 110*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(40*d) + x*(3*A*a**3*b/2 + 2*A*a*b**3 + 3*B*a**4/8 + 3*B*a**2*b**2 + B*b**4) + (8*A*a**4 + 60*A*a**2*b**2 + 15*A*b**4 + 40*B*a**3*b + 60*B*a*b**3)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_306():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*cos(c + d*x)**6
    F = A*a*(a + b*sec(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) + a**2*(25*A*a**2 + 48*A*b**2 + 72*B*a*b)*sin(c + d*x)*cos(c + d*x)**3/(120*d) + a*(a + b*sec(c + d*x))**2*(3*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)**4/(10*d) - a*(16*A*a**2*b + 13*A*b**3 + 4*B*a**3 + 27*B*a*b**2)*sin(c + d*x)**3/(15*d) + x*(5*A*a**4/16 + 9*A*a**2*b**2/4 + A*b**4/2 + 3*B*a**3*b/2 + 2*B*a*b**3) + (5*A*a**4 + 36*A*a**2*b**2 + 8*A*b**4 + 24*B*a**3*b + 32*B*a*b**3)*sin(c + d*x)*cos(c + d*x)/(16*d) + (48*A*a**3*b + 53*A*a*b**3 + 12*B*a**4 + 87*B*a**2*b**2 + 15*B*b**4)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_307():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a + b*sec(c + d*x))
    F = B*tan(c + d*x)*sec(c + d*x)**2/(3*b*d) - 2*a**3*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*sqrt(a - b)*sqrt(a + b)) + (A*b - B*a)*tan(c + d*x)*sec(c + d*x)/(2*b**2*d) - (3*A*a*b - 3*B*a**2 - 2*B*b**2)*tan(c + d*x)/(3*b**3*d) + (2*a**2 + b**2)*(A*b - B*a)*atanh(sin(c + d*x))/(2*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_308():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))
    F = B*tan(c + d*x)*sec(c + d*x)/(2*b*d) + 2*a**2*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*sqrt(a - b)*sqrt(a + b)) + (A*b - B*a)*tan(c + d*x)/(b**2*d) - (2*A*a*b - 2*B*a**2 - B*b**2)*atanh(sin(c + d*x))/(2*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_309():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))
    F = B*tan(c + d*x)/(b*d) - 2*a*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*sqrt(a - b)*sqrt(a + b)) + (A*b - B*a)*atanh(sin(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_310():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))
    F = B*atanh(sin(c + d*x))/(b*d) + (2*A*b - 2*B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b*d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_311():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))
    F = A*x/a - (2*A*b - 2*B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_312():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))
    F = A*sin(c + d*x)/(a*d) + 2*b*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - x*(A*b - B*a)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_313():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))
    F = A*sin(c + d*x)*cos(c + d*x)/(2*a*d) - (A*b - B*a)*sin(c + d*x)/(a**2*d) - 2*b**2*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*sqrt(a - b)*sqrt(a + b)) + x*(A*a**2 + 2*A*b**2 - 2*B*a*b)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_314():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a + b*sec(c + d*x))
    F = A*sin(c + d*x)*cos(c + d*x)**2/(3*a*d) - (A*b - B*a)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + (2*A*a**2 + 3*A*b**2 - 3*B*a*b)*sin(c + d*x)/(3*a**3*d) + 2*b**3*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*sqrt(a - b)*sqrt(a + b)) - x*(a**2 + 2*b**2)*(A*b - B*a)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_315():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**4/(a + b*sec(c + d*x))
    F = A*sin(c + d*x)*cos(c + d*x)**3/(4*a*d) - (A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d) + (3*A*a**2 + 4*A*b**2 - 4*B*a*b)*sin(c + d*x)*cos(c + d*x)/(8*a**3*d) - (2*a**2 + 3*b**2)*(A*b - B*a)*sin(c + d*x)/(3*a**4*d) - 2*b**4*(A*b - B*a)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*sqrt(a - b)*sqrt(a + b)) + x*(3*A*a**4 + 4*A*a**2*b**2 + 8*A*b**4 - 4*B*a**3*b - 8*B*a*b**3)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_316():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = 2*a**2*(2*A*a**2*b - 3*A*b**3 - 3*B*a**3 + 4*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**2/(b*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - (2*A*a*b - 3*B*a**2 + B*b**2)*tan(c + d*x)*sec(c + d*x)/(2*b**2*d*(a**2 - b**2)) + (2*A*a**2*b - A*b**3 - 3*B*a**3 + 2*B*a*b**2)*tan(c + d*x)/(b**3*d*(a**2 - b**2)) - (4*A*a*b - 6*B*a**2 - B*b**2)*atanh(sin(c + d*x))/(2*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_317():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = B*tan(c + d*x)/(b**2*d) - a**2*(A*b - B*a)*tan(c + d*x)/(b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*a*(A*a**2*b - 2*A*b**3 - 2*B*a**3 + 3*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (A*b - 2*B*a)*atanh(sin(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_318():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = B*atanh(sin(c + d*x))/(b**2*d) + a*(A*b - B*a)*tan(c + d*x)/(b*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - (2*A*b**3 + 2*B*a**3 - 4*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_319():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))**2
    F = -(A*b - B*a)*tan(c + d*x)/(d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (2*A*a - 2*B*b)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_320():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**2
    F = A*x/a**2 + b*(A*b - B*a)*tan(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - (4*A*a**2*b - 2*A*b**3 - 2*B*a**3)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_321():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))**2
    F = b*(A*b - B*a)*sin(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (A*a**2 - 2*A*b**2 + B*a*b)*sin(c + d*x)/(a**2*d*(a**2 - b**2)) + 2*b*(3*A*a**2*b - 2*A*b**3 - 2*B*a**3 + B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - x*(2*A*b - B*a)/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_322():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (A*a**2 - 3*A*b**2 + 2*B*a*b)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d*(a**2 - b**2)) - (2*A*a**2*b - 3*A*b**3 - B*a**3 + 2*B*a*b**2)*sin(c + d*x)/(a**3*d*(a**2 - b**2)) - 2*b**2*(4*A*a**2*b - 3*A*b**3 - 3*B*a**3 + 2*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x*(A*a**2 + 6*A*b**2 - 4*B*a*b)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_323():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (A*a**2 - 4*A*b**2 + 3*B*a*b)*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(a**2 - b**2)) - (2*A*a**2*b - 4*A*b**3 - B*a**3 + 3*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d*(a**2 - b**2)) + (2*A*a**4 + 7*A*a**2*b**2 - 12*A*b**4 - 6*B*a**3*b + 9*B*a*b**3)*sin(c + d*x)/(3*a**4*d*(a**2 - b**2)) + 2*b**3*(5*A*a**2*b - 4*A*b**3 - 4*B*a**3 + 3*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - x*(2*A*a**2*b + 8*A*b**3 - B*a**3 - 6*B*a*b**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_324():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**5/(a + b*sec(c + d*x))**3
    F = a**2*(6*A*a**4*b - 15*A*a**2*b**3 + 12*A*b**5 - 12*B*a**5 + 29*B*a**3*b**2 - 20*B*a*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**3/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + a*(2*A*a**2*b - 5*A*b**3 - 4*B*a**3 + 7*B*a*b**2)*tan(c + d*x)*sec(c + d*x)**2/(2*b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - (3*A*a**3*b - 6*A*a*b**3 - 6*B*a**4 + 10*B*a**2*b**2 - B*b**4)*tan(c + d*x)*sec(c + d*x)/(2*b**3*d*(a**2 - b**2)**2) + (6*A*a**4*b - 11*A*a**2*b**3 + 2*A*b**5 - 12*B*a**5 + 21*B*a**3*b**2 - 6*B*a*b**4)*tan(c + d*x)/(2*b**4*d*(a**2 - b**2)**2) - (6*A*a*b - 12*B*a**2 - B*b**2)*atanh(sin(c + d*x))/(2*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_325():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a + b*sec(c + d*x))**3
    F = -a**2*(A*a**2*b - 4*A*b**3 - 3*B*a**3 + 6*B*a*b**2)*tan(c + d*x)/(2*b**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**2/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) - a*(2*A*a**4*b - 5*A*a**2*b**3 + 6*A*b**5 - 6*B*a**5 + 15*B*a**3*b**2 - 12*B*a*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - (A*a*b - 3*B*a**2 + 2*B*b**2)*tan(c + d*x)/(2*b**3*d*(a**2 - b**2)) + (A*b - 3*B*a)*atanh(sin(c + d*x))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_326():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))**3
    F = B*atanh(sin(c + d*x))/(b**3*d) - a**2*(A*b - B*a)*tan(c + d*x)/(2*b**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + a*(A*a**2*b - 4*A*b**3 - 3*B*a**3 + 6*B*a*b**2)*tan(c + d*x)/(2*b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (A*a**2*b**3 + 2*A*b**5 - 2*B*a**5 + 5*B*a**3*b**2 - 6*B*a*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_327():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = a*(A*b - B*a)*tan(c + d*x)/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) - (3*A*a*b - B*a**2 - 2*B*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (A*a**2*b + 2*A*b**3 + B*a**3 - 4*B*a*b**2)*tan(c + d*x)/(2*b*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_328():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))**3
    F = -(3*A*a*b - B*a**2 - 2*B*b**2)*tan(c + d*x)/(2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - (A*b - B*a)*tan(c + d*x)/(d*(a + b*sec(c + d*x))**2*(2*a**2 - 2*b**2)) + (2*A*a**2 + A*b**2 - 3*B*a*b)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_329():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**3
    F = A*x/a**3 + b*(A*b - B*a)*tan(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + b*(5*A*a**2*b - 2*A*b**3 - 3*B*a**3)*tan(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - (6*A*a**4*b - 5*A*a**2*b**3 + 2*A*b**5 - 2*B*a**5 - B*a**3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_330():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))**3
    F = b*(A*b - B*a)*sin(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + b*(6*A*a**2*b - 3*A*b**3 - 4*B*a**3 + B*a*b**2)*sin(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (2*A*a**4 - 11*A*a**2*b**2 + 6*A*b**4 + 5*B*a**3*b - 2*B*a*b**3)*sin(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) + b*(12*A*a**4*b - 15*A*a**2*b**3 + 6*A*b**5 - 6*B*a**5 + 5*B*a**3*b**2 - 2*B*a*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - x*(3*A*b - B*a)/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_331():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + b*(7*A*a**2*b - 4*A*b**3 - 5*B*a**3 + 2*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (A*a**4 - 10*A*a**2*b**2 + 6*A*b**4 + 6*B*a**3*b - 3*B*a*b**3)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) - (6*A*a**4*b - 21*A*a**2*b**3 + 12*A*b**5 - 2*B*a**5 + 11*B*a**3*b**2 - 6*B*a*b**4)*sin(c + d*x)/(2*a**4*d*(a**2 - b**2)**2) - b**2*(20*A*a**4*b - 29*A*a**2*b**3 + 12*A*b**5 - 12*B*a**5 + 15*B*a**3*b**2 - 6*B*a*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + x*(A*a**2 + 12*A*b**2 - 6*B*a*b)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_332():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**5/(a + b*sec(c + d*x))**4
    F = -a**2*(A*a**4*b - 2*A*a**2*b**3 + 6*A*b**5 - 4*B*a**5 + 11*B*a**3*b**2 - 12*B*a*b**4)*tan(c + d*x)/(2*b**4*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**3/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + a*(A*a**2*b - 6*A*b**3 - 4*B*a**3 + 9*B*a*b**2)*tan(c + d*x)*sec(c + d*x)**2/(6*b**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) - a*(2*A*a**6*b - 7*A*a**4*b**3 + 8*A*a**2*b**5 - 8*A*b**7 - 8*B*a**7 + 28*B*a**5*b**2 - 35*B*a**3*b**4 + 20*B*a*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - (3*A*a**3*b - 8*A*a*b**3 - 12*B*a**4 + 23*B*a**2*b**2 - 6*B*b**4)*tan(c + d*x)/(6*b**4*d*(a**2 - b**2)**2) + (A*b - 4*B*a)*atanh(sin(c + d*x))/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_333():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a + b*sec(c + d*x))**4
    F = B*atanh(sin(c + d*x))/(b**4*d) + a**2*(5*A*b**3 + 3*B*a**3 - 8*B*a*b**2)*tan(c + d*x)/(6*b**3*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**2/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) - a*(A*a**2*b**3 - 16*A*b**5 + 9*B*a**5 - 28*B*a**3*b**2 + 34*B*a*b**4)*tan(c + d*x)/(6*b**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - (3*A*a**2*b**5 + 2*A*b**7 + 2*B*a**7 - 7*B*a**5*b**2 + 8*B*a**3*b**4 - 8*B*a*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_334():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))**4
    F = -a**2*(A*b - B*a)*tan(c + d*x)/(3*b**2*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + a*(A*a**2*b - 6*A*b**3 - 4*B*a**3 + 9*B*a*b**2)*tan(c + d*x)/(6*b**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + (A*a**3 + 4*A*a*b**2 - 3*B*a**2*b - 2*B*b**3)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (A*a**4*b - 10*A*a**2*b**3 - 6*A*b**5 + 2*B*a**5 - 5*B*a**3*b**2 + 18*B*a*b**4)*tan(c + d*x)/(6*b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_335():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))**4
    F = a*(A*b - B*a)*tan(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) - (4*A*a**2*b + A*b**3 - B*a**3 - 4*B*a*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (2*A*a**3*b + 13*A*a*b**3 + B*a**4 - 10*B*a**2*b**2 - 6*B*b**4)*tan(c + d*x)/(6*b*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + (2*A*a**2*b + 3*A*b**3 + B*a**3 - 6*B*a*b**2)*tan(c + d*x)/(6*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_336():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))**4
    F = -(11*A*a**2*b + 4*A*b**3 - 2*B*a**3 - 13*B*a*b**2)*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - (5*A*a*b - 2*B*a**2 - 3*B*b**2)*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) - (A*b - B*a)*tan(c + d*x)/(d*(a + b*sec(c + d*x))**3*(3*a**2 - 3*b**2)) + (2*A*a**3 + 3*A*a*b**2 - 4*B*a**2*b - B*b**3)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_337():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**4
    F = A*x/a**4 + b*(A*b - B*a)*tan(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + b*(8*A*a**2*b - 3*A*b**3 - 5*B*a**3)*tan(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b*(26*A*a**4*b - 17*A*a**2*b**3 + 6*A*b**5 - 11*B*a**5 - 4*B*a**3*b**2)*tan(c + d*x)/(6*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - (8*A*a**6*b - 8*A*a**4*b**3 + 7*A*a**2*b**5 - 2*A*b**7 - 2*B*a**7 - 3*B*a**5*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_338():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))**4
    F = b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + b*(9*A*a**2*b - 4*A*b**3 - 6*B*a**3 + B*a*b**2)*sin(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b*(12*A*a**4*b - 11*A*a**2*b**3 + 4*A*b**5 - 6*B*a**5 + 2*B*a**3*b**2 - B*a*b**4)*sin(c + d*x)/(2*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + (6*A*a**6 - 65*A*a**4*b**2 + 68*A*a**2*b**4 - 24*A*b**6 + 26*B*a**5*b - 17*B*a**3*b**3 + 6*B*a*b**5)*sin(c + d*x)/(6*a**4*d*(a**2 - b**2)**3) + b*(20*A*a**6*b - 35*A*a**4*b**3 + 28*A*a**2*b**5 - 8*A*b**7 - 8*B*a**7 + 8*B*a**5*b**2 - 7*B*a**3*b**4 + 2*B*a*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - x*(4*A*b - B*a)/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_339():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))**4
    F = b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + b*(10*A*a**2*b - 5*A*b**3 - 7*B*a**3 + 2*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b*(48*A*a**4*b - 53*A*a**2*b**3 + 20*A*b**5 - 27*B*a**5 + 20*B*a**3*b**2 - 8*B*a*b**4)*sin(c + d*x)*cos(c + d*x)/(6*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + (A*a**6 - 23*A*a**4*b**2 + 27*A*a**2*b**4 - 10*A*b**6 + 12*B*a**5*b - 11*B*a**3*b**3 + 4*B*a*b**5)*sin(c + d*x)*cos(c + d*x)/(2*a**4*d*(a**2 - b**2)**3) - (24*A*a**6*b - 146*A*a**4*b**3 + 167*A*a**2*b**5 - 60*A*b**7 - 6*B*a**7 + 65*B*a**5*b**2 - 68*B*a**3*b**4 + 24*B*a*b**6)*sin(c + d*x)/(6*a**5*d*(a**2 - b**2)**3) - b**2*(40*A*a**6*b - 84*A*a**4*b**3 + 69*A*a**2*b**5 - 20*A*b**7 - 20*B*a**7 + 35*B*a**5*b**2 - 28*B*a**3*b**4 + 8*B*a*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**6*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + x*(A*a**2 + 20*A*b**2 - 8*B*a*b)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_340():
    f = (B*sec(c + d*x) + B*b/a)/(a + b*sec(c + d*x))
    F = B*b*x/a**2 + 2*B*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_341():
    f = (B*a/b + B*sec(c + d*x))/(a + b*sec(c + d*x))
    F = B*x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_342():
    f = (a + b*sec(c + d*x))/(a*sec(c + d*x) + b)**2
    F = -a*tan(c + d*x)/(b*d*(a*sec(c + d*x) + b)) + a*x/b**2 - 2*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_343():
    f = (sec(c + d*x) + 3)/(2 - sec(c + d*x))
    F = 3*x/2 - 5*sqrt(3)*log(-sqrt(3)*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(6*d) + 5*sqrt(3)*log(sqrt(3)*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_344():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sec(c + d*x)**4
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)**3/(9*d) + sqrt(a + b*sec(c + d*x))*(18*A*b + 2*B*a)*tan(c + d*x)*sec(c + d*x)**2/(63*b*d) + sqrt(a + b*sec(c + d*x))*(18*A*a*b - 12*B*a**2 + 98*B*b**2)*tan(c + d*x)*sec(c + d*x)/(315*b**2*d) - sqrt(a + b*sec(c + d*x))*(24*A*a**2*b - 150*A*b**3 - 16*B*a**3 - 26*B*a*b**2)*tan(c + d*x)/(315*b**3*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(-16*B*a**3 + 12*a**2*b*(2*A - B) + 18*a*b**2*(A - 2*B) + 3*b**3*(25*A - 49*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**4*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(24*A*a**3*b + 57*A*a*b**3 - 16*B*a**4 - 24*B*a**2*b**2 + 147*B*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_345():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sec(c + d*x)**3
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)**2/(7*d) + sqrt(a + b*sec(c + d*x))*(14*A*b + 2*B*a)*tan(c + d*x)*sec(c + d*x)/(35*b*d) + sqrt(a + b*sec(c + d*x))*(14*A*a*b - 8*B*a**2 + 50*B*b**2)*tan(c + d*x)/(105*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(-8*B*a**2 + 2*a*b*(7*A - 3*B) + b**2*(63*A - 25*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(14*A*a**2*b - 63*A*b**3 - 8*B*a**3 - 19*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_346():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sec(c + d*x)**2
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*b*d) + sqrt(a + b*sec(c + d*x))*(10*A*b - 4*B*a)*tan(c + d*x)/(15*b*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(5*A*b - 2*B*a - 9*B*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(5*A*a*b - 2*B*a**2 + 9*B*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_347():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sec(c + d*x)
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(3*A - B)*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(3*A*b + B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_348():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))
    F = ((Integer(-2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + ((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_349():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)
    F = ((Symbol('A') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('A') + (Integer(2) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_350():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)**2
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * (Symbol('A') + (Integer(2) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_351():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)**3
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(6) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a')) + Symbol('b')) * ((Integer(8) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * Symbol('b'))) + (Integer(6) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(6) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(6) * Symbol('a') * Symbol('B'))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_352():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)*sec(c + d*x)/(9*b*d) + (a + b*sec(c + d*x))**(sympy.S(5)/2)*(18*A*b - 8*B*a)*tan(c + d*x)/(63*b**2*d) - (a + b*sec(c + d*x))**(sympy.S(3)/2)*(36*A*a*b - 16*B*a**2 - 98*B*b**2)*tan(c + d*x)/(315*b**2*d) - sqrt(a + b*sec(c + d*x))*(36*A*a**2*b - 150*A*b**3 - 16*B*a**3 - 78*B*a*b**2)*tan(c + d*x)/(315*b**2*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*B*a**3 - 6*a**2*b*(3*A - B) - 3*a*b**2*(57*A - 13*B) + 3*b**3*(25*A - 49*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(18*A*a**3*b - 246*A*a*b**3 - 8*B*a**4 - 33*B*a**2*b**2 - 147*B*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_353():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*b*d) + (a + b*sec(c + d*x))**(sympy.S(3)/2)*(14*A*b - 4*B*a)*tan(c + d*x)/(35*b*d) + sqrt(a + b*sec(c + d*x))*(42*A*a*b - 12*B*a**2 + 50*B*b**2)*tan(c + d*x)/(105*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(6*B*a**2 - a*(21*A*b - 57*B*b) + b**2*(63*A - 25*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(21*A*a**2*b + 63*A*b**3 - 6*B*a**3 + 82*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_354():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) + sqrt(a + b*sec(c + d*x))*(10*A*b + 6*B*a)*tan(c + d*x)/(15*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(15*A*a - 5*A*b - 3*B*a + 9*B*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(20*A*a*b + 3*B*a**2 + 9*B*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_355():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(-2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('A')) + (Integer(-1) * Symbol('B')))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * ((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('B'))))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_356():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('b') * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Symbol('a') * (Symbol('A') + (Integer(4) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_357():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(5) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B')) + (Integer(8) * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_358():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(6) * Symbol('a') * Symbol('B'))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_359():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)*sec(c + d*x)/(11*b*d) + (a + b*sec(c + d*x))**(sympy.S(7)/2)*(22*A*b - 8*B*a)*tan(c + d*x)/(99*b**2*d) - (a + b*sec(c + d*x))**(sympy.S(5)/2)*(44*A*a*b - 16*B*a**2 - 162*B*b**2)*tan(c + d*x)/(693*b**2*d) - (a + b*sec(c + d*x))**(sympy.S(3)/2)*(220*A*a**2*b - 1078*A*b**3 - 80*B*a**3 - 670*B*a*b**2)*tan(c + d*x)/(3465*b**2*d) - sqrt(a + b*sec(c + d*x))*(220*A*a**3*b - 2508*A*a*b**3 - 80*B*a**4 - 570*B*a**2*b**2 - 1350*B*b**4)*tan(c + d*x)/(3465*b**2*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(40*B*a**4 - a**3*(110*A*b - 30*B*b) - 15*a**2*b**2*(121*A - 19*B) + 6*a*b**3*(209*A - 505*B) - 3*b**4*(539*A - 225*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3465*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(110*A*a**4*b - 3069*A*a**2*b**3 - 1617*A*b**5 - 40*B*a**5 - 255*B*a**3*b**2 - 3705*B*a*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3465*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_360():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(9*b*d) + (a + b*sec(c + d*x))**(sympy.S(5)/2)*(18*A*b - 4*B*a)*tan(c + d*x)/(63*b*d) + (a + b*sec(c + d*x))**(sympy.S(3)/2)*(90*A*a*b - 20*B*a**2 + 98*B*b**2)*tan(c + d*x)/(315*b*d) + sqrt(a + b*sec(c + d*x))*(90*A*a**2*b + 150*A*b**3 - 20*B*a**3 + 228*B*a*b**2)*tan(c + d*x)/(315*b*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(-10*B*a**3 + 15*a**2*b*(3*A - 11*B) - 6*a*b**2*(60*A - 19*B) + 3*b**3*(25*A - 49*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(45*A*a**3*b + 435*A*a*b**3 - 10*B*a**4 + 279*B*a**2*b**2 + 147*B*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_361():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = 2*B*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*d) + (a + b*sec(c + d*x))**(sympy.S(3)/2)*(14*A*b + 10*B*a)*tan(c + d*x)/(35*d) + sqrt(a + b*sec(c + d*x))*(112*A*a*b + 30*B*a**2 + 50*B*b**2)*tan(c + d*x)/(105*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(15*a**2*(7*A - B) - 8*a*b*(7*A - 15*B) + b**2*(63*A - 25*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(161*A*a**2*b + 63*A*b**3 + 15*B*a**3 + 145*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_362():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(-2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(35) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(23) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(9) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(45) * Symbol('A')) + (Integer(-1) * (Integer(23) * Symbol('B'))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(35) * Symbol('A')) + (Integer(-1) * (Integer(17) * Symbol('B')))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(5) * Symbol('A')) + (Integer(-1) * (Integer(9) * Symbol('B'))))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_363():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(6) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('b') * ((Integer(9) * Symbol('A')) + (Integer(-1) * (Integer(7) * Symbol('B'))))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('A')) + (Integer(-1) * Symbol('B'))))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(6) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_364():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(2) * Symbol('B')))) + (Integer(3) * Symbol('a') * Symbol('b') * ((Integer(3) * Symbol('A')) + (Integer(8) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(20) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_365():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(33) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(54) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(26) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(54) * Symbol('a') * Symbol('b') * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(20) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(30) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(33) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(54) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_366():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**4
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(284) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(128) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(264) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(192) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * ((Integer(9) * Symbol('A')) + (Integer(16) * Symbol('B')))) + (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(71) * Symbol('A')) + (Integer(52) * Symbol('B')))) + (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(59) * Symbol('A')) + (Integer(132) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(48) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(120) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(160) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(40) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(64) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(284) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(128) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(264) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((((Integer(36) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(59) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(104) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Integer(11) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('a') * Symbol('B'))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_367():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/sqrt(a + b*sec(c + d*x))
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)/(5*b*d) + sqrt(a + b*sec(c + d*x))*(10*A*b - 8*B*a)*tan(c + d*x)/(15*b**2*d) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(-8*B*a**2 + 2*a*b*(5*A + B) + b**2*(5*A - 9*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(10*A*a*b - 8*B*a**2 - 9*B*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_368():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*b*d) - 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(3*A*b - B*(2*a + b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(3*A*b - 2*B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_369():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = B*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(A - B)*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_370():
    f = (A + B*sec(c + d*x))/sqrt(a + b*sec(c + d*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_371():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = ((Symbol('A') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_372():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = ((Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('A') + (Integer(2) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_373():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/sqrt(a + b*sec(c + d*x))
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(18) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(10) * Symbol('a') * Symbol('A') * Symbol('b'))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(18) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(18) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('B')))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_374():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*b**2*d) - 2*a**2*(A*b - B*a)*tan(c + d*x)/(b**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*a + 2*b)*(3*A*b - B*(4*a + b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**3*d*sqrt(a + b)) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(12*A*a**2*b - 6*A*b**3 - 16*B*a**3 + 10*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**4*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_375():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(A*b - B*a)*tan(c + d*x)/(b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*A*b - 2*B*(2*a + b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*sqrt(a + b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*A*a*b - 4*B*a**2 + 2*B*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**3*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_376():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -(2*A*b - 2*B*a)*tan(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*A + 2*B)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*sqrt(a + b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*A*b + 2*B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_377():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_378():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = (((((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * Symbol('A') * Symbol('b')) + (Symbol('a') * (Symbol('A') + (Integer(-1) * (Integer(2) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_379():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(-1) * (((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Symbol('a') * Symbol('b') * ((Integer(5) * Symbol('A')) + (Integer(-1) * (Integer(12) * Symbol('B'))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(2) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_380():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((((Integer(16) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(41) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(-1) * (Integer(42) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Integer(90) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((((Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(7) * Symbol('A')) + (Integer(-1) * (Integer(18) * Symbol('B'))))) + (Integer(4) * (Symbol('a'))**(Integer(3)) * ((Integer(4) * Symbol('A')) + (Integer(3) * Symbol('B')))) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('A') + (Integer(5) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * (Symbol('a'))**(Integer(4)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(30) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('B')))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(16) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(41) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(-1) * (Integer(42) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Integer(90) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_381():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**4/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(3*A*a**2*b - 7*A*b**3 - 6*B*a**3 + 10*B*a*b**2)*tan(c + d*x)/(3*b**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + 2*a*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)**2/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - sqrt(a + b*sec(c + d*x))*(2*A*a*b - 4*B*a**2 + 2*B*b**2)*tan(c + d*x)/(3*b**3*d*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(32*B*a**4 - 2*a**3*(8*A*b - 12*B*b) - 4*a**2*b**2*(3*A + 8*B) + 18*a*b**3*(A - B) + 2*b**4*(3*A - B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**4*d*sqrt(a + b)*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-16*A*a**4*b + 30*A*a**2*b**3 - 6*A*b**5 + 32*B*a**5 - 56*B*a**3*b**2 + 16*B*a*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**5*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_382():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(A*b - B*a)*tan(c + d*x)/(3*b**2*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*a*(2*A*a**2*b - 6*A*b**3 - 5*B*a**3 + 9*B*a*b**2)*tan(c + d*x)/(3*b**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-16*B*a**3 + 4*a**2*b*(A - 3*B) + 6*a*b**2*(A + 3*B) - 6*b**3*(A - B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**3*d*sqrt(a + b)*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*A*a**3*b - 12*A*a*b**3 - 16*B*a**4 + 30*B*a**2*b**2 - 6*B*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**4*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_383():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*tan(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + (2*A*a**2*b + 6*A*b**3 + 4*B*a**3 - 12*B*a*b**2)*tan(c + d*x)/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*B*a**2 + 2*a*b*(A + 3*B) - 6*b**2*(A + B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**2*d*sqrt(a + b)*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*A*a**2*b + 6*A*b**3 + 4*B*a**3 - 12*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**3*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_384():
    f = (A + B*sec(c + d*x))*sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -(8*A*a*b - 2*B*a**2 - 6*B*b**2)*tan(c + d*x)/(3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - (2*A*b - 2*B*a)*tan(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(6*A*a - 2*A*b + 2*B*a - 6*B*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-8*A*a*b + 2*B*a**2 + 6*B*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_385():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_386():
    f = (A + B*sec(c + d*x))*cos(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(14) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(5) * Symbol('A')) + (Integer(-1) * (Integer(6) * Symbol('B'))))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('A') + (Integer(-1) * (Integer(4) * Symbol('B')))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(21) * Symbol('A')) + (Integer(2) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(14) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_387():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(-1) * (((Integer(33) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(104) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(60) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(7) * Symbol('A')) + (Integer(-1) * (Integer(12) * Symbol('B'))))) + (Integer(6) * (Symbol('a'))**(Integer(4)) * (Symbol('A') + (Integer(2) * Symbol('B')))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(27) * Symbol('A')) + (Integer(4) * Symbol('B'))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(27) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(84) * Symbol('b') * Symbol('B'))))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(20) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('A') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(27) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(33) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(104) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(60) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_388():
    f = (A*sec(e + f*x) + A)*sec(e + f*x)/sqrt(a + b*sec(e + f*x))
    F = -2*A*sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*cot(e + f*x)*elliptic_e(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_389():
    f = (-A*sec(e + f*x) + A)*sec(e + f*x)/sqrt(a + b*sec(e + f*x))
    F = 2*A*sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a - b)*(a + b)*cot(e + f*x)*elliptic_e(asin(sqrt(a + b*sec(e + f*x))/sqrt(a - b)), (a - b)/(a + b))/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_390():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + (10*A*a + 6*B*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (10*A*a + 6*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_391():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sqrt(sec(c + d*x))
    F = 2*B*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (6*A*a + 2*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (2*A*b + 2*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/d - (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_392():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/sqrt(sec(c + d*x))
    F = 2*B*b*sin(c + d*x)*sqrt(sec(c + d*x))/d + (2*A*a - 2*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_393():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (2*A*a + 6*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_394():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (6*A*a + 10*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_395():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + (10*A*a + 14*B*b)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (10*A*a + 14*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (6*A*b + 6*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_396():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*b*(7*A*b + 9*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + (10*B*b**2 + 14*a*(2*A*b + B*a))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (10*B*b**2 + 14*a*(2*A*b + B*a))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_397():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sqrt(sec(c + d*x))
    F = 2*B*b*(a + b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*b*(5*A*b + 7*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + (6*B*b**2 + 10*a*(2*A*b + B*a))*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (6*B*b**2 + 10*a*(2*A*b + B*a))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (6*A*a**2 + 2*A*b**2 + 4*B*a*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_398():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sqrt(sec(c + d*x))
    F = 2*B*b*(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*b*(3*A*b + 5*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (12*A*a*b + 6*B*a**2 + 2*B*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_399():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*B*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/d - (2*B*b**2 - 2*a*(2*A*b + B*a))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*A*a**2 + 6*A*b**2 + 12*B*a*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_400():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(2*A*b + B*a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (6*A*a**2 + 10*A*b**2 + 20*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (4*A*a*b + 2*B*a**2 + 6*B*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_401():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(2*A*b + B*a)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (10*A*a**2 + 14*A*b**2 + 28*B*a*b)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (10*A*a**2 + 14*A*b**2 + 28*B*a*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (12*A*a*b + 6*B*a**2 + 10*B*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_402():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a*(2*A*b + B*a)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + (14*B*b**2 + 10*a*(2*A*b + B*a))*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (14*B*b**2 + 10*a*(2*A*b + B*a))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (14*A*a**2 + 18*A*b**2 + 36*B*a*b)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + (14*A*a**2 + 18*A*b**2 + 36*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_403():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*sec(c + d*x))**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(9*d) + 2*b**2*(9*A*b + 13*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d) + 2*b*(27*A*a*b + 22*B*a**2 + 7*B*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(45*d) + (30*A*a**3 + 54*A*a*b**2 + 54*B*a**2*b + 14*B*b**3)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) - (30*A*a**3 + 54*A*a*b**2 + 54*B*a**2*b + 14*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_404():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*sqrt(sec(c + d*x))
    F = 2*B*b*(a + b*sec(c + d*x))**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*b**2*(7*A*b + 11*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*b*(21*A*a*b + 18*B*a**2 + 5*B*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (42*A*a**3 + 42*A*a*b**2 + 42*B*a**2*b + 10*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (30*A*a**2*b + 6*A*b**3 + 10*B*a**3 + 18*B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (30*A*a**2*b + 6*A*b**3 + 10*B*a**3 + 18*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_405():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sqrt(sec(c + d*x))
    F = 2*B*b*(a + b*sec(c + d*x))**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(5*A*b + 9*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + 2*b*(15*A*a*b + 14*B*a**2 + 3*B*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + (10*A*a**3 - 30*A*a*b**2 - 30*B*a**2*b - 6*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (18*A*a**2*b + 2*A*b**3 + 6*B*a**3 + 6*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_406():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) - 2*b**2*(A*a - B*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) - 2*b*(2*A*a**2 - 3*A*b**2 - 9*B*a*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + (2*A*a**3 + 18*A*a*b**2 + 18*B*a**2*b + 2*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_407():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(9*A*b + 5*B*a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) - 2*b**2*(A*a - 5*B*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + (6*A*a**3 + 30*A*a*b**2 + 30*B*a**2*b - 10*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (6*A*a**2*b + 6*A*b**3 + 2*B*a**3 + 18*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_408():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(11*A*b + 7*B*a)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(5*A*a**2 + 18*A*b**2 + 21*B*a*b)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (10*A*a**3 + 42*A*a*b**2 + 42*B*a**2*b + 42*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (18*A*a**2*b + 10*A*b**3 + 6*B*a**3 + 30*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_409():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**2*(13*A*b + 9*B*a)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(7*A*a**2 + 22*A*b**2 + 27*B*a*b)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + (14*A*a**3 + 54*A*a*b**2 + 54*B*a**2*b + 30*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + (30*A*a**2*b + 14*A*b**3 + 10*B*a**3 + 42*B*a*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (30*A*a**2*b + 14*A*b**3 + 10*B*a**3 + 42*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_410():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**2*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 2*a**2*(15*A*b + 11*B*a)*sin(c + d*x)/(99*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a*(9*A*a**2 + 26*A*b**2 + 33*B*a*b)*sin(c + d*x)/(77*d*sec(c + d*x)**(sympy.S(5)/2)) + (90*A*a**3 + 330*A*a*b**2 + 330*B*a**2*b + 154*B*b**3)*sin(c + d*x)/(231*d*sqrt(sec(c + d*x))) + (90*A*a**3 + 330*A*a*b**2 + 330*B*a**2*b + 154*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(231*d) + (42*A*a**2*b + 18*A*b**3 + 14*B*a**3 + 54*B*a*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + (42*A*a**2*b + 18*A*b**3 + 14*B*a**3 + 54*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_411():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))
    F = ((Integer(2) * ((Integer(5) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_412():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_413():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_414():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a + b*sec(c + d*x))
    F = ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_415():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_416():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_417():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('a') * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_418():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_419():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**2
    F = ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_420():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_421():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**2
    F = ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_422():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*sqrt(sec(c + d*x)))
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_423():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_424():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)/(a + b*sec(c + d*x))**3
    F = (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(65) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(86) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(63) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(65) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(13) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_425():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**3
    F = ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(35) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_426():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**3
    F = (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_427():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**3
    F = (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(5)))) + ((Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_428():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**3
    F = ((((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_429():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*sqrt(sec(c + d*x)))
    F = ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_430():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(6)) * Symbol('A')) + (Integer(128) * (Symbol('a'))**(Integer(4)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(223) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(6))) + (Integer(-1) * (Integer(72) * (Symbol('a'))**(Integer(5)) * Symbol('b') * Symbol('B'))) + (Integer(99) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(45) * Symbol('a') * (Symbol('b'))**(Integer(5)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(12) * (Symbol('a'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(63) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(13) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_431():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_432():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sqrt(sec(c + d*x))
    F = ((((Integer(2) * Symbol('a') * Symbol('A')) + (Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('B') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_433():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/sqrt(sec(c + d*x))
    F = ((Integer(2) * Symbol('a') * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_434():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*A*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*b + 6*B*a)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_435():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(2*A*b + 10*B*a)*sin(c + d*x)/(15*a*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*a**2 + 2*b**2)*(2*A*b - 5*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 - 4*A*b**2 + 10*B*a*b)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_436():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(2*A*b + 14*B*a)*sin(c + d*x)/(35*a*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 - 8*A*b**2 + 14*B*a*b)*sin(c + d*x)/(105*a**2*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 + 8*A*b**2 - 14*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(105*a**3*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(38*A*a**2*b + 16*A*b**3 + 126*B*a**3 - 28*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_437():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(42) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(17) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(8) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_438():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(7) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(12) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_439():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_440():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_441():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(12*A*b + 10*B*a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(3*A*b + 5*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 + 6*A*b**2 + 40*B*a*b)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_442():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(16*A*b + 14*B*a)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 + 6*A*b**2 + 84*B*a*b)*sin(c + d*x)/(105*a*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 - 6*A*b**2 + 21*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(105*a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(164*A*a**2*b - 12*A*b**3 + 126*B*a**3 + 42*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_443():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*sec(c + d*x))*(20*A*b + 18*B*a)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(98*A*a**2 + 6*A*b**2 + 144*B*a*b)*sin(c + d*x)/(315*a*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(176*A*a**2*b - 8*A*b**3 + 150*B*a**3 + 18*B*a*b**2)*sin(c + d*x)/(315*a**2*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(39*A*a**2*b + 8*A*b**3 + 75*B*a**3 - 18*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(315*a**3*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*A*a**4 + 66*A*a**2*b**2 + 16*A*b**4 + 492*B*a**3*b - 36*B*a*b**3)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_444():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(472) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(133) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(356) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(192) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(40) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(160) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(64) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((((Integer(104) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(59) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(36) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(11) * Symbol('a') * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_445():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = ((((Integer(48) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(66) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(59) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(30) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(8) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_446():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(20) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_447():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(4) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_448():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(10) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(15) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(23) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(35) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(15) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_449():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*sqrt(a + b*sec(c + d*x))*(10*A*b + 7*B*a)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 + 90*A*b**2 + 154*B*a*b)*sin(c + d*x)/(105*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 + 15*A*b**2 + 56*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(105*a*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(290*A*a**2*b + 30*A*b**3 + 126*B*a**3 + 322*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_450():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a*sqrt(a + b*sec(c + d*x))*(4*A*b + 3*B*a)*sin(c + d*x)/(21*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(98*A*a**2 + 150*A*b**2 + 270*B*a*b)*sin(c + d*x)/(315*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(326*A*a**2*b + 10*A*b**3 + 150*B*a**3 + 270*B*a*b**2)*sin(c + d*x)/(315*a*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(114*A*a**2*b - 10*A*b**3 + 75*B*a**3 + 45*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(315*a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*A*a**4 + 558*A*a**2*b**2 - 20*A*b**4 + 870*B*a**3*b + 90*B*a*b**3)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_451():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 2*a*sqrt(a + b*sec(c + d*x))*(14*A*b + 11*B*a)*sin(c + d*x)/(99*d*sec(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*sec(c + d*x))*(162*A*a**2 + 226*A*b**2 + 418*B*a*b)*sin(c + d*x)/(693*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(2290*A*a**2*b + 30*A*b**3 + 1078*B*a**3 + 1650*B*a*b**2)*sin(c + d*x)/(3465*a*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(1350*A*a**4 + 2050*A*a**2*b**2 - 40*A*b**4 + 3586*B*a**3*b + 110*B*a*b**3)*sin(c + d*x)/(3465*a**2*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(675*A*a**4 + 285*A*a**2*b**2 + 40*A*b**4 + 1254*B*a**3*b - 110*B*a*b**3)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3465*a**3*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(7410*A*a**4*b + 510*A*a**2*b**3 + 80*A*b**5 + 3234*B*a**5 + 6138*B*a**3*b**2 - 220*B*a*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3465*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_452():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*sec(c + d*x))
    F = ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('B') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_453():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*sec(c + d*x))
    F = ((Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('B') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_454():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/sqrt(a + b*sec(c + d*x))
    F = ((Integer(2) * Symbol('A') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_455():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = 2*A*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_456():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**2 + 4*A*b**2 - 6*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))) - sqrt(a + b*sec(c + d*x))*(4*A*b - 6*B*a)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_457():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*sec(c + d*x))*(8*A*b - 10*B*a)*sin(c + d*x)/(15*a**2*d*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-14*A*a**2*b - 16*A*b**3 + 10*B*a**3 + 20*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**3*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 + 16*A*b**2 - 20*B*a*b)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_458():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_459():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_460():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))) - (2*A*b - 2*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(2*A*b - 2*B*a)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_461():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-4*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 4*A*b**2 + 2*B*a*b)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_462():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 8*A*b**2 + 6*B*a*b)*sin(c + d*x)/(3*a**2*d*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**2 + 16*A*b**2 - 12*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**3*d*sqrt(a + b*sec(c + d*x))) - sqrt(a + b*sec(c + d*x))*(10*A*a**2*b - 16*A*b**3 - 6*B*a**3 + 12*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_463():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 12*A*b**2 + 10*B*a*b)*sin(c + d*x)/(5*a**2*d*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*sec(c + d*x))*(18*A*a**2*b - 48*A*b**3 - 10*B*a**3 + 40*B*a*b**2)*sin(c + d*x)/(15*a**3*d*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-24*A*a**2*b - 96*A*b**3 + 10*B*a**3 + 80*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**4*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**4 + 48*A*a**2*b**2 - 96*A*b**4 - 50*B*a**3*b + 80*B*a*b**3)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_464():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_465():
    f = (A + B*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + (4*A*a**2*b + 4*A*b**3 + 2*B*a**3 - 10*B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - sqrt(a + b*sec(c + d*x))*(6*A*a**2 + 2*A*b**2 - 8*B*a*b)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_466():
    f = (A + B*sec(c + d*x))*sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -(2*A*b - 2*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) - (10*A*a**2*b - 2*A*b**3 - 4*B*a**3 - 4*B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(6*A*a**2 - 4*A*b**2 - 2*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_467():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*b*(8*A*a**2*b - 4*A*b**3 - 5*B*a**3 + B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-18*A*a**2*b + 16*A*b**3 + 6*B*a**3 - 4*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(6*A*a**4 - 30*A*a**2*b**2 + 16*A*b**4 + 12*B*a**3*b - 4*B*a*b**3)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_468():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(sec(c + d*x))) + 2*b*(10*A*a**2*b - 6*A*b**3 - 7*B*a**3 + 3*B*a*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*a**4 - 26*A*a**2*b**2 + 16*A*b**4 + 16*B*a**3*b - 8*B*a*b**3)*sin(c + d*x)/(3*a**3*d*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**4 + 32*A*a**2*b**2 - 32*A*b**4 - 18*B*a**3*b + 16*B*a*b**3)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**4*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - sqrt(a + b*sec(c + d*x))*(16*A*a**4*b - 56*A*a**2*b**3 + 32*A*b**5 - 6*B*a**5 + 30*B*a**3*b**2 - 16*B*a*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_469():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*(12*A*a**2*b - 8*A*b**3 - 9*B*a**3 + 5*B*a*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(6*A*a**4 - 142*A*a**2*b**2 + 96*A*b**4 + 100*B*a**3*b - 60*B*a*b**3)*sin(c + d*x)/(15*a**3*d*(a**2 - b**2)**2*sec(c + d*x)**(sympy.S(3)/2)) - sqrt(a + b*sec(c + d*x))*(28*A*a**4*b - 196*A*a**2*b**3 + 128*A*b**5 - 10*B*a**5 + 130*B*a**3*b**2 - 80*B*a*b**4)*sin(c + d*x)/(15*a**4*d*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-34*A*a**4*b - 232*A*a**2*b**3 + 256*A*b**5 + 10*B*a**5 + 160*B*a**3*b**2 - 160*B*a*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**5*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(18*A*a**6 + 110*A*a**4*b**2 - 424*A*a**2*b**4 + 256*A*b**6 - 80*B*a**5*b + 280*B*a**3*b**3 - 160*B*a*b**5)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**5*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_470():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(2)/3)
    F = ((sympy.sqrt(Integer(2)) * Symbol('B') * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Symbol('A') * sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_471():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = ((sympy.sqrt(Integer(2)) * Symbol('B') * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Symbol('A') * sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_472():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = ((sympy.sqrt(Integer(2)) * Symbol('B') * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(3))**(Integer(-1)), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**((Integer(3))**(Integer(-1))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Symbol('A') * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))))**(Integer(-1)), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_473():
    f = (A + B*sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(2)/3)
    F = ((sympy.sqrt(Integer(2)) * Symbol('B') * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(2) * (Integer(3))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Symbol('A') * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1)), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_474():
    f = (c*sec(e + f*x))**n*(A + B*sec(e + f*x))*(a + b*sec(e + f*x))**m
    F = sympy.Function('Unintegrable')((((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (Symbol('A') + (Symbol('B') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_475():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**4*sec(c + d*x)**m
    F = B*b*(a + b*sec(c + d*x))**3*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 4)) + b**2*(2*A*a*b*(m + 4)**2 + B*a**2*(m**2 + 9*m + 26) + B*b**2*(m + 3)**2)*sin(c + d*x)*sec(c + d*x)**(m + 2)/(d*(m + 2)*(m + 3)*(m + 4)) + b*(a + b*sec(c + d*x))**2*(A*b*(m + 4) + B*a*(m + 7))*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 3)*(m + 4)) + b*(A*a**2*b*(5*m**2 + 37*m + 68) + A*b**3*(m**2 + 6*m + 8) + 2*B*a**3*(m**2 + 8*m + 19) + 4*B*a*b**2*(m**2 + 6*m + 8))*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 1)*(m + 3)*(m + 4)) - (A*a**4*(m**2 + 4*m + 3) + 6*A*a**2*b**2*m*(m + 3) + A*b**4*m*(m + 2) + 4*B*a**3*b*m*(m + 3) + 4*B*a*b**3*m*(m + 2))*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - m)*(m + 1)*(m + 3)*sqrt(sin(c + d*x)**2)) + (4*A*a**3*b*(m**2 + 6*m + 8) + 4*A*a*b**3*(m**2 + 5*m + 4) + B*a**4*(m**2 + 6*m + 8) + 6*B*a**2*b**2*(m**2 + 5*m + 4) + B*b**4*(m**2 + 4*m + 3))*sin(c + d*x)*hyper((sympy.S.Half, -m/2), (1 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*m*(m + 2)*(m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_476():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**3*sec(c + d*x)**m
    F = B*b*(a + b*sec(c + d*x))**2*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 3)) + b**2*(A*b*(m + 3) + B*a*(m + 5))*sin(c + d*x)*sec(c + d*x)**(m + 2)/(d*(m + 2)*(m + 3)) + b*(3*A*a*b*(m + 3) + 2*B*a**2*(m + 4) + B*b**2*(m + 2))*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 1)*(m + 3)) - (A*a**3*(m**2 + 4*m + 3) + 3*A*a*b**2*m*(m + 3) + 3*B*a**2*b*m*(m + 3) + B*b**3*m*(m + 2))*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - m**2)*(m + 3)*sqrt(sin(c + d*x)**2)) + (3*A*a**2*b*(m + 2) + A*b**3*(m + 1) + B*a**3*(m + 2) + 3*B*a*b**2*(m + 1))*sin(c + d*x)*hyper((sympy.S.Half, -m/2), (1 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*m*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_477():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sec(c + d*x)**m
    F = B*b*(a + b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 2)) + b*(A*b*(m + 2) + B*a*(m + 3))*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 1)*(m + 2)) - (A*a**2*(m + 1) + A*b**2*m + 2*B*a*b*m)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - m**2)*sqrt(sin(c + d*x)**2)) + (B*b**2*(m + 1) + a*(m + 2)*(2*A*b + B*a))*sin(c + d*x)*hyper((sympy.S.Half, -m/2), (1 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*m*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_478():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sec(c + d*x)**m
    F = B*b*sin(c + d*x)*sec(c + d*x)**(m + 1)/(d*(m + 1)) - (A*a*(m + 1) + B*b*m)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - m**2)*sqrt(sin(c + d*x)**2)) + (A*b + B*a)*sin(c + d*x)*hyper((sympy.S.Half, -m/2), (1 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*m*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_479():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*(A + B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*(5*A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*a*(5*A + 7*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_480():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*(A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a*(3*A + 5*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_481():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_482():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*B*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a*(A - B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_483():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*B*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(3*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_484():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*(A + B)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a*(5*A + 3*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 2*a*(5*A + 3*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_485():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(9*d) + 4*a**2*(5*A + 6*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 4*a**2*(5*A + 6*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**2*(8*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 4*a**2*(8*A + 9*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 2*a**2*(11*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_486():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 4*a**2*(3*A + 4*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*(6*A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 4*a**2*(6*A + 7*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*a**2*(9*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_487():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 4*a**2*(A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 4*a**2*(4*A + 5*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a**2*(7*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_488():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 4*A*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**2*(A - 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a**2*(2*A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_489():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))
    F = -4*B*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(3*A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a**2*(3*A + 5*B)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_490():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/sqrt(cos(c + d*x))
    F = 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(2*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 4*a**2*(5*A + 4*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 4*a**2*(5*A + 4*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a**2*(5*A + 7*B)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_491():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 4*a**2*(4*A + 3*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 4*a**2*(4*A + 3*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*(7*A + 6*B)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(7*A + 6*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*a**2*(7*A + 9*B)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_492():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(d*(a*cos(c + d*x) + a)) - (5*A - 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) - (5*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d) + (7*A - 5*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) + (21*A - 15*B)*elliptic_e(c/2 + d*x/2, 2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_493():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(a*cos(c + d*x) + a)) - (3*A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (5*A - 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) + (5*A - 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_494():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) - (A - B)*elliptic_f(c/2 + d*x/2, 2)/(a*d) + (3*A - B)*elliptic_e(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_495():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) - (A - B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (A + B)*elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_496():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - (A - 3*B)*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) + (A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (A - B)*elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_497():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (3*A - 5*B)*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) - (3*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d) + (3*A - 3*B)*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) - (3*A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_498():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(3*d*(a*cos(c + d*x) + a)**2) - (3*A - 2*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(a**2*d*(cos(c + d*x) + 1)) - (15*A - 10*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) - (15*A - 10*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + (56*A - 35*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a**2*d) + (56*A - 35*B)*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_499():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*(a*cos(c + d*x) + a)**2) - (7*A - 4*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - (7*A - 4*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(cos(c + d*x) + 1)) + (10*A - 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) + (10*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_500():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d*(a*cos(c + d*x) + a)**2) + (4*A - B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - (5*A - 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - (5*A - 2*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_501():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x)))
    F = -A*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + A*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1)) - (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + (2*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_502():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2))
    F = B*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - B*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1)) + (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + (A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_503():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) - (A - 4*B)*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) + (A - 4*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + (2*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + (2*A - 5*B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_504():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) + (4*A - 7*B)*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) - (4*A - 7*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + (4*A - 7*B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2)) - (5*A - 10*B)*sin(c + d*x)/(3*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) - (5*A - 10*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_505():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d*(a*cos(c + d*x) + a)**3) - (119*A - 49*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(30*d*(a**3*cos(c + d*x) + a**3)) - (2*A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*a*d*(a*cos(c + d*x) + a)**2) + (33*A - 13*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*a**3*d) + (33*A - 13*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (119*A - 49*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_506():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(5*d*(a*cos(c + d*x) + a)**3) - (13*A - 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a**3*cos(c + d*x) + a**3)) - (8*A - 3*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a*d*(a*cos(c + d*x) + a)**2) - (13*A - 3*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) + (49*A - 9*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_507():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*(a*cos(c + d*x) + a)**3) + (9*A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - (6*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + (3*A + B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (9*A + B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_508():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) + (4*A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) - (A - B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + (A + B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_509():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) - (A + 9*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) + (A - 6*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + (A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) + (A + 9*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_510():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x))) + (3*A - 13*B)*sin(c + d*x)/(6*d*(a**3*cos(c + d*x) + a**3)*sqrt(cos(c + d*x))) + (3*A - 8*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) + (3*A - 13*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (9*A - 49*B)*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) + (9*A - 49*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_511():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(8*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d*sqrt(a*sec(c + d*x) + a)) + 4*a*(8*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*sec(c + d*x) + a)) + 16*a*(8*A + 9*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 32*a*(8*A + 9*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_512():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(6*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*sec(c + d*x) + a)) + 8*a*(6*A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 16*a*(6*A + 7*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_513():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(a*sec(c + d*x) + a)) + 2*a*(4*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 4*a*(4*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_514():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_515():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*A*a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*B*sqrt(a)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_516():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/sqrt(cos(c + d*x))
    F = B*a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a)*(2*A + B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_517():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*sin(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a)*(4*A + 3*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + a*(4*A + 3*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_518():
    f = (A + B*sec(c + d*x))*sqrt(a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = B*a*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a)*(6*A + 5*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a*(6*A + 5*B)*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a*(6*A + 5*B)*sin(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_519():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(11*d) + 2*a**2*(12*A + 11*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(99*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(168*A + 187*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(693*d*sqrt(a*sec(c + d*x) + a)) + 4*a**2*(168*A + 187*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(1155*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*(168*A + 187*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3465*d*sqrt(a*sec(c + d*x) + a)) + 32*a**2*(168*A + 187*B)*sin(c + d*x)/(3465*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_520():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*a**2*(10*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(34*A + 39*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*sec(c + d*x) + a)) + 8*a**2*(34*A + 39*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*(34*A + 39*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_521():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a**2*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*(52*A + 63*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 4*a**2*(52*A + 63*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_522():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 8*a**2*(3*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*(3*A + 5*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_523():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*B*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**2*(4*A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_524():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + a**(sympy.S(3)/2)*(2*A + 3*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + a**2*(2*A - B)*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_525():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(3)/2)*(12*A + 7*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + a**2*(4*A + 5*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_526():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(3)/2)*(14*A + 11*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a**2*(6*A + 7*B)*sin(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + a**2*(14*A + 11*B)*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_527():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = B*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(7)/2)) + a**(sympy.S(3)/2)*(88*A + 75*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + a**2*(8*A + 9*B)*sin(c + d*x)/(24*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + a**2*(88*A + 75*B)*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**2*(88*A + 75*B)*sin(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_528():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(11*d) + 2*a**3*(194*A + 209*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(693*d*sqrt(a*sec(c + d*x) + a)) + 2*a**3*(710*A + 803*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(1155*d*sqrt(a*sec(c + d*x) + a)) + 8*a**3*(710*A + 803*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3465*d*sqrt(a*sec(c + d*x) + a)) + 16*a**3*(710*A + 803*B)*sin(c + d*x)/(3465*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(14*A + 11*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(99*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_529():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*a**3*(124*A + 135*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*sec(c + d*x) + a)) + 2*a**3*(292*A + 345*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 4*a**3*(292*A + 345*B)*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(4*A + 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_530():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 64*a**3*(5*A + 7*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 16*a**2*(5*A + 7*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d) + 2*a*(5*A + 7*B)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_531():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*B*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**3*(32*A + 35*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(8*A + 5*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_532():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + a**(sympy.S(5)/2)*(2*A + 5*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + a**3*(14*A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - a**2*(2*A - 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_533():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(2*d*sqrt(cos(c + d*x))) + a**(sympy.S(5)/2)*(20*A + 19*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + a**3*(4*A - 9*B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + a**2*(4*A + 7*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_534():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(5)/2)*(38*A + 25*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a**3*(54*A + 49*B)*sin(c + d*x)/(24*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**2*(2*A + 3*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_535():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(5)/2)*(200*A + 163*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + a**3*(104*A + 95*B)*sin(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + a**3*(200*A + 163*B)*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**2*(8*A + 11*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(24*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_536():
    f = (A + B*sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = B*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(7)/2)) + a**(sympy.S(5)/2)*(326*A + 283*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(128*d) + a**3*(170*A + 157*B)*sin(c + d*x)/(240*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + a**3*(326*A + 283*B)*sin(c + d*x)/(128*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**3*(326*A + 283*B)*sin(c + d*x)/(192*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + a**2*(10*A + 13*B)*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(40*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_537():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*A*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d*sqrt(a*sec(c + d*x) + a)) - (2*A - 14*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*sec(c + d*x) + a)) + (62*A - 14*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) - (86*A - 182*B)*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_538():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*A*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(a*sec(c + d*x) + a)) - (2*A - 10*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + (26*A - 10*B)*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_539():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) - (2*A - 6*B)*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_540():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/sqrt(a*sec(c + d*x) + a)
    F = 2*A*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_541():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = 2*B*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_542():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = B*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d) + (2*A - B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_543():
    f = (A + B*sec(c + d*x))/(sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = B*sin(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + (4*A - B)*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d) - (4*A - 7*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_544():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (9*A - 5*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(10*a*d*sqrt(a*sec(c + d*x) + a)) - (39*A - 35*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(30*a*d*sqrt(a*sec(c + d*x) + a)) + (147*A - 95*B)*sin(c + d*x)/(30*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(15*A - 11*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_545():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (7*A - 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*a*d*sqrt(a*sec(c + d*x) + a)) - (19*A - 15*B)*sin(c + d*x)/(6*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(11*A - 7*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_546():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + (5*A - B)*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(7*A - 3*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_547():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(3*A + B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_548():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*B*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) + (A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(A - 5*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_549():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)) - (A - 3*B)*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + (2*A - 3*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) - sqrt(2)*(5*A - 9*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_550():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2)) - (A - 2*B)*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + (6*A - 7*B)*sin(c + d*x)/(4*a*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(9*A - 13*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d) - (12*A - 19*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_551():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (21*A - 13*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (157*A - 85*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(80*a**2*d*sqrt(a*sec(c + d*x) + a)) - (787*A - 475*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(240*a**2*d*sqrt(a*sec(c + d*x) + a)) + (2671*A - 1495*B)*sin(c + d*x)/(240*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(283*A - 163*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_552():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - (17*A - 9*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + (95*A - 39*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(48*a**2*d*sqrt(a*sec(c + d*x) + a)) - (299*A - 147*B)*sin(c + d*x)/(48*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(163*A - 75*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_553():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - (13*A - 5*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + (49*A - 9*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(75*A - 19*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_554():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = (A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) + (5*A + 3*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) - (9*A - B)*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(19*A + 5*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_555():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)) + (5*A + 3*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(5*A + 3*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_556():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*B*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) + (A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)) + (3*A - 11*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 43*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_557():
    f = (A + B*sec(c + d*x))/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)) + (7*A - 15*B)*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)) - (11*A - 35*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + (2*A - 5*B)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) - sqrt(2)*(43*A - 115*B)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_558():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + (10*A*a + 14*B*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (10*A*a + 14*B*b)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (2*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + (6*A*b + 6*B*a)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_559():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + (6*A*a + 10*B*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_560():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (2*A*a + 6*B*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (2*A*b + 2*B*a)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_561():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))*sqrt(cos(c + d*x))
    F = 2*B*b*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + (2*A*a - 2*B*b)*elliptic_e(c/2 + d*x/2, 2)/d + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_562():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/sqrt(cos(c + d*x))
    F = 2*B*b*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (6*A*a + 2*B*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (-2*A*b - 2*B*a)*elliptic_e(c/2 + d*x/2, 2)/d + (2*A*b + 2*B*a)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_563():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + (-10*A*a - 6*B*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (10*A*a + 6*B*b)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + (2*A*b + 2*B*a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_564():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*cos(c + d*x) + b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*a*(9*A*b + 7*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + (10*A*a**2 + 14*b*(A*b + 2*B*a))*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (10*A*a**2 + 14*b*(A*b + 2*B*a))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (12*A*a*b + 6*B*a**2 + 10*B*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_565():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*cos(c + d*x) + b)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*a*(7*A*b + 5*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + (6*A*a**2 + 10*b*(A*b + 2*B*a))*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (4*A*a*b + 2*B*a**2 + 6*B*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_566():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*B*b**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + (2*A*a**2 + 6*A*b**2 + 12*B*a*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (4*A*a*b + 2*B*a**2 - 2*B*b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_567():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2*sqrt(cos(c + d*x))
    F = 2*B*b**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*b*(A*b + 2*B*a)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/d + (12*A*a*b + 6*B*a**2 + 2*B*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_568():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/sqrt(cos(c + d*x))
    F = 2*B*b**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*b*(A*b + 2*B*a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (6*A*a**2 + 2*A*b**2 + 4*B*a*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (20*A*a*b + 10*B*a**2 + 6*B*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (20*A*a*b + 10*B*a**2 + 6*B*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_569():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b**2*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*b*(A*b + 2*B*a)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (28*A*a*b + 14*B*a**2 + 10*B*b**2)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + (28*A*a*b + 14*B*a**2 + 10*B*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_570():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_571():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = ((Integer(-2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_572():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a + b*sec(c + d*x))
    F = ((Integer(2) * Symbol('A') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_573():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    F = ((Integer(2) * Symbol('A') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_574():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(-2) * Symbol('B') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_575():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(-2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_576():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(2) * ((Integer(5) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_577():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_578():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**2
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_579():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*sqrt(cos(c + d*x)))
    F = ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_580():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_581():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2))
    F = ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_582():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(7)/2))
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_583():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**3
    F = (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(6)) * Symbol('A')) + (Integer(128) * (Symbol('a'))**(Integer(4)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(223) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(6))) + (Integer(-1) * (Integer(72) * (Symbol('a'))**(Integer(5)) * Symbol('b') * Symbol('B'))) + (Integer(99) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(45) * Symbol('a') * (Symbol('b'))**(Integer(5)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('a'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(63) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(13) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_584():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**3
    F = ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_585():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*sqrt(cos(c + d*x)))
    F = ((((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_586():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(5)))) + ((Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_587():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(5)/2))
    F = (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_588():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(7)/2))
    F = ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(35) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_589():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(9)/2))
    F = (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(65) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(86) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(63) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(65) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(13) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_590():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + sqrt(a + b*sec(c + d*x))*(2*A*b + 14*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*a*d) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 - 8*A*b**2 + 14*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*a**2*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 + 8*A*b**2 - 14*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(105*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(38*A*a**2*b + 16*A*b**3 + 126*B*a**3 - 28*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_591():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + sqrt(a + b*sec(c + d*x))*(2*A*b + 10*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*a**2 + 2*b**2)*(2*A*b - 5*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 - 4*A*b**2 + 10*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_592():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*A*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*b + 6*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_593():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))
    F = ((Integer(2) * Symbol('a') * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_594():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/sqrt(cos(c + d*x))
    F = ((((Integer(2) * Symbol('a') * Symbol('A')) + (Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_595():
    f = (A + B*sec(c + d*x))*sqrt(a + b*sec(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_596():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + sqrt(a + b*sec(c + d*x))*(20*A*b + 18*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + sqrt(a + b*sec(c + d*x))*(98*A*a**2 + 6*A*b**2 + 144*B*a*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(315*a*d) + sqrt(a + b*sec(c + d*x))*(176*A*a**2*b - 8*A*b**3 + 150*B*a**3 + 18*B*a*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*a**2*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(39*A*a**2*b + 8*A*b**3 + 75*B*a**3 - 18*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(315*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*A*a**4 + 66*A*a**2*b**2 + 16*A*b**4 + 492*B*a**3*b - 36*B*a*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_597():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + sqrt(a + b*sec(c + d*x))*(16*A*b + 14*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 + 6*A*b**2 + 84*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*a*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 - 6*A*b**2 + 21*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(164*A*a**2*b - 12*A*b**3 + 126*B*a**3 + 42*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_598():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + sqrt(a + b*sec(c + d*x))*(12*A*b + 10*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(3*A*b + 5*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 + 6*A*b**2 + 40*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_599():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_600():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_601():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(7) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(12) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_602():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(42) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(17) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_603():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(11*d) + 2*a*sqrt(a + b*sec(c + d*x))*(14*A*b + 11*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(99*d) + sqrt(a + b*sec(c + d*x))*(162*A*a**2 + 226*A*b**2 + 418*B*a*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(693*d) + sqrt(a + b*sec(c + d*x))*(2290*A*a**2*b + 30*A*b**3 + 1078*B*a**3 + 1650*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3465*a*d) + sqrt(a + b*sec(c + d*x))*(1350*A*a**4 + 2050*A*a**2*b**2 - 40*A*b**4 + 3586*B*a**3*b + 110*B*a*b**3)*sin(c + d*x)*sqrt(cos(c + d*x))/(3465*a**2*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(675*A*a**4 + 285*A*a**2*b**2 + 40*A*b**4 + 1254*B*a**3*b - 110*B*a*b**3)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3465*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(7410*A*a**4*b + 510*A*a**2*b**3 + 80*A*b**5 + 3234*B*a**5 + 6138*B*a**3*b**2 - 220*B*a*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3465*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_604():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*a*sqrt(a + b*sec(c + d*x))*(4*A*b + 3*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(21*d) + sqrt(a + b*sec(c + d*x))*(98*A*a**2 + 150*A*b**2 + 270*B*a*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(315*d) + sqrt(a + b*sec(c + d*x))*(326*A*a**2*b + 10*A*b**3 + 150*B*a**3 + 270*B*a*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*a*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(114*A*a**2*b - 10*A*b**3 + 75*B*a**3 + 45*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*A*a**4 + 558*A*a**2*b**2 - 20*A*b**4 + 870*B*a**3*b + 90*B*a*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_605():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*sqrt(a + b*sec(c + d*x))*(10*A*b + 7*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + sqrt(a + b*sec(c + d*x))*(50*A*a**2 + 90*A*b**2 + 154*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*(25*A*a**2 + 15*A*b**2 + 56*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(105*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(290*A*a**2*b + 30*A*b**3 + 126*B*a**3 + 322*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_606():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(10) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(23) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(35) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(15) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_607():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(4) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_608():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(20) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_609():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = ((((Integer(48) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(66) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(59) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(30) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_610():
    f = (A + B*sec(c + d*x))*(a + b*sec(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = ((((Integer(472) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(133) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(356) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(192) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(40) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(160) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(64) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(11) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(104) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(59) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(36) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_611():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*sec(c + d*x))
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) - sqrt(a + b*sec(c + d*x))*(8*A*b - 10*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a**2*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-14*A*a**2*b - 16*A*b**3 + 10*B*a**3 + 20*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**2 + 16*A*b**2 - 20*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_612():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*sec(c + d*x))
    F = 2*A*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**2 + 4*A*b**2 - 6*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) - sqrt(a + b*sec(c + d*x))*(4*A*b - 6*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_613():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/sqrt(a + b*sec(c + d*x))
    F = 2*A*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_614():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    F = ((Integer(2) * Symbol('A') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_615():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = ((Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_616():
    f = (A + B*sec(c + d*x))/(sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_617():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 12*A*b**2 + 10*B*a*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a**2*d*(a**2 - b**2)) - sqrt(a + b*sec(c + d*x))*(18*A*a**2*b - 48*A*b**3 - 10*B*a**3 + 40*B*a*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a**3*d*(a**2 - b**2)) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-24*A*a**2*b - 96*A*b**3 + 10*B*a**3 + 80*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**4*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**4 + 48*A*a**2*b**2 - 96*A*b**4 - 50*B*a**3*b + 80*B*a*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_618():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 8*A*b**2 + 6*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(a**2 - b**2)) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**2 + 16*A*b**2 - 12*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) - sqrt(a + b*sec(c + d*x))*(10*A*a**2*b - 16*A*b**3 - 6*B*a**3 + 12*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_619():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-4*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*a**2 - 4*A*b**2 + 2*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_620():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = 2*A*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) - (2*A*b - 2*B*a)*sin(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*A*b - 2*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_621():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_622():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_623():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(12) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_624():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*b*(12*A*a**2*b - 8*A*b**3 - 9*B*a**3 + 5*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*sec(c + d*x))*(6*A*a**4 - 142*A*a**2*b**2 + 96*A*b**4 + 100*B*a**3*b - 60*B*a*b**3)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a**3*d*(a**2 - b**2)**2) - sqrt(a + b*sec(c + d*x))*(28*A*a**4*b - 196*A*a**2*b**3 + 128*A*b**5 - 10*B*a**5 + 130*B*a**3*b**2 - 80*B*a*b**4)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a**4*d*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-34*A*a**4*b - 232*A*a**2*b**3 + 256*A*b**5 + 10*B*a**5 + 160*B*a**3*b**2 - 160*B*a*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**5*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*A*a**6 + 110*A*a**4*b**2 - 424*A*a**2*b**4 + 256*A*b**6 - 80*B*a**5*b + 280*B*a**3*b**3 - 160*B*a*b**5)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**5*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_625():
    f = (A + B*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*b*(10*A*a**2*b - 6*A*b**3 - 7*B*a**3 + 3*B*a*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*sec(c + d*x))*(2*A*a**4 - 26*A*a**2*b**2 + 16*A*b**4 + 16*B*a**3*b - 8*B*a*b**3)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**3*d*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*A*a**4 + 32*A*a**2*b**2 - 32*A*b**4 - 18*B*a**3*b + 16*B*a*b**3)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - sqrt(a + b*sec(c + d*x))*(16*A*a**4*b - 56*A*a**2*b**3 + 32*A*b**5 - 6*B*a**5 + 30*B*a**3*b**2 - 16*B*a*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_626():
    f = (A + B*sec(c + d*x))*sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + 2*b*(8*A*a**2*b - 4*A*b**3 - 5*B*a**3 + B*a*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-18*A*a**2*b + 16*A*b**3 + 6*B*a**3 - 4*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*A*a**4 - 30*A*a**2*b**2 + 16*A*b**4 + 12*B*a**3*b - 4*B*a*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_627():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = -(2*A*b - 2*B*a)*sin(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(cos(c + d*x))) - (10*A*a**2*b - 2*A*b**3 - 4*B*a**3 - 4*B*a*b**2)*sin(c + d*x)/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(6*A*a**2 - 4*A*b**2 - 2*B*a*b)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_628():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*a*(A*b - B*a)*sin(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + (4*A*a**2*b + 4*A*b**3 + 2*B*a**3 - 10*B*a*b**2)*sin(c + d*x)/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(-2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - sqrt(a + b*sec(c + d*x))*(6*A*a**2 + 2*A*b**2 - 8*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_629():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_3_1_a_plus_b_sec_pow_m_d_sec_pow_n_A_plus_B_sec_630():
    f = (A + B*sec(c + d*x))/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(-1) * (((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F

