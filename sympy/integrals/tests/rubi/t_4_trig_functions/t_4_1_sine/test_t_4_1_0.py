"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.0 (a sin)^m (b trg)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_1():
    f = sin(a + b*x)
    F = -cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_2():
    f = sin(a + b*x)**2
    F = x/2 - sin(a + b*x)*cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_3():
    f = sin(a + b*x)**3
    F = cos(a + b*x)**3/(3*b) - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_4():
    f = sin(a + b*x)**4
    F = 3*x/8 - sin(a + b*x)**3*cos(a + b*x)/(4*b) - 3*sin(a + b*x)*cos(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_5():
    f = sin(a + b*x)**5
    F = -cos(a + b*x)**5/(5*b) + 2*cos(a + b*x)**3/(3*b) - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_6():
    f = sin(a + b*x)**6
    F = 5*x/16 - sin(a + b*x)**5*cos(a + b*x)/(6*b) - 5*sin(a + b*x)**3*cos(a + b*x)/(24*b) - 5*sin(a + b*x)*cos(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_7():
    f = sin(a + b*x)**7
    F = cos(a + b*x)**7/(7*b) - 3*cos(a + b*x)**5/(5*b) + cos(a + b*x)**3/b - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_8():
    f = sin(a + b*x)**8
    F = 35*x/128 - sin(a + b*x)**7*cos(a + b*x)/(8*b) - 7*sin(a + b*x)**5*cos(a + b*x)/(48*b) - 35*sin(a + b*x)**3*cos(a + b*x)/(192*b) - 35*sin(a + b*x)*cos(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_9():
    f = sin(b*x)**(sympy.S(7)/2)
    F = -2*sin(b*x)**(sympy.S(5)/2)*cos(b*x)/(7*b) - 10*sqrt(sin(b*x))*cos(b*x)/(21*b) + 10*elliptic_f(b*x/2 - pi/4, 2)/(21*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_10():
    f = sin(b*x)**(sympy.S(5)/2)
    F = -2*sin(b*x)**(sympy.S(3)/2)*cos(b*x)/(5*b) + 6*elliptic_e(b*x/2 - pi/4, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_11():
    f = sin(b*x)**(sympy.S(3)/2)
    F = -2*sqrt(sin(b*x))*cos(b*x)/(3*b) + 2*elliptic_f(b*x/2 - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_12():
    f = sqrt(sin(b*x))
    F = 2*elliptic_e(b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_13():
    f = 1/sqrt(sin(b*x))
    F = 2*elliptic_f(b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_14():
    f = sin(b*x)**(sympy.S(-3)/2)
    F = -2*elliptic_e(b*x/2 - pi/4, 2)/b - 2*cos(b*x)/(b*sqrt(sin(b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_15():
    f = sin(b*x)**(sympy.S(-5)/2)
    F = 2*elliptic_f(b*x/2 - pi/4, 2)/(3*b) - 2*cos(b*x)/(3*b*sin(b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_16():
    f = sin(b*x)**(sympy.S(-7)/2)
    F = -6*elliptic_e(b*x/2 - pi/4, 2)/(5*b) - 6*cos(b*x)/(5*b*sqrt(sin(b*x))) - 2*cos(b*x)/(5*b*sin(b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_17():
    f = sin(a + b*x)**(sympy.S(7)/2)
    F = -2*sin(a + b*x)**(sympy.S(5)/2)*cos(a + b*x)/(7*b) - 10*sqrt(sin(a + b*x))*cos(a + b*x)/(21*b) + 10*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(21*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_18():
    f = sin(a + b*x)**(sympy.S(5)/2)
    F = -2*sin(a + b*x)**(sympy.S(3)/2)*cos(a + b*x)/(5*b) + 6*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_19():
    f = sin(a + b*x)**(sympy.S(3)/2)
    F = -2*sqrt(sin(a + b*x))*cos(a + b*x)/(3*b) + 2*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_20():
    f = sqrt(sin(a + b*x))
    F = 2*elliptic_e(a/2 + b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_21():
    f = 1/sqrt(sin(a + b*x))
    F = 2*elliptic_f(a/2 + b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_22():
    f = sin(a + b*x)**(sympy.S(-3)/2)
    F = -2*elliptic_e(a/2 + b*x/2 - pi/4, 2)/b - 2*cos(a + b*x)/(b*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_23():
    f = sin(a + b*x)**(sympy.S(-5)/2)
    F = 2*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b) - 2*cos(a + b*x)/(3*b*sin(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_24():
    f = sin(a + b*x)**(sympy.S(-7)/2)
    F = -6*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b) - 6*cos(a + b*x)/(5*b*sqrt(sin(a + b*x))) - 2*cos(a + b*x)/(5*b*sin(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_25():
    f = (c*sin(a + b*x))**(sympy.S(7)/2)
    F = 10*c**4*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(21*b*sqrt(c*sin(a + b*x))) - 10*c**3*sqrt(c*sin(a + b*x))*cos(a + b*x)/(21*b) - 2*c*(c*sin(a + b*x))**(sympy.S(5)/2)*cos(a + b*x)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_26():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)
    F = 6*c**2*sqrt(c*sin(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b*sqrt(sin(a + b*x))) - 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_27():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)
    F = 2*c**2*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b*sqrt(c*sin(a + b*x))) - 2*c*sqrt(c*sin(a + b*x))*cos(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_28():
    f = sqrt(c*sin(a + b*x))
    F = 2*sqrt(c*sin(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(b*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_29():
    f = 1/sqrt(c*sin(a + b*x))
    F = 2*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(b*sqrt(c*sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_30():
    f = (c*sin(a + b*x))**(sympy.S(-3)/2)
    F = -2*cos(a + b*x)/(b*c*sqrt(c*sin(a + b*x))) - 2*sqrt(c*sin(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(b*c**2*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_31():
    f = (c*sin(a + b*x))**(sympy.S(-5)/2)
    F = -2*cos(a + b*x)/(3*b*c*(c*sin(a + b*x))**(sympy.S(3)/2)) + 2*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b*c**2*sqrt(c*sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_32():
    f = (c*sin(a + b*x))**(sympy.S(-7)/2)
    F = -2*cos(a + b*x)/(5*b*c*(c*sin(a + b*x))**(sympy.S(5)/2)) - 6*cos(a + b*x)/(5*b*c**3*sqrt(c*sin(a + b*x))) - 6*sqrt(c*sin(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b*c**4*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_33():
    f = (c*sin(a + b*x))**(sympy.S(4)/3)
    F = 3*(c*sin(a + b*x))**(sympy.S(7)/3)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), sin(a + b*x)**2)/(7*b*c*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_34():
    f = (c*sin(a + b*x))**(sympy.S(2)/3)
    F = 3*(c*sin(a + b*x))**(sympy.S(5)/3)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), sin(a + b*x)**2)/(5*b*c*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_35():
    f = (c*sin(a + b*x))**(sympy.S(-4)/3)
    F = -3*cos(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), sin(a + b*x)**2)/(b*c*(c*sin(a + b*x))**(sympy.S(1)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_36():
    f = sin(a + b*x)**n
    F = sin(a + b*x)**(n + 1)*cos(a + b*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*(n + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_37():
    f = (c*sin(a + b*x))**n
    F = (c*sin(a + b*x))**(n + 1)*cos(a + b*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(n + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_38():
    f = (a*sin(e + f*x))**m*(b*sin(e + f*x))**n
    F = (a*sin(e + f*x))**(m + 1)*(b*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(a*f*(m + n + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_39():
    f = sin(a + b*x)*cos(a + b*x)**3
    F = -cos(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_40():
    f = sin(a + b*x)*cos(a + b*x)**2
    F = -cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_41():
    f = sin(a + b*x)*cos(a + b*x)
    F = sin(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_42():
    f = sin(a + b*x)*sec(a + b*x)
    F = -log(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_43():
    f = sin(a + b*x)*sec(a + b*x)**2
    F = sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_44():
    f = sin(a + b*x)*sec(a + b*x)**3
    F = sec(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_45():
    f = sin(a + b*x)*sec(a + b*x)**4
    F = sec(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_46():
    f = sin(a + b*x)**2*cos(a + b*x)**7
    F = -sin(a + b*x)**9/(9*b) + 3*sin(a + b*x)**7/(7*b) - 3*sin(a + b*x)**5/(5*b) + sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_47():
    f = sin(a + b*x)**2*cos(a + b*x)**5
    F = sin(a + b*x)**7/(7*b) - 2*sin(a + b*x)**5/(5*b) + sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_48():
    f = sin(a + b*x)**2*cos(a + b*x)**3
    F = -sin(a + b*x)**5/(5*b) + sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_49():
    f = sin(a + b*x)**2*cos(a + b*x)
    F = sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_50():
    f = sin(a + b*x)**2*sec(a + b*x)**2
    F = -x + tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_51():
    f = sin(a + b*x)**2*sec(a + b*x)**4
    F = tan(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_52():
    f = sin(a + b*x)**2*sec(a + b*x)**6
    F = tan(a + b*x)**5/(5*b) + tan(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_53():
    f = sin(a + b*x)**2*sec(a + b*x)**8
    F = tan(a + b*x)**7/(7*b) + 2*tan(a + b*x)**5/(5*b) + tan(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_54():
    f = sin(a + b*x)**2*sec(a + b*x)**10
    F = tan(a + b*x)**9/(9*b) + 3*tan(a + b*x)**7/(7*b) + 3*tan(a + b*x)**5/(5*b) + tan(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_55():
    f = sin(a + b*x)**2*cos(a + b*x)**6
    F = 5*x/128 - sin(a + b*x)*cos(a + b*x)**7/(8*b) + sin(a + b*x)*cos(a + b*x)**5/(48*b) + 5*sin(a + b*x)*cos(a + b*x)**3/(192*b) + 5*sin(a + b*x)*cos(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_56():
    f = sin(a + b*x)**2*cos(a + b*x)**4
    F = x/16 - sin(a + b*x)*cos(a + b*x)**5/(6*b) + sin(a + b*x)*cos(a + b*x)**3/(24*b) + sin(a + b*x)*cos(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_57():
    f = sin(a + b*x)**2*cos(a + b*x)**2
    F = x/8 - sin(a + b*x)*cos(a + b*x)**3/(4*b) + sin(a + b*x)*cos(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_58():
    f = sin(a + b*x)**2
    F = x/2 - sin(a + b*x)*cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_59():
    f = sin(a + b*x)**2*sec(a + b*x)
    F = -sin(a + b*x)/b + atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_60():
    f = sin(a + b*x)**2*sec(a + b*x)**3
    F = tan(a + b*x)*sec(a + b*x)/(2*b) - atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_61():
    f = sin(a + b*x)**2*sec(a + b*x)**5
    F = tan(a + b*x)*sec(a + b*x)**3/(4*b) - tan(a + b*x)*sec(a + b*x)/(8*b) - atanh(sin(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_62():
    f = sin(a + b*x)**2*sec(a + b*x)**7
    F = tan(a + b*x)*sec(a + b*x)**5/(6*b) - tan(a + b*x)*sec(a + b*x)**3/(24*b) - tan(a + b*x)*sec(a + b*x)/(16*b) - atanh(sin(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_63():
    f = sin(a + b*x)**3*cos(a + b*x)**5
    F = cos(a + b*x)**8/(8*b) - cos(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_64():
    f = sin(a + b*x)**3*cos(a + b*x)**4
    F = cos(a + b*x)**7/(7*b) - cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_65():
    f = sin(a + b*x)**3*cos(a + b*x)**3
    F = -sin(a + b*x)**6/(6*b) + sin(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_66():
    f = sin(a + b*x)**3*cos(a + b*x)**2
    F = cos(a + b*x)**5/(5*b) - cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_67():
    f = sin(a + b*x)**3*cos(a + b*x)
    F = sin(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_68():
    f = sin(a + b*x)**3*sec(a + b*x)
    F = -log(cos(a + b*x))/b + cos(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_69():
    f = sin(a + b*x)**3*sec(a + b*x)**2
    F = cos(a + b*x)/b + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_70():
    f = sin(a + b*x)**3*sec(a + b*x)**3
    F = log(cos(a + b*x))/b + tan(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_71():
    f = sin(a + b*x)**3*sec(a + b*x)**4
    F = sec(a + b*x)**3/(3*b) - sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_72():
    f = sin(a + b*x)**3*sec(a + b*x)**5
    F = tan(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_73():
    f = sin(a + b*x)**3*sec(a + b*x)**6
    F = sec(a + b*x)**5/(5*b) - sec(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_74():
    f = sin(a + b*x)**3*sec(a + b*x)**7
    F = sec(a + b*x)**6/(6*b) - sec(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_75():
    f = sin(a + b*x)**3*sec(a + b*x)**8
    F = sec(a + b*x)**7/(7*b) - sec(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_76():
    f = sin(a + b*x)**3*sec(a + b*x)**9
    F = sec(a + b*x)**8/(8*b) - sec(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_77():
    f = sin(a + b*x)**4*cos(a + b*x)**7
    F = -sin(a + b*x)**11/(11*b) + sin(a + b*x)**9/(3*b) - 3*sin(a + b*x)**7/(7*b) + sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_78():
    f = sin(a + b*x)**4*cos(a + b*x)**5
    F = sin(a + b*x)**9/(9*b) - 2*sin(a + b*x)**7/(7*b) + sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_79():
    f = sin(a + b*x)**4*cos(a + b*x)**3
    F = -sin(a + b*x)**7/(7*b) + sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_80():
    f = sin(a + b*x)**4*cos(a + b*x)
    F = sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_81():
    f = sin(a + b*x)**4*sec(a + b*x)**2
    F = -3*x/2 - sin(a + b*x)**2*tan(a + b*x)/(2*b) + 3*tan(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_82():
    f = sin(a + b*x)**4*sec(a + b*x)**4
    F = x + tan(a + b*x)**3/(3*b) - tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_83():
    f = sin(a + b*x)**4*sec(a + b*x)**6
    F = tan(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_84():
    f = sin(a + b*x)**4*sec(a + b*x)**8
    F = tan(a + b*x)**7/(7*b) + tan(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_85():
    f = sin(a + b*x)**4*sec(a + b*x)**10
    F = tan(a + b*x)**9/(9*b) + 2*tan(a + b*x)**7/(7*b) + tan(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_86():
    f = sin(a + b*x)**4*cos(a + b*x)**6
    F = 3*x/256 - sin(a + b*x)**3*cos(a + b*x)**7/(10*b) - 3*sin(a + b*x)*cos(a + b*x)**7/(80*b) + sin(a + b*x)*cos(a + b*x)**5/(160*b) + sin(a + b*x)*cos(a + b*x)**3/(128*b) + 3*sin(a + b*x)*cos(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_87():
    f = sin(a + b*x)**4*cos(a + b*x)**4
    F = 3*x/128 - sin(a + b*x)**3*cos(a + b*x)**5/(8*b) - sin(a + b*x)*cos(a + b*x)**5/(16*b) + sin(a + b*x)*cos(a + b*x)**3/(64*b) + 3*sin(a + b*x)*cos(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_88():
    f = sin(a + b*x)**4*cos(a + b*x)**2
    F = x/16 - sin(a + b*x)**3*cos(a + b*x)**3/(6*b) - sin(a + b*x)*cos(a + b*x)**3/(8*b) + sin(a + b*x)*cos(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_89():
    f = sin(a + b*x)**4
    F = 3*x/8 - sin(a + b*x)**3*cos(a + b*x)/(4*b) - 3*sin(a + b*x)*cos(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_90():
    f = sin(a + b*x)**4*sec(a + b*x)
    F = -sin(a + b*x)**3/(3*b) - sin(a + b*x)/b + atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_91():
    f = sin(a + b*x)**4*sec(a + b*x)**3
    F = sin(a + b*x)*tan(a + b*x)**2/(2*b) + 3*sin(a + b*x)/(2*b) - 3*atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_92():
    f = sin(a + b*x)**4*sec(a + b*x)**5
    F = tan(a + b*x)**3*sec(a + b*x)/(4*b) - 3*tan(a + b*x)*sec(a + b*x)/(8*b) + 3*atanh(sin(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_93():
    f = sin(a + b*x)**4*sec(a + b*x)**7
    F = tan(a + b*x)**3*sec(a + b*x)**3/(6*b) - tan(a + b*x)*sec(a + b*x)**3/(8*b) + tan(a + b*x)*sec(a + b*x)/(16*b) + atanh(sin(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_94():
    f = sin(a + b*x)**4*sec(a + b*x)**9
    F = tan(a + b*x)**3*sec(a + b*x)**5/(8*b) - tan(a + b*x)*sec(a + b*x)**5/(16*b) + tan(a + b*x)*sec(a + b*x)**3/(64*b) + 3*tan(a + b*x)*sec(a + b*x)/(128*b) + 3*atanh(sin(a + b*x))/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_95():
    f = sin(a + b*x)**5*cos(a + b*x)**7
    F = -cos(a + b*x)**12/(12*b) + cos(a + b*x)**10/(5*b) - cos(a + b*x)**8/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_96():
    f = sin(a + b*x)**5*cos(a + b*x)**6
    F = -cos(a + b*x)**11/(11*b) + 2*cos(a + b*x)**9/(9*b) - cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_97():
    f = sin(a + b*x)**5*cos(a + b*x)**5
    F = sin(a + b*x)**10/(10*b) - sin(a + b*x)**8/(4*b) + sin(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_98():
    f = sin(a + b*x)**5*cos(a + b*x)**4
    F = -cos(a + b*x)**9/(9*b) + 2*cos(a + b*x)**7/(7*b) - cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_99():
    f = sin(a + b*x)**5*cos(a + b*x)**3
    F = -sin(a + b*x)**8/(8*b) + sin(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_100():
    f = sin(a + b*x)**5*cos(a + b*x)**2
    F = -cos(a + b*x)**7/(7*b) + 2*cos(a + b*x)**5/(5*b) - cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_101():
    f = sin(a + b*x)**5*cos(a + b*x)
    F = sin(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_102():
    f = sin(a + b*x)**5*sec(a + b*x)
    F = -log(cos(a + b*x))/b - cos(a + b*x)**4/(4*b) + cos(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_103():
    f = sin(a + b*x)**5*sec(a + b*x)**2
    F = -cos(a + b*x)**3/(3*b) + 2*cos(a + b*x)/b + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_104():
    f = sin(a + b*x)**5*sec(a + b*x)**3
    F = 2*log(cos(a + b*x))/b - cos(a + b*x)**2/(2*b) + sec(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_105():
    f = sin(a + b*x)**5*sec(a + b*x)**4
    F = -cos(a + b*x)/b + sec(a + b*x)**3/(3*b) - 2*sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_106():
    f = sin(a + b*x)**5*sec(a + b*x)**5
    F = -log(cos(a + b*x))/b + tan(a + b*x)**4/(4*b) - tan(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_107():
    f = sin(a + b*x)**5*sec(a + b*x)**6
    F = sec(a + b*x)**5/(5*b) - 2*sec(a + b*x)**3/(3*b) + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_108():
    f = sin(a + b*x)**5*sec(a + b*x)**7
    F = tan(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_109():
    f = sin(a + b*x)**5*sec(a + b*x)**8
    F = sec(a + b*x)**7/(7*b) - 2*sec(a + b*x)**5/(5*b) + sec(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_110():
    f = sin(a + b*x)**5*sec(a + b*x)**9
    F = tan(a + b*x)**8/(8*b) + tan(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_111():
    f = sin(a + b*x)**5*sec(a + b*x)**10
    F = sec(a + b*x)**9/(9*b) - 2*sec(a + b*x)**7/(7*b) + sec(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_112():
    f = sin(a + b*x)**5*sec(a + b*x)**11
    F = sec(a + b*x)**10/(10*b) - sec(a + b*x)**8/(4*b) + sec(a + b*x)**6/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_113():
    f = sin(a + b*x)**5*sec(a + b*x)**12
    F = sec(a + b*x)**11/(11*b) - 2*sec(a + b*x)**9/(9*b) + sec(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_114():
    f = sin(a + b*x)**5*sec(a + b*x)**13
    F = sec(a + b*x)**12/(12*b) - sec(a + b*x)**10/(5*b) + sec(a + b*x)**8/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_115():
    f = sin(a + b*x)**6*sec(a + b*x)**3
    F = sin(a + b*x)**3*tan(a + b*x)**2/(2*b) + 5*sin(a + b*x)**3/(6*b) + 5*sin(a + b*x)/(2*b) - 5*atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_116():
    f = sin(a + b*x)**7*sec(a + b*x)**6
    F = cos(a + b*x)/b + sec(a + b*x)**5/(5*b) - sec(a + b*x)**3/b + 3*sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_117():
    f = cos(a + b*x)**6/sin(a + b*x)
    F = cos(a + b*x)**5/(5*b) + cos(a + b*x)**3/(3*b) + cos(a + b*x)/b - atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_118():
    f = cos(a + b*x)**5/sin(a + b*x)
    F = log(sin(a + b*x))/b + sin(a + b*x)**4/(4*b) - sin(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_119():
    f = cos(a + b*x)**4/sin(a + b*x)
    F = cos(a + b*x)**3/(3*b) + cos(a + b*x)/b - atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_120():
    f = cos(a + b*x)**3/sin(a + b*x)
    F = log(sin(a + b*x))/b - sin(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_121():
    f = cos(a + b*x)**2/sin(a + b*x)
    F = cos(a + b*x)/b - atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_122():
    f = cos(a + b*x)/sin(a + b*x)
    F = log(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_123():
    f = sec(a + b*x)/sin(a + b*x)
    F = log(tan(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_124():
    f = sec(a + b*x)**2/sin(a + b*x)
    F = -atanh(cos(a + b*x))/b + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_125():
    f = sec(a + b*x)**3/sin(a + b*x)
    F = log(tan(a + b*x))/b + tan(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_126():
    f = sec(a + b*x)**4/sin(a + b*x)
    F = -atanh(cos(a + b*x))/b + sec(a + b*x)**3/(3*b) + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_127():
    f = sec(a + b*x)**5/sin(a + b*x)
    F = log(tan(a + b*x))/b + tan(a + b*x)**4/(4*b) + tan(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_128():
    f = sec(a + b*x)**6/sin(a + b*x)
    F = -atanh(cos(a + b*x))/b + sec(a + b*x)**5/(5*b) + sec(a + b*x)**3/(3*b) + sec(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_129():
    f = sec(a + b*x)**7/sin(a + b*x)
    F = log(tan(a + b*x))/b + tan(a + b*x)**6/(6*b) + 3*tan(a + b*x)**4/(4*b) + 3*tan(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_130():
    f = cos(a + b*x)**7/sin(a + b*x)**2
    F = -sin(a + b*x)**5/(5*b) + sin(a + b*x)**3/b - 3*sin(a + b*x)/b - csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_131():
    f = cos(a + b*x)**6/sin(a + b*x)**2
    F = -15*x/8 + cos(a + b*x)**4*cot(a + b*x)/(4*b) + 5*cos(a + b*x)**2*cot(a + b*x)/(8*b) - 15*cot(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_132():
    f = cos(a + b*x)**5/sin(a + b*x)**2
    F = sin(a + b*x)**3/(3*b) - 2*sin(a + b*x)/b - csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_133():
    f = cos(a + b*x)**4/sin(a + b*x)**2
    F = -3*x/2 + cos(a + b*x)**2*cot(a + b*x)/(2*b) - 3*cot(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_134():
    f = cos(a + b*x)**3/sin(a + b*x)**2
    F = -sin(a + b*x)/b - csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_135():
    f = cos(a + b*x)**2/sin(a + b*x)**2
    F = -x - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_136():
    f = cos(a + b*x)/sin(a + b*x)**2
    F = -csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_137():
    f = sec(a + b*x)/sin(a + b*x)**2
    F = atanh(sin(a + b*x))/b - csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_138():
    f = sec(a + b*x)**2/sin(a + b*x)**2
    F = tan(a + b*x)/b - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_139():
    f = sec(a + b*x)**3/sin(a + b*x)**2
    F = 3*atanh(sin(a + b*x))/(2*b) + csc(a + b*x)*sec(a + b*x)**2/(2*b) - 3*csc(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_140():
    f = sec(a + b*x)**4/sin(a + b*x)**2
    F = tan(a + b*x)**3/(3*b) + 2*tan(a + b*x)/b - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_141():
    f = sec(a + b*x)**5/sin(a + b*x)**2
    F = 15*atanh(sin(a + b*x))/(8*b) + csc(a + b*x)*sec(a + b*x)**4/(4*b) + 5*csc(a + b*x)*sec(a + b*x)**2/(8*b) - 15*csc(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_142():
    f = cos(a + b*x)**7/sin(a + b*x)**3
    F = -3*log(sin(a + b*x))/b - sin(a + b*x)**4/(4*b) + 3*sin(a + b*x)**2/(2*b) - csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_143():
    f = cos(a + b*x)**6/sin(a + b*x)**3
    F = -cos(a + b*x)**3*cot(a + b*x)**2/(2*b) - 5*cos(a + b*x)**3/(6*b) - 5*cos(a + b*x)/(2*b) + 5*atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_144():
    f = cos(a + b*x)**5/sin(a + b*x)**3
    F = -2*log(sin(a + b*x))/b + sin(a + b*x)**2/(2*b) - csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_145():
    f = cos(a + b*x)**4/sin(a + b*x)**3
    F = -cos(a + b*x)*cot(a + b*x)**2/(2*b) - 3*cos(a + b*x)/(2*b) + 3*atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_146():
    f = cos(a + b*x)**3/sin(a + b*x)**3
    F = -log(sin(a + b*x))/b - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_147():
    f = cos(a + b*x)**2/sin(a + b*x)**3
    F = -cot(a + b*x)*csc(a + b*x)/(2*b) + atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_148():
    f = cos(a + b*x)/sin(a + b*x)**3
    F = -csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_149():
    f = sec(a + b*x)/sin(a + b*x)**3
    F = log(tan(a + b*x))/b - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_150():
    f = sec(a + b*x)**2/sin(a + b*x)**3
    F = -3*atanh(cos(a + b*x))/(2*b) - csc(a + b*x)**2*sec(a + b*x)/(2*b) + 3*sec(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_151():
    f = sec(a + b*x)**3/sin(a + b*x)**3
    F = 2*log(tan(a + b*x))/b + tan(a + b*x)**2/(2*b) - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_152():
    f = sec(a + b*x)**4/sin(a + b*x)**3
    F = -5*atanh(cos(a + b*x))/(2*b) - csc(a + b*x)**2*sec(a + b*x)**3/(2*b) + 5*sec(a + b*x)**3/(6*b) + 5*sec(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_153():
    f = sec(a + b*x)**5/sin(a + b*x)**3
    F = 3*log(tan(a + b*x))/b + tan(a + b*x)**4/(4*b) + 3*tan(a + b*x)**2/(2*b) - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_154():
    f = cos(a + b*x)**9/sin(a + b*x)**4
    F = sin(a + b*x)**5/(5*b) - 4*sin(a + b*x)**3/(3*b) + 6*sin(a + b*x)/b - csc(a + b*x)**3/(3*b) + 4*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_155():
    f = cos(a + b*x)**8/sin(a + b*x)**4
    F = 35*x/8 + cos(a + b*x)**4*cot(a + b*x)**3/(4*b) + 7*cos(a + b*x)**2*cot(a + b*x)**3/(8*b) - 35*cot(a + b*x)**3/(24*b) + 35*cot(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_156():
    f = cos(a + b*x)**7/sin(a + b*x)**4
    F = -sin(a + b*x)**3/(3*b) + 3*sin(a + b*x)/b - csc(a + b*x)**3/(3*b) + 3*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_157():
    f = cos(a + b*x)**6/sin(a + b*x)**4
    F = 5*x/2 + cos(a + b*x)**2*cot(a + b*x)**3/(2*b) - 5*cot(a + b*x)**3/(6*b) + 5*cot(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_158():
    f = cos(a + b*x)**5/sin(a + b*x)**4
    F = sin(a + b*x)/b - csc(a + b*x)**3/(3*b) + 2*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_159():
    f = cos(a + b*x)**4/sin(a + b*x)**4
    F = x - cot(a + b*x)**3/(3*b) + cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_160():
    f = cos(a + b*x)**3/sin(a + b*x)**4
    F = -csc(a + b*x)**3/(3*b) + csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_161():
    f = cos(a + b*x)**2/sin(a + b*x)**4
    F = -cot(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_162():
    f = cos(a + b*x)/sin(a + b*x)**4
    F = -csc(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_163():
    f = sec(a + b*x)/sin(a + b*x)**4
    F = atanh(sin(a + b*x))/b - csc(a + b*x)**3/(3*b) - csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_164():
    f = sec(a + b*x)**2/sin(a + b*x)**4
    F = tan(a + b*x)/b - cot(a + b*x)**3/(3*b) - 2*cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_165():
    f = sec(a + b*x)**3/sin(a + b*x)**4
    F = 5*atanh(sin(a + b*x))/(2*b) + csc(a + b*x)**3*sec(a + b*x)**2/(2*b) - 5*csc(a + b*x)**3/(6*b) - 5*csc(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_166():
    f = sec(a + b*x)**4/sin(a + b*x)**4
    F = tan(a + b*x)**3/(3*b) + 3*tan(a + b*x)/b - cot(a + b*x)**3/(3*b) - 3*cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_167():
    f = sec(a + b*x)**5/sin(a + b*x)**4
    F = 35*atanh(sin(a + b*x))/(8*b) + csc(a + b*x)**3*sec(a + b*x)**4/(4*b) + 7*csc(a + b*x)**3*sec(a + b*x)**2/(8*b) - 35*csc(a + b*x)**3/(24*b) - 35*csc(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_168():
    f = cos(a + b*x)**9/sin(a + b*x)**5
    F = 6*log(sin(a + b*x))/b + sin(a + b*x)**4/(4*b) - 2*sin(a + b*x)**2/b - csc(a + b*x)**4/(4*b) + 2*csc(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_169():
    f = cos(a + b*x)**8/sin(a + b*x)**5
    F = -cos(a + b*x)**3*cot(a + b*x)**4/(4*b) + 7*cos(a + b*x)**3*cot(a + b*x)**2/(8*b) + 35*cos(a + b*x)**3/(24*b) + 35*cos(a + b*x)/(8*b) - 35*atanh(cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_170():
    f = cos(a + b*x)**7/sin(a + b*x)**5
    F = 3*log(sin(a + b*x))/b - sin(a + b*x)**2/(2*b) - csc(a + b*x)**4/(4*b) + 3*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_171():
    f = cos(a + b*x)**6/sin(a + b*x)**5
    F = -cos(a + b*x)*cot(a + b*x)**4/(4*b) + 5*cos(a + b*x)*cot(a + b*x)**2/(8*b) + 15*cos(a + b*x)/(8*b) - 15*atanh(cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_172():
    f = cos(a + b*x)**5/sin(a + b*x)**5
    F = log(sin(a + b*x))/b - cot(a + b*x)**4/(4*b) + cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_173():
    f = cos(a + b*x)**4/sin(a + b*x)**5
    F = -cot(a + b*x)**3*csc(a + b*x)/(4*b) + 3*cot(a + b*x)*csc(a + b*x)/(8*b) - 3*atanh(cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_174():
    f = cos(a + b*x)**3/sin(a + b*x)**5
    F = -cot(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_175():
    f = cos(a + b*x)**2/sin(a + b*x)**5
    F = -cot(a + b*x)*csc(a + b*x)**3/(4*b) + cot(a + b*x)*csc(a + b*x)/(8*b) + atanh(cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_176():
    f = cos(a + b*x)/sin(a + b*x)**5
    F = -csc(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_177():
    f = sec(a + b*x)/sin(a + b*x)**5
    F = log(tan(a + b*x))/b - cot(a + b*x)**4/(4*b) - cot(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_178():
    f = sec(a + b*x)**2/sin(a + b*x)**5
    F = -15*atanh(cos(a + b*x))/(8*b) - csc(a + b*x)**4*sec(a + b*x)/(4*b) - 5*csc(a + b*x)**2*sec(a + b*x)/(8*b) + 15*sec(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_179():
    f = sec(a + b*x)**3/sin(a + b*x)**5
    F = 3*log(tan(a + b*x))/b + tan(a + b*x)**2/(2*b) - cot(a + b*x)**4/(4*b) - 3*cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_180():
    f = sec(a + b*x)**4/sin(a + b*x)**5
    F = -35*atanh(cos(a + b*x))/(8*b) - csc(a + b*x)**4*sec(a + b*x)**3/(4*b) - 7*csc(a + b*x)**2*sec(a + b*x)**3/(8*b) + 35*sec(a + b*x)**3/(24*b) + 35*sec(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_181():
    f = sec(a + b*x)**5/sin(a + b*x)**5
    F = 6*log(tan(a + b*x))/b + tan(a + b*x)**4/(4*b) + 2*tan(a + b*x)**2/b - cot(a + b*x)**4/(4*b) - 2*cot(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_182():
    f = cos(x)**2/sin(x)**6
    F = -cot(x)**5/5 - cot(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_183():
    f = cos(x)**3/sin(x)**7
    F = -csc(x)**6/6 + csc(x)**4/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_184():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)
    F = -2*(d*cos(a + b*x))**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_185():
    f = sqrt(d*cos(a + b*x))*sin(a + b*x)
    F = -2*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_186():
    f = sin(a + b*x)/sqrt(d*cos(a + b*x))
    F = -2*sqrt(d*cos(a + b*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_187():
    f = sin(a + b*x)/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 2/(b*d*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_188():
    f = sin(a + b*x)/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_189():
    f = sin(a + b*x)/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_190():
    f = sin(a + b*x)/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_191():
    f = (d*cos(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)**2
    F = 28*d**4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(195*b*sqrt(cos(a + b*x))) + 28*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(585*b) + 4*d*(d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)/(117*b) - 2*(d*cos(a + b*x))**(sympy.S(11)/2)*sin(a + b*x)/(13*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_192():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)**2
    F = 20*d**4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(231*b*sqrt(d*cos(a + b*x))) + 20*d**3*sqrt(d*cos(a + b*x))*sin(a + b*x)/(231*b) + 4*d*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(77*b) - 2*(d*cos(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)/(11*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_193():
    f = (d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**2
    F = 4*d**2*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(15*b*sqrt(cos(a + b*x))) + 4*d*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(45*b) - 2*(d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)/(9*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_194():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**2
    F = 4*d**2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(21*b*sqrt(d*cos(a + b*x))) + 4*d*sqrt(d*cos(a + b*x))*sin(a + b*x)/(21*b) - 2*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_195():
    f = sqrt(d*cos(a + b*x))*sin(a + b*x)**2
    F = 4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*sqrt(cos(a + b*x))) - 2*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_196():
    f = sin(a + b*x)**2/sqrt(d*cos(a + b*x))
    F = 4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*sqrt(d*cos(a + b*x))) - 2*sqrt(d*cos(a + b*x))*sin(a + b*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_197():
    f = sin(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 2*sin(a + b*x)/(b*d*sqrt(d*cos(a + b*x))) - 4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*d**2*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_198():
    f = sin(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2*sin(a + b*x)/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) - 4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*d**2*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_199():
    f = sin(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2*sin(a + b*x)/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) - 4*sin(a + b*x)/(5*b*d**3*sqrt(d*cos(a + b*x))) + 4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*d**4*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_200():
    f = sin(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2*sin(a + b*x)/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2)) - 4*sin(a + b*x)/(21*b*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)) - 4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(21*b*d**4*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_201():
    f = sqrt(d*cos(a + b*x))*sin(a + b*x)**3
    F = -2*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b*d) + 2*(d*cos(a + b*x))**(sympy.S(7)/2)/(7*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_202():
    f = sin(a + b*x)**3/sqrt(d*cos(a + b*x))
    F = -2*sqrt(d*cos(a + b*x))/(b*d) + 2*(d*cos(a + b*x))**(sympy.S(5)/2)/(5*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_203():
    f = sin(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 2/(b*d*sqrt(d*cos(a + b*x))) + 2*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_204():
    f = sin(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) + 2*sqrt(d*cos(a + b*x))/(b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_205():
    f = sin(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) - 2/(b*d**3*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_206():
    f = sin(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2)) - 2/(3*b*d**3*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_207():
    f = sin(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(11)/2)
    F = 2/(9*b*d*(d*cos(a + b*x))**(sympy.S(9)/2)) - 2/(5*b*d**3*(d*cos(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_208():
    f = (d*cos(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)**4
    F = 56*d**4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(1105*b*sqrt(cos(a + b*x))) + 56*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(3315*b) + 8*d*(d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)/(663*b) - 2*(d*cos(a + b*x))**(sympy.S(11)/2)*sin(a + b*x)**3/(17*b*d) - 12*(d*cos(a + b*x))**(sympy.S(11)/2)*sin(a + b*x)/(221*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_209():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)**4
    F = 8*d**4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(231*b*sqrt(d*cos(a + b*x))) + 8*d**3*sqrt(d*cos(a + b*x))*sin(a + b*x)/(231*b) + 8*d*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(385*b) - 2*(d*cos(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)**3/(15*b*d) - 4*(d*cos(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)/(55*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_210():
    f = (d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**4
    F = 8*d**2*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(65*b*sqrt(cos(a + b*x))) + 8*d*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(195*b) - 2*(d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)**3/(13*b*d) - 4*(d*cos(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)/(39*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_211():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**4
    F = 8*d**2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(77*b*sqrt(d*cos(a + b*x))) + 8*d*sqrt(d*cos(a + b*x))*sin(a + b*x)/(77*b) - 2*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**3/(11*b*d) - 12*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(77*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_212():
    f = sqrt(d*cos(a + b*x))*sin(a + b*x)**4
    F = 8*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(15*b*sqrt(cos(a + b*x))) - 2*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**3/(9*b*d) - 4*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(15*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_213():
    f = sin(a + b*x)**4/sqrt(d*cos(a + b*x))
    F = 8*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(7*b*sqrt(d*cos(a + b*x))) - 2*sqrt(d*cos(a + b*x))*sin(a + b*x)**3/(7*b*d) - 4*sqrt(d*cos(a + b*x))*sin(a + b*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_214():
    f = sin(a + b*x)**4/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 2*sin(a + b*x)**3/(b*d*sqrt(d*cos(a + b*x))) - 24*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*d**2*sqrt(cos(a + b*x))) + 12*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(5*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_215():
    f = sin(a + b*x)**4/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2*sin(a + b*x)**3/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) - 8*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*d**2*sqrt(d*cos(a + b*x))) + 4*sqrt(d*cos(a + b*x))*sin(a + b*x)/(3*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_216():
    f = sin(a + b*x)**4/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2*sin(a + b*x)**3/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) - 12*sin(a + b*x)/(5*b*d**3*sqrt(d*cos(a + b*x))) + 24*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*d**4*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_217():
    f = sin(a + b*x)**4/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2*sin(a + b*x)**3/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2)) - 4*sin(a + b*x)/(7*b*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)) + 8*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(7*b*d**4*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_218():
    f = sin(a + b*x)**5*cos(a + b*x)**(sympy.S(3)/2)
    F = -2*cos(a + b*x)**(sympy.S(13)/2)/(13*b) + 4*cos(a + b*x)**(sympy.S(9)/2)/(9*b) - 2*cos(a + b*x)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_219():
    f = (d*cos(a + b*x))**(sympy.S(9)/2)*csc(a + b*x)
    F = d**(sympy.S(9)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/b - d**(sympy.S(9)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/b + 2*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b) + 2*d*(d*cos(a + b*x))**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_220():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)*csc(a + b*x)
    F = -d**(sympy.S(7)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/b - d**(sympy.S(7)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/b + 2*d**3*sqrt(d*cos(a + b*x))/b + 2*d*(d*cos(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_221():
    f = (d*cos(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)
    F = d**(sympy.S(5)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/b - d**(sympy.S(5)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/b + 2*d*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_222():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)
    F = -d**(sympy.S(3)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/b - d**(sympy.S(3)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/b + 2*d*sqrt(d*cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_223():
    f = sqrt(d*cos(a + b*x))*csc(a + b*x)
    F = sqrt(d)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/b - sqrt(d)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_224():
    f = csc(a + b*x)/sqrt(d*cos(a + b*x))
    F = -atan(sqrt(d*cos(a + b*x))/sqrt(d))/(b*sqrt(d)) - atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_225():
    f = csc(a + b*x)/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 2/(b*d*sqrt(d*cos(a + b*x))) + atan(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(3)/2)) - atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_226():
    f = csc(a + b*x)/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) - atan(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(5)/2)) - atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_227():
    f = csc(a + b*x)/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 2/(b*d**3*sqrt(d*cos(a + b*x))) + atan(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(7)/2)) - atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_228():
    f = csc(a + b*x)/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2)) + 2/(3*b*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)) - atan(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(9)/2)) - atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(b*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_229():
    f = (d*cos(a + b*x))**(sympy.S(11)/2)*csc(a + b*x)**2
    F = -15*d**6*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(7*b*sqrt(d*cos(a + b*x))) - 15*d**5*sqrt(d*cos(a + b*x))*sin(a + b*x)/(7*b) - 9*d**3*(d*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(7*b) - d*(d*cos(a + b*x))**(sympy.S(9)/2)*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_230():
    f = (d*cos(a + b*x))**(sympy.S(9)/2)*csc(a + b*x)**2
    F = -21*d**4*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*sqrt(cos(a + b*x))) - 7*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(5*b) - d*(d*cos(a + b*x))**(sympy.S(7)/2)*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_231():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)*csc(a + b*x)**2
    F = -5*d**4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*sqrt(d*cos(a + b*x))) - 5*d**3*sqrt(d*cos(a + b*x))*sin(a + b*x)/(3*b) - d*(d*cos(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_232():
    f = (d*cos(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**2
    F = -3*d**2*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*sqrt(cos(a + b*x))) - d*(d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_233():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**2
    F = -d**2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(b*sqrt(d*cos(a + b*x))) - d*sqrt(d*cos(a + b*x))*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_234():
    f = sqrt(d*cos(a + b*x))*csc(a + b*x)**2
    F = -sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*sqrt(cos(a + b*x))) - (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_235():
    f = csc(a + b*x)**2/sqrt(d*cos(a + b*x))
    F = sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(b*sqrt(d*cos(a + b*x))) - sqrt(d*cos(a + b*x))*csc(a + b*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_236():
    f = csc(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = 3*sin(a + b*x)/(b*d*sqrt(d*cos(a + b*x))) - csc(a + b*x)/(b*d*sqrt(d*cos(a + b*x))) - 3*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*d**2*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_237():
    f = csc(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 5*sin(a + b*x)/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) - csc(a + b*x)/(b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) + 5*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*d**2*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_238():
    f = csc(a + b*x)**2/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 7*sin(a + b*x)/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) - csc(a + b*x)/(b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 21*sin(a + b*x)/(5*b*d**3*sqrt(d*cos(a + b*x))) - 21*sqrt(d*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*d**4*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_239():
    f = (d*cos(a + b*x))**(sympy.S(11)/2)*csc(a + b*x)**3
    F = 9*d**(sympy.S(11)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) + 9*d**(sympy.S(11)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - 9*d**5*sqrt(d*cos(a + b*x))/(2*b) - 9*d**3*(d*cos(a + b*x))**(sympy.S(5)/2)/(10*b) - d*(d*cos(a + b*x))**(sympy.S(9)/2)*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_240():
    f = (d*cos(a + b*x))**(sympy.S(9)/2)*csc(a + b*x)**3
    F = -7*d**(sympy.S(9)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) + 7*d**(sympy.S(9)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - 7*d**3*(d*cos(a + b*x))**(sympy.S(3)/2)/(6*b) - d*(d*cos(a + b*x))**(sympy.S(7)/2)*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_241():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)*csc(a + b*x)**3
    F = 5*d**(sympy.S(7)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) + 5*d**(sympy.S(7)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - 5*d**3*sqrt(d*cos(a + b*x))/(2*b) - d*(d*cos(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_242():
    f = (d*cos(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**3
    F = -3*d**(sympy.S(5)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) + 3*d**(sympy.S(5)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - d*(d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_243():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**3
    F = d**(sympy.S(3)/2)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) + d**(sympy.S(3)/2)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - d*sqrt(d*cos(a + b*x))*csc(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_244():
    f = sqrt(d*cos(a + b*x))*csc(a + b*x)**3
    F = sqrt(d)*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - sqrt(d)*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b) - (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_245():
    f = csc(a + b*x)**3/sqrt(d*cos(a + b*x))
    F = -sqrt(d*cos(a + b*x))*csc(a + b*x)**2/(2*b*d) - 3*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*sqrt(d)) - 3*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_246():
    f = csc(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = -csc(a + b*x)**2/(2*b*d*sqrt(d*cos(a + b*x))) + 5/(2*b*d*sqrt(d*cos(a + b*x))) + 5*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(3)/2)) - 5*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_247():
    f = csc(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = -csc(a + b*x)**2/(2*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) + 7/(6*b*d*(d*cos(a + b*x))**(sympy.S(3)/2)) - 7*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(5)/2)) - 7*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_248():
    f = csc(a + b*x)**3/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = -csc(a + b*x)**2/(2*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 9/(10*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 9/(2*b*d**3*sqrt(d*cos(a + b*x))) + 9*atan(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(7)/2)) - 9*atanh(sqrt(d*cos(a + b*x))/sqrt(d))/(4*b*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_249():
    f = (d*cos(a + b*x))**(sympy.S(1)/5)*sin(a + b*x)
    F = -5*(d*cos(a + b*x))**(sympy.S(6)/5)/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_250():
    f = sqrt(sin(x))*cos(x)**3
    F = -2*sin(x)**(sympy.S(7)/2)/7 + 2*sin(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_251():
    f = sin(x)**(sympy.S(3)/2)*cos(x)**3
    F = -2*sin(x)**(sympy.S(9)/2)/9 + 2*sin(x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_252():
    f = sin(x)**(sympy.S(5)/2)*cos(x)**3
    F = -2*sin(x)**(sympy.S(11)/2)/11 + 2*sin(x)**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_253():
    f = cos(x)**3/sqrt(sin(x))
    F = -2*sin(x)**(sympy.S(5)/2)/5 + 2*sqrt(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_254():
    f = sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 7*d**4*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(20*b*sqrt(sin(2*a + 2*b*x))) + 7*d**3*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)/(30*b*c) + d*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(7)/2)/(5*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_255():
    f = sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(5)/2)
    F = d**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(2*b*sqrt(sin(2*a + 2*b*x))) + d*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_256():
    f = sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))
    F = sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_257():
    f = sqrt(c*sin(a + b*x))/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = -2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(sin(2*a + 2*b*x))) + 2*(c*sin(a + b*x))**(sympy.S(3)/2)/(b*c*d*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_258():
    f = sqrt(c*sin(a + b*x))/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = -4*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(5*b*d**4*sqrt(sin(2*a + 2*b*x))) + 2*(c*sin(a + b*x))**(sympy.S(3)/2)/(5*b*c*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 4*(c*sin(a + b*x))**(sympy.S(3)/2)/(5*b*c*d**3*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_259():
    f = sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(3)/2)
    F = sqrt(2)*sqrt(c)*d**(sympy.S(3)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(16*b) - sqrt(2)*sqrt(c)*d**(sympy.S(3)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(16*b) - sqrt(2)*sqrt(c)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(8*b) + sqrt(2)*sqrt(c)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(8*b) + d*(c*sin(a + b*x))**(sympy.S(3)/2)*sqrt(d*cos(a + b*x))/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_260():
    f = sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x))
    F = sqrt(2)*sqrt(c)*log(sqrt(c)*tan(a + b*x) + sqrt(c) - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(4*b*sqrt(d)) - sqrt(2)*sqrt(c)*log(sqrt(c)*tan(a + b*x) + sqrt(c) + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(4*b*sqrt(d)) - sqrt(2)*sqrt(c)*atan(1 - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(2*b*sqrt(d)) + sqrt(2)*sqrt(c)*atan(1 + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(2*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_261():
    f = sqrt(c*sin(a + b*x))/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 2*(c*sin(a + b*x))**(sympy.S(3)/2)/(3*b*c*d*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_262():
    f = sqrt(c*sin(a + b*x))/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2*(c*sin(a + b*x))**(sympy.S(3)/2)/(7*b*c*d*(d*cos(a + b*x))**(sympy.S(7)/2)) + 8*(c*sin(a + b*x))**(sympy.S(3)/2)/(21*b*c*d**3*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_263():
    f = sqrt(c*sin(a + b*x))/(d*cos(a + b*x))**(sympy.S(13)/2)
    F = 2*(c*sin(a + b*x))**(sympy.S(3)/2)/(11*b*c*d*(d*cos(a + b*x))**(sympy.S(11)/2)) + 16*(c*sin(a + b*x))**(sympy.S(3)/2)/(77*b*c*d**3*(d*cos(a + b*x))**(sympy.S(7)/2)) + 64*(c*sin(a + b*x))**(sympy.S(3)/2)/(231*b*c*d**5*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_264():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)
    F = c**2*d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(12*b*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + c*d*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))/(6*b) - c*sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(5)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_265():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/sqrt(d*cos(a + b*x))
    F = c**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(2*b*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) - c*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_266():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = -c**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b*d**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + 2*c*sqrt(c*sin(a + b*x))/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_267():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = -2*c**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(21*b*d**4*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + 2*c*sqrt(c*sin(a + b*x))/(7*b*d*(d*cos(a + b*x))**(sympy.S(7)/2)) - 2*c*sqrt(c*sin(a + b*x))/(21*b*d**3*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_268():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)*sqrt(d*cos(a + b*x))
    F = -sqrt(2)*c**(sympy.S(3)/2)*sqrt(d)*log(-sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(16*b) + sqrt(2)*c**(sympy.S(3)/2)*sqrt(d)*log(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(16*b) - sqrt(2)*c**(sympy.S(3)/2)*sqrt(d)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) - 1)/(8*b) - sqrt(2)*c**(sympy.S(3)/2)*sqrt(d)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) + 1)/(8*b) - c*sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(3)/2)/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_269():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = sqrt(2)*c**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(4*b*d**(sympy.S(3)/2)) - sqrt(2)*c**(sympy.S(3)/2)*log(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(4*b*d**(sympy.S(3)/2)) + sqrt(2)*c**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) - 1)/(2*b*d**(sympy.S(3)/2)) + sqrt(2)*c**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) + 1)/(2*b*d**(sympy.S(3)/2)) + 2*c*sqrt(c*sin(a + b*x))/(b*d*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_270():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 2*(c*sin(a + b*x))**(sympy.S(5)/2)/(5*b*c*d*(d*cos(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_271():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(11)/2)
    F = 2*c*sqrt(c*sin(a + b*x))/(9*b*d*(d*cos(a + b*x))**(sympy.S(9)/2)) - 2*c*sqrt(c*sin(a + b*x))/(45*b*d**3*(d*cos(a + b*x))**(sympy.S(5)/2)) - 8*c*sqrt(c*sin(a + b*x))/(45*b*d**5*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_272():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)/(d*cos(a + b*x))**(sympy.S(15)/2)
    F = 2*c*sqrt(c*sin(a + b*x))/(13*b*d*(d*cos(a + b*x))**(sympy.S(13)/2)) - 2*c*sqrt(c*sin(a + b*x))/(117*b*d**3*(d*cos(a + b*x))**(sympy.S(9)/2)) - 16*c*sqrt(c*sin(a + b*x))/(585*b*d**5*(d*cos(a + b*x))**(sympy.S(5)/2)) - 64*c*sqrt(c*sin(a + b*x))/(585*b*d**7*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_273():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)*(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 3*c**2*d**4*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(40*b*sqrt(sin(2*a + 2*b*x))) + c*d**3*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)/(20*b) + 3*c*d*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(7)/2)/(70*b) - c*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(11)/2)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_274():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)*(d*cos(a + b*x))**(sympy.S(5)/2)
    F = 3*c**2*d**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(20*b*sqrt(sin(2*a + 2*b*x))) + c*d*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)/(10*b) - c*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(7)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_275():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)*sqrt(d*cos(a + b*x))
    F = c**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(2*b*sqrt(sin(2*a + 2*b*x))) - c*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_276():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = -3*c**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(sin(2*a + 2*b*x))) + 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(b*d*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_277():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(7)/2)
    F = 6*c**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(5*b*d**4*sqrt(sin(2*a + 2*b*x))) + 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(5*b*d*(d*cos(a + b*x))**(sympy.S(5)/2)) - 6*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(5*b*d**3*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_278():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(11)/2)
    F = 4*c**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))*elliptic_e(a + b*x - pi/4, 2)/(15*b*d**6*sqrt(sin(2*a + 2*b*x))) + 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(9*b*d*(d*cos(a + b*x))**(sympy.S(9)/2)) - 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(15*b*d**3*(d*cos(a + b*x))**(sympy.S(5)/2)) - 4*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(15*b*d**5*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_279():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/sqrt(d*cos(a + b*x))
    F = 3*sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(16*b*sqrt(d)) - 3*sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(16*b*sqrt(d)) - 3*sqrt(2)*c**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(8*b*sqrt(d)) + 3*sqrt(2)*c**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(8*b*sqrt(d)) - c*(c*sin(a + b*x))**(sympy.S(3)/2)*sqrt(d*cos(a + b*x))/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_280():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = -sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(4*b*d**(sympy.S(5)/2)) + sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*tan(a + b*x) + sqrt(c) + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/sqrt(d*cos(a + b*x)))/(4*b*d**(sympy.S(5)/2)) + sqrt(2)*c**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(2*b*d**(sympy.S(5)/2)) - sqrt(2)*c**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d)*sqrt(c*sin(a + b*x))/(sqrt(c)*sqrt(d*cos(a + b*x))))/(2*b*d**(sympy.S(5)/2)) + 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(3*b*d*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_281():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(9)/2)
    F = 2*(c*sin(a + b*x))**(sympy.S(7)/2)/(7*b*c*d*(d*cos(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_282():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(13)/2)
    F = 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(11*b*d*(d*cos(a + b*x))**(sympy.S(11)/2)) - 6*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(77*b*d**3*(d*cos(a + b*x))**(sympy.S(7)/2)) - 8*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(77*b*d**5*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_283():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)/(d*cos(a + b*x))**(sympy.S(17)/2)
    F = 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(15*b*d*(d*cos(a + b*x))**(sympy.S(15)/2)) - 2*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(55*b*d**3*(d*cos(a + b*x))**(sympy.S(11)/2)) - 16*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(385*b*d**5*(d*cos(a + b*x))**(sympy.S(7)/2)) - 64*c*(c*sin(a + b*x))**(sympy.S(3)/2)/(1155*b*d**7*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_284():
    f = sin(a + b*x)**(sympy.S(7)/2)/cos(a + b*x)**(sympy.S(7)/2)
    F = -sqrt(2)*log(cot(a + b*x) + 1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) + sqrt(2)*log(cot(a + b*x) + 1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) + 2*sin(a + b*x)**(sympy.S(5)/2)/(5*b*cos(a + b*x)**(sympy.S(5)/2)) - 2*sqrt(sin(a + b*x))/(b*sqrt(cos(a + b*x))) + sqrt(2)*atan(1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b) - sqrt(2)*atan(1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_285():
    f = sin(x)**(sympy.S(3)/2)/cos(x)**(sympy.S(7)/2)
    F = 2*sin(x)**(sympy.S(5)/2)/(5*cos(x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_286():
    f = sqrt(sin(x))/sqrt(cos(x))
    F = sqrt(2)*log(-sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + tan(x) + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + tan(x) + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + 1)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_287():
    f = sin(x)**(sympy.S(5)/2)/sqrt(cos(x))
    F = 3*sqrt(2)*log(-sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + tan(x) + 1)/16 - 3*sqrt(2)*log(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + tan(x) + 1)/16 - sin(x)**(sympy.S(3)/2)*sqrt(cos(x))/2 + 3*sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) - 1)/8 + 3*sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + 1)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_288():
    f = (d*cos(a + b*x))**(sympy.S(7)/2)/sqrt(c*sin(a + b*x))
    F = 5*d**4*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(12*b*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + 5*d**3*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))/(6*b*c) + d*sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(5)/2)/(3*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_289():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)/sqrt(c*sin(a + b*x))
    F = d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(2*b*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + d*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_290():
    f = 1/(sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x)))
    F = sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(b*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_291():
    f = 1/(sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(5)/2))
    F = 2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b*d**2*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + 2*sqrt(c*sin(a + b*x))/(3*b*c*d*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_292():
    f = 1/(sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(9)/2))
    F = 4*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(7*b*d**4*sqrt(c*sin(a + b*x))*sqrt(d*cos(a + b*x))) + 2*sqrt(c*sin(a + b*x))/(7*b*c*d*(d*cos(a + b*x))**(sympy.S(7)/2)) + 4*sqrt(c*sin(a + b*x))/(7*b*c*d**3*(d*cos(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_293():
    f = sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x))
    F = -sqrt(2)*sqrt(d)*log(-sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(4*b*sqrt(c)) + sqrt(2)*sqrt(d)*log(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/sqrt(c*sin(a + b*x)) + sqrt(d)*cot(a + b*x) + sqrt(d))/(4*b*sqrt(c)) - sqrt(2)*sqrt(d)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) - 1)/(2*b*sqrt(c)) - sqrt(2)*sqrt(d)*atan(sqrt(2)*sqrt(c)*sqrt(d*cos(a + b*x))/(sqrt(d)*sqrt(c*sin(a + b*x))) + 1)/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_294():
    f = 1/(sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(3)/2))
    F = 2*sqrt(c*sin(a + b*x))/(b*c*d*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_295():
    f = 1/(sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(7)/2))
    F = 2*sqrt(c*sin(a + b*x))/(5*b*c*d*(d*cos(a + b*x))**(sympy.S(5)/2)) + 8*sqrt(c*sin(a + b*x))/(5*b*c*d**3*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_296():
    f = 1/(sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(sympy.S(11)/2))
    F = 2*sqrt(c*sin(a + b*x))/(9*b*c*d*(d*cos(a + b*x))**(sympy.S(9)/2)) + 16*sqrt(c*sin(a + b*x))/(45*b*c*d**3*(d*cos(a + b*x))**(sympy.S(5)/2)) + 64*sqrt(c*sin(a + b*x))/(45*b*c*d**5*sqrt(d*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_297():
    f = sqrt(cos(a + b*x))/sqrt(sin(a + b*x))
    F = -sqrt(2)*log(cot(a + b*x) + 1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) + sqrt(2)*log(cot(a + b*x) + 1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) + sqrt(2)*atan(1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b) - sqrt(2)*atan(1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_298():
    f = cos(a + b*x)**(sympy.S(3)/2)/sin(a + b*x)**(sympy.S(3)/2)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + tan(a + b*x) + 1)/(4*b) + sqrt(2)*log(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + tan(a + b*x) + 1)/(4*b) - sqrt(2)*atan(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) - 1)/(2*b) - sqrt(2)*atan(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + 1)/(2*b) - 2*sqrt(cos(a + b*x))/(b*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_299():
    f = cos(a + b*x)**(sympy.S(5)/2)/sin(a + b*x)**(sympy.S(5)/2)
    F = sqrt(2)*log(cot(a + b*x) + 1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) - sqrt(2)*log(cot(a + b*x) + 1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(4*b) - sqrt(2)*atan(1 - sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b) + sqrt(2)*atan(1 + sqrt(2)*sqrt(cos(a + b*x))/sqrt(sin(a + b*x)))/(2*b) - 2*cos(a + b*x)**(sympy.S(3)/2)/(3*b*sin(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_300():
    f = cos(a + b*x)**(sympy.S(7)/2)/sin(a + b*x)**(sympy.S(7)/2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + tan(a + b*x) + 1)/(4*b) - sqrt(2)*log(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + tan(a + b*x) + 1)/(4*b) + sqrt(2)*atan(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) - 1)/(2*b) + sqrt(2)*atan(sqrt(2)*sqrt(sin(a + b*x))/sqrt(cos(a + b*x)) + 1)/(2*b) + 2*sqrt(cos(a + b*x))/(b*sqrt(sin(a + b*x))) - 2*cos(a + b*x)**(sympy.S(5)/2)/(5*b*sin(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_301():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*cos(e + f*x)**4
    F = 3*(b*sin(e + f*x))**(sympy.S(4)/3)*cos(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(2)/3), (sympy.S(5)/3,), sin(e + f*x)**2)/(4*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_302():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*cos(e + f*x)**2
    F = 3*(b*sin(e + f*x))**(sympy.S(4)/3)*cos(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(2)/3), (sympy.S(5)/3,), sin(e + f*x)**2)/(4*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_303():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(4)/3)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), sin(e + f*x)**2)/(4*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_304():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*sec(e + f*x)**2
    F = 3*(b*sin(e + f*x))**(sympy.S(4)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(2)/3, sympy.S(3)/2), (sympy.S(5)/3,), sin(e + f*x)**2)*sec(e + f*x)/(4*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_305():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*sec(e + f*x)**4
    F = 3*(b*sin(e + f*x))**(sympy.S(4)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(2)/3, sympy.S(5)/2), (sympy.S(5)/3,), sin(e + f*x)**2)*sec(e + f*x)/(4*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_306():
    f = (b*sin(e + f*x))**(sympy.S(5)/3)*cos(e + f*x)**4
    F = 3*(b*sin(e + f*x))**(sympy.S(8)/3)*cos(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(4)/3), (sympy.S(7)/3,), sin(e + f*x)**2)/(8*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_307():
    f = (b*sin(e + f*x))**(sympy.S(5)/3)*cos(e + f*x)**2
    F = 3*(b*sin(e + f*x))**(sympy.S(8)/3)*cos(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(4)/3), (sympy.S(7)/3,), sin(e + f*x)**2)/(8*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_308():
    f = (b*sin(e + f*x))**(sympy.S(5)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(8)/3)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), sin(e + f*x)**2)/(8*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_309():
    f = (b*sin(e + f*x))**(sympy.S(5)/3)*sec(e + f*x)**2
    F = 3*(b*sin(e + f*x))**(sympy.S(8)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(4)/3, sympy.S(3)/2), (sympy.S(7)/3,), sin(e + f*x)**2)*sec(e + f*x)/(8*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_310():
    f = (b*sin(e + f*x))**(sympy.S(5)/3)*sec(e + f*x)**4
    F = 3*(b*sin(e + f*x))**(sympy.S(8)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(4)/3, sympy.S(5)/2), (sympy.S(7)/3,), sin(e + f*x)**2)*sec(e + f*x)/(8*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_311():
    f = cos(e + f*x)**4/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(2)/3)*cos(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(1)/3), (sympy.S(4)/3,), sin(e + f*x)**2)/(2*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_312():
    f = cos(e + f*x)**2/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(2)/3)*cos(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(1)/3), (sympy.S(4)/3,), sin(e + f*x)**2)/(2*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_313():
    f = (b*sin(e + f*x))**(sympy.S(-1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(2)/3)*cos(e + f*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), sin(e + f*x)**2)/(2*b*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_314():
    f = sec(e + f*x)**2/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(2)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(1)/3, sympy.S(3)/2), (sympy.S(4)/3,), sin(e + f*x)**2)*sec(e + f*x)/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_315():
    f = sec(e + f*x)**4/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sin(e + f*x))**(sympy.S(2)/3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(1)/3, sympy.S(5)/2), (sympy.S(4)/3,), sin(e + f*x)**2)*sec(e + f*x)/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_316():
    f = cos(e + f*x)**4/(b*sin(e + f*x))**(sympy.S(5)/3)
    F = -3*cos(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(-1)/3), (sympy.S(2)/3,), sin(e + f*x)**2)/(2*b*f*(b*sin(e + f*x))**(sympy.S(2)/3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_317():
    f = cos(e + f*x)**2/(b*sin(e + f*x))**(sympy.S(5)/3)
    F = -3*cos(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(-1)/3), (sympy.S(2)/3,), sin(e + f*x)**2)/(2*b*f*(b*sin(e + f*x))**(sympy.S(2)/3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_318():
    f = (b*sin(e + f*x))**(sympy.S(-5)/3)
    F = -3*cos(e + f*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), sin(e + f*x)**2)/(2*b*f*(b*sin(e + f*x))**(sympy.S(2)/3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_319():
    f = sec(e + f*x)**2/(b*sin(e + f*x))**(sympy.S(5)/3)
    F = -3*sqrt(cos(e + f*x)**2)*hyper((sympy.S(-1)/3, sympy.S(3)/2), (sympy.S(2)/3,), sin(e + f*x)**2)*sec(e + f*x)/(2*b*f*(b*sin(e + f*x))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_320():
    f = sec(e + f*x)**4/(b*sin(e + f*x))**(sympy.S(5)/3)
    F = -3*sqrt(cos(e + f*x)**2)*hyper((sympy.S(-1)/3, sympy.S(5)/2), (sympy.S(2)/3,), sin(e + f*x)**2)*sec(e + f*x)/(2*b*f*(b*sin(e + f*x))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_321():
    f = sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3)
    F = -log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) + log(sin(a + b*x)**(sympy.S(4)/3)/cos(a + b*x)**(sympy.S(4)/3) - sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) - sqrt(3)*atan(sqrt(3)*(-2*sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_322():
    f = sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3)
    F = sqrt(3)*log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) - sqrt(3)*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) - sqrt(3)*log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + sqrt(3)*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + atan(sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3))/b + atan(2*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) - sqrt(3))/(2*b) + atan(2*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + sqrt(3))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_323():
    f = sin(a + b*x)**(sympy.S(4)/3)/cos(a + b*x)**(sympy.S(4)/3)
    F = sqrt(3)*log(1 - sqrt(3)*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3) + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(4*b) - sqrt(3)*log(1 + sqrt(3)*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3) + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(4*b) + 3*sin(a + b*x)**(sympy.S(1)/3)/(b*cos(a + b*x)**(sympy.S(1)/3)) + atan(cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/b - atan(sqrt(3) - 2*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/(2*b) + atan(sqrt(3) + 2*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_324():
    f = sin(a + b*x)**(sympy.S(5)/3)/cos(a + b*x)**(sympy.S(5)/3)
    F = -log(1 + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(2*b) + log(1 - cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3) + cos(a + b*x)**(sympy.S(4)/3)/sin(a + b*x)**(sympy.S(4)/3))/(4*b) + 3*sin(a + b*x)**(sympy.S(2)/3)/(2*b*cos(a + b*x)**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(1 - 2*cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_325():
    f = sin(a + b*x)**(sympy.S(7)/3)/cos(a + b*x)**(sympy.S(7)/3)
    F = log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) - log(sin(a + b*x)**(sympy.S(4)/3)/cos(a + b*x)**(sympy.S(4)/3) - sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) + 3*sin(a + b*x)**(sympy.S(4)/3)/(4*b*cos(a + b*x)**(sympy.S(4)/3)) + sqrt(3)*atan(sqrt(3)*(-2*sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_326():
    f = cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3)
    F = log(1 + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(2*b) - log(1 - cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3) + cos(a + b*x)**(sympy.S(4)/3)/sin(a + b*x)**(sympy.S(4)/3))/(4*b) + sqrt(3)*atan(sqrt(3)*(1 - 2*cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_327():
    f = cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3)
    F = -sqrt(3)*log(1 - sqrt(3)*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3) + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(4*b) + sqrt(3)*log(1 + sqrt(3)*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3) + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(4*b) - atan(cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/b + atan(sqrt(3) - 2*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/(2*b) - atan(sqrt(3) + 2*cos(a + b*x)**(sympy.S(1)/3)/sin(a + b*x)**(sympy.S(1)/3))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_328():
    f = cos(a + b*x)**(sympy.S(4)/3)/sin(a + b*x)**(sympy.S(4)/3)
    F = -sqrt(3)*log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) - sqrt(3)*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + sqrt(3)*log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + sqrt(3)*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) - atan(sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3))/b - atan(2*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) - sqrt(3))/(2*b) - atan(2*sin(a + b*x)**(sympy.S(1)/3)/cos(a + b*x)**(sympy.S(1)/3) + sqrt(3))/(2*b) - 3*cos(a + b*x)**(sympy.S(1)/3)/(b*sin(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_329():
    f = cos(a + b*x)**(sympy.S(5)/3)/sin(a + b*x)**(sympy.S(5)/3)
    F = log(sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) - log(sin(a + b*x)**(sympy.S(4)/3)/cos(a + b*x)**(sympy.S(4)/3) - sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) + sqrt(3)*atan(sqrt(3)*(-2*sin(a + b*x)**(sympy.S(2)/3)/cos(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b) - 3*cos(a + b*x)**(sympy.S(2)/3)/(2*b*sin(a + b*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_330():
    f = cos(a + b*x)**(sympy.S(7)/3)/sin(a + b*x)**(sympy.S(7)/3)
    F = -log(1 + cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/(2*b) + log(1 - cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3) + cos(a + b*x)**(sympy.S(4)/3)/sin(a + b*x)**(sympy.S(4)/3))/(4*b) - sqrt(3)*atan(sqrt(3)*(1 - 2*cos(a + b*x)**(sympy.S(2)/3)/sin(a + b*x)**(sympy.S(2)/3))/3)/(2*b) - 3*cos(a + b*x)**(sympy.S(4)/3)/(4*b*sin(a + b*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_331():
    f = cos(x)**(sympy.S(2)/3)/sin(x)**(sympy.S(8)/3)
    F = -3*cos(x)**(sympy.S(5)/3)/(5*sin(x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_332():
    f = sin(x)**(sympy.S(2)/3)/cos(x)**(sympy.S(8)/3)
    F = 3*sin(x)**(sympy.S(5)/3)/(5*cos(x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_333():
    f = sin(e + f*x)**m*cos(e + f*x)**n
    F = -(sin(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)**(m - 1)*cos(e + f*x)**(n + 1)*hyper((sympy.S.Half - m/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_334():
    f = (d*cos(e + f*x))**n*sin(e + f*x)**m
    F = -(d*cos(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)**(m - 1)*hyper((sympy.S.Half - m/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_335():
    f = (b*sin(e + f*x))**m*cos(e + f*x)**n
    F = -b*(b*sin(e + f*x))**(m - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)**(n + 1)*hyper((sympy.S.Half - m/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_336():
    f = (b*sin(e + f*x))**m*(d*cos(e + f*x))**n
    F = -b*(b*sin(e + f*x))**(m - 1)*(d*cos(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*hyper((sympy.S.Half - m/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_337():
    f = (c*sin(a + b*x))**m*cos(a + b*x)**5
    F = (c*sin(a + b*x))**(m + 1)/(b*c*(m + 1)) - 2*(c*sin(a + b*x))**(m + 3)/(b*c**3*(m + 3)) + (c*sin(a + b*x))**(m + 5)/(b*c**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_338():
    f = (c*sin(a + b*x))**m*cos(a + b*x)**3
    F = (c*sin(a + b*x))**(m + 1)/(b*c*(m + 1)) - (c*sin(a + b*x))**(m + 3)/(b*c**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_339():
    f = (c*sin(a + b*x))**m*cos(a + b*x)
    F = (c*sin(a + b*x))**(m + 1)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_340():
    f = (c*sin(a + b*x))**m*sec(a + b*x)
    F = (c*sin(a + b*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_341():
    f = (c*sin(a + b*x))**m*sec(a + b*x)**3
    F = (c*sin(a + b*x))**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_342():
    f = (c*sin(a + b*x))**m*cos(a + b*x)**4
    F = (c*sin(a + b*x))**(m + 1)*cos(a + b*x)*hyper((sympy.S(-3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_343():
    f = (c*sin(a + b*x))**m*cos(a + b*x)**2
    F = (c*sin(a + b*x))**(m + 1)*cos(a + b*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_344():
    f = (c*sin(a + b*x))**m
    F = (c*sin(a + b*x))**(m + 1)*cos(a + b*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_345():
    f = (c*sin(a + b*x))**m*sec(a + b*x)**2
    F = (c*sin(a + b*x))**(m + 1)*sqrt(cos(a + b*x)**2)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)*sec(a + b*x)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_346():
    f = (c*sin(a + b*x))**m*sec(a + b*x)**4
    F = (c*sin(a + b*x))**(m + 1)*sqrt(cos(a + b*x)**2)*hyper((sympy.S(5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)*sec(a + b*x)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_347():
    f = (c*sin(a + b*x))**m*(d*cos(a + b*x))**(sympy.S(3)/2)
    F = d*(c*sin(a + b*x))**(m + 1)*sqrt(d*cos(a + b*x))*hyper((sympy.S(-1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1)*(cos(a + b*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_348():
    f = (c*sin(a + b*x))**m*sqrt(d*cos(a + b*x))
    F = d*(c*sin(a + b*x))**(m + 1)*(cos(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*sqrt(d*cos(a + b*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_349():
    f = (c*sin(a + b*x))**m/sqrt(d*cos(a + b*x))
    F = d*(c*sin(a + b*x))**(m + 1)*(cos(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(d*cos(a + b*x))**(sympy.S(3)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_350():
    f = (c*sin(a + b*x))**m/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = (c*sin(a + b*x))**(m + 1)*(cos(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(5)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*d*sqrt(d*cos(a + b*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_351():
    f = (c*sin(a + b*x))**m/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = (c*sin(a + b*x))**(m + 1)*(cos(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(7)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*d*(d*cos(a + b*x))**(sympy.S(3)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_352():
    f = (d*cos(a + b*x))**n*sin(a + b*x)**5
    F = -(d*cos(a + b*x))**(n + 1)/(b*d*(n + 1)) + 2*(d*cos(a + b*x))**(n + 3)/(b*d**3*(n + 3)) - (d*cos(a + b*x))**(n + 5)/(b*d**5*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_353():
    f = (d*cos(a + b*x))**n*sin(a + b*x)**3
    F = -(d*cos(a + b*x))**(n + 1)/(b*d*(n + 1)) + (d*cos(a + b*x))**(n + 3)/(b*d**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_354():
    f = (d*cos(a + b*x))**n*sin(a + b*x)
    F = -(d*cos(a + b*x))**(n + 1)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_355():
    f = (d*cos(a + b*x))**n*csc(a + b*x)
    F = -(d*cos(a + b*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_356():
    f = (d*cos(a + b*x))**n*csc(a + b*x)**3
    F = -(d*cos(a + b*x))**(n + 1)*hyper((2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_357():
    f = (d*cos(a + b*x))**n*csc(a + b*x)**5
    F = -(d*cos(a + b*x))**(n + 1)*hyper((3, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_358():
    f = (d*cos(a + b*x))**n*sin(a + b*x)**4
    F = -(d*cos(a + b*x))**(n + 1)*sin(a + b*x)*hyper((sympy.S(-3)/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_359():
    f = (d*cos(a + b*x))**n*sin(a + b*x)**2
    F = -(d*cos(a + b*x))**(n + 1)*sin(a + b*x)*hyper((sympy.S(-1)/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_360():
    f = (d*cos(a + b*x))**n
    F = -(d*cos(a + b*x))**(n + 1)*sin(a + b*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_361():
    f = (d*cos(a + b*x))**n*csc(a + b*x)**2
    F = -(d*cos(a + b*x))**(n + 1)*sqrt(sin(a + b*x)**2)*csc(a + b*x)*hyper((sympy.S(3)/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_362():
    f = (d*cos(a + b*x))**n*csc(a + b*x)**4
    F = -(d*cos(a + b*x))**(n + 1)*sqrt(sin(a + b*x)**2)*csc(a + b*x)*hyper((sympy.S(5)/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_363():
    f = (c*sin(a + b*x))**(sympy.S(5)/2)*(d*cos(a + b*x))**n
    F = -c*(c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**(n + 1)*hyper((sympy.S(-3)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1)*(sin(a + b*x)**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_364():
    f = (c*sin(a + b*x))**(sympy.S(3)/2)*(d*cos(a + b*x))**n
    F = -c*sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**(n + 1)*hyper((sympy.S(-1)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(n + 1)*(sin(a + b*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_365():
    f = sqrt(c*sin(a + b*x))*(d*cos(a + b*x))**n
    F = -c*(d*cos(a + b*x))**(n + 1)*(sin(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*sqrt(c*sin(a + b*x))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_366():
    f = (d*cos(a + b*x))**n/sqrt(c*sin(a + b*x))
    F = -c*(d*cos(a + b*x))**(n + 1)*(sin(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*d*(c*sin(a + b*x))**(sympy.S(3)/2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_367():
    f = (d*cos(a + b*x))**n/(c*sin(a + b*x))**(sympy.S(3)/2)
    F = -(d*cos(a + b*x))**(n + 1)*(sin(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(5)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*c*d*sqrt(c*sin(a + b*x))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_368():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**7
    F = 2*b**7/(13*f*(b*sec(e + f*x))**(sympy.S(13)/2)) - 2*b**5/(3*f*(b*sec(e + f*x))**(sympy.S(9)/2)) + 6*b**3/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2)) - 2*b/(f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_369():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**5
    F = -2*b**5/(9*f*(b*sec(e + f*x))**(sympy.S(9)/2)) + 4*b**3/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2)) - 2*b/(f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_370():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**3
    F = 2*b**3/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2)) - 2*b/(f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_371():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)
    F = -2*b/(f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_372():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)
    F = sqrt(b)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/f - sqrt(b)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_373():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)**3
    F = 3*sqrt(b)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) - 3*sqrt(b)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) - (b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**2/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_374():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)**5
    F = 21*sqrt(b)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(32*f) - 21*sqrt(b)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(32*f) - 7*(b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**2/(16*b*f) - (b*sec(e + f*x))**(sympy.S(7)/2)*cot(e + f*x)**4/(4*b**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_375():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**6
    F = -2*b*sin(e + f*x)**5/(11*f*sqrt(b*sec(e + f*x))) - 20*b*sin(e + f*x)**3/(77*f*sqrt(b*sec(e + f*x))) - 40*b*sin(e + f*x)/(77*f*sqrt(b*sec(e + f*x))) + 80*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(77*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_376():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**4
    F = -2*b*sin(e + f*x)**3/(7*f*sqrt(b*sec(e + f*x))) - 4*b*sin(e + f*x)/(7*f*sqrt(b*sec(e + f*x))) + 8*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_377():
    f = sqrt(b*sec(e + f*x))*sin(e + f*x)**2
    F = -2*b*sin(e + f*x)/(3*f*sqrt(b*sec(e + f*x))) + 4*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_378():
    f = sqrt(b*sec(e + f*x))
    F = 2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_379():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)**2
    F = -b*csc(e + f*x)/(f*sqrt(b*sec(e + f*x))) + sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_380():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)**4
    F = -b*csc(e + f*x)**3/(3*f*sqrt(b*sec(e + f*x))) - 5*b*csc(e + f*x)/(6*f*sqrt(b*sec(e + f*x))) + 5*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_381():
    f = sqrt(b*sec(e + f*x))*csc(e + f*x)**6
    F = -b*csc(e + f*x)**5/(5*f*sqrt(b*sec(e + f*x))) - 3*b*csc(e + f*x)**3/(10*f*sqrt(b*sec(e + f*x))) - 3*b*csc(e + f*x)/(4*f*sqrt(b*sec(e + f*x))) + 3*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_382():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**7
    F = 2*b**7/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) - 6*b**5/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2)) + 2*b**3/(f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 2*b*sqrt(b*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_383():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**5
    F = -2*b**5/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2)) + 4*b**3/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 2*b*sqrt(b*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_384():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**3
    F = 2*b**3/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 2*b*sqrt(b*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_385():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)
    F = 2*b*sqrt(b*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_386():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)
    F = -b**(sympy.S(3)/2)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/f - b**(sympy.S(3)/2)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/f + 2*b*sqrt(b*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_387():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**3
    F = -5*b**(sympy.S(3)/2)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) - 5*b**(sympy.S(3)/2)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) + 5*b*sqrt(b*sec(e + f*x))/(2*f) - (b*sec(e + f*x))**(sympy.S(5)/2)*cot(e + f*x)**2/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_388():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**6
    F = 20*b**3*sin(e + f*x)**3/(9*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 8*b**3*sin(e + f*x)/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 16*b**2*elliptic_e(e/2 + f*x/2, 2)/(3*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(b*sec(e + f*x))*sin(e + f*x)**5/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_389():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**4
    F = 12*b**3*sin(e + f*x)/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 24*b**2*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(b*sec(e + f*x))*sin(e + f*x)**3/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_390():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**2
    F = -4*b**2*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(b*sec(e + f*x))*sin(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_391():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(b*sec(e + f*x))*sin(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_392():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**2
    F = -3*b**2*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 3*b*sqrt(b*sec(e + f*x))*sin(e + f*x)/f - b*sqrt(b*sec(e + f*x))*csc(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_393():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**4
    F = -7*b**2*elliptic_e(e/2 + f*x/2, 2)/(2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))) + 7*b*sqrt(b*sec(e + f*x))*sin(e + f*x)/(2*f) - b*sqrt(b*sec(e + f*x))*csc(e + f*x)**3/(3*f) - 7*b*sqrt(b*sec(e + f*x))*csc(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_394():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**7
    F = 2*b**7/(9*f*(b*sec(e + f*x))**(sympy.S(9)/2)) - 6*b**5/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2)) + 6*b**3/(f*sqrt(b*sec(e + f*x))) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_395():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**5
    F = -2*b**5/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2)) + 4*b**3/(f*sqrt(b*sec(e + f*x))) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_396():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**3
    F = 2*b**3/(f*sqrt(b*sec(e + f*x))) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_397():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)
    F = 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_398():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*csc(e + f*x)
    F = b**(sympy.S(5)/2)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/f - b**(sympy.S(5)/2)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/f + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_399():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*csc(e + f*x)**3
    F = 7*b**(sympy.S(5)/2)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) - 7*b**(sympy.S(5)/2)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*f) + 7*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(6*f) - (b*sec(e + f*x))**(sympy.S(7)/2)*cot(e + f*x)**2/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_400():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*csc(e + f*x)**5
    F = 77*b**(sympy.S(5)/2)*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(32*f) - 77*b**(sympy.S(5)/2)*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(32*f) + 77*b*(b*sec(e + f*x))**(sympy.S(3)/2)/(48*f) - 11*(b*sec(e + f*x))**(sympy.S(7)/2)*cot(e + f*x)**2/(16*b*f) - (b*sec(e + f*x))**(sympy.S(11)/2)*cot(e + f*x)**4/(4*b**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_401():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**6
    F = 20*b**3*sin(e + f*x)**3/(21*f*sqrt(b*sec(e + f*x))) + 40*b**3*sin(e + f*x)/(21*f*sqrt(b*sec(e + f*x))) - 80*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*f) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**5/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_402():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**4
    F = 4*b**3*sin(e + f*x)/(3*f*sqrt(b*sec(e + f*x))) - 8*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_403():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)**2
    F = -4*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_404():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)
    F = 2*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_405():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*csc(e + f*x)**2
    F = -5*b**3*csc(e + f*x)/(3*f*sqrt(b*sec(e + f*x))) + 5*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f) + 2*b*(b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_406():
    f = (b*sec(e + f*x))**(sympy.S(5)/2)*csc(e + f*x)**4
    F = -5*b**3*csc(e + f*x)/(2*f*sqrt(b*sec(e + f*x))) + 5*b**2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(2*f) - b*(b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**3/(3*f) + b*(b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_407():
    f = sin(e + f*x)**7/sqrt(b*sec(e + f*x))
    F = 2*b**7/(15*f*(b*sec(e + f*x))**(sympy.S(15)/2)) - 6*b**5/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) + 6*b**3/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2)) - 2*b/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_408():
    f = sin(e + f*x)**5/sqrt(b*sec(e + f*x))
    F = -2*b**5/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) + 4*b**3/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2)) - 2*b/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_409():
    f = sin(e + f*x)**3/sqrt(b*sec(e + f*x))
    F = 2*b**3/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2)) - 2*b/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_410():
    f = sin(e + f*x)/sqrt(b*sec(e + f*x))
    F = -2*b/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_411():
    f = csc(e + f*x)/sqrt(b*sec(e + f*x))
    F = -atan(sqrt(b*sec(e + f*x))/sqrt(b))/(sqrt(b)*f) - atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_412():
    f = csc(e + f*x)**3/sqrt(b*sec(e + f*x))
    F = -sqrt(b*sec(e + f*x))*cot(e + f*x)**2/(2*b*f) - atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*sqrt(b)*f) - atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_413():
    f = csc(e + f*x)**5/sqrt(b*sec(e + f*x))
    F = -5*sqrt(b*sec(e + f*x))*cot(e + f*x)**2/(16*b*f) - (b*sec(e + f*x))**(sympy.S(5)/2)*cot(e + f*x)**4/(4*b**3*f) - 5*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(32*sqrt(b)*f) - 5*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(32*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_414():
    f = sin(e + f*x)**6/sqrt(b*sec(e + f*x))
    F = -2*b*sin(e + f*x)**5/(13*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 20*b*sin(e + f*x)**3/(117*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 8*b*sin(e + f*x)/(39*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 16*elliptic_e(e/2 + f*x/2, 2)/(39*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_415():
    f = sin(e + f*x)**4/sqrt(b*sec(e + f*x))
    F = -2*b*sin(e + f*x)**3/(9*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 4*b*sin(e + f*x)/(15*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 8*elliptic_e(e/2 + f*x/2, 2)/(15*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_416():
    f = sin(e + f*x)**2/sqrt(b*sec(e + f*x))
    F = -2*b*sin(e + f*x)/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 4*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_417():
    f = 1/sqrt(b*sec(e + f*x))
    F = 2*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_418():
    f = csc(e + f*x)**2/sqrt(b*sec(e + f*x))
    F = -b*csc(e + f*x)/(f*(b*sec(e + f*x))**(sympy.S(3)/2)) - elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_419():
    f = csc(e + f*x)**4/sqrt(b*sec(e + f*x))
    F = -b*csc(e + f*x)**3/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - b*csc(e + f*x)/(2*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - elliptic_e(e/2 + f*x/2, 2)/(2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_420():
    f = csc(e + f*x)**6/sqrt(b*sec(e + f*x))
    F = -b*csc(e + f*x)**5/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 7*b*csc(e + f*x)**3/(30*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 7*b*csc(e + f*x)/(20*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 7*elliptic_e(e/2 + f*x/2, 2)/(20*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_421():
    f = sin(e + f*x)**7/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*b**7/(17*f*(b*sec(e + f*x))**(sympy.S(17)/2)) - 6*b**5/(13*f*(b*sec(e + f*x))**(sympy.S(13)/2)) + 2*b**3/(3*f*(b*sec(e + f*x))**(sympy.S(9)/2)) - 2*b/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_422():
    f = sin(e + f*x)**5/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -2*b**5/(13*f*(b*sec(e + f*x))**(sympy.S(13)/2)) + 4*b**3/(9*f*(b*sec(e + f*x))**(sympy.S(9)/2)) - 2*b/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_423():
    f = sin(e + f*x)**3/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*b**3/(9*f*(b*sec(e + f*x))**(sympy.S(9)/2)) - 2*b/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_424():
    f = sin(e + f*x)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -2*b/(5*f*(b*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_425():
    f = csc(e + f*x)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = 2/(b*f*sqrt(b*sec(e + f*x))) + atan(sqrt(b*sec(e + f*x))/sqrt(b))/(b**(sympy.S(3)/2)*f) - atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_426():
    f = csc(e + f*x)**3/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -(b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**2/(2*b**3*f) - atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*b**(sympy.S(3)/2)*f) + atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_427():
    f = csc(e + f*x)**5/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -(b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**4/(4*b**3*f) - 3*(b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**2/(16*b**3*f) - 3*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(32*b**(sympy.S(3)/2)*f) + 3*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(32*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_428():
    f = sin(e + f*x)**4/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -2*b*sin(e + f*x)**3/(11*f*(b*sec(e + f*x))**(sympy.S(5)/2)) - 12*b*sin(e + f*x)/(77*f*(b*sec(e + f*x))**(sympy.S(5)/2)) + 8*sin(e + f*x)/(77*b*f*sqrt(b*sec(e + f*x))) + 8*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(77*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_429():
    f = sin(e + f*x)**2/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -2*b*sin(e + f*x)/(7*f*(b*sec(e + f*x))**(sympy.S(5)/2)) + 4*sin(e + f*x)/(21*b*f*sqrt(b*sec(e + f*x))) + 4*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_430():
    f = (b*sec(e + f*x))**(sympy.S(-3)/2)
    F = 2*sin(e + f*x)/(3*b*f*sqrt(b*sec(e + f*x))) + 2*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_431():
    f = csc(e + f*x)**2/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -csc(e + f*x)/(b*f*sqrt(b*sec(e + f*x))) - sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_432():
    f = csc(e + f*x)**4/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -csc(e + f*x)**3/(3*b*f*sqrt(b*sec(e + f*x))) + csc(e + f*x)/(6*b*f*sqrt(b*sec(e + f*x))) - sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(6*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_433():
    f = csc(e + f*x)**6/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = -csc(e + f*x)**5/(5*b*f*sqrt(b*sec(e + f*x))) + csc(e + f*x)**3/(30*b*f*sqrt(b*sec(e + f*x))) + csc(e + f*x)/(12*b*f*sqrt(b*sec(e + f*x))) - sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(12*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_434():
    f = sin(e + f*x)**7/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = 2*b**7/(19*f*(b*sec(e + f*x))**(sympy.S(19)/2)) - 2*b**5/(5*f*(b*sec(e + f*x))**(sympy.S(15)/2)) + 6*b**3/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) - 2*b/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_435():
    f = sin(e + f*x)**5/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -2*b**5/(15*f*(b*sec(e + f*x))**(sympy.S(15)/2)) + 4*b**3/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) - 2*b/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_436():
    f = sin(e + f*x)**3/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = 2*b**3/(11*f*(b*sec(e + f*x))**(sympy.S(11)/2)) - 2*b/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_437():
    f = sin(e + f*x)/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -2*b/(7*f*(b*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_438():
    f = csc(e + f*x)/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = 2/(3*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - atan(sqrt(b*sec(e + f*x))/sqrt(b))/(b**(sympy.S(5)/2)*f) - atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_439():
    f = csc(e + f*x)**3/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(b*sec(e + f*x))*cot(e + f*x)**2/(2*b**3*f) + 3*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(4*b**(sympy.S(5)/2)*f) + 3*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(4*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_440():
    f = csc(e + f*x)**5/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(b*sec(e + f*x))*cot(e + f*x)**4/(4*b**3*f) - sqrt(b*sec(e + f*x))*cot(e + f*x)**2/(16*b**3*f) + 3*atan(sqrt(b*sec(e + f*x))/sqrt(b))/(32*b**(sympy.S(5)/2)*f) + 3*atanh(sqrt(b*sec(e + f*x))/sqrt(b))/(32*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_441():
    f = sin(e + f*x)**4/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -2*b*sin(e + f*x)**3/(13*f*(b*sec(e + f*x))**(sympy.S(7)/2)) - 4*b*sin(e + f*x)/(39*f*(b*sec(e + f*x))**(sympy.S(7)/2)) + 8*sin(e + f*x)/(195*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 8*elliptic_e(e/2 + f*x/2, 2)/(65*b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_442():
    f = sin(e + f*x)**2/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -2*b*sin(e + f*x)/(9*f*(b*sec(e + f*x))**(sympy.S(7)/2)) + 4*sin(e + f*x)/(45*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 4*elliptic_e(e/2 + f*x/2, 2)/(15*b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_443():
    f = (b*sec(e + f*x))**(sympy.S(-5)/2)
    F = 2*sin(e + f*x)/(5*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 6*elliptic_e(e/2 + f*x/2, 2)/(5*b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_444():
    f = csc(e + f*x)**2/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -csc(e + f*x)/(b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 3*elliptic_e(e/2 + f*x/2, 2)/(b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_445():
    f = csc(e + f*x)**4/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -csc(e + f*x)**3/(3*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + csc(e + f*x)/(2*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + elliptic_e(e/2 + f*x/2, 2)/(2*b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_446():
    f = csc(e + f*x)**6/(b*sec(e + f*x))**(sympy.S(5)/2)
    F = -csc(e + f*x)**5/(5*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + csc(e + f*x)**3/(10*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 3*csc(e + f*x)/(20*b*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 3*elliptic_e(e/2 + f*x/2, 2)/(20*b**2*f*sqrt(b*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_447():
    f = (a*sin(e + f*x))**(sympy.S(9)/2)*sqrt(b*sec(e + f*x))
    F = 21*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(128*sqrt(b)*f) - 21*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(128*sqrt(b)*f) - 21*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(64*sqrt(b)*f) + 21*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(64*sqrt(b)*f) - 7*a**3*b*(a*sin(e + f*x))**(sympy.S(3)/2)/(16*f*sqrt(b*sec(e + f*x))) - a*b*(a*sin(e + f*x))**(sympy.S(7)/2)/(4*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_448():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*sec(e + f*x))
    F = 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(16*sqrt(b)*f) - 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(16*sqrt(b)*f) - 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(8*sqrt(b)*f) + 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(8*sqrt(b)*f) - a*b*(a*sin(e + f*x))**(sympy.S(3)/2)/(2*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_449():
    f = sqrt(a*sin(e + f*x))*sqrt(b*sec(e + f*x))
    F = sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(4*sqrt(b)*f) - sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(4*sqrt(b)*f) - sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(2*sqrt(b)*f) + sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_450():
    f = sqrt(b*sec(e + f*x))/(a*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*b/(a*f*sqrt(a*sin(e + f*x))*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_451():
    f = sqrt(b*sec(e + f*x))/(a*sin(e + f*x))**(sympy.S(7)/2)
    F = -2*b/(5*a*f*(a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*sec(e + f*x))) - 8*b/(5*a**3*f*sqrt(a*sin(e + f*x))*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_452():
    f = sqrt(b*sec(e + f*x))/(a*sin(e + f*x))**(sympy.S(11)/2)
    F = -2*b/(9*a*f*(a*sin(e + f*x))**(sympy.S(9)/2)*sqrt(b*sec(e + f*x))) - 16*b/(45*a**3*f*(a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*sec(e + f*x))) - 64*b/(45*a**5*f*sqrt(a*sin(e + f*x))*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_453():
    f = (a*sin(e + f*x))**(sympy.S(7)/2)*sqrt(b*sec(e + f*x))
    F = 5*a**4*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(12*f*sqrt(a*sin(e + f*x))) - 5*a**3*b*sqrt(a*sin(e + f*x))/(6*f*sqrt(b*sec(e + f*x))) - a*b*(a*sin(e + f*x))**(sympy.S(5)/2)/(3*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_454():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))
    F = a**2*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(2*f*sqrt(a*sin(e + f*x))) - a*b*sqrt(a*sin(e + f*x))/(f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_455():
    f = sqrt(b*sec(e + f*x))/sqrt(a*sin(e + f*x))
    F = sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_456():
    f = sqrt(b*sec(e + f*x))/(a*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*b/(3*a*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))) + 2*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(3*a**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_457():
    f = sqrt(b*sec(e + f*x))/(a*sin(e + f*x))**(sympy.S(9)/2)
    F = -2*b/(7*a*f*(a*sin(e + f*x))**(sympy.S(7)/2)*sqrt(b*sec(e + f*x))) - 4*b/(7*a**3*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))) + 4*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(7*a**4*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_458():
    f = sin(e + f*x)**(sympy.S(9)/2)/sqrt(b*sec(e + f*x))
    F = -b*sin(e + f*x)**(sympy.S(7)/2)/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)) - 7*b*sin(e + f*x)**(sympy.S(3)/2)/(30*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + 7*sqrt(sin(e + f*x))*elliptic_e(e + f*x - pi/4, 2)/(20*f*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_459():
    f = sin(e + f*x)**(sympy.S(5)/2)/sqrt(b*sec(e + f*x))
    F = -b*sin(e + f*x)**(sympy.S(3)/2)/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)) + sqrt(sin(e + f*x))*elliptic_e(e + f*x - pi/4, 2)/(2*f*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_460():
    f = sqrt(sin(e + f*x))/sqrt(b*sec(e + f*x))
    F = sqrt(sin(e + f*x))*elliptic_e(e + f*x - pi/4, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_461():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(3)/2))
    F = -2*b/(f*(b*sec(e + f*x))**(sympy.S(3)/2)*sqrt(sin(e + f*x))) - 2*sqrt(sin(e + f*x))*elliptic_e(e + f*x - pi/4, 2)/(f*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_462():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(7)/2))
    F = -4*b/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sqrt(sin(e + f*x))) - 2*b/(5*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(5)/2)) - 4*sqrt(sin(e + f*x))*elliptic_e(e + f*x - pi/4, 2)/(5*f*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_463():
    f = sin(e + f*x)**(sympy.S(3)/2)/sqrt(b*sec(e + f*x))
    F = -sqrt(2)*sqrt(b)*log(sqrt(b)*cot(e + f*x) + sqrt(b) - sqrt(2)*sqrt(b*cos(e + f*x))/sqrt(sin(e + f*x)))/(16*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) + sqrt(2)*sqrt(b)*log(sqrt(b)*cot(e + f*x) + sqrt(b) + sqrt(2)*sqrt(b*cos(e + f*x))/sqrt(sin(e + f*x)))/(16*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) + sqrt(2)*sqrt(b)*atan(1 - sqrt(2)*sqrt(b*cos(e + f*x))/(sqrt(b)*sqrt(sin(e + f*x))))/(8*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) - sqrt(2)*sqrt(b)*atan(1 + sqrt(2)*sqrt(b*cos(e + f*x))/(sqrt(b)*sqrt(sin(e + f*x))))/(8*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) - b*sqrt(sin(e + f*x))/(2*f*(b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_464():
    f = 1/(sqrt(b*sec(e + f*x))*sqrt(sin(e + f*x)))
    F = -sqrt(2)*sqrt(b)*log(sqrt(b)*cot(e + f*x) + sqrt(b) - sqrt(2)*sqrt(b*cos(e + f*x))/sqrt(sin(e + f*x)))/(4*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) + sqrt(2)*sqrt(b)*log(sqrt(b)*cot(e + f*x) + sqrt(b) + sqrt(2)*sqrt(b*cos(e + f*x))/sqrt(sin(e + f*x)))/(4*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) + sqrt(2)*sqrt(b)*atan(1 - sqrt(2)*sqrt(b*cos(e + f*x))/(sqrt(b)*sqrt(sin(e + f*x))))/(2*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))) - sqrt(2)*sqrt(b)*atan(1 + sqrt(2)*sqrt(b*cos(e + f*x))/(sqrt(b)*sqrt(sin(e + f*x))))/(2*f*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_465():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(5)/2))
    F = -2*b/(3*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_466():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(9)/2))
    F = -8*b/(21*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(3)/2)) - 2*b/(7*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_467():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(13)/2))
    F = -64*b/(231*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(3)/2)) - 16*b/(77*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(7)/2)) - 2*b/(11*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_468():
    f = 1/(sqrt(b*sec(e + f*x))*sin(e + f*x)**(sympy.S(17)/2))
    F = -256*b/(1155*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(3)/2)) - 64*b/(385*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(7)/2)) - 8*b/(55*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(11)/2)) - 2*b/(15*f*(b*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_469():
    f = (a*sin(e + f*x))**(sympy.S(9)/2)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = 7*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(512*b**(sympy.S(5)/2)*f) - 7*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(512*b**(sympy.S(5)/2)*f) - 7*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(256*b**(sympy.S(5)/2)*f) + 7*sqrt(2)*a**(sympy.S(9)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(256*b**(sympy.S(5)/2)*f) - 7*a**3*(a*sin(e + f*x))**(sympy.S(3)/2)/(192*b*f*sqrt(b*sec(e + f*x))) - a*(a*sin(e + f*x))**(sympy.S(7)/2)/(48*b*f*sqrt(b*sec(e + f*x))) + (a*sin(e + f*x))**(sympy.S(11)/2)/(6*a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_470():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(128*b**(sympy.S(5)/2)*f) - 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(128*b**(sympy.S(5)/2)*f) - 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(64*b**(sympy.S(5)/2)*f) + 3*sqrt(2)*a**(sympy.S(5)/2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(64*b**(sympy.S(5)/2)*f) - a*(a*sin(e + f*x))**(sympy.S(3)/2)/(16*b*f*sqrt(b*sec(e + f*x))) + (a*sin(e + f*x))**(sympy.S(7)/2)/(4*a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_471():
    f = sqrt(a*sin(e + f*x))/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(16*b**(sympy.S(5)/2)*f) - sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(16*b**(sympy.S(5)/2)*f) - sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(8*b**(sympy.S(5)/2)*f) + sqrt(2)*sqrt(a)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(8*b**(sympy.S(5)/2)*f) + (a*sin(e + f*x))**(sympy.S(3)/2)/(2*a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_472():
    f = 1/((a*sin(e + f*x))**(sympy.S(3)/2)*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = -2/(a*b*f*sqrt(a*sin(e + f*x))*sqrt(b*sec(e + f*x))) - sqrt(2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(4*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*f) + sqrt(2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*log(sqrt(a)*tan(e + f*x) + sqrt(a) + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/sqrt(b*cos(e + f*x)))/(4*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*f) + sqrt(2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 - sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*f) - sqrt(2)*sqrt(b*cos(e + f*x))*sqrt(b*sec(e + f*x))*atan(1 + sqrt(2)*sqrt(b)*sqrt(a*sin(e + f*x))/(sqrt(a)*sqrt(b*cos(e + f*x))))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_473():
    f = 1/((a*sin(e + f*x))**(sympy.S(7)/2)*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = -2*b/(5*a*f*(a*sin(e + f*x))**(sympy.S(5)/2)*(b*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_474():
    f = (a*sin(e + f*x))**(sympy.S(7)/2)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = a**4*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(24*b**2*f*sqrt(a*sin(e + f*x))) - a**3*sqrt(a*sin(e + f*x))/(12*b*f*sqrt(b*sec(e + f*x))) - a*(a*sin(e + f*x))**(sympy.S(5)/2)/(30*b*f*sqrt(b*sec(e + f*x))) + (a*sin(e + f*x))**(sympy.S(9)/2)/(5*a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_475():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)/(b*sec(e + f*x))**(sympy.S(3)/2)
    F = a**2*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(12*b**2*f*sqrt(a*sin(e + f*x))) - a*sqrt(a*sin(e + f*x))/(6*b*f*sqrt(b*sec(e + f*x))) + (a*sin(e + f*x))**(sympy.S(5)/2)/(3*a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_476():
    f = 1/(sqrt(a*sin(e + f*x))*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(2*b**2*f*sqrt(a*sin(e + f*x))) + sqrt(a*sin(e + f*x))/(a*b*f*sqrt(b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_477():
    f = 1/((a*sin(e + f*x))**(sympy.S(5)/2)*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = -2/(3*a*b*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))) - sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(3*a**2*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_478():
    f = 1/((a*sin(e + f*x))**(sympy.S(9)/2)*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = -2/(7*a*b*f*(a*sin(e + f*x))**(sympy.S(7)/2)*sqrt(b*sec(e + f*x))) + 2/(21*a**3*b*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))) - 2*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(21*a**4*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_479():
    f = 1/((a*sin(e + f*x))**(sympy.S(13)/2)*(b*sec(e + f*x))**(sympy.S(3)/2))
    F = -2/(11*a*b*f*(a*sin(e + f*x))**(sympy.S(11)/2)*sqrt(b*sec(e + f*x))) + 2/(77*a**3*b*f*(a*sin(e + f*x))**(sympy.S(7)/2)*sqrt(b*sec(e + f*x))) + 4/(77*a**5*b*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*sec(e + f*x))) - 4*sqrt(b*sec(e + f*x))*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)/(77*a**6*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_480():
    f = (c*sin(a + b*x))**m*(d*sec(a + b*x))**(sympy.S(5)/2)
    F = d*(c*sin(a + b*x))**(m + 1)*(d*sec(a + b*x))**(sympy.S(3)/2)*(cos(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(7)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_481():
    f = (c*sin(a + b*x))**m*(d*sec(a + b*x))**(sympy.S(3)/2)
    F = d*(c*sin(a + b*x))**(m + 1)*sqrt(d*sec(a + b*x))*(cos(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(5)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_482():
    f = (c*sin(a + b*x))**m*sqrt(d*sec(a + b*x))
    F = (c*sin(a + b*x))**(m + 1)*(d*sec(a + b*x))**(sympy.S(3)/2)*(cos(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_483():
    f = (c*sin(a + b*x))**m/sqrt(d*sec(a + b*x))
    F = (c*sin(a + b*x))**(m + 1)*sqrt(d*sec(a + b*x))*(cos(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_484():
    f = (c*sin(a + b*x))**m/(d*sec(a + b*x))**(sympy.S(3)/2)
    F = (c*sin(a + b*x))**(m + 1)*hyper((sympy.S(-1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*c*d*sqrt(d*sec(a + b*x))*(m + 1)*(cos(a + b*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_485():
    f = sin(e + f*x)**m*sec(e + f*x)**n
    F = -(sin(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)**(m - 1)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_486():
    f = (a*sin(e + f*x))**m*sec(e + f*x)**n
    F = -a*(a*sin(e + f*x))**(m - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_487():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**m
    F = -b*(b*sec(e + f*x))**(n - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)**(m - 1)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_488():
    f = (a*sin(e + f*x))**m*(b*sec(e + f*x))**n
    F = -a*b*(a*sin(e + f*x))**(m - 1)*(b*sec(e + f*x))**(n - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_489():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**5
    F = -b**5*(b*sec(e + f*x))**(n - 5)/(f*(5 - n)) + 2*b**3*(b*sec(e + f*x))**(n - 3)/(f*(3 - n)) - b*(b*sec(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_490():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**3
    F = b**3*(b*sec(e + f*x))**(n - 3)/(f*(3 - n)) - b*(b*sec(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_491():
    f = (b*sec(e + f*x))**n*sin(e + f*x)
    F = -b*(b*sec(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_492():
    f = (b*sec(e + f*x))**n*csc(e + f*x)
    F = -(b*sec(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sec(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_493():
    f = (b*sec(e + f*x))**n*csc(e + f*x)**3
    F = (b*sec(e + f*x))**(n + 3)*hyper((2, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), sec(e + f*x)**2)/(b**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_494():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**6
    F = -b*(b*sec(e + f*x))**(n - 1)*sin(e + f*x)*hyper((sympy.S(-5)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_495():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**4
    F = -b*(b*sec(e + f*x))**(n - 1)*sin(e + f*x)*hyper((sympy.S(-3)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_496():
    f = (b*sec(e + f*x))**n*sin(e + f*x)**2
    F = -b*(b*sec(e + f*x))**(n - 1)*sin(e + f*x)*hyper((sympy.S(-1)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_497():
    f = (b*sec(e + f*x))**n
    F = -b*(b*sec(e + f*x))**(n - 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_498():
    f = (b*sec(e + f*x))**n*csc(e + f*x)**2
    F = -b*(b*sec(e + f*x))**(n - 1)*sqrt(sin(e + f*x)**2)*csc(e + f*x)*hyper((sympy.S(3)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_499():
    f = (b*sec(e + f*x))**n*csc(e + f*x)**4
    F = -b*(b*sec(e + f*x))**(n - 1)*sqrt(sin(e + f*x)**2)*csc(e + f*x)*hyper((sympy.S(5)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_500():
    f = (b*sec(a + b*x))**n*(c*sin(a + b*x))**(sympy.S(3)/2)
    F = -c*(b*sec(a + b*x))**(n - 1)*sqrt(c*sin(a + b*x))*hyper((sympy.S(-1)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)/((1 - n)*(sin(a + b*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_501():
    f = (b*sec(a + b*x))**n*sqrt(c*sin(a + b*x))
    F = -c*(b*sec(a + b*x))**(n - 1)*(sin(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)/(sqrt(c*sin(a + b*x))*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_502():
    f = (b*sec(a + b*x))**n/sqrt(c*sin(a + b*x))
    F = -c*(b*sec(a + b*x))**(n - 1)*(sin(a + b*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)/((c*sin(a + b*x))**(sympy.S(3)/2)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_503():
    f = (b*sec(a + b*x))**n/(c*sin(a + b*x))**(sympy.S(3)/2)
    F = -(b*sec(a + b*x))**(n - 1)*(sin(a + b*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(5)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)/(c*sqrt(c*sin(a + b*x))*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_504():
    f = sqrt(d*csc(e + f*x))*sin(e + f*x)**4
    F = -2*d**3*cos(e + f*x)/(7*f*(d*csc(e + f*x))**(sympy.S(5)/2)) - 10*d*cos(e + f*x)/(21*f*sqrt(d*csc(e + f*x))) + 10*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(21*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_505():
    f = sqrt(d*csc(e + f*x))*sin(e + f*x)**3
    F = -2*d**2*cos(e + f*x)/(5*f*(d*csc(e + f*x))**(sympy.S(3)/2)) + 6*d*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_506():
    f = sqrt(d*csc(e + f*x))*sin(e + f*x)**2
    F = -2*d*cos(e + f*x)/(3*f*sqrt(d*csc(e + f*x))) + 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_507():
    f = sqrt(d*csc(e + f*x))*sin(e + f*x)
    F = 2*d*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_508():
    f = sqrt(d*csc(e + f*x))
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_509():
    f = sqrt(d*csc(e + f*x))*csc(e + f*x)
    F = -2*d*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 2*sqrt(d*csc(e + f*x))*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_510():
    f = sqrt(d*csc(e + f*x))*csc(e + f*x)**2
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*f) - 2*(d*csc(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_511():
    f = sqrt(d*csc(e + f*x))*csc(e + f*x)**3
    F = -6*d*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 6*sqrt(d*csc(e + f*x))*cos(e + f*x)/(5*f) - 2*(d*csc(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(5*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_512():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**5
    F = -2*d**4*cos(e + f*x)/(7*f*(d*csc(e + f*x))**(sympy.S(5)/2)) - 10*d**2*cos(e + f*x)/(21*f*sqrt(d*csc(e + f*x))) + 10*d*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(21*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_513():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**4
    F = -2*d**3*cos(e + f*x)/(5*f*(d*csc(e + f*x))**(sympy.S(3)/2)) + 6*d**2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_514():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**3
    F = -2*d**2*cos(e + f*x)/(3*f*sqrt(d*csc(e + f*x))) + 2*d*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_515():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)**2
    F = 2*d**2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_516():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)
    F = 2*d*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_517():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)
    F = -2*d**2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 2*d*sqrt(d*csc(e + f*x))*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_518():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)
    F = 2*d*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*f) - 2*(d*csc(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_519():
    f = (d*csc(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**2
    F = -6*d**2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 6*d*sqrt(d*csc(e + f*x))*cos(e + f*x)/(5*f) - 2*(d*csc(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_520():
    f = sin(e + f*x)**3/sqrt(d*csc(e + f*x))
    F = -2*d**2*cos(e + f*x)/(7*f*(d*csc(e + f*x))**(sympy.S(5)/2)) - 10*cos(e + f*x)/(21*f*sqrt(d*csc(e + f*x))) + 10*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(21*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_521():
    f = sin(e + f*x)**2/sqrt(d*csc(e + f*x))
    F = -2*d*cos(e + f*x)/(5*f*(d*csc(e + f*x))**(sympy.S(3)/2)) + 6*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_522():
    f = sin(e + f*x)/sqrt(d*csc(e + f*x))
    F = -2*cos(e + f*x)/(3*f*sqrt(d*csc(e + f*x))) + 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_523():
    f = 1/sqrt(d*csc(e + f*x))
    F = 2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_524():
    f = csc(e + f*x)/sqrt(d*csc(e + f*x))
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_525():
    f = csc(e + f*x)**2/sqrt(d*csc(e + f*x))
    F = -2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 2*sqrt(d*csc(e + f*x))*cos(e + f*x)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_526():
    f = csc(e + f*x)**3/sqrt(d*csc(e + f*x))
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d*f) - 2*(d*csc(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(3*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_527():
    f = sin(e + f*x)**2/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = -2*d*cos(e + f*x)/(7*f*(d*csc(e + f*x))**(sympy.S(5)/2)) - 10*cos(e + f*x)/(21*d*f*sqrt(d*csc(e + f*x))) + 10*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(21*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_528():
    f = sin(e + f*x)/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = -2*cos(e + f*x)/(5*f*(d*csc(e + f*x))**(sympy.S(3)/2)) + 6*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*d*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_529():
    f = (d*csc(e + f*x))**(sympy.S(-3)/2)
    F = -2*cos(e + f*x)/(3*d*f*sqrt(d*csc(e + f*x))) + 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_530():
    f = csc(e + f*x)/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = 2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(d*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_531():
    f = csc(e + f*x)**2/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_532():
    f = csc(e + f*x)**3/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = -2*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(d*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 2*sqrt(d*csc(e + f*x))*cos(e + f*x)/(d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_533():
    f = csc(e + f*x)**4/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = 2*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d**2*f) - 2*(d*csc(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(3*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_534():
    f = csc(e + f*x)**5/(d*csc(e + f*x))**(sympy.S(3)/2)
    F = -6*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*d*f*sqrt(d*csc(e + f*x))*sqrt(sin(e + f*x))) - 6*sqrt(d*csc(e + f*x))*cos(e + f*x)/(5*d**2*f) - 2*(d*csc(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(5*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_0_a_sin_pow_m_b_trg_pow_n_535():
    f = (a*sin(e + f*x))**m*(b*csc(e + f*x))**n
    F = (a*sin(e + f*x))**(m + 1)*(b*csc(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, m/2 - n/2 + sympy.S.Half), (m/2 - n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(a*f*(m - n + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F

