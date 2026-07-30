"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.0 (a sec)^m (b trg)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_1():
    f = sec(a + b*x)
    F = atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_2():
    f = sec(a + b*x)**2
    F = tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_3():
    f = sec(a + b*x)**3
    F = tan(a + b*x)*sec(a + b*x)/(2*b) + atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_4():
    f = sec(a + b*x)**4
    F = tan(a + b*x)**3/(3*b) + tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_5():
    f = sec(a + b*x)**5
    F = tan(a + b*x)*sec(a + b*x)**3/(4*b) + 3*tan(a + b*x)*sec(a + b*x)/(8*b) + 3*atanh(sin(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_6():
    f = sec(a + b*x)**6
    F = tan(a + b*x)**5/(5*b) + 2*tan(a + b*x)**3/(3*b) + tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_7():
    f = sec(a + b*x)**7
    F = tan(a + b*x)*sec(a + b*x)**5/(6*b) + 5*tan(a + b*x)*sec(a + b*x)**3/(24*b) + 5*tan(a + b*x)*sec(a + b*x)/(16*b) + 5*atanh(sin(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_8():
    f = sec(a + b*x)**8
    F = tan(a + b*x)**7/(7*b) + 3*tan(a + b*x)**5/(5*b) + tan(a + b*x)**3/b + tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_9():
    f = sec(a + b*x)**(sympy.S(7)/2)
    F = 2*sin(a + b*x)*sec(a + b*x)**(sympy.S(5)/2)/(5*b) + 6*sin(a + b*x)*sqrt(sec(a + b*x))/(5*b) - 6*sqrt(cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_10():
    f = sec(a + b*x)**(sympy.S(5)/2)
    F = 2*sin(a + b*x)*sec(a + b*x)**(sympy.S(3)/2)/(3*b) + 2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_11():
    f = sec(a + b*x)**(sympy.S(3)/2)
    F = 2*sin(a + b*x)*sqrt(sec(a + b*x))/b - 2*sqrt(cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_12():
    f = sqrt(sec(a + b*x))
    F = 2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_13():
    f = 1/sqrt(sec(a + b*x))
    F = 2*sqrt(cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_14():
    f = sec(a + b*x)**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)/(3*b*sqrt(sec(a + b*x))) + 2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_15():
    f = sec(a + b*x)**(sympy.S(-5)/2)
    F = 2*sin(a + b*x)/(5*b*sec(a + b*x)**(sympy.S(3)/2)) + 6*sqrt(cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_16():
    f = sec(a + b*x)**(sympy.S(-7)/2)
    F = 10*sin(a + b*x)/(21*b*sqrt(sec(a + b*x))) + 2*sin(a + b*x)/(7*b*sec(a + b*x)**(sympy.S(5)/2)) + 10*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)*sqrt(sec(a + b*x))/(21*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_17():
    f = (c*sec(a + b*x))**(sympy.S(7)/2)
    F = -6*c**4*elliptic_e(a/2 + b*x/2, 2)/(5*b*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))) + 6*c**3*sqrt(c*sec(a + b*x))*sin(a + b*x)/(5*b) + 2*c*(c*sec(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_18():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)
    F = 2*c**2*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b) + 2*c*(c*sec(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_19():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)
    F = -2*c**2*elliptic_e(a/2 + b*x/2, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))) + 2*c*sqrt(c*sec(a + b*x))*sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_20():
    f = sqrt(c*sec(a + b*x))
    F = 2*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_21():
    f = 1/sqrt(c*sec(a + b*x))
    F = 2*elliptic_e(a/2 + b*x/2, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_22():
    f = (c*sec(a + b*x))**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)/(3*b*c*sqrt(c*sec(a + b*x))) + 2*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_23():
    f = (c*sec(a + b*x))**(sympy.S(-5)/2)
    F = 2*sin(a + b*x)/(5*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)) + 6*elliptic_e(a/2 + b*x/2, 2)/(5*b*c**2*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_24():
    f = (c*sec(a + b*x))**(sympy.S(-7)/2)
    F = 2*sin(a + b*x)/(7*b*c*(c*sec(a + b*x))**(sympy.S(5)/2)) + 10*sin(a + b*x)/(21*b*c**3*sqrt(c*sec(a + b*x))) + 10*sqrt(c*sec(a + b*x))*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(21*b*c**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_25():
    f = sec(a + b*x)**(sympy.S(4)/3)
    F = 3*sin(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(a + b*x)**2)*sec(a + b*x)**(sympy.S(1)/3)/(b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_26():
    f = sec(a + b*x)**(sympy.S(2)/3)
    F = -3*sin(a + b*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(a + b*x)**2)/(b*sqrt(sin(a + b*x)**2)*sec(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_27():
    f = sec(a + b*x)**(sympy.S(1)/3)
    F = -3*sin(a + b*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(a + b*x)**2)/(2*b*sqrt(sin(a + b*x)**2)*sec(a + b*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_28():
    f = sec(a + b*x)**(sympy.S(-1)/3)
    F = -3*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(a + b*x)**2)/(4*b*sqrt(sin(a + b*x)**2)*sec(a + b*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_29():
    f = sec(a + b*x)**(sympy.S(-2)/3)
    F = -3*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(a + b*x)**2)/(5*b*sqrt(sin(a + b*x)**2)*sec(a + b*x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_30():
    f = sec(a + b*x)**(sympy.S(-4)/3)
    F = -3*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(a + b*x)**2)/(7*b*sqrt(sin(a + b*x)**2)*sec(a + b*x)**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_31():
    f = (c*sec(a + b*x))**(sympy.S(4)/3)
    F = 3*c*(c*sec(a + b*x))**(sympy.S(1)/3)*sin(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(a + b*x)**2)/(b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_32():
    f = (c*sec(a + b*x))**(sympy.S(2)/3)
    F = -3*c*sin(a + b*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(a + b*x)**2)/(b*(c*sec(a + b*x))**(sympy.S(1)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_33():
    f = (c*sec(a + b*x))**(sympy.S(1)/3)
    F = -3*c*sin(a + b*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(a + b*x)**2)/(2*b*(c*sec(a + b*x))**(sympy.S(2)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_34():
    f = (c*sec(a + b*x))**(sympy.S(-1)/3)
    F = -3*c*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(a + b*x)**2)/(4*b*(c*sec(a + b*x))**(sympy.S(4)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_35():
    f = (c*sec(a + b*x))**(sympy.S(-2)/3)
    F = -3*c*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(a + b*x)**2)/(5*b*(c*sec(a + b*x))**(sympy.S(5)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_36():
    f = (c*sec(a + b*x))**(sympy.S(-4)/3)
    F = -3*c*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(a + b*x)**2)/(7*b*(c*sec(a + b*x))**(sympy.S(7)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_37():
    f = sec(a + b*x)**n
    F = -sin(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)*sec(a + b*x)**(n - 1)/(b*(1 - n)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_38():
    f = (c*sec(a + b*x))**n
    F = -c*(c*sec(a + b*x))**(n - 1)*sin(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(a + b*x)**2)/(b*(1 - n)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_39():
    f = (sec(x)**2)**(sympy.S(7)/2)
    F = (sec(x)**2)**(sympy.S(5)/2)*tan(x)/6 + 5*(sec(x)**2)**(sympy.S(3)/2)*tan(x)/24 + 5*sqrt(sec(x)**2)*tan(x)/16 + 5*asinh(tan(x))/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_40():
    f = (sec(x)**2)**(sympy.S(5)/2)
    F = (sec(x)**2)**(sympy.S(3)/2)*tan(x)/4 + 3*sqrt(sec(x)**2)*tan(x)/8 + 3*asinh(tan(x))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_41():
    f = (sec(x)**2)**(sympy.S(3)/2)
    F = sqrt(sec(x)**2)*tan(x)/2 + asinh(tan(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_42():
    f = sqrt(sec(x)**2)
    F = asinh(tan(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_43():
    f = 1/sqrt(sec(x)**2)
    F = tan(x)/sqrt(sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_44():
    f = (sec(x)**2)**(sympy.S(-3)/2)
    F = 2*tan(x)/(3*sqrt(sec(x)**2)) + tan(x)/(3*(sec(x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_45():
    f = (sec(x)**2)**(sympy.S(-5)/2)
    F = 8*tan(x)/(15*sqrt(sec(x)**2)) + 4*tan(x)/(15*(sec(x)**2)**(sympy.S(3)/2)) + tan(x)/(5*(sec(x)**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_46():
    f = (sec(x)**2)**(sympy.S(-7)/2)
    F = 16*tan(x)/(35*sqrt(sec(x)**2)) + 8*tan(x)/(35*(sec(x)**2)**(sympy.S(3)/2)) + 6*tan(x)/(35*(sec(x)**2)**(sympy.S(5)/2)) + tan(x)/(7*(sec(x)**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_47():
    f = (a*sec(x)**2)**(sympy.S(7)/2)
    F = 5*a**(sympy.S(7)/2)*atanh(sqrt(a)*tan(x)/sqrt(a*sec(x)**2))/16 + 5*a**3*sqrt(a*sec(x)**2)*tan(x)/16 + 5*a**2*(a*sec(x)**2)**(sympy.S(3)/2)*tan(x)/24 + a*(a*sec(x)**2)**(sympy.S(5)/2)*tan(x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_48():
    f = (a*sec(x)**2)**(sympy.S(5)/2)
    F = 3*a**(sympy.S(5)/2)*atanh(sqrt(a)*tan(x)/sqrt(a*sec(x)**2))/8 + 3*a**2*sqrt(a*sec(x)**2)*tan(x)/8 + a*(a*sec(x)**2)**(sympy.S(3)/2)*tan(x)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_49():
    f = (a*sec(x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tan(x)/sqrt(a*sec(x)**2))/2 + a*sqrt(a*sec(x)**2)*tan(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_50():
    f = sqrt(a*sec(x)**2)
    F = sqrt(a)*atanh(sqrt(a)*tan(x)/sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_51():
    f = 1/sqrt(a*sec(x)**2)
    F = tan(x)/sqrt(a*sec(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_52():
    f = (a*sec(x)**2)**(sympy.S(-3)/2)
    F = tan(x)/(3*(a*sec(x)**2)**(sympy.S(3)/2)) + 2*tan(x)/(3*a*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_53():
    f = (a*sec(x)**2)**(sympy.S(-5)/2)
    F = tan(x)/(5*(a*sec(x)**2)**(sympy.S(5)/2)) + 4*tan(x)/(15*a*(a*sec(x)**2)**(sympy.S(3)/2)) + 8*tan(x)/(15*a**2*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_54():
    f = (a*sec(x)**2)**(sympy.S(-7)/2)
    F = tan(x)/(7*(a*sec(x)**2)**(sympy.S(7)/2)) + 6*tan(x)/(35*a*(a*sec(x)**2)**(sympy.S(5)/2)) + 8*tan(x)/(35*a**2*(a*sec(x)**2)**(sympy.S(3)/2)) + 16*tan(x)/(35*a**3*sqrt(a*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_55():
    f = (a*sec(x)**3)**(sympy.S(5)/2)
    F = 154*a**2*sqrt(a*sec(x)**3)*sin(x)*cos(x)/195 - 154*a**2*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2)*elliptic_e(x/2, 2)/195 + 2*a**2*sqrt(a*sec(x)**3)*tan(x)*sec(x)**4/13 + 22*a**2*sqrt(a*sec(x)**3)*tan(x)*sec(x)**2/117 + 154*a**2*sqrt(a*sec(x)**3)*tan(x)/585
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_56():
    f = (a*sec(x)**3)**(sympy.S(3)/2)
    F = 10*a*sqrt(a*sec(x)**3)*sin(x)/21 + 10*a*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2)*elliptic_f(x/2, 2)/21 + 2*a*sqrt(a*sec(x)**3)*tan(x)*sec(x)/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_57():
    f = sqrt(a*sec(x)**3)
    F = 2*sqrt(a*sec(x)**3)*sin(x)*cos(x) - 2*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2)*elliptic_e(x/2, 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_58():
    f = 1/sqrt(a*sec(x)**3)
    F = 2*tan(x)/(3*sqrt(a*sec(x)**3)) + 2*elliptic_f(x/2, 2)/(3*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_59():
    f = (a*sec(x)**3)**(sympy.S(-3)/2)
    F = 2*sin(x)*cos(x)**2/(9*a*sqrt(a*sec(x)**3)) + 14*sin(x)/(45*a*sqrt(a*sec(x)**3)) + 14*elliptic_e(x/2, 2)/(15*a*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_60():
    f = (a*sec(x)**3)**(sympy.S(-5)/2)
    F = 2*sin(x)*cos(x)**5/(15*a**2*sqrt(a*sec(x)**3)) + 26*sin(x)*cos(x)**3/(165*a**2*sqrt(a*sec(x)**3)) + 78*sin(x)*cos(x)/(385*a**2*sqrt(a*sec(x)**3)) + 26*tan(x)/(77*a**2*sqrt(a*sec(x)**3)) + 26*elliptic_f(x/2, 2)/(77*a**2*sqrt(a*sec(x)**3)*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_61():
    f = (a*sec(x)**4)**(sympy.S(7)/2)
    F = a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**11/13 + 6*a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**9/11 + 5*a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**7/3 + 20*a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**5/7 + 3*a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**3 + 2*a**3*sqrt(a*sec(x)**4)*sin(x)**2*tan(x) + a**3*sqrt(a*sec(x)**4)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_62():
    f = (a*sec(x)**4)**(sympy.S(5)/2)
    F = a**2*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**7/9 + 4*a**2*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**5/7 + 6*a**2*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**3/5 + 4*a**2*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)/3 + a**2*sqrt(a*sec(x)**4)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_63():
    f = (a*sec(x)**4)**(sympy.S(3)/2)
    F = a*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)**3/5 + 2*a*sqrt(a*sec(x)**4)*sin(x)**2*tan(x)/3 + a*sqrt(a*sec(x)**4)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_64():
    f = sqrt(a*sec(x)**4)
    F = sqrt(a*sec(x)**4)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_65():
    f = 1/sqrt(a*sec(x)**4)
    F = x*sec(x)**2/(2*sqrt(a*sec(x)**4)) + tan(x)/(2*sqrt(a*sec(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_66():
    f = (a*sec(x)**4)**(sympy.S(-3)/2)
    F = 5*x*sec(x)**2/(16*a*sqrt(a*sec(x)**4)) + sin(x)*cos(x)**3/(6*a*sqrt(a*sec(x)**4)) + 5*sin(x)*cos(x)/(24*a*sqrt(a*sec(x)**4)) + 5*tan(x)/(16*a*sqrt(a*sec(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_67():
    f = (a*sec(x)**4)**(sympy.S(-5)/2)
    F = 63*x*sec(x)**2/(256*a**2*sqrt(a*sec(x)**4)) + sin(x)*cos(x)**7/(10*a**2*sqrt(a*sec(x)**4)) + 9*sin(x)*cos(x)**5/(80*a**2*sqrt(a*sec(x)**4)) + 21*sin(x)*cos(x)**3/(160*a**2*sqrt(a*sec(x)**4)) + 21*sin(x)*cos(x)/(128*a**2*sqrt(a*sec(x)**4)) + 63*tan(x)/(256*a**2*sqrt(a*sec(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_68():
    f = ((b*sec(c + d*x))**p)**n
    F = -((b*sec(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(-n*p + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_69():
    f = (a*(b*sec(c + d*x))**p)**n
    F = -(a*(b*sec(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(-n*p + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_70():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**4
    F = 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 10*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*b*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_71():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**3
    F = -6*b*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_72():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**2
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_73():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)
    F = -2*b*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_74():
    f = sqrt(b*sec(c + d*x))
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_75():
    f = sqrt(b*sec(c + d*x))*cos(c + d*x)
    F = 2*b*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_76():
    f = sqrt(b*sec(c + d*x))*cos(c + d*x)**2
    F = 2*b*sin(c + d*x)/(3*d*sqrt(b*sec(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_77():
    f = sqrt(b*sec(c + d*x))*cos(c + d*x)**3
    F = 2*b**2*sin(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*b*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_78():
    f = sqrt(b*sec(c + d*x))*cos(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*b*sin(c + d*x)/(21*d*sqrt(b*sec(c + d*x))) + 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_79():
    f = sqrt(b*sec(c + d*x))*cos(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*b**2*sin(c + d*x)/(45*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*b*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_80():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 10*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 10*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_81():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = -6*b**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_82():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_83():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*b**2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*b*sqrt(b*sec(c + d*x))*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_84():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_85():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 2*b**2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_86():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 2*b**2*sin(c + d*x)/(3*d*sqrt(b*sec(c + d*x))) + 2*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_87():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*b**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_88():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*b**2*sin(c + d*x)/(21*d*sqrt(b*sec(c + d*x))) + 10*b*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_89():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**6
    F = 2*b**5*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*b**3*sin(c + d*x)/(45*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*b**2*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_90():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 10*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 10*b*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_91():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = -6*b**3*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_92():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*b*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_93():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = -2*b**3*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_94():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 2*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_95():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 2*b**3*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_96():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(3*d*sqrt(b*sec(c + d*x))) + 2*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_97():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*b**3*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_98():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**6
    F = 2*b**5*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*b**3*sin(c + d*x)/(21*d*sqrt(b*sec(c + d*x))) + 10*b**2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_99():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**7
    F = 2*b**6*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*b**4*sin(c + d*x)/(45*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*b**3*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_100():
    f = (b*sec(c + d*x))**(sympy.S(7)/2)
    F = -6*b**4*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*b**3*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*b*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_101():
    f = sec(c + d*x)**5/sqrt(b*sec(c + d*x))
    F = 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d) + 10*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*b**2*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_102():
    f = sec(c + d*x)**4/sqrt(b*sec(c + d*x))
    F = -6*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*b*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_103():
    f = sec(c + d*x)**3/sqrt(b*sec(c + d*x))
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d) + 2*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_104():
    f = sec(c + d*x)**2/sqrt(b*sec(c + d*x))
    F = -2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_105():
    f = sec(c + d*x)/sqrt(b*sec(c + d*x))
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_106():
    f = 1/sqrt(b*sec(c + d*x))
    F = 2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_107():
    f = cos(c + d*x)/sqrt(b*sec(c + d*x))
    F = 2*sin(c + d*x)/(3*d*sqrt(b*sec(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_108():
    f = cos(c + d*x)**2/sqrt(b*sec(c + d*x))
    F = 2*b*sin(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_109():
    f = cos(c + d*x)**3/sqrt(b*sec(c + d*x))
    F = 2*b**2*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*sin(c + d*x)/(21*d*sqrt(b*sec(c + d*x))) + 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_110():
    f = cos(c + d*x)**4/sqrt(b*sec(c + d*x))
    F = 2*b**3*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*b*sin(c + d*x)/(45*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_111():
    f = sec(c + d*x)**6/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d) + 10*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*b**3*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_112():
    f = sec(c + d*x)**5/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = -6*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*b**2*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_113():
    f = sec(c + d*x)**4/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d) + 2*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_114():
    f = sec(c + d*x)**3/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_115():
    f = sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_116():
    f = sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_117():
    f = (b*sec(c + d*x))**(sympy.S(-3)/2)
    F = 2*sin(c + d*x)/(3*b*d*sqrt(b*sec(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_118():
    f = cos(c + d*x)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*sin(c + d*x)/(5*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_119():
    f = cos(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*sin(c + d*x)/(21*b*d*sqrt(b*sec(c + d*x))) + 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_120():
    f = cos(c + d*x)**3/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*sin(c + d*x)/(45*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_121():
    f = sec(c + d*x)**7/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**3*d) + 10*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*b**4*d) + 2*(b*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(7*b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_122():
    f = sec(c + d*x)**6/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = -6*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*sqrt(b*sec(c + d*x))*sin(c + d*x)/(5*b**3*d) + 2*(b*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_123():
    f = sec(c + d*x)**5/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d) + 2*(b*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_124():
    f = sec(c + d*x)**4/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_125():
    f = sec(c + d*x)**3/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_126():
    f = sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_127():
    f = sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*sin(c + d*x)/(3*b**2*d*sqrt(b*sec(c + d*x))) + 2*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_128():
    f = (b*sec(c + d*x))**(sympy.S(-5)/2)
    F = 2*sin(c + d*x)/(5*b*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 6*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_129():
    f = cos(c + d*x)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*sin(c + d*x)/(7*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*sin(c + d*x)/(21*b**2*d*sqrt(b*sec(c + d*x))) + 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_130():
    f = cos(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b*sin(c + d*x)/(9*d*(b*sec(c + d*x))**(sympy.S(7)/2)) + 14*sin(c + d*x)/(45*b*d*(b*sec(c + d*x))**(sympy.S(3)/2)) + 14*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_131():
    f = (b*sec(c + d*x))**(sympy.S(-7)/2)
    F = 2*sin(c + d*x)/(7*b*d*(b*sec(c + d*x))**(sympy.S(5)/2)) + 10*sin(c + d*x)/(21*b**3*d*sqrt(b*sec(c + d*x))) + 10*sqrt(b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_132():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)
    F = sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d) + 3*sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d) + 3*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_133():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)
    F = sqrt(b*sec(c + d*x))*sin(c + d*x)**3*sec(c + d*x)**(sympy.S(5)/2)/(3*d) + sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_134():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d) + sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_135():
    f = sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_136():
    f = sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x))
    F = sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_137():
    f = sqrt(b*sec(c + d*x))/sqrt(sec(c + d*x))
    F = x*sqrt(b*sec(c + d*x))/sqrt(sec(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_138():
    f = sqrt(b*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_139():
    f = sqrt(b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = x*sqrt(b*sec(c + d*x))/(2*sqrt(sec(c + d*x))) + sqrt(b*sec(c + d*x))*sin(c + d*x)/(2*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_140():
    f = sqrt(b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/2)
    F = -sqrt(b*sec(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(sec(c + d*x))) + sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_141():
    f = sqrt(b*sec(c + d*x))/sec(c + d*x)**(sympy.S(9)/2)
    F = 3*x*sqrt(b*sec(c + d*x))/(8*sqrt(sec(c + d*x))) + 3*sqrt(b*sec(c + d*x))*sin(c + d*x)/(8*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(b*sec(c + d*x))*sin(c + d*x)/(4*d*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_142():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = b*sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d) + 3*b*sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d) + 3*b*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_143():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = b*sqrt(b*sec(c + d*x))*sin(c + d*x)**3*sec(c + d*x)**(sympy.S(5)/2)/(3*d) + b*sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_144():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = b*sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d) + b*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_145():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = b*sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_146():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = b*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_147():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = b*x*sqrt(b*sec(c + d*x))/sqrt(sec(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_148():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_149():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = b*x*sqrt(b*sec(c + d*x))/(2*sqrt(sec(c + d*x))) + b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(2*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_150():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = -b*sqrt(b*sec(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(sec(c + d*x))) + b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_151():
    f = (b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(11)/2)
    F = 3*b*x*sqrt(b*sec(c + d*x))/(8*sqrt(sec(c + d*x))) + 3*b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(8*d*sec(c + d*x)**(sympy.S(3)/2)) + b*sqrt(b*sec(c + d*x))*sin(c + d*x)/(4*d*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_152():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)**5*sec(c + d*x)**(sympy.S(9)/2)/(5*d) + 2*b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)**3*sec(c + d*x)**(sympy.S(5)/2)/(3*d) + b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_153():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)**3*sec(c + d*x)**(sympy.S(5)/2)/(3*d) + b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_154():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d) + b**2*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_155():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_156():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = b**2*sqrt(b*sec(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_157():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = b**2*x*sqrt(b*sec(c + d*x))/sqrt(sec(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_158():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_159():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = b**2*x*sqrt(b*sec(c + d*x))/(2*sqrt(sec(c + d*x))) + b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(2*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_160():
    f = (b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(11)/2)
    F = -b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(sec(c + d*x))) + b**2*sqrt(b*sec(c + d*x))*sin(c + d*x)/(d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_161():
    f = sec(c + d*x)**(sympy.S(7)/2)/sqrt(b*sec(c + d*x))
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(b*sec(c + d*x))) + atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_162():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(b*sec(c + d*x))
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_163():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(b*sec(c + d*x))
    F = atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_164():
    f = sqrt(sec(c + d*x))/sqrt(b*sec(c + d*x))
    F = x*sqrt(sec(c + d*x))/sqrt(b*sec(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_165():
    f = 1/(sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_166():
    f = 1/(sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = x*sqrt(sec(c + d*x))/(2*sqrt(b*sec(c + d*x))) + sin(c + d*x)/(2*d*sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_167():
    f = 1/(sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)**3*sqrt(sec(c + d*x))/(3*d*sqrt(b*sec(c + d*x))) + sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_168():
    f = sec(c + d*x)**(sympy.S(9)/2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*b*d*sqrt(b*sec(c + d*x))) + atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(2*b*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_169():
    f = sec(c + d*x)**(sympy.S(7)/2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(b*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_170():
    f = sec(c + d*x)**(sympy.S(5)/2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(b*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_171():
    f = sec(c + d*x)**(sympy.S(3)/2)/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = x*sqrt(sec(c + d*x))/(b*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_172():
    f = sqrt(sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(3)/2)
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(b*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_173():
    f = 1/((b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = x*sqrt(sec(c + d*x))/(2*b*sqrt(b*sec(c + d*x))) + sin(c + d*x)/(2*b*d*sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_174():
    f = 1/((b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)**3*sqrt(sec(c + d*x))/(3*b*d*sqrt(b*sec(c + d*x))) + sin(c + d*x)*sqrt(sec(c + d*x))/(b*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_175():
    f = 1/((b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 3*x*sqrt(sec(c + d*x))/(8*b*sqrt(b*sec(c + d*x))) + 3*sin(c + d*x)/(8*b*d*sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x))) + sin(c + d*x)/(4*b*d*sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_176():
    f = sec(c + d*x)**(sympy.S(11)/2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*b**2*d*sqrt(b*sec(c + d*x))) + atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(2*b**2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_177():
    f = sec(c + d*x)**(sympy.S(9)/2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(b**2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_178():
    f = sec(c + d*x)**(sympy.S(7)/2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = atanh(sin(c + d*x))*sqrt(sec(c + d*x))/(b**2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_179():
    f = sec(c + d*x)**(sympy.S(5)/2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = x*sqrt(sec(c + d*x))/(b**2*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_180():
    f = sec(c + d*x)**(sympy.S(3)/2)/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(b**2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_181():
    f = sqrt(sec(c + d*x))/(b*sec(c + d*x))**(sympy.S(5)/2)
    F = x*sqrt(sec(c + d*x))/(2*b**2*sqrt(b*sec(c + d*x))) + sin(c + d*x)/(2*b**2*d*sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_182():
    f = 1/((b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = -sin(c + d*x)**3*sqrt(sec(c + d*x))/(3*b**2*d*sqrt(b*sec(c + d*x))) + sin(c + d*x)*sqrt(sec(c + d*x))/(b**2*d*sqrt(b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_183():
    f = 1/((b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 3*x*sqrt(sec(c + d*x))/(8*b**2*sqrt(b*sec(c + d*x))) + 3*sin(c + d*x)/(8*b**2*d*sqrt(b*sec(c + d*x))*sqrt(sec(c + d*x))) + sin(c + d*x)/(4*b**2*d*sqrt(b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_184():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)**2
    F = 3*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_185():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)
    F = 3*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_186():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*b*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_187():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)
    F = -3*b**2*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_188():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)**2
    F = -3*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*d*(b*sec(c + d*x))**(sympy.S(8)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_189():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)**2
    F = 3*(b*sec(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_190():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)
    F = 3*(b*sec(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_191():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)
    F = 3*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_192():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*cos(c + d*x)
    F = -3*b**2*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*(b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_193():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*cos(c + d*x)**2
    F = -3*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*d*(b*sec(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_194():
    f = sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = 3*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_195():
    f = sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_196():
    f = (b*sec(c + d*x))**(sympy.S(-1)/3)
    F = -3*b*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_197():
    f = cos(c + d*x)/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*b**2*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_198():
    f = cos(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*d*(b*sec(c + d*x))**(sympy.S(10)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_199():
    f = sec(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_200():
    f = sec(c + d*x)/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*(b*sec(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_201():
    f = (b*sec(c + d*x))**(sympy.S(-4)/3)
    F = -3*b*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*d*(b*sec(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_202():
    f = cos(c + d*x)/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*b**2*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*d*(b*sec(c + d*x))**(sympy.S(10)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_203():
    f = cos(c + d*x)**2/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*b**3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*d*(b*sec(c + d*x))**(sympy.S(13)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_204():
    f = (b*sec(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)**m
    F = 3*b*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/6), (sympy.S(5)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**m/(d*(3*m + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_205():
    f = (b*sec(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**m
    F = -3*(b*sec(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/6 - m/2), (sympy.S(7)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_206():
    f = (b*sec(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)**m
    F = -3*(b*sec(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/3 - m/2), (sympy.S(4)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(2 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_207():
    f = sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3 - m/2), (sympy.S(5)/3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(b*sec(c + d*x))**(sympy.S(1)/3)*(4 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_208():
    f = sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(2)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6 - m/2), (sympy.S(11)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(b*sec(c + d*x))**(sympy.S(2)/3)*(5 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_209():
    f = sec(c + d*x)**m/(b*sec(c + d*x))**(sympy.S(4)/3)
    F = -3*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6 - m/2), (sympy.S(13)/6 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 2)/(b*d*(b*sec(c + d*x))**(sympy.S(1)/3)*(7 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_210():
    f = (b*sec(c + d*x))**n*sec(c + d*x)**m
    F = -(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -m/2 - n/2 + sympy.S.Half), (-m/2 - n/2 + sympy.S(3)/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(-m - n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_211():
    f = (b*sec(c + d*x))**n*sec(c + d*x)**2
    F = (b*sec(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_212():
    f = (b*sec(c + d*x))**n*sec(c + d*x)
    F = (b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_213():
    f = (b*sec(c + d*x))**n
    F = -b*(b*sec(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_214():
    f = (b*sec(c + d*x))**n*cos(c + d*x)
    F = -b**2*(b*sec(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_215():
    f = (b*sec(c + d*x))**n*cos(c + d*x)**2
    F = -b**3*(b*sec(c + d*x))**(n - 3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), cos(c + d*x)**2)/(d*(3 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_216():
    f = (b*sec(c + d*x))**n*cos(c + d*x)**3
    F = -b**4*(b*sec(c + d*x))**(n - 4)*sin(c + d*x)*hyper((sympy.S.Half, 2 - n/2), (3 - n/2,), cos(c + d*x)**2)/(d*(4 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_217():
    f = (b*sec(c + d*x))**n*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-3)/4), (sympy.S(1)/4 - n/2,), cos(c + d*x)**2)*sec(c + d*x)**(sympy.S(3)/2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_218():
    f = (b*sec(c + d*x))**n*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/4), (sympy.S(3)/4 - n/2,), cos(c + d*x)**2)*sqrt(sec(c + d*x))/(d*(2*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_219():
    f = (b*sec(c + d*x))**n*sqrt(sec(c + d*x))
    F = -2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(1)/4 - n/2), (sympy.S(5)/4 - n/2,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_220():
    f = (b*sec(c + d*x))**n/sqrt(sec(c + d*x))
    F = -2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/4 - n/2), (sympy.S(7)/4 - n/2,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_221():
    f = (b*sec(c + d*x))**n/sec(c + d*x)**(sympy.S(3)/2)
    F = -2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/4 - n/2), (sympy.S(9)/4 - n/2,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_222():
    f = (b*sec(c + d*x))**n/sec(c + d*x)**(sympy.S(5)/2)
    F = -2*(b*sec(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/4 - n/2), (sympy.S(11)/4 - n/2,), cos(c + d*x)**2)/(d*(7 - 2*n)*sqrt(sin(c + d*x)**2)*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_223():
    f = (d*sec(a + b*x))**(sympy.S(7)/2)*sin(a + b*x)
    F = 2*d*(d*sec(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_224():
    f = (d*sec(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)
    F = 2*d*(d*sec(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_225():
    f = (d*sec(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)
    F = 2*d*sqrt(d*sec(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_226():
    f = sqrt(d*sec(a + b*x))*sin(a + b*x)
    F = -2*d/(b*sqrt(d*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_227():
    f = sin(a + b*x)/sqrt(d*sec(a + b*x))
    F = -2*d/(3*b*(d*sec(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_228():
    f = (d*sec(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**3
    F = 2*d**3/(b*sqrt(d*sec(a + b*x))) + 2*d*(d*sec(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_229():
    f = (d*sec(a + b*x))**(sympy.S(9)/2)*sin(a + b*x)**3
    F = -2*d**3*(d*sec(a + b*x))**(sympy.S(3)/2)/(3*b) + 2*d*(d*sec(a + b*x))**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_230():
    f = sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(9)/2)
    F = -4*c*d**3*(d*csc(a + b*x))**(sympy.S(3)/2)/(7*b*sqrt(c*sec(a + b*x))) - 2*c*d*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b*sqrt(c*sec(a + b*x))) + 4*d**4*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_231():
    f = sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(7)/2)
    F = -8*c*d**3*sqrt(d*csc(a + b*x))/(5*b*sqrt(c*sec(a + b*x))) - 2*c*d*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_232():
    f = sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(5)/2)
    F = -2*c*d*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b*sqrt(c*sec(a + b*x))) + 2*d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_233():
    f = sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)
    F = -2*c*d*sqrt(d*csc(a + b*x))/(b*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_234():
    f = sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))
    F = sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_235():
    f = sqrt(c*sec(a + b*x))/sqrt(d*csc(a + b*x))
    F = sqrt(2)*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_236():
    f = sqrt(c*sec(a + b*x))/(d*csc(a + b*x))**(sympy.S(3)/2)
    F = -c/(b*d*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))) + sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(2*b*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_237():
    f = sqrt(c*sec(a + b*x))/(d*csc(a + b*x))**(sympy.S(5)/2)
    F = -c/(2*b*d*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(16*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - 3*sqrt(2)*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(16*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(8*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(8*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_238():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(9)/2)
    F = 64*c*d**5*sqrt(c*sec(a + b*x))/(21*b*sqrt(d*csc(a + b*x))) - 16*c*d**3*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)/(21*b) - 2*c*d*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_239():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(7)/2)
    F = -24*c**2*d**4*elliptic_e(a + b*x - pi/4, 2)/(5*b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 24*c*d**5*sqrt(c*sec(a + b*x))/(5*b*(d*csc(a + b*x))**(sympy.S(3)/2)) - 12*c*d**3*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))/(5*b) - 2*c*d*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_240():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(5)/2)
    F = 8*c*d**3*sqrt(c*sec(a + b*x))/(3*b*sqrt(d*csc(a + b*x))) - 2*c*d*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_241():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)
    F = -4*c**2*d**2*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 4*c*d**3*sqrt(c*sec(a + b*x))/(b*(d*csc(a + b*x))**(sympy.S(3)/2)) - 2*c*d*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_242():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))
    F = 2*c*d*sqrt(c*sec(a + b*x))/(b*sqrt(d*csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_243():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)/sqrt(d*csc(a + b*x))
    F = -2*c**2*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 2*c*d*sqrt(c*sec(a + b*x))/(b*(d*csc(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_244():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)/(d*csc(a + b*x))**(sympy.S(3)/2)
    F = sqrt(2)*c**2*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*d**2*sqrt(c*sec(a + b*x))) - sqrt(2)*c**2*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*d**2*sqrt(c*sec(a + b*x))) - sqrt(2)*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*d**2*sqrt(c*sec(a + b*x))) - sqrt(2)*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*d**2*sqrt(c*sec(a + b*x))) + 2*c*sqrt(c*sec(a + b*x))/(b*d*sqrt(d*csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_245():
    f = (c*sec(a + b*x))**(sympy.S(3)/2)/(d*csc(a + b*x))**(sympy.S(5)/2)
    F = -3*c**2*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 2*c*sqrt(c*sec(a + b*x))/(b*d*(d*csc(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_246():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(9)/2)
    F = 40*c**2*d**4*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(21*b) + 40*c*d**5*(c*sec(a + b*x))**(sympy.S(3)/2)/(21*b*sqrt(d*csc(a + b*x))) - 20*c*d**3*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)/(21*b) - 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_247():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(7)/2)
    F = -64*c**3*d**3*sqrt(d*csc(a + b*x))/(15*b*sqrt(c*sec(a + b*x))) + 16*c*d**3*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))/(15*b) - 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_248():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(5)/2)
    F = 4*c**2*d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b) + 4*c*d**3*(c*sec(a + b*x))**(sympy.S(3)/2)/(3*b*sqrt(d*csc(a + b*x))) - 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_249():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)
    F = -8*c**3*d*sqrt(d*csc(a + b*x))/(3*b*sqrt(c*sec(a + b*x))) + 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_250():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)*sqrt(d*csc(a + b*x))
    F = 2*c**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b) + 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)/(3*b*sqrt(d*csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_251():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)/sqrt(d*csc(a + b*x))
    F = 2*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)/(3*b*(d*csc(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_252():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)/(d*csc(a + b*x))**(sympy.S(3)/2)
    F = -c**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b*d**2) + 2*c*(c*sec(a + b*x))**(sympy.S(3)/2)/(3*b*d*sqrt(d*csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_253():
    f = (c*sec(a + b*x))**(sympy.S(5)/2)/(d*csc(a + b*x))**(sympy.S(5)/2)
    F = -sqrt(2)*c**2*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*c**2*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*c**2*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*c**2*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + 2*c*(c*sec(a + b*x))**(sympy.S(3)/2)/(3*b*d*(d*csc(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_254():
    f = (d*csc(a + b*x))**(sympy.S(9)/2)/sqrt(c*sec(a + b*x))
    F = -8*c*d**3*(d*csc(a + b*x))**(sympy.S(3)/2)/(21*b*(c*sec(a + b*x))**(sympy.S(3)/2)) - 2*c*d*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b*(c*sec(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_255():
    f = (d*csc(a + b*x))**(sympy.S(7)/2)/sqrt(c*sec(a + b*x))
    F = -4*c*d**3*sqrt(d*csc(a + b*x))/(5*b*(c*sec(a + b*x))**(sympy.S(3)/2)) - 2*c*d*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b*(c*sec(a + b*x))**(sympy.S(3)/2)) - 4*d**4*elliptic_e(a + b*x - pi/4, 2)/(5*b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_256():
    f = (d*csc(a + b*x))**(sympy.S(5)/2)/sqrt(c*sec(a + b*x))
    F = -2*c*d*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b*(c*sec(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_257():
    f = (d*csc(a + b*x))**(sympy.S(3)/2)/sqrt(c*sec(a + b*x))
    F = -2*c*d*sqrt(d*csc(a + b*x))/(b*(c*sec(a + b*x))**(sympy.S(3)/2)) - 2*d**2*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_258():
    f = sqrt(d*csc(a + b*x))/sqrt(c*sec(a + b*x))
    F = -sqrt(2)*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_259():
    f = 1/(sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x)))
    F = elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_260():
    f = 1/(sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2))
    F = -c/(2*b*d*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))) - sqrt(2)*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(16*b*d**2*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(16*b*d**2*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(8*b*d**2*sqrt(c*sec(a + b*x))) + sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(8*b*d**2*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_261():
    f = 1/(sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(5)/2))
    F = -c/(3*b*d*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)) + elliptic_e(a + b*x - pi/4, 2)/(2*b*d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_262():
    f = (d*csc(a + b*x))**(sympy.S(11)/2)/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = 8*d**5*sqrt(d*csc(a + b*x))/(45*b*c*sqrt(c*sec(a + b*x))) + 2*d**3*(d*csc(a + b*x))**(sympy.S(5)/2)/(45*b*c*sqrt(c*sec(a + b*x))) - 2*d*(d*csc(a + b*x))**(sympy.S(9)/2)/(9*b*c*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_263():
    f = (d*csc(a + b*x))**(sympy.S(9)/2)/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = 2*d**3*(d*csc(a + b*x))**(sympy.S(3)/2)/(21*b*c*sqrt(c*sec(a + b*x))) - 2*d*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b*c*sqrt(c*sec(a + b*x))) - 2*d**4*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(21*b*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_264():
    f = (d*csc(a + b*x))**(sympy.S(7)/2)/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = -2*c*d*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b*(c*sec(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_265():
    f = (d*csc(a + b*x))**(sympy.S(5)/2)/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = -2*d*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b*c*sqrt(c*sec(a + b*x))) - d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(3*b*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_266():
    f = (d*csc(a + b*x))**(sympy.S(3)/2)/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = -2*d*sqrt(d*csc(a + b*x))/(b*c*sqrt(c*sec(a + b*x))) - sqrt(2)*d**2*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*d**2*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(4*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*d**2*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*d**2*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_267():
    f = sqrt(d*csc(a + b*x))/(c*sec(a + b*x))**(sympy.S(3)/2)
    F = d/(b*c*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))) + sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(2*b*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_268():
    f = 1/((c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x)))
    F = d/(2*b*c*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)) + sqrt(2)*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(16*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - sqrt(2)*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(16*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(8*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(8*b*c**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_269():
    f = 1/((c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2))
    F = -c/(3*b*d*(c*sec(a + b*x))**(sympy.S(5)/2)*sqrt(d*csc(a + b*x))) + 1/(6*b*c*d*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))) + sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)/(12*b*c**2*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_270():
    f = 1/((c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(5)/2))
    F = -c/(4*b*d*(c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)) + 3/(16*b*c*d*sqrt(c*sec(a + b*x))*(d*csc(a + b*x))**(sympy.S(3)/2)) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(128*b*c**2*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) - 3*sqrt(2)*sqrt(c*sec(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)/(128*b*c**2*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(64*b*c**2*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))) + 3*sqrt(2)*sqrt(c*sec(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(64*b*c**2*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_271():
    f = (d*csc(a + b*x))**(sympy.S(9)/2)/(c*sec(a + b*x))**(sympy.S(5)/2)
    F = -2*c*d*(d*csc(a + b*x))**(sympy.S(7)/2)/(7*b*(c*sec(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_272():
    f = (d*csc(a + b*x))**(sympy.S(7)/2)/(c*sec(a + b*x))**(sympy.S(5)/2)
    F = 6*d**3*sqrt(d*csc(a + b*x))/(5*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)) - 2*d*(d*csc(a + b*x))**(sympy.S(5)/2)/(5*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)) + 6*d**4*elliptic_e(a + b*x - pi/4, 2)/(5*b*c**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_273():
    f = (d*csc(a + b*x))**(sympy.S(5)/2)/(c*sec(a + b*x))**(sympy.S(5)/2)
    F = -2*d*(d*csc(a + b*x))**(sympy.S(3)/2)/(3*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)) + sqrt(2)*d**2*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*c**2*sqrt(c*sec(a + b*x))) - sqrt(2)*d**2*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(4*b*c**2*sqrt(c*sec(a + b*x))) - sqrt(2)*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(2*b*c**2*sqrt(c*sec(a + b*x))) - sqrt(2)*d**2*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(2*b*c**2*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_274():
    f = (d*csc(a + b*x))**(sympy.S(3)/2)/(c*sec(a + b*x))**(sympy.S(5)/2)
    F = -2*d*sqrt(d*csc(a + b*x))/(b*c*(c*sec(a + b*x))**(sympy.S(3)/2)) - 3*d**2*elliptic_e(a + b*x - pi/4, 2)/(b*c**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_275():
    f = sqrt(d*csc(a + b*x))/(c*sec(a + b*x))**(sympy.S(5)/2)
    F = d/(2*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))) - 3*sqrt(2)*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(16*b*c**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(16*b*c**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(8*b*c**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(8*b*c**2*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_276():
    f = 1/((c*sec(a + b*x))**(sympy.S(5)/2)*sqrt(d*csc(a + b*x)))
    F = d/(3*b*c*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)) + elliptic_e(a + b*x - pi/4, 2)/(2*b*c**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_277():
    f = 1/((c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(3)/2))
    F = -c/(4*b*d*(c*sec(a + b*x))**(sympy.S(7)/2)*sqrt(d*csc(a + b*x))) + 1/(16*b*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))) - 3*sqrt(2)*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(128*b*c**2*d**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(128*b*c**2*d**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(64*b*c**2*d**2*sqrt(c*sec(a + b*x))) + 3*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(64*b*c**2*d**2*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_278():
    f = 1/((c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(5)/2))
    F = -c/(5*b*d*(c*sec(a + b*x))**(sympy.S(7)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)) + 1/(10*b*c*d*(c*sec(a + b*x))**(sympy.S(3)/2)*(d*csc(a + b*x))**(sympy.S(3)/2)) + 3*elliptic_e(a + b*x - pi/4, 2)/(20*b*c**2*d**2*sqrt(c*sec(a + b*x))*sqrt(d*csc(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_279():
    f = 1/((c*sec(a + b*x))**(sympy.S(5)/2)*(d*csc(a + b*x))**(sympy.S(7)/2))
    F = -c/(6*b*d*(c*sec(a + b*x))**(sympy.S(7)/2)*(d*csc(a + b*x))**(sympy.S(5)/2)) - 5*c/(48*b*d**3*(c*sec(a + b*x))**(sympy.S(7)/2)*sqrt(d*csc(a + b*x))) + 5/(192*b*c*d**3*(c*sec(a + b*x))**(sympy.S(3)/2)*sqrt(d*csc(a + b*x))) - 5*sqrt(2)*sqrt(d*csc(a + b*x))*log(-sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(512*b*c**2*d**4*sqrt(c*sec(a + b*x))) + 5*sqrt(2)*sqrt(d*csc(a + b*x))*log(sqrt(2)*sqrt(tan(a + b*x)) + tan(a + b*x) + 1)*sqrt(tan(a + b*x))/(512*b*c**2*d**4*sqrt(c*sec(a + b*x))) + 5*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) - 1)/(256*b*c**2*d**4*sqrt(c*sec(a + b*x))) + 5*sqrt(2)*sqrt(d*csc(a + b*x))*sqrt(tan(a + b*x))*atan(sqrt(2)*sqrt(tan(a + b*x)) + 1)/(256*b*c**2*d**4*sqrt(c*sec(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_280():
    f = csc(e + f*x)**n*sec(e + f*x)**m
    F = (cos(e + f*x)**2)**(m/2 + sympy.S.Half)*csc(e + f*x)**(n - 1)*hyper((sympy.S.Half - n/2, m/2 + sympy.S.Half), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)*sec(e + f*x)**(m + 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_281():
    f = (a*sec(e + f*x))**m*csc(e + f*x)**n
    F = (a*sec(e + f*x))**(m + 1)*(cos(e + f*x)**2)**(m/2 + sympy.S.Half)*csc(e + f*x)**(n - 1)*hyper((sympy.S.Half - n/2, m/2 + sympy.S.Half), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(a*f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_282():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**m
    F = b*(b*csc(e + f*x))**(n - 1)*(cos(e + f*x)**2)**(m/2 + sympy.S.Half)*hyper((sympy.S.Half - n/2, m/2 + sympy.S.Half), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)*sec(e + f*x)**(m + 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_283():
    f = (a*sec(e + f*x))**m*(b*csc(e + f*x))**n
    F = b*(a*sec(e + f*x))**(m + 1)*(b*csc(e + f*x))**(n - 1)*(cos(e + f*x)**2)**(m/2 + sympy.S.Half)*hyper((sympy.S.Half - n/2, m/2 + sympy.S.Half), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(a*f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_284():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**5
    F = (b*csc(e + f*x))**(n + 5)*hyper((3, n/2 + sympy.S(5)/2), (n/2 + sympy.S(7)/2,), csc(e + f*x)**2)/(b**5*f*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_285():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**3
    F = -(b*csc(e + f*x))**(n + 3)*hyper((2, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), csc(e + f*x)**2)/(b**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_286():
    f = (b*csc(e + f*x))**n*sec(e + f*x)
    F = (b*csc(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), csc(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_287():
    f = (b*csc(e + f*x))**n*cos(e + f*x)
    F = b*(b*csc(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_288():
    f = (b*csc(e + f*x))**n*cos(e + f*x)**3
    F = -b**3*(b*csc(e + f*x))**(n - 3)/(f*(3 - n)) + b*(b*csc(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_289():
    f = (b*csc(e + f*x))**n*cos(e + f*x)**5
    F = b**5*(b*csc(e + f*x))**(n - 5)/(f*(5 - n)) - 2*b**3*(b*csc(e + f*x))**(n - 3)/(f*(3 - n)) + b*(b*csc(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_290():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**6
    F = b*(b*csc(e + f*x))**(n - 1)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(7)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)*sec(e + f*x)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_291():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**4
    F = b*(b*csc(e + f*x))**(n - 1)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(5)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)*sec(e + f*x)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_292():
    f = (b*csc(e + f*x))**n*sec(e + f*x)**2
    F = b*(b*csc(e + f*x))**(n - 1)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(3)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)*sec(e + f*x)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_293():
    f = (b*csc(e + f*x))**n
    F = b*(b*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_294():
    f = (b*csc(e + f*x))**n*cos(e + f*x)**2
    F = b*(b*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S(-1)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_295():
    f = (b*csc(e + f*x))**n*cos(e + f*x)**4
    F = b*(b*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S(-3)/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_296():
    f = (b*csc(e + f*x))**n*(c*sec(e + f*x))**(sympy.S(3)/2)
    F = b*(b*csc(e + f*x))**(n - 1)*(c*sec(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(5)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(c*f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_297():
    f = (b*csc(e + f*x))**n*sqrt(c*sec(e + f*x))
    F = b*(b*csc(e + f*x))**(n - 1)*(c*sec(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(c*f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_298():
    f = (b*csc(e + f*x))**n/sqrt(c*sec(e + f*x))
    F = b*(b*csc(e + f*x))**(n - 1)*sqrt(c*sec(e + f*x))*(cos(e + f*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(c*f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_0_a_sec_pow_m_b_trg_pow_n_299():
    f = (b*csc(e + f*x))**n/(c*sec(e + f*x))**(sympy.S(3)/2)
    F = b*(b*csc(e + f*x))**(n - 1)*hyper((sympy.S(-1)/4, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(c*f*sqrt(c*sec(e + f*x))*(1 - n)*(cos(e + f*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F

