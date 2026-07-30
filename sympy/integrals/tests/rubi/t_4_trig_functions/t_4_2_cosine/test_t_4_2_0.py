"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.0 (a cos)^m (b trg)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_1():
    f = cos(a + b*x)
    F = sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_2():
    f = cos(a + b*x)**2
    F = x/2 + sin(a + b*x)*cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_3():
    f = cos(a + b*x)**3
    F = -sin(a + b*x)**3/(3*b) + sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_4():
    f = cos(a + b*x)**4
    F = 3*x/8 + sin(a + b*x)*cos(a + b*x)**3/(4*b) + 3*sin(a + b*x)*cos(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_5():
    f = cos(a + b*x)**5
    F = sin(a + b*x)**5/(5*b) - 2*sin(a + b*x)**3/(3*b) + sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_6():
    f = cos(a + b*x)**6
    F = 5*x/16 + sin(a + b*x)*cos(a + b*x)**5/(6*b) + 5*sin(a + b*x)*cos(a + b*x)**3/(24*b) + 5*sin(a + b*x)*cos(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_7():
    f = cos(a + b*x)**7
    F = -sin(a + b*x)**7/(7*b) + 3*sin(a + b*x)**5/(5*b) - sin(a + b*x)**3/b + sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_8():
    f = cos(a + b*x)**8
    F = 35*x/128 + sin(a + b*x)*cos(a + b*x)**7/(8*b) + 7*sin(a + b*x)*cos(a + b*x)**5/(48*b) + 35*sin(a + b*x)*cos(a + b*x)**3/(192*b) + 35*sin(a + b*x)*cos(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_9():
    f = cos(a + b*x)**(sympy.S(7)/2)
    F = 2*sin(a + b*x)*cos(a + b*x)**(sympy.S(5)/2)/(7*b) + 10*sin(a + b*x)*sqrt(cos(a + b*x))/(21*b) + 10*elliptic_f(a/2 + b*x/2, 2)/(21*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_10():
    f = cos(a + b*x)**(sympy.S(5)/2)
    F = 2*sin(a + b*x)*cos(a + b*x)**(sympy.S(3)/2)/(5*b) + 6*elliptic_e(a/2 + b*x/2, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_11():
    f = cos(a + b*x)**(sympy.S(3)/2)
    F = 2*sin(a + b*x)*sqrt(cos(a + b*x))/(3*b) + 2*elliptic_f(a/2 + b*x/2, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_12():
    f = sqrt(cos(a + b*x))
    F = 2*elliptic_e(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_13():
    f = 1/sqrt(cos(a + b*x))
    F = 2*elliptic_f(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_14():
    f = cos(a + b*x)**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)/(b*sqrt(cos(a + b*x))) - 2*elliptic_e(a/2 + b*x/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_15():
    f = cos(a + b*x)**(sympy.S(-5)/2)
    F = 2*sin(a + b*x)/(3*b*cos(a + b*x)**(sympy.S(3)/2)) + 2*elliptic_f(a/2 + b*x/2, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_16():
    f = cos(a + b*x)**(sympy.S(-7)/2)
    F = 6*sin(a + b*x)/(5*b*sqrt(cos(a + b*x))) + 2*sin(a + b*x)/(5*b*cos(a + b*x)**(sympy.S(5)/2)) - 6*elliptic_e(a/2 + b*x/2, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_17():
    f = (c*cos(a + b*x))**(sympy.S(7)/2)
    F = 10*c**4*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(21*b*sqrt(c*cos(a + b*x))) + 10*c**3*sqrt(c*cos(a + b*x))*sin(a + b*x)/(21*b) + 2*c*(c*cos(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_18():
    f = (c*cos(a + b*x))**(sympy.S(5)/2)
    F = 6*c**2*sqrt(c*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*sqrt(cos(a + b*x))) + 2*c*(c*cos(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_19():
    f = (c*cos(a + b*x))**(sympy.S(3)/2)
    F = 2*c**2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*sqrt(c*cos(a + b*x))) + 2*c*sqrt(c*cos(a + b*x))*sin(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_20():
    f = sqrt(c*cos(a + b*x))
    F = 2*sqrt(c*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_21():
    f = 1/sqrt(c*cos(a + b*x))
    F = 2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(b*sqrt(c*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_22():
    f = (c*cos(a + b*x))**(sympy.S(-3)/2)
    F = 2*sin(a + b*x)/(b*c*sqrt(c*cos(a + b*x))) - 2*sqrt(c*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(b*c**2*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_23():
    f = (c*cos(a + b*x))**(sympy.S(-5)/2)
    F = 2*sin(a + b*x)/(3*b*c*(c*cos(a + b*x))**(sympy.S(3)/2)) + 2*sqrt(cos(a + b*x))*elliptic_f(a/2 + b*x/2, 2)/(3*b*c**2*sqrt(c*cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_24():
    f = (c*cos(a + b*x))**(sympy.S(-7)/2)
    F = 2*sin(a + b*x)/(5*b*c*(c*cos(a + b*x))**(sympy.S(5)/2)) + 6*sin(a + b*x)/(5*b*c**3*sqrt(c*cos(a + b*x))) - 6*sqrt(c*cos(a + b*x))*elliptic_e(a/2 + b*x/2, 2)/(5*b*c**4*sqrt(cos(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_25():
    f = cos(a + b*x)**(sympy.S(4)/3)
    F = -3*sin(a + b*x)*cos(a + b*x)**(sympy.S(7)/3)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(a + b*x)**2)/(7*b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_26():
    f = cos(a + b*x)**(sympy.S(2)/3)
    F = -3*sin(a + b*x)*cos(a + b*x)**(sympy.S(5)/3)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(a + b*x)**2)/(5*b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_27():
    f = cos(a + b*x)**(sympy.S(1)/3)
    F = -3*sin(a + b*x)*cos(a + b*x)**(sympy.S(4)/3)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(a + b*x)**2)/(4*b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_28():
    f = cos(a + b*x)**(sympy.S(-1)/3)
    F = -3*sin(a + b*x)*cos(a + b*x)**(sympy.S(2)/3)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(a + b*x)**2)/(2*b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_29():
    f = cos(a + b*x)**(sympy.S(-2)/3)
    F = -3*sin(a + b*x)*cos(a + b*x)**(sympy.S(1)/3)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(a + b*x)**2)/(b*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_30():
    f = cos(a + b*x)**(sympy.S(-4)/3)
    F = 3*sin(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(a + b*x)**2)/(b*sqrt(sin(a + b*x)**2)*cos(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_31():
    f = (c*cos(a + b*x))**(sympy.S(4)/3)
    F = -3*(c*cos(a + b*x))**(sympy.S(7)/3)*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(a + b*x)**2)/(7*b*c*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_32():
    f = (c*cos(a + b*x))**(sympy.S(2)/3)
    F = -3*(c*cos(a + b*x))**(sympy.S(5)/3)*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(a + b*x)**2)/(5*b*c*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_33():
    f = (c*cos(a + b*x))**(sympy.S(1)/3)
    F = -3*(c*cos(a + b*x))**(sympy.S(4)/3)*sin(a + b*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(a + b*x)**2)/(4*b*c*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_34():
    f = (c*cos(a + b*x))**(sympy.S(-1)/3)
    F = -3*(c*cos(a + b*x))**(sympy.S(2)/3)*sin(a + b*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(a + b*x)**2)/(2*b*c*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_35():
    f = (c*cos(a + b*x))**(sympy.S(-2)/3)
    F = -3*(c*cos(a + b*x))**(sympy.S(1)/3)*sin(a + b*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(a + b*x)**2)/(b*c*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_36():
    f = (c*cos(a + b*x))**(sympy.S(-4)/3)
    F = 3*sin(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(a + b*x)**2)/(b*c*(c*cos(a + b*x))**(sympy.S(1)/3)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_37():
    f = cos(a + b*x)**n
    F = -sin(a + b*x)*cos(a + b*x)**(n + 1)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*(n + 1)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_38():
    f = (c*cos(a + b*x))**n
    F = -(c*cos(a + b*x))**(n + 1)*sin(a + b*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*c*(n + 1)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_39():
    f = (a*cos(x)**2)**(sympy.S(5)/2)
    F = 8*a**2*sqrt(a*cos(x)**2)*tan(x)/15 + 4*a*(a*cos(x)**2)**(sympy.S(3)/2)*tan(x)/15 + (a*cos(x)**2)**(sympy.S(5)/2)*tan(x)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_40():
    f = (a*cos(x)**2)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*cos(x)**2)*tan(x)/3 + (a*cos(x)**2)**(sympy.S(3)/2)*tan(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_41():
    f = sqrt(a*cos(x)**2)
    F = sqrt(a*cos(x)**2)*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_42():
    f = 1/sqrt(a*cos(x)**2)
    F = cos(x)*atanh(sin(x))/sqrt(a*cos(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_43():
    f = (a*cos(x)**2)**(sympy.S(-3)/2)
    F = cos(x)*atanh(sin(x))/(2*a*sqrt(a*cos(x)**2)) + tan(x)/(2*a*sqrt(a*cos(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_44():
    f = (a*cos(x)**2)**(sympy.S(-5)/2)
    F = tan(x)/(4*a*(a*cos(x)**2)**(sympy.S(3)/2)) + 3*cos(x)*atanh(sin(x))/(8*a**2*sqrt(a*cos(x)**2)) + 3*tan(x)/(8*a**2*sqrt(a*cos(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_45():
    f = (a*cos(x)**3)**(sympy.S(5)/2)
    F = 2*a**2*sqrt(a*cos(x)**3)*sin(x)*cos(x)**5/15 + 26*a**2*sqrt(a*cos(x)**3)*sin(x)*cos(x)**3/165 + 78*a**2*sqrt(a*cos(x)**3)*sin(x)*cos(x)/385 + 26*a**2*sqrt(a*cos(x)**3)*tan(x)/77 + 26*a**2*sqrt(a*cos(x)**3)*elliptic_f(x/2, 2)/(77*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_46():
    f = (a*cos(x)**3)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*cos(x)**3)*sin(x)*cos(x)**2/9 + 14*a*sqrt(a*cos(x)**3)*sin(x)/45 + 14*a*sqrt(a*cos(x)**3)*elliptic_e(x/2, 2)/(15*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_47():
    f = sqrt(a*cos(x)**3)
    F = 2*sqrt(a*cos(x)**3)*tan(x)/3 + 2*sqrt(a*cos(x)**3)*elliptic_f(x/2, 2)/(3*cos(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_48():
    f = 1/sqrt(a*cos(x)**3)
    F = 2*sin(x)*cos(x)/sqrt(a*cos(x)**3) - 2*cos(x)**(sympy.S(3)/2)*elliptic_e(x/2, 2)/sqrt(a*cos(x)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_49():
    f = (a*cos(x)**3)**(sympy.S(-3)/2)
    F = 10*sin(x)/(21*a*sqrt(a*cos(x)**3)) + 10*cos(x)**(sympy.S(3)/2)*elliptic_f(x/2, 2)/(21*a*sqrt(a*cos(x)**3)) + 2*tan(x)*sec(x)/(7*a*sqrt(a*cos(x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_50():
    f = (a*cos(x)**3)**(sympy.S(-5)/2)
    F = 154*sin(x)*cos(x)/(195*a**2*sqrt(a*cos(x)**3)) - 154*cos(x)**(sympy.S(3)/2)*elliptic_e(x/2, 2)/(195*a**2*sqrt(a*cos(x)**3)) + 2*tan(x)*sec(x)**4/(13*a**2*sqrt(a*cos(x)**3)) + 22*tan(x)*sec(x)**2/(117*a**2*sqrt(a*cos(x)**3)) + 154*tan(x)/(585*a**2*sqrt(a*cos(x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_51():
    f = (a*cos(x)**4)**(sympy.S(5)/2)
    F = 63*a**2*x*sqrt(a*cos(x)**4)*sec(x)**2/256 + a**2*sqrt(a*cos(x)**4)*sin(x)*cos(x)**7/10 + 9*a**2*sqrt(a*cos(x)**4)*sin(x)*cos(x)**5/80 + 21*a**2*sqrt(a*cos(x)**4)*sin(x)*cos(x)**3/160 + 21*a**2*sqrt(a*cos(x)**4)*sin(x)*cos(x)/128 + 63*a**2*sqrt(a*cos(x)**4)*tan(x)/256
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_52():
    f = (a*cos(x)**4)**(sympy.S(3)/2)
    F = 5*a*x*sqrt(a*cos(x)**4)*sec(x)**2/16 + a*sqrt(a*cos(x)**4)*sin(x)*cos(x)**3/6 + 5*a*sqrt(a*cos(x)**4)*sin(x)*cos(x)/24 + 5*a*sqrt(a*cos(x)**4)*tan(x)/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_53():
    f = sqrt(a*cos(x)**4)
    F = x*sqrt(a*cos(x)**4)*sec(x)**2/2 + sqrt(a*cos(x)**4)*tan(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_54():
    f = 1/sqrt(a*cos(x)**4)
    F = sin(x)*cos(x)/sqrt(a*cos(x)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_55():
    f = (a*cos(x)**4)**(sympy.S(-3)/2)
    F = sin(x)**2*tan(x)**3/(5*a*sqrt(a*cos(x)**4)) + 2*sin(x)**2*tan(x)/(3*a*sqrt(a*cos(x)**4)) + sin(x)*cos(x)/(a*sqrt(a*cos(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_56():
    f = (a*cos(x)**4)**(sympy.S(-5)/2)
    F = sin(x)**2*tan(x)**7/(9*a**2*sqrt(a*cos(x)**4)) + 4*sin(x)**2*tan(x)**5/(7*a**2*sqrt(a*cos(x)**4)) + 6*sin(x)**2*tan(x)**3/(5*a**2*sqrt(a*cos(x)**4)) + 4*sin(x)**2*tan(x)/(3*a**2*sqrt(a*cos(x)**4)) + sin(x)*cos(x)/(a**2*sqrt(a*cos(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_57():
    f = (b*cos(c + d*x)**m)**n
    F = -(b*cos(c + d*x)**m)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, m*n/2 + sympy.S.Half), (m*n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_58():
    f = (c*cos(a + b*x)**m)**(sympy.S(5)/2)
    F = -2*c**2*sqrt(c*cos(a + b*x)**m)*sin(a + b*x)*cos(a + b*x)**(2*m + 1)*hyper((sympy.S.Half, 5*m/4 + sympy.S.Half), (5*m/4 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*(5*m + 2)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_59():
    f = (c*cos(a + b*x)**m)**(sympy.S(3)/2)
    F = -2*c*sqrt(c*cos(a + b*x)**m)*sin(a + b*x)*cos(a + b*x)**(m + 1)*hyper((sympy.S.Half, 3*m/4 + sympy.S.Half), (3*m/4 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*(3*m + 2)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_60():
    f = sqrt(c*cos(a + b*x)**m)
    F = -2*sqrt(c*cos(a + b*x)**m)*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, m/4 + sympy.S.Half), (m/4 + sympy.S(3)/2,), cos(a + b*x)**2)/(b*(m + 2)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_61():
    f = 1/sqrt(c*cos(a + b*x)**m)
    F = -2*sin(a + b*x)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - m/4), (sympy.S(3)/2 - m/4,), cos(a + b*x)**2)/(b*sqrt(c*cos(a + b*x)**m)*(2 - m)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_62():
    f = (c*cos(a + b*x)**m)**(sympy.S(-3)/2)
    F = -2*sin(a + b*x)*cos(a + b*x)**(1 - m)*hyper((sympy.S.Half, sympy.S.Half - 3*m/4), (sympy.S(3)/2 - 3*m/4,), cos(a + b*x)**2)/(b*c*sqrt(c*cos(a + b*x)**m)*(2 - 3*m)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_63():
    f = (c*cos(a + b*x)**m)**(sympy.S(-5)/2)
    F = -2*sin(a + b*x)*cos(a + b*x)**(1 - 2*m)*hyper((sympy.S.Half, sympy.S.Half - 5*m/4), (sympy.S(3)/2 - 5*m/4,), cos(a + b*x)**2)/(b*c**2*sqrt(c*cos(a + b*x)**m)*(2 - 5*m)*sqrt(sin(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_64():
    f = (c*cos(a + b*x)**m)**(1/m)
    F = (c*cos(a + b*x)**m)**(1/m)*tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_65():
    f = (a*(b*cos(c + d*x))**p)**n
    F = -(a*(b*cos(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(n*p + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_66():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**5
    F = 30*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*sqrt(b*cos(c + d*x))) + 30*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*b**2*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_67():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**4
    F = 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 14*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*b*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_68():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**3
    F = 10*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_69():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**2
    F = 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_70():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)
    F = 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_71():
    f = sqrt(b*cos(c + d*x))
    F = 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_72():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)
    F = 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_73():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)**2
    F = 2*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_74():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)**3
    F = 2*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_75():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*b*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_76():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*b**2*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_77():
    f = sqrt(b*cos(c + d*x))*sec(c + d*x)**6
    F = 2*b**5*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*b**3*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*b*sin(c + d*x)/(15*d*sqrt(b*cos(c + d*x))) - 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_78():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**4
    F = 30*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*sqrt(b*cos(c + d*x))) + 30*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*b*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_79():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 14*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 14*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_80():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 10*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_81():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = 6*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_82():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_83():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_84():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 2*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_85():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 2*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_86():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_87():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*b**2*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_88():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**6
    F = 2*b**5*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*b**3*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_89():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**7
    F = 2*b**6*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*b**4*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*b**2*sin(c + d*x)/(15*d*sqrt(b*cos(c + d*x))) - 14*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_90():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 30*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*sqrt(b*cos(c + d*x))) + 30*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_91():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 14*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 14*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_92():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = 10*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_93():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)
    F = 6*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_94():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = 2*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_95():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 2*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_96():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 2*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_97():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**4
    F = 2*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_98():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**5
    F = 2*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_99():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**6
    F = 2*b**5*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*b**3*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_100():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**7
    F = 2*b**6*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*b**4*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_101():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**8
    F = 2*b**7*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*b**5*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*b**3*sin(c + d*x)/(15*d*sqrt(b*cos(c + d*x))) - 14*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_102():
    f = (b*cos(c + d*x))**(sympy.S(7)/2)
    F = 10*b**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*b**3*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*b*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_103():
    f = cos(c + d*x)**6/sqrt(b*cos(c + d*x))
    F = 30*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*sqrt(b*cos(c + d*x))) + 30*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*b*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*b**3*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_104():
    f = cos(c + d*x)**5/sqrt(b*cos(c + d*x))
    F = 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(cos(c + d*x))) + 14*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*b**2*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_105():
    f = cos(c + d*x)**4/sqrt(b*cos(c + d*x))
    F = 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_106():
    f = cos(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_107():
    f = cos(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_108():
    f = cos(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_109():
    f = 1/sqrt(b*cos(c + d*x))
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_110():
    f = sec(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_111():
    f = sec(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_112():
    f = sec(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 2*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_113():
    f = sec(c + d*x)**4/sqrt(b*cos(c + d*x))
    F = 2*b**3*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*b*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_114():
    f = sec(c + d*x)**5/sqrt(b*cos(c + d*x))
    F = 2*b**4*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*b**2*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*sin(c + d*x)/(15*d*sqrt(b*cos(c + d*x))) - 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_115():
    f = cos(c + d*x)**7/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 30*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*b*d*sqrt(b*cos(c + d*x))) + 30*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*b**2*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*b**4*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_116():
    f = cos(c + d*x)**6/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(cos(c + d*x))) + 14*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*b**3*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_117():
    f = cos(c + d*x)**5/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x))) + 10*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**2*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_118():
    f = cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + 2*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_119():
    f = cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + 2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_120():
    f = cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_121():
    f = cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_122():
    f = (b*cos(c + d*x))**(sympy.S(-3)/2)
    F = 2*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_123():
    f = sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_124():
    f = sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*sin(c + d*x)/(5*b*d*sqrt(b*cos(c + d*x))) - 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_125():
    f = sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*sin(c + d*x)/(21*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_126():
    f = sec(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b**3*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*b*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*sin(c + d*x)/(15*b*d*sqrt(b*cos(c + d*x))) - 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_127():
    f = cos(c + d*x)**8/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 30*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*b**2*d*sqrt(b*cos(c + d*x))) + 30*sqrt(b*cos(c + d*x))*sin(c + d*x)/(77*b**3*d) + 18*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*b**5*d) + 2*(b*cos(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(11*b**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_128():
    f = cos(c + d*x)**7/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b**3*d*sqrt(cos(c + d*x))) + 14*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*b**4*d) + 2*(b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_129():
    f = cos(c + d*x)**6/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x))) + 10*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**3*d) + 2*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_130():
    f = cos(c + d*x)**5/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + 2*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_131():
    f = cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + 2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_132():
    f = cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_133():
    f = cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_134():
    f = cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) - 2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_135():
    f = (b*cos(c + d*x))**(sympy.S(-5)/2)
    F = 2*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_136():
    f = sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*sin(c + d*x)/(5*b**2*d*sqrt(b*cos(c + d*x))) - 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_137():
    f = sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b*sin(c + d*x)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/2)) + 10*sin(c + d*x)/(21*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_138():
    f = sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)/(9*d*(b*cos(c + d*x))**(sympy.S(9)/2)) + 14*sin(c + d*x)/(45*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 14*sin(c + d*x)/(15*b**2*d*sqrt(b*cos(c + d*x))) - 14*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_139():
    f = (b*cos(c + d*x))**(sympy.S(-7)/2)
    F = 2*sin(c + d*x)/(5*b*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*sin(c + d*x)/(5*b**3*d*sqrt(b*cos(c + d*x))) - 6*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_140():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)
    F = 3*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_141():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = -sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_142():
    f = sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_143():
    f = sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))
    F = sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_144():
    f = sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    F = x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_145():
    f = sqrt(b*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_146():
    f = sqrt(b*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_147():
    f = sqrt(b*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_148():
    f = sqrt(b*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_149():
    f = sqrt(b*cos(c + d*x))/cos(c + d*x)**(sympy.S(11)/2)
    F = 3*sqrt(b*cos(c + d*x))*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + 3*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_150():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 3*b*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_151():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = -b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_152():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_153():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_154():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_155():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_156():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_157():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_158():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_159():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(13)/2)
    F = 3*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + 3*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_160():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**5/(5*d*sqrt(cos(c + d*x))) - 2*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_161():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 3*b**2*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_162():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = -b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_163():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_164():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_165():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_166():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_167():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_168():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_169():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(13)/2)
    F = b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_170():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(15)/2)
    F = 3*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(8*d*cos(c + d*x)**(sympy.S(5)/2)) + b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(9)/2)) + 3*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_171():
    f = cos(c + d*x)**(sympy.S(11)/2)/sqrt(b*cos(c + d*x))
    F = sin(c + d*x)**5*sqrt(cos(c + d*x))/(5*d*sqrt(b*cos(c + d*x))) - 2*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x))) + sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_172():
    f = cos(c + d*x)**(sympy.S(9)/2)/sqrt(b*cos(c + d*x))
    F = 3*x*sqrt(cos(c + d*x))/(8*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*d*sqrt(b*cos(c + d*x))) + 3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_173():
    f = cos(c + d*x)**(sympy.S(7)/2)/sqrt(b*cos(c + d*x))
    F = -sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x))) + sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_174():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(b*cos(c + d*x))
    F = x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_175():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(b*cos(c + d*x))
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_176():
    f = sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    F = x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_177():
    f = 1/(sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_178():
    f = 1/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_179():
    f = 1/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_180():
    f = 1/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = sin(c + d*x)**3/(3*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_181():
    f = 1/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(9)/2))
    F = 3*sin(c + d*x)/(8*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sin(c + d*x)/(4*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + 3*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_182():
    f = cos(c + d*x)**(sympy.S(11)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 3*x*sqrt(cos(c + d*x))/(8*b*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b*d*sqrt(b*cos(c + d*x))) + 3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_183():
    f = cos(c + d*x)**(sympy.S(9)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = -sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b*d*sqrt(b*cos(c + d*x))) + sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_184():
    f = cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_185():
    f = cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_186():
    f = cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_187():
    f = sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_188():
    f = 1/((b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_189():
    f = 1/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_190():
    f = 1/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)**3/(3*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_191():
    f = 1/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = 3*sin(c + d*x)/(8*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sin(c + d*x)/(4*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + 3*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_192():
    f = cos(c + d*x)**(sympy.S(13)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 3*x*sqrt(cos(c + d*x))/(8*b**2*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*b**2*d*sqrt(b*cos(c + d*x))) + 3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_193():
    f = cos(c + d*x)**(sympy.S(11)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = -sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b**2*d*sqrt(b*cos(c + d*x))) + sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_194():
    f = cos(c + d*x)**(sympy.S(9)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_195():
    f = cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_196():
    f = cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_197():
    f = cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_198():
    f = sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_199():
    f = 1/((b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_200():
    f = 1/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)**3/(3*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_201():
    f = 1/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 3*sin(c + d*x)/(8*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + sin(c + d*x)/(4*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)) + 3*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(8*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_202():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)**m
    F = -3*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(3*m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_203():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)**2
    F = -3*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_204():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_205():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_206():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_207():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)**2
    F = 3*b*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_208():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)**3
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_209():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*cos(c + d*x)**m
    F = -3*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(3*m + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_210():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*cos(c + d*x)**2
    F = -3*(b*cos(c + d*x))**(sympy.S(11)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(17)/6,), cos(c + d*x)**2)/(11*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_211():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*cos(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_212():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_213():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_214():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**2
    F = 3*b*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_215():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**3
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_216():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*cos(c + d*x)**m
    F = -3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_217():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*cos(c + d*x)**2
    F = -3*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_218():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*cos(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_219():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_220():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)
    F = -3*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_221():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)**2
    F = -3*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_222():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)**3
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_223():
    f = cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_224():
    f = cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_225():
    f = cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_226():
    f = (b*cos(c + d*x))**(sympy.S(-1)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_227():
    f = sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_228():
    f = sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*b*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_229():
    f = sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_230():
    f = cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/6), (m/2 + sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_231():
    f = cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_232():
    f = cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_233():
    f = (b*cos(c + d*x))**(sympy.S(-2)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_234():
    f = sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_235():
    f = sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*b*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_236():
    f = sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-4)/3, sympy.S.Half), (sympy.S(-1)/3,), cos(c + d*x)**2)/(8*d*(b*cos(c + d*x))**(sympy.S(8)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_237():
    f = cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2 + sympy.S(-1)/6), (m/2 + sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(1 - 3*m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_238():
    f = cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_239():
    f = cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_240():
    f = (b*cos(c + d*x))**(sympy.S(-4)/3)
    F = 3*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_241():
    f = sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_242():
    f = sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*b*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_243():
    f = sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*b**2*sin(c + d*x)*hyper((sympy.S(-5)/3, sympy.S.Half), (sympy.S(-2)/3,), cos(c + d*x)**2)/(10*d*(b*cos(c + d*x))**(sympy.S(10)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_244():
    f = (a*cos(e + f*x))**m*(b*cos(e + f*x))**n
    F = -(a*cos(e + f*x))**(m + 1)*(b*cos(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*(m + n + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_245():
    f = (b*cos(c + d*x))**n*cos(c + d*x)**2
    F = -(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_246():
    f = (b*cos(c + d*x))**n*cos(c + d*x)
    F = -(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_247():
    f = (b*cos(c + d*x))**n
    F = -(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_248():
    f = (b*cos(c + d*x))**n*sec(c + d*x)
    F = -(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_249():
    f = (b*cos(c + d*x))**n*sec(c + d*x)**2
    F = b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_250():
    f = (b*cos(c + d*x))**n*sec(c + d*x)**3
    F = b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_251():
    f = (b*cos(c + d*x))**n*sec(c + d*x)**4
    F = b**3*(b*cos(c + d*x))**(n - 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/2), (n/2 + sympy.S(-1)/2,), cos(c + d*x)**2)/(d*(3 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_252():
    f = (b*cos(c + d*x))**n*cos(c + d*x)**(sympy.S(5)/2)
    F = -2*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_253():
    f = (b*cos(c + d*x))**n*cos(c + d*x)**(sympy.S(3)/2)
    F = -2*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_254():
    f = (b*cos(c + d*x))**n*sqrt(cos(c + d*x))
    F = -2*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_255():
    f = (b*cos(c + d*x))**n/sqrt(cos(c + d*x))
    F = -2*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_256():
    f = (b*cos(c + d*x))**n/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_257():
    f = (b*cos(c + d*x))**n/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_258():
    f = (b*cos(c + d*x))**n/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_259():
    f = (b*cos(c + d*x))**n/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-7)/4), (n/2 + sympy.S(-3)/4,), cos(c + d*x)**2)/(d*(7 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_260():
    f = (a*cos(e + f*x))**m*(b*sec(e + f*x))**n
    F = -(a*cos(e + f*x))**(m + 1)*(b*sec(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, m/2 - n/2 + sympy.S.Half), (m/2 - n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*(m - n + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_261():
    f = cos(a + b*x)*sqrt(csc(a + b*x))
    F = 2/(b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_262():
    f = cos(a + b*x)/sqrt(csc(a + b*x))
    F = 2/(3*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_263():
    f = cos(a + b*x)**2*sqrt(csc(a + b*x))
    F = 4*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b) + 2*cos(a + b*x)/(3*b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_264():
    f = cos(a + b*x)**2/sqrt(csc(a + b*x))
    F = 4*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b) + 2*cos(a + b*x)/(5*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_265():
    f = cos(x)**3*csc(x)**(sympy.S(9)/2)
    F = -2*csc(x)**(sympy.S(7)/2)/7 + 2*csc(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_266():
    f = cos(a + b*x)**3*sqrt(csc(a + b*x))
    F = 2/(b*sqrt(csc(a + b*x))) - 2/(5*b*csc(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_267():
    f = cos(a + b*x)**3/sqrt(csc(a + b*x))
    F = 2/(3*b*csc(a + b*x)**(sympy.S(3)/2)) - 2/(7*b*csc(a + b*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_268():
    f = cos(a + b*x)**4*sqrt(csc(a + b*x))
    F = 8*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(7*b) + 2*cos(a + b*x)**3/(7*b*sqrt(csc(a + b*x))) + 4*cos(a + b*x)/(7*b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_269():
    f = cos(a + b*x)**4/sqrt(csc(a + b*x))
    F = 8*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(15*b) + 2*cos(a + b*x)**3/(9*b*csc(a + b*x)**(sympy.S(3)/2)) + 4*cos(a + b*x)/(15*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_270():
    f = cos(x)*csc(x)**(sympy.S(7)/3)
    F = -3*csc(x)**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_271():
    f = sqrt(csc(a + b*x))*sec(a + b*x)
    F = -atan(sqrt(csc(a + b*x)))/b + atanh(sqrt(csc(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_272():
    f = sec(a + b*x)/sqrt(csc(a + b*x))
    F = atan(sqrt(csc(a + b*x)))/b + atanh(sqrt(csc(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_273():
    f = sqrt(csc(a + b*x))*sec(a + b*x)**2
    F = sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/b + sec(a + b*x)/(b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_274():
    f = sec(a + b*x)**2/sqrt(csc(a + b*x))
    F = -sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/b + sec(a + b*x)/(b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_275():
    f = sqrt(csc(a + b*x))*sec(a + b*x)**3
    F = -3*atan(sqrt(csc(a + b*x)))/(4*b) + 3*atanh(sqrt(csc(a + b*x)))/(4*b) + sec(a + b*x)**2/(2*b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_276():
    f = sec(a + b*x)**3/sqrt(csc(a + b*x))
    F = atan(sqrt(csc(a + b*x)))/(4*b) + atanh(sqrt(csc(a + b*x)))/(4*b) + sec(a + b*x)**2/(2*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_277():
    f = sqrt(csc(a + b*x))*sec(a + b*x)**4
    F = 5*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(6*b) + sec(a + b*x)**3/(3*b*sqrt(csc(a + b*x))) + 5*sec(a + b*x)/(6*b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_278():
    f = sec(a + b*x)**4/sqrt(csc(a + b*x))
    F = -sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(2*b) + sec(a + b*x)**3/(3*b*csc(a + b*x)**(sympy.S(3)/2)) + sec(a + b*x)/(2*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_279():
    f = (d*cos(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**p
    F = d*sqrt(d*cos(a + b*x))*csc(a + b*x)**(p - 1)*hyper((sympy.S(-1)/4, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), sin(a + b*x)**2)/(b*(1 - p)*(cos(a + b*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_280():
    f = sqrt(d*cos(a + b*x))*csc(a + b*x)**p
    F = d*(cos(a + b*x)**2)**(sympy.S(1)/4)*csc(a + b*x)**(p - 1)*hyper((sympy.S(1)/4, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), sin(a + b*x)**2)/(b*sqrt(d*cos(a + b*x))*(1 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_281():
    f = csc(a + b*x)**p/sqrt(d*cos(a + b*x))
    F = d*(cos(a + b*x)**2)**(sympy.S(3)/4)*csc(a + b*x)**(p - 1)*hyper((sympy.S(3)/4, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), sin(a + b*x)**2)/(b*(d*cos(a + b*x))**(sympy.S(3)/2)*(1 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_282():
    f = csc(a + b*x)**p/(d*cos(a + b*x))**(sympy.S(3)/2)
    F = (cos(a + b*x)**2)**(sympy.S(1)/4)*csc(a + b*x)**(p - 1)*hyper((sympy.S(5)/4, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), sin(a + b*x)**2)/(b*d*sqrt(d*cos(a + b*x))*(1 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_283():
    f = csc(a + b*x)**p/(d*cos(a + b*x))**(sympy.S(5)/2)
    F = (cos(a + b*x)**2)**(sympy.S(3)/4)*csc(a + b*x)**(p - 1)*hyper((sympy.S(7)/4, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), sin(a + b*x)**2)/(b*d*(d*cos(a + b*x))**(sympy.S(3)/2)*(1 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_284():
    f = cos(e + f*x)**m*csc(e + f*x)**n
    F = (cos(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)**(m - 1)*csc(e + f*x)**(n - 1)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_285():
    f = (a*cos(e + f*x))**m*csc(e + f*x)**n
    F = a*(a*cos(e + f*x))**(m - 1)*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*csc(e + f*x)**(n - 1)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_286():
    f = (b*csc(e + f*x))**n*cos(e + f*x)**m
    F = b*(b*csc(e + f*x))**(n - 1)*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)**(m - 1)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_287():
    f = (a*cos(e + f*x))**m*(b*csc(e + f*x))**n
    F = a*b*(a*cos(e + f*x))**(m - 1)*(b*csc(e + f*x))**(n - 1)*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*hyper((sympy.S.Half - m/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_288():
    f = (a*cos(e + f*x))**m*(b*csc(e + f*x))**(sympy.S(5)/2)
    F = -b*(a*cos(e + f*x))**(m + 1)*(b*csc(e + f*x))**(sympy.S(3)/2)*(sin(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(7)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_289():
    f = (a*cos(e + f*x))**m*(b*csc(e + f*x))**(sympy.S(3)/2)
    F = -b*(a*cos(e + f*x))**(m + 1)*sqrt(b*csc(e + f*x))*(sin(e + f*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(5)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_290():
    f = (a*cos(e + f*x))**m*sqrt(b*csc(e + f*x))
    F = -(a*cos(e + f*x))**(m + 1)*(b*csc(e + f*x))**(sympy.S(3)/2)*(sin(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*b*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_291():
    f = (a*cos(e + f*x))**m/sqrt(b*csc(e + f*x))
    F = -(a*cos(e + f*x))**(m + 1)*sqrt(b*csc(e + f*x))*(sin(e + f*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*b*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_0_a_cos_pow_m_b_trg_pow_n_292():
    f = (a*cos(e + f*x))**m/(b*csc(e + f*x))**(sympy.S(3)/2)
    F = -(a*cos(e + f*x))**(m + 1)*hyper((sympy.S(-1)/4, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*b*f*sqrt(b*csc(e + f*x))*(m + 1)*(sin(e + f*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F

