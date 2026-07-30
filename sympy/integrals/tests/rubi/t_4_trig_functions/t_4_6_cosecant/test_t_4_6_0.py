"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.0 (a csc)^m (b trg)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_1():
    f = csc(a + b*x)
    F = -atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_2():
    f = csc(a + b*x)**2
    F = -cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_3():
    f = csc(a + b*x)**3
    F = -cot(a + b*x)*csc(a + b*x)/(2*b) - atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_4():
    f = csc(a + b*x)**4
    F = -cot(a + b*x)**3/(3*b) - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_5():
    f = csc(a + b*x)**5
    F = -cot(a + b*x)*csc(a + b*x)**3/(4*b) - 3*cot(a + b*x)*csc(a + b*x)/(8*b) - 3*atanh(cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_6():
    f = csc(a + b*x)**6
    F = -cot(a + b*x)**5/(5*b) - 2*cot(a + b*x)**3/(3*b) - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_7():
    f = csc(a + b*x)**7
    F = -cot(a + b*x)*csc(a + b*x)**5/(6*b) - 5*cot(a + b*x)*csc(a + b*x)**3/(24*b) - 5*cot(a + b*x)*csc(a + b*x)/(16*b) - 5*atanh(cos(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_8():
    f = csc(a + b*x)**8
    F = -cot(a + b*x)**7/(7*b) - 3*cot(a + b*x)**5/(5*b) - cot(a + b*x)**3/b - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_9():
    f = csc(a + b*x)**(sympy.S(7)/2)
    F = -6*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b) - 2*cos(a + b*x)*csc(a + b*x)**(sympy.S(5)/2)/(5*b) - 6*cos(a + b*x)*sqrt(csc(a + b*x))/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_10():
    f = csc(a + b*x)**(sympy.S(5)/2)
    F = 2*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b) - 2*cos(a + b*x)*csc(a + b*x)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_11():
    f = csc(a + b*x)**(sympy.S(3)/2)
    F = -2*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/b - 2*cos(a + b*x)*sqrt(csc(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_12():
    f = sqrt(csc(a + b*x))
    F = 2*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_13():
    f = 1/sqrt(csc(a + b*x))
    F = 2*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_14():
    f = csc(a + b*x)**(sympy.S(-3)/2)
    F = 2*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b) - 2*cos(a + b*x)/(3*b*sqrt(csc(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_15():
    f = csc(a + b*x)**(sympy.S(-5)/2)
    F = 6*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b) - 2*cos(a + b*x)/(5*b*csc(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_16():
    f = csc(a + b*x)**(sympy.S(-7)/2)
    F = 10*sqrt(sin(a + b*x))*sqrt(csc(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(21*b) - 10*cos(a + b*x)/(21*b*sqrt(csc(a + b*x))) - 2*cos(a + b*x)/(7*b*csc(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_17():
    f = (c*csc(a + b*x))**(sympy.S(7)/2)
    F = -6*c**4*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))) - 6*c**3*sqrt(c*csc(a + b*x))*cos(a + b*x)/(5*b) - 2*c*(c*csc(a + b*x))**(sympy.S(5)/2)*cos(a + b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_18():
    f = (c*csc(a + b*x))**(sympy.S(5)/2)
    F = 2*c**2*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b) - 2*c*(c*csc(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_19():
    f = (c*csc(a + b*x))**(sympy.S(3)/2)
    F = -2*c**2*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(b*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))) - 2*c*sqrt(c*csc(a + b*x))*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_20():
    f = sqrt(c*csc(a + b*x))
    F = 2*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_21():
    f = 1/sqrt(c*csc(a + b*x))
    F = 2*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(b*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_22():
    f = (c*csc(a + b*x))**(sympy.S(-3)/2)
    F = -2*cos(a + b*x)/(3*b*c*sqrt(c*csc(a + b*x))) + 2*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(3*b*c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_23():
    f = (c*csc(a + b*x))**(sympy.S(-5)/2)
    F = -2*cos(a + b*x)/(5*b*c*(c*csc(a + b*x))**(sympy.S(3)/2)) + 6*elliptic_e(a/2 + b*x/2 - pi/4, 2)/(5*b*c**2*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_24():
    f = (c*csc(a + b*x))**(sympy.S(-7)/2)
    F = -2*cos(a + b*x)/(7*b*c*(c*csc(a + b*x))**(sympy.S(5)/2)) - 10*cos(a + b*x)/(21*b*c**3*sqrt(c*csc(a + b*x))) + 10*sqrt(c*csc(a + b*x))*sqrt(sin(a + b*x))*elliptic_f(a/2 + b*x/2 - pi/4, 2)/(21*b*c**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_25():
    f = csc(a + b*x)**(sympy.S(4)/3)
    F = -3*cos(a + b*x)*csc(a + b*x)**(sympy.S(1)/3)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), sin(a + b*x)**2)/(b*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_26():
    f = csc(a + b*x)**(sympy.S(2)/3)
    F = 3*cos(a + b*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), sin(a + b*x)**2)/(b*sqrt(cos(a + b*x)**2)*csc(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_27():
    f = csc(a + b*x)**(sympy.S(1)/3)
    F = 3*cos(a + b*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), sin(a + b*x)**2)/(2*b*sqrt(cos(a + b*x)**2)*csc(a + b*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_28():
    f = csc(a + b*x)**(sympy.S(-1)/3)
    F = 3*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), sin(a + b*x)**2)/(4*b*sqrt(cos(a + b*x)**2)*csc(a + b*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_29():
    f = csc(a + b*x)**(sympy.S(-2)/3)
    F = 3*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), sin(a + b*x)**2)/(5*b*sqrt(cos(a + b*x)**2)*csc(a + b*x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_30():
    f = csc(a + b*x)**(sympy.S(-4)/3)
    F = 3*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), sin(a + b*x)**2)/(7*b*sqrt(cos(a + b*x)**2)*csc(a + b*x)**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_31():
    f = (c*csc(a + b*x))**(sympy.S(4)/3)
    F = -3*c*(c*csc(a + b*x))**(sympy.S(1)/3)*cos(a + b*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), sin(a + b*x)**2)/(b*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_32():
    f = (c*csc(a + b*x))**(sympy.S(2)/3)
    F = 3*c*cos(a + b*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), sin(a + b*x)**2)/(b*(c*csc(a + b*x))**(sympy.S(1)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_33():
    f = (c*csc(a + b*x))**(sympy.S(1)/3)
    F = 3*c*cos(a + b*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), sin(a + b*x)**2)/(2*b*(c*csc(a + b*x))**(sympy.S(2)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_34():
    f = (c*csc(a + b*x))**(sympy.S(-1)/3)
    F = 3*c*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), sin(a + b*x)**2)/(4*b*(c*csc(a + b*x))**(sympy.S(4)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_35():
    f = (c*csc(a + b*x))**(sympy.S(-2)/3)
    F = 3*c*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), sin(a + b*x)**2)/(5*b*(c*csc(a + b*x))**(sympy.S(5)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_36():
    f = (c*csc(a + b*x))**(sympy.S(-4)/3)
    F = 3*c*cos(a + b*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), sin(a + b*x)**2)/(7*b*(c*csc(a + b*x))**(sympy.S(7)/3)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_37():
    f = csc(a + b*x)**n
    F = cos(a + b*x)*csc(a + b*x)**(n - 1)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(a + b*x)**2)/(b*(1 - n)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_38():
    f = (c*csc(a + b*x))**n
    F = c*(c*csc(a + b*x))**(n - 1)*cos(a + b*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(a + b*x)**2)/(b*(1 - n)*sqrt(cos(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_39():
    f = (csc(x)**2)**(sympy.S(7)/2)
    F = -(csc(x)**2)**(sympy.S(5)/2)*cot(x)/6 - 5*(csc(x)**2)**(sympy.S(3)/2)*cot(x)/24 - 5*sqrt(csc(x)**2)*cot(x)/16 - 5*asinh(cot(x))/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_40():
    f = (csc(x)**2)**(sympy.S(5)/2)
    F = -(csc(x)**2)**(sympy.S(3)/2)*cot(x)/4 - 3*sqrt(csc(x)**2)*cot(x)/8 - 3*asinh(cot(x))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_41():
    f = (csc(x)**2)**(sympy.S(3)/2)
    F = -sqrt(csc(x)**2)*cot(x)/2 - asinh(cot(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_42():
    f = sqrt(csc(x)**2)
    F = -asinh(cot(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_43():
    f = 1/sqrt(csc(x)**2)
    F = -cot(x)/sqrt(csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_44():
    f = (csc(x)**2)**(sympy.S(-3)/2)
    F = -2*cot(x)/(3*sqrt(csc(x)**2)) - cot(x)/(3*(csc(x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_45():
    f = (csc(x)**2)**(sympy.S(-5)/2)
    F = -8*cot(x)/(15*sqrt(csc(x)**2)) - 4*cot(x)/(15*(csc(x)**2)**(sympy.S(3)/2)) - cot(x)/(5*(csc(x)**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_46():
    f = (csc(x)**2)**(sympy.S(-7)/2)
    F = -16*cot(x)/(35*sqrt(csc(x)**2)) - 8*cot(x)/(35*(csc(x)**2)**(sympy.S(3)/2)) - 6*cot(x)/(35*(csc(x)**2)**(sympy.S(5)/2)) - cot(x)/(7*(csc(x)**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_47():
    f = (a*csc(x)**2)**(sympy.S(7)/2)
    F = -5*a**(sympy.S(7)/2)*atanh(sqrt(a)*cot(x)/sqrt(a*csc(x)**2))/16 - 5*a**3*sqrt(a*csc(x)**2)*cot(x)/16 - 5*a**2*(a*csc(x)**2)**(sympy.S(3)/2)*cot(x)/24 - a*(a*csc(x)**2)**(sympy.S(5)/2)*cot(x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_48():
    f = (a*csc(x)**2)**(sympy.S(5)/2)
    F = -3*a**(sympy.S(5)/2)*atanh(sqrt(a)*cot(x)/sqrt(a*csc(x)**2))/8 - 3*a**2*sqrt(a*csc(x)**2)*cot(x)/8 - a*(a*csc(x)**2)**(sympy.S(3)/2)*cot(x)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_49():
    f = (a*csc(x)**2)**(sympy.S(3)/2)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a)*cot(x)/sqrt(a*csc(x)**2))/2 - a*sqrt(a*csc(x)**2)*cot(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_50():
    f = sqrt(a*csc(x)**2)
    F = -sqrt(a)*atanh(sqrt(a)*cot(x)/sqrt(a*csc(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_51():
    f = 1/sqrt(a*csc(x)**2)
    F = -cot(x)/sqrt(a*csc(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_52():
    f = (a*csc(x)**2)**(sympy.S(-3)/2)
    F = -cot(x)/(3*(a*csc(x)**2)**(sympy.S(3)/2)) - 2*cot(x)/(3*a*sqrt(a*csc(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_53():
    f = (a*csc(x)**2)**(sympy.S(-5)/2)
    F = -cot(x)/(5*(a*csc(x)**2)**(sympy.S(5)/2)) - 4*cot(x)/(15*a*(a*csc(x)**2)**(sympy.S(3)/2)) - 8*cot(x)/(15*a**2*sqrt(a*csc(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_54():
    f = (a*csc(x)**2)**(sympy.S(-7)/2)
    F = -cot(x)/(7*(a*csc(x)**2)**(sympy.S(7)/2)) - 6*cot(x)/(35*a*(a*csc(x)**2)**(sympy.S(5)/2)) - 8*cot(x)/(35*a**2*(a*csc(x)**2)**(sympy.S(3)/2)) - 16*cot(x)/(35*a**3*sqrt(a*csc(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_55():
    f = (a*csc(x)**3)**(sympy.S(5)/2)
    F = -154*a**2*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2)*elliptic_e(x/2 - pi/4, 2)/195 - 154*a**2*sqrt(a*csc(x)**3)*sin(x)*cos(x)/195 - 2*a**2*sqrt(a*csc(x)**3)*cot(x)*csc(x)**4/13 - 22*a**2*sqrt(a*csc(x)**3)*cot(x)*csc(x)**2/117 - 154*a**2*sqrt(a*csc(x)**3)*cot(x)/585
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_56():
    f = (a*csc(x)**3)**(sympy.S(3)/2)
    F = 10*a*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2)*elliptic_f(x/2 - pi/4, 2)/21 - 10*a*sqrt(a*csc(x)**3)*cos(x)/21 - 2*a*sqrt(a*csc(x)**3)*cot(x)*csc(x)/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_57():
    f = sqrt(a*csc(x)**3)
    F = -2*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2)*elliptic_e(x/2 - pi/4, 2) - 2*sqrt(a*csc(x)**3)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_58():
    f = 1/sqrt(a*csc(x)**3)
    F = -2*cot(x)/(3*sqrt(a*csc(x)**3)) + 2*elliptic_f(x/2 - pi/4, 2)/(3*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_59():
    f = (a*csc(x)**3)**(sympy.S(-3)/2)
    F = -2*sin(x)**2*cos(x)/(9*a*sqrt(a*csc(x)**3)) - 14*cos(x)/(45*a*sqrt(a*csc(x)**3)) + 14*elliptic_e(x/2 - pi/4, 2)/(15*a*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_60():
    f = (a*csc(x)**3)**(sympy.S(-5)/2)
    F = -2*sin(x)**5*cos(x)/(15*a**2*sqrt(a*csc(x)**3)) - 26*sin(x)**3*cos(x)/(165*a**2*sqrt(a*csc(x)**3)) - 78*sin(x)*cos(x)/(385*a**2*sqrt(a*csc(x)**3)) - 26*cot(x)/(77*a**2*sqrt(a*csc(x)**3)) + 26*elliptic_f(x/2 - pi/4, 2)/(77*a**2*sqrt(a*csc(x)**3)*sin(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_61():
    f = (a*csc(x)**4)**(sympy.S(7)/2)
    F = -a**3*sqrt(a*csc(x)**4)*sin(x)*cos(x) - a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**11/13 - 6*a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**9/11 - 5*a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**7/3 - 20*a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**5/7 - 3*a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**3 - 2*a**3*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_62():
    f = (a*csc(x)**4)**(sympy.S(5)/2)
    F = -a**2*sqrt(a*csc(x)**4)*sin(x)*cos(x) - a**2*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**7/9 - 4*a**2*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**5/7 - 6*a**2*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**3/5 - 4*a**2*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_63():
    f = (a*csc(x)**4)**(sympy.S(3)/2)
    F = -a*sqrt(a*csc(x)**4)*sin(x)*cos(x) - a*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)**3/5 - 2*a*sqrt(a*csc(x)**4)*cos(x)**2*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_64():
    f = sqrt(a*csc(x)**4)
    F = -sqrt(a*csc(x)**4)*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_65():
    f = 1/sqrt(a*csc(x)**4)
    F = x*csc(x)**2/(2*sqrt(a*csc(x)**4)) - cot(x)/(2*sqrt(a*csc(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_66():
    f = (a*csc(x)**4)**(sympy.S(-3)/2)
    F = 5*x*csc(x)**2/(16*a*sqrt(a*csc(x)**4)) - sin(x)**3*cos(x)/(6*a*sqrt(a*csc(x)**4)) - 5*sin(x)*cos(x)/(24*a*sqrt(a*csc(x)**4)) - 5*cot(x)/(16*a*sqrt(a*csc(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_67():
    f = (a*csc(x)**4)**(sympy.S(-5)/2)
    F = 63*x*csc(x)**2/(256*a**2*sqrt(a*csc(x)**4)) - sin(x)**7*cos(x)/(10*a**2*sqrt(a*csc(x)**4)) - 9*sin(x)**5*cos(x)/(80*a**2*sqrt(a*csc(x)**4)) - 21*sin(x)**3*cos(x)/(160*a**2*sqrt(a*csc(x)**4)) - 21*sin(x)*cos(x)/(128*a**2*sqrt(a*csc(x)**4)) - 63*cot(x)/(256*a**2*sqrt(a*csc(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_68():
    f = ((b*csc(c + d*x))**p)**n
    F = ((b*csc(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(-n*p + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_69():
    f = (a*(b*csc(c + d*x))**p)**n
    F = (a*(b*csc(c + d*x))**p)**n*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(-n*p + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_0_a_csc_pow_m_b_trg_pow_n_70():
    f = (a*csc(e + f*x))**m*(b*csc(e + f*x))**n
    F = a*(a*csc(e + f*x))**(m - 1)*(b*csc(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -m/2 - n/2 + sympy.S.Half), (-m/2 - n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(-m - n + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F

