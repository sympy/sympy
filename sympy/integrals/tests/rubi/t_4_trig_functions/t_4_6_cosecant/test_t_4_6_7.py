"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.7 (d trig)^m (a+b (c csc)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_1():
    f = (a + b*csc(c + d*x)**2)**4
    F = a**4*x - b**4*cot(c + d*x)**7/(7*d) - b**3*(4*a + 3*b)*cot(c + d*x)**5/(5*d) - b**2*(6*a**2 + 8*a*b + 3*b**2)*cot(c + d*x)**3/(3*d) - b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_2():
    f = (a + b*csc(c + d*x)**2)**3
    F = a**3*x - b**3*cot(c + d*x)**5/(5*d) - b**2*(3*a + 2*b)*cot(c + d*x)**3/(3*d) - b*(3*a**2 + 3*a*b + b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_3():
    f = (a + b*csc(c + d*x)**2)**2
    F = a**2*x - b**2*cot(c + d*x)**3/(3*d) - b*(2*a + b)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_4():
    f = a + b*csc(c + d*x)**2
    F = a*x - b*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_5():
    f = 1/(a + b*csc(c + d*x)**2)
    F = -sqrt(b)*atan(sqrt(a + b)*tan(c + d*x)/sqrt(b))/(a*d*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_6():
    f = (a + b*csc(c + d*x)**2)**(-2)
    F = b*cot(c + d*x)/(2*a*d*(a + b)*(a + b*cot(c + d*x)**2 + b)) + sqrt(b)*(3*a + 2*b)*atan(sqrt(b)*cot(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_7():
    f = (a + b*csc(c + d*x)**2)**(-3)
    F = b*cot(c + d*x)/(4*a*d*(a + b)*(a + b*cot(c + d*x)**2 + b)**2) + b*(7*a + 4*b)*cot(c + d*x)/(8*a**2*d*(a + b)**2*(a + b*cot(c + d*x)**2 + b)) + sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atan(sqrt(b)*cot(c + d*x)/sqrt(a + b))/(8*a**3*d*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_8():
    f = (a + b*csc(c + d*x)**2)**(-4)
    F = b*cot(c + d*x)/(6*a*d*(a + b)*(a + b*cot(c + d*x)**2 + b)**3) + b*(11*a + 6*b)*cot(c + d*x)/(24*a**2*d*(a + b)**2*(a + b*cot(c + d*x)**2 + b)**2) + b*(19*a**2 + 22*a*b + 8*b**2)*cot(c + d*x)/(16*a**3*d*(a + b)**3*(a + b*cot(c + d*x)**2 + b)) + sqrt(b)*(35*a**3 + 70*a**2*b + 56*a*b**2 + 16*b**3)*atan(sqrt(b)*cot(c + d*x)/sqrt(a + b))/(16*a**4*d*(a + b)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_9():
    f = (a + b*csc(c + d*x)**2)**(sympy.S(5)/2)
    F = -a**(sympy.S(5)/2)*atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/d - sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/(8*d) - b*(7*a + 3*b)*sqrt(a + b*cot(c + d*x)**2 + b)*cot(c + d*x)/(8*d) - b*(a + b*cot(c + d*x)**2 + b)**(sympy.S(3)/2)*cot(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_10():
    f = (a + b*csc(c + d*x)**2)**(sympy.S(3)/2)
    F = -a**(sympy.S(3)/2)*atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/d - sqrt(b)*(3*a + b)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/(2*d) - b*sqrt(a + b*cot(c + d*x)**2 + b)*cot(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_11():
    f = sqrt(a + b*csc(c + d*x)**2)
    F = -sqrt(a)*atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/d - sqrt(b)*atanh(sqrt(b)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_12():
    f = (a + b*csc(c + d*x)**2)**(sympy.S(-3)/2)
    F = b*cot(c + d*x)/(a*d*(a + b)*sqrt(a + b*cot(c + d*x)**2 + b)) - atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_13():
    f = (a + b*csc(c + d*x)**2)**(sympy.S(-5)/2)
    F = b*cot(c + d*x)/(3*a*d*(a + b)*(a + b*cot(c + d*x)**2 + b)**(sympy.S(3)/2)) + b*(5*a + 3*b)*cot(c + d*x)/(3*a**2*d*(a + b)**2*sqrt(a + b*cot(c + d*x)**2 + b)) - atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_14():
    f = (a + b*csc(c + d*x)**2)**(sympy.S(-7)/2)
    F = b*cot(c + d*x)/(5*a*d*(a + b)*(a + b*cot(c + d*x)**2 + b)**(sympy.S(5)/2)) + b*(9*a + 5*b)*cot(c + d*x)/(15*a**2*d*(a + b)**2*(a + b*cot(c + d*x)**2 + b)**(sympy.S(3)/2)) + b*(33*a**2 + 40*a*b + 15*b**2)*cot(c + d*x)/(15*a**3*d*(a + b)**3*sqrt(a + b*cot(c + d*x)**2 + b)) - atan(sqrt(a)*cot(c + d*x)/sqrt(a + b*cot(c + d*x)**2 + b))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_15():
    f = (csc(x)**2 + 1)**(sympy.S(3)/2)
    F = -sqrt(cot(x)**2 + 2)*cot(x)/2 - 2*asinh(sqrt(2)*cot(x)/2) - atan(cot(x)/sqrt(cot(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_16():
    f = sqrt(csc(x)**2 + 1)
    F = -asinh(sqrt(2)*cot(x)/2) - atan(cot(x)/sqrt(cot(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_17():
    f = 1/sqrt(csc(x)**2 + 1)
    F = -atan(cot(x)/sqrt(cot(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_18():
    f = (1 - csc(x)**2)**(sympy.S(3)/2)
    F = sqrt(-cot(x)**2)*log(sin(x))*tan(x) + sqrt(-cot(x)**2)*cot(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_19():
    f = sqrt(1 - csc(x)**2)
    F = sqrt(-cot(x)**2)*log(sin(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_20():
    f = 1/sqrt(1 - csc(x)**2)
    F = -log(cos(x))*cot(x)/sqrt(-cot(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_21():
    f = (csc(x)**2 - 1)**(sympy.S(3)/2)
    F = -(cot(x)**2)**(sympy.S(3)/2)*tan(x)/2 - sqrt(cot(x)**2)*log(sin(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_22():
    f = sqrt(csc(x)**2 - 1)
    F = sqrt(cot(x)**2)*log(sin(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_23():
    f = 1/sqrt(csc(x)**2 - 1)
    F = -log(cos(x))*cot(x)/sqrt(cot(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_24():
    f = (-csc(x)**2 - 1)**(sympy.S(3)/2)
    F = sqrt(-cot(x)**2 - 2)*cot(x)/2 - 2*atan(cot(x)/sqrt(-cot(x)**2 - 2)) - atanh(cot(x)/sqrt(-cot(x)**2 - 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_25():
    f = sqrt(-csc(x)**2 - 1)
    F = atan(cot(x)/sqrt(-cot(x)**2 - 2)) + atanh(cot(x)/sqrt(-cot(x)**2 - 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_7_d_trig_pow_m_a_plus_b_c_csc_pow_n_pow_p_26():
    f = 1/sqrt(-csc(x)**2 - 1)
    F = -atanh(cot(x)/sqrt(-cot(x)**2 - 2))
    assert integrate(f, x) == F

