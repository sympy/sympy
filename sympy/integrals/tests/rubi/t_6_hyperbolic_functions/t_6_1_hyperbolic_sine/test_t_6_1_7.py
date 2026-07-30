"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.1 Hyperbolic sine/6.1.7 hyper^m (a+b sinh^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_1():
    f = (a + b*sinh(c + d*x)**2)*sinh(c + d*x)**4
    F = b*sinh(c + d*x)**5*cosh(c + d*x)/(6*d) + x*(6*a - 5*b)/16 + (6*a - 5*b)*sinh(c + d*x)**3*cosh(c + d*x)/(24*d) - (6*a - 5*b)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_2():
    f = (a + b*sinh(c + d*x)**2)*sinh(c + d*x)**3
    F = b*cosh(c + d*x)**5/(5*d) + (a - 2*b)*cosh(c + d*x)**3/(3*d) - (a - b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_3():
    f = (a + b*sinh(c + d*x)**2)*sinh(c + d*x)**2
    F = b*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - x*(4*a - 3*b)/8 + (4*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_4():
    f = (a + b*sinh(c + d*x)**2)*sinh(c + d*x)
    F = b*cosh(c + d*x)**3/(3*d) + (a - b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_5():
    f = a + b*sinh(c + d*x)**2
    F = a*x - b*x/2 + b*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_6():
    f = (a + b*sinh(c + d*x)**2)*csch(c + d*x)
    F = -a*atanh(cosh(c + d*x))/d + b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_7():
    f = (a + b*sinh(c + d*x)**2)*csch(c + d*x)**2
    F = -a*coth(c + d*x)/d + b*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_8():
    f = (a + b*sinh(c + d*x)**2)*csch(c + d*x)**3
    F = -a*coth(c + d*x)*csch(c + d*x)/(2*d) + (a - 2*b)*atanh(cosh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_9():
    f = (a + b*sinh(c + d*x)**2)*csch(c + d*x)**4
    F = -a*coth(c + d*x)*csch(c + d*x)**2/(3*d) + (2*a - 3*b)*coth(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_10():
    f = (a + b*sinh(c + d*x)**2)**2*sinh(c + d*x)**4
    F = b**2*sinh(c + d*x)**5*cosh(c + d*x)**3/(8*d) + b*(16*a - 13*b)*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) + x*(3*a**2/8 - 5*a*b/8 + 35*b**2/128) + (48*a**2 - 208*a*b + 139*b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) - (80*a**2 - 176*a*b + 93*b**2)*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_11():
    f = (a + b*sinh(c + d*x)**2)**2*sinh(c + d*x)**3
    F = b**2*cosh(c + d*x)**7/(7*d) + b*(2*a - 3*b)*cosh(c + d*x)**5/(5*d) + (a - 3*b)*(a - b)*cosh(c + d*x)**3/(3*d) - (a - b)**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_12():
    f = (a + b*sinh(c + d*x)**2)**2*sinh(c + d*x)
    F = b**2*cosh(c + d*x)**5/(5*d) + b*(2*a - 2*b)*cosh(c + d*x)**3/(3*d) + (a - b)**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_13():
    f = (a + b*sinh(c + d*x)**2)**2
    F = b**2*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) + b*(8*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) + x*(a**2 - a*b + 3*b**2/8)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_14():
    f = (a + b*sinh(c + d*x)**2)**2*csch(c + d*x)
    F = -a**2*atanh(cosh(c + d*x))/d + b**2*cosh(c + d*x)**3/(3*d) + b*(2*a - b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_15():
    f = (a + b*sinh(c + d*x)**2)**2*csch(c + d*x)**3
    F = -a**2*coth(c + d*x)*csch(c + d*x)/(2*d) + a*(a - 4*b)*atanh(cosh(c + d*x))/(2*d) + b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_16():
    f = (a + b*sinh(c + d*x)**2)**2*csch(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) + a*(a - 2*b)*coth(c + d*x)/d + b**2*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_17():
    f = (a + b*sinh(c + d*x)**2)**3*sinh(c + d*x)**4
    F = -b*(a*(14*a - 9*b) - (a - b)*(22*a - 21*b)*tanh(c + d*x)**2)*sinh(c + d*x)**3*cosh(c + d*x)**3/(160*d) + x*(3*a/64 - 9*b/256)*(8*a**2 - 14*a*b + 7*b**2) + (a - (a - b)*tanh(c + d*x)**2)**3*sinh(c + d*x)**3*cosh(c + d*x)**7/(10*d) + (a - (a - b)*tanh(c + d*x)**2)**2*(6*a - 9*b)*sinh(c + d*x)**3*cosh(c + d*x)**5/(80*d) + (48*a**3 - 272*a**2*b + 314*a*b**2 - 105*b**3)*sinh(c + d*x)*cosh(c + d*x)**3/(640*d) - (576*a**3 - 1744*a**2*b + 1678*a*b**2 - 525*b**3)*sinh(c + d*x)*cosh(c + d*x)/(1280*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_18():
    f = (a + b*sinh(c + d*x)**2)**3*sinh(c + d*x)**3
    F = b**3*cosh(c + d*x)**9/(9*d) + b**2*(3*a - 4*b)*cosh(c + d*x)**7/(7*d) + b*(a - b)*(3*a - 6*b)*cosh(c + d*x)**5/(5*d) + (a - 4*b)*(a - b)**2*cosh(c + d*x)**3/(3*d) - (a - b)**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_19():
    f = (a + b*sinh(c + d*x)**2)**3*sinh(c + d*x)**2
    F = b*(24*a**2 - 64*a*b + 35*b**2)*sinh(c + d*x)**3*cosh(c + d*x)/(192*d) + x*(-a**3/2 + 9*a**2*b/8 - 15*a*b**2/16 + 35*b**3/128) + (a + b*sinh(c + d*x)**2)**3*sinh(c + d*x)*cosh(c + d*x)/(8*d) + (a + b*sinh(c + d*x)**2)**2*(6*a - 7*b)*sinh(c + d*x)*cosh(c + d*x)/(48*d) + (96*a**3 - 376*a**2*b + 360*a*b**2 - 105*b**3)*sinh(c + d*x)*cosh(c + d*x)/(384*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_20():
    f = (a + b*sinh(c + d*x)**2)**3*sinh(c + d*x)
    F = b**3*cosh(c + d*x)**7/(7*d) + b**2*(3*a - 3*b)*cosh(c + d*x)**5/(5*d) + b*(a - b)**2*cosh(c + d*x)**3/d + (a - b)**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_21():
    f = (a + b*sinh(c + d*x)**2)**3
    F = b**2*(10*a - 5*b)*sinh(c + d*x)**3*cosh(c + d*x)/(24*d) + b*(a + b*sinh(c + d*x)**2)**2*sinh(c + d*x)*cosh(c + d*x)/(6*d) + b*(64*a**2 - 54*a*b + 15*b**2)*sinh(c + d*x)*cosh(c + d*x)/(48*d) + x*(a/8 - b/16)*(8*a**2 - 8*a*b + 5*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_22():
    f = (a + b*sinh(c + d*x)**2)**3*csch(c + d*x)
    F = -a**3*atanh(cosh(c + d*x))/d + b**3*cosh(c + d*x)**5/(5*d) + b**2*(3*a - 2*b)*cosh(c + d*x)**3/(3*d) + b*(3*a**2 - 3*a*b + b**2)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_23():
    f = (a + b*sinh(c + d*x)**2)**3*csch(c + d*x)**2
    F = -a*(2*a + b)*(4*a + b)*coth(c + d*x)/(8*d) + 3*b*x*(8*a**2 - 4*a*b + b**2)/8 + b*(a - (a - b)*tanh(c + d*x)**2)**2*cosh(c + d*x)**4*coth(c + d*x)/(4*d) + b*(a*(4*a + b) - (a - b)*(4*a - 3*b)*tanh(c + d*x)**2)*cosh(c + d*x)**2*coth(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_24():
    f = (a + b*sinh(c + d*x)**2)**3*csch(c + d*x)**3
    F = -a**3*coth(c + d*x)*csch(c + d*x)/(2*d) + a**2*(a - 6*b)*atanh(cosh(c + d*x))/(2*d) + b**3*cosh(c + d*x)**3/(3*d) + b**2*(3*a - b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_25():
    f = (a + b*sinh(c + d*x)**2)**3*csch(c + d*x)**4
    F = -a**2*(2*a + 3*b)*coth(c + d*x)**3/(6*d) + a*(2*a**2 - 5*a*b - 2*b**2)*coth(c + d*x)/(2*d) + b**2*x*(3*a - b/2) + b*(a - (a - b)*tanh(c + d*x)**2)**2*cosh(c + d*x)**2*coth(c + d*x)**3/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_26():
    f = sinh(c + d*x)**7/(a + b*sinh(c + d*x)**2)
    F = -a**3*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(b**(sympy.S(7)/2)*d*sqrt(a - b)) + cosh(c + d*x)**5/(5*b*d) - (a + 2*b)*cosh(c + d*x)**3/(3*b**2*d) + (a**2 + a*b + b**2)*cosh(c + d*x)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_27():
    f = sinh(c + d*x)**6/(a + b*sinh(c + d*x)**2)
    F = -a**(sympy.S(5)/2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(b**3*d*sqrt(a - b)) + sinh(c + d*x)**3*cosh(c + d*x)/(4*b*d) - (4*a + 3*b)*sinh(c + d*x)*cosh(c + d*x)/(8*b**2*d) + x*(8*a**2 + 4*a*b + 3*b**2)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_28():
    f = sinh(c + d*x)**5/(a + b*sinh(c + d*x)**2)
    F = a**2*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(b**(sympy.S(5)/2)*d*sqrt(a - b)) + cosh(c + d*x)**3/(3*b*d) - (a + b)*cosh(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_29():
    f = sinh(c + d*x)**4/(a + b*sinh(c + d*x)**2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(b**2*d*sqrt(a - b)) + sinh(c + d*x)*cosh(c + d*x)/(2*b*d) - x*(2*a + b)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_30():
    f = sinh(c + d*x)**3/(a + b*sinh(c + d*x)**2)
    F = -a*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(b**(sympy.S(3)/2)*d*sqrt(a - b)) + cosh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_31():
    f = sinh(c + d*x)**2/(a + b*sinh(c + d*x)**2)
    F = -sqrt(a)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(b*d*sqrt(a - b)) + x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_32():
    f = sinh(c + d*x)/(a + b*sinh(c + d*x)**2)
    F = atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(sqrt(b)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_33():
    f = 1/(a + b*sinh(c + d*x)**2)
    F = atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_34():
    f = csch(c + d*x)/(a + b*sinh(c + d*x)**2)
    F = -sqrt(b)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(a*d*sqrt(a - b)) - atanh(cosh(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_35():
    f = csch(c + d*x)**2/(a + b*sinh(c + d*x)**2)
    F = -coth(c + d*x)/(a*d) - b*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(3)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_36():
    f = csch(c + d*x)**3/(a + b*sinh(c + d*x)**2)
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d) + b**(sympy.S(3)/2)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(a**2*d*sqrt(a - b)) + (a + 2*b)*atanh(cosh(c + d*x))/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_37():
    f = csch(c + d*x)**4/(a + b*sinh(c + d*x)**2)
    F = -coth(c + d*x)**3/(3*a*d) + (a + b)*coth(c + d*x)/(a**2*d) + b**2*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(5)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_38():
    f = csch(c + d*x)**5/(a + b*sinh(c + d*x)**2)
    F = -coth(c + d*x)*csch(c + d*x)**3/(4*a*d) + (3*a + 4*b)*coth(c + d*x)*csch(c + d*x)/(8*a**2*d) - b**(sympy.S(5)/2)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(a**3*d*sqrt(a - b)) - (3*a**2 + 4*a*b + 8*b**2)*atanh(cosh(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_39():
    f = csch(c + d*x)**6/(a + b*sinh(c + d*x)**2)
    F = -coth(c + d*x)**5/(5*a*d) + (2*a + b)*coth(c + d*x)**3/(3*a**2*d) - (a**2 + a*b + b**2)*coth(c + d*x)/(a**3*d) - b**3*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(7)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_40():
    f = sinh(c + d*x)**4/(a + b*sinh(c + d*x)**2)**2
    F = -sqrt(a)*(2*a - 3*b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*b**2*d*(a - b)**(sympy.S(3)/2)) - a*tanh(c + d*x)/(b*d*(a - (a - b)*tanh(c + d*x)**2)*(2*a - 2*b)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_41():
    f = sinh(c + d*x)**3/(a + b*sinh(c + d*x)**2)**2
    F = -a*cosh(c + d*x)/(b*d*(2*a - 2*b)*(a + b*cosh(c + d*x)**2 - b)) + (a - 2*b)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(2*b**(sympy.S(3)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_42():
    f = sinh(c + d*x)**2/(a + b*sinh(c + d*x)**2)**2
    F = sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*sinh(c + d*x)**2)*(2*a - 2*b)) - atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*sqrt(a)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_43():
    f = sinh(c + d*x)/(a + b*sinh(c + d*x)**2)**2
    F = cosh(c + d*x)/(d*(2*a - 2*b)*(a + b*cosh(c + d*x)**2 - b)) + atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(2*sqrt(b)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_44():
    f = (a + b*sinh(c + d*x)**2)**(-2)
    F = -b*sinh(c + d*x)*cosh(c + d*x)/(2*a*d*(a - b)*(a + b*sinh(c + d*x)**2)) + (2*a - b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_45():
    f = csch(c + d*x)/(a + b*sinh(c + d*x)**2)**2
    F = -b*cosh(c + d*x)/(2*a*d*(a - b)*(a + b*cosh(c + d*x)**2 - b)) - sqrt(b)*(3*a - 2*b)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(2*a**2*d*(a - b)**(sympy.S(3)/2)) - atanh(cosh(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_46():
    f = csch(c + d*x)**2/(a + b*sinh(c + d*x)**2)**2
    F = -coth(c + d*x)/(a*d*(a - (a - b)*tanh(c + d*x)**2)) + (2*a**2 - 4*a*b + 3*b**2)*tanh(c + d*x)/(2*a**2*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)) - b*(4*a - 3*b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(5)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_47():
    f = csch(c + d*x)**3/(a + b*sinh(c + d*x)**2)**2
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d*(a + b*cosh(c + d*x)**2 - b)) - b*(a - 2*b)*cosh(c + d*x)/(2*a**2*d*(a - b)*(a + b*cosh(c + d*x)**2 - b)) + b**(sympy.S(3)/2)*(5*a - 4*b)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(2*a**3*d*(a - b)**(sympy.S(3)/2)) + (a + 4*b)*atanh(cosh(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_48():
    f = csch(c + d*x)**4/(a + b*sinh(c + d*x)**2)**2
    F = -b*csch(c + d*x)**3*sech(c + d*x)/(2*a*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)) - (2*a - 5*b)*coth(c + d*x)**3/(6*a**2*d*(a - b)) + (2*a**2 + a*b - 5*b**2)*coth(c + d*x)/(2*a**3*d*(a - b)) + b**2*(6*a - 5*b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_49():
    f = sinh(c + d*x)**4/(a + b*sinh(c + d*x)**2)**3
    F = tanh(c + d*x)**3/(d*(a - (a - b)*tanh(c + d*x)**2)**2*(4*a - 4*b)) - 3*tanh(c + d*x)/(8*d*(a - b)**2*(a - (a - b)*tanh(c + d*x)**2)) + 3*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*sqrt(a)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_50():
    f = sinh(c + d*x)**3/(a + b*sinh(c + d*x)**2)**3
    F = -a*cosh(c + d*x)/(b*d*(4*a - 4*b)*(a + b*cosh(c + d*x)**2 - b)**2) + (a - 4*b)*cosh(c + d*x)/(8*b*d*(a - b)**2*(a + b*cosh(c + d*x)**2 - b)) + (a - 4*b)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(8*b**(sympy.S(3)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_51():
    f = sinh(c + d*x)**2/(a + b*sinh(c + d*x)**2)**3
    F = sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*sinh(c + d*x)**2)**2*(4*a - 4*b)) + (2*a + b)*sinh(c + d*x)*cosh(c + d*x)/(8*a*d*(a - b)**2*(a + b*sinh(c + d*x)**2)) - (4*a - b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_52():
    f = sinh(c + d*x)/(a + b*sinh(c + d*x)**2)**3
    F = cosh(c + d*x)/(d*(4*a - 4*b)*(a + b*cosh(c + d*x)**2 - b)**2) + 3*cosh(c + d*x)/(8*d*(a - b)**2*(a + b*cosh(c + d*x)**2 - b)) + 3*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(8*sqrt(b)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_53():
    f = (a + b*sinh(c + d*x)**2)**(-3)
    F = -b*sinh(c + d*x)*cosh(c + d*x)/(4*a*d*(a - b)*(a + b*sinh(c + d*x)**2)**2) - b*(6*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*sinh(c + d*x)**2)) + (8*a**2 - 8*a*b + 3*b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_54():
    f = csch(c + d*x)/(a + b*sinh(c + d*x)**2)**3
    F = -b*cosh(c + d*x)/(4*a*d*(a - b)*(a + b*cosh(c + d*x)**2 - b)**2) - b*(7*a - 4*b)*cosh(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*cosh(c + d*x)**2 - b)) - sqrt(b)*(15*a**2 - 20*a*b + 8*b**2)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(8*a**3*d*(a - b)**(sympy.S(5)/2)) - atanh(cosh(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_55():
    f = csch(c + d*x)**2/(a + b*sinh(c + d*x)**2)**3
    F = -b*csch(c + d*x)*sech(c + d*x)**3/(4*a*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)**2) - b*(4*a - 5*b - (4*a - b)*tanh(c + d*x)**2)*coth(c + d*x)/(8*a**2*d*(a - b)**2*(a - (a - b)*tanh(c + d*x)**2)) - (2*a - 3*b)*(4*a - 5*b)*coth(c + d*x)/(8*a**3*d*(a - b)**2) - 3*b*(8*a**2 - 12*a*b + 5*b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_56():
    f = csch(c + d*x)**3/(a + b*sinh(c + d*x)**2)**3
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d*(a + b*cosh(c + d*x)**2 - b)**2) - b*(2*a - 3*b)*cosh(c + d*x)/(4*a**2*d*(a - b)*(a + b*cosh(c + d*x)**2 - b)**2) - b*(a - 4*b)*(4*a - 3*b)*cosh(c + d*x)/(8*a**3*d*(a - b)**2*(a + b*cosh(c + d*x)**2 - b)) + b**(sympy.S(3)/2)*(35*a**2 - 56*a*b + 24*b**2)*atan(sqrt(b)*cosh(c + d*x)/sqrt(a - b))/(8*a**4*d*(a - b)**(sympy.S(5)/2)) + (a + 6*b)*atanh(cosh(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_57():
    f = csch(c + d*x)**4/(a + b*sinh(c + d*x)**2)**3
    F = -b*csch(c + d*x)**3*sech(c + d*x)**3/(4*a*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)**2) - b*(10*a - 7*b)*csch(c + d*x)**3*sech(c + d*x)/(8*a**2*d*(a - b)**2*(a - (a - b)*tanh(c + d*x)**2)) - (8*a**2 - 52*a*b + 35*b**2)*coth(c + d*x)**3/(24*a**3*d*(a - b)**2) + (8*a**3 - 4*a**2*b - 45*a*b**2 + 35*b**3)*coth(c + d*x)/(8*a**4*d*(a - b)**2) + b**2*(48*a**2 - 80*a*b + 35*b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(9)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_58():
    f = 1/(sinh(x)**2 + 1)
    F = tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_59():
    f = (sinh(x)**2 + 1)**(-2)
    F = -tanh(x)**3/3 + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_60():
    f = (sinh(x)**2 + 1)**(-3)
    F = tanh(x)**5/5 - 2*tanh(x)**3/3 + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_61():
    f = 1/(1 - sinh(x)**2)
    F = sqrt(2)*atanh(sqrt(2)*tanh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_62():
    f = (1 - sinh(x)**2)**(-2)
    F = 3*sqrt(2)*atanh(sqrt(2)*tanh(x))/8 + sinh(x)*cosh(x)/(4 - 4*sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_63():
    f = (1 - sinh(x)**2)**(-3)
    F = 19*sqrt(2)*atanh(sqrt(2)*tanh(x))/64 + 9*sinh(x)*cosh(x)/(32 - 32*sinh(x)**2) + sinh(x)*cosh(x)/(8*(1 - sinh(x)**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_64():
    f = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**3
    F = -(a + 3*b)*sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(8*b*f) + (a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*cosh(e + f*x)/(4*b*f) - (a - b)*(a + 3*b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(8*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_65():
    f = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)
    F = sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(2*f) + (a - b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_66():
    f = sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/f + sqrt(b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_67():
    f = sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**3
    F = -sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(2*f) + (a - b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_68():
    f = sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**5
    F = (3*a + b)*sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(8*a*f) - (a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*coth(e + f*x)*csch(e + f*x)**3/(4*a*f) - (a - b)*(3*a + b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_69():
    f = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**4
    F = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**3*cosh(e + f*x)/(5*f) + (a - 4*b)*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(15*b*f) - (a - 4*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a**2 + 3*a*b - 8*b**2)*tanh(e + f*x)/(15*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a**2 + 3*a*b - 8*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_70():
    f = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**2
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)*elliptic_f(I*e + I*f*x, b/a)/(3*b*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - I*(a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(3*b*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_71():
    f = sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_72():
    f = sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**2
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/f - sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_73():
    f = sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**4
    F = -sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)*csch(e + f*x)**2/(3*f) - b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*tanh(e + f*x)/(3*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*coth(e + f*x)/(3*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_74():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)**3
    F = -(a - b)*(a + 5*b)*sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(16*b*f) - (a + 5*b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*cosh(e + f*x)/(24*b*f) + (a + b*cosh(e + f*x)**2 - b)**(sympy.S(5)/2)*cosh(e + f*x)/(6*b*f) - (a - b)**2*(a + 5*b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(16*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_75():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)
    F = (3*a - 3*b)*sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(8*f) + (a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*cosh(e + f*x)/(4*f) + 3*(a - b)**2*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(8*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_76():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/f + sqrt(b)*(3*a - b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*f) + b*sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_77():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**3
    F = sqrt(a)*(a - 3*b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*f) - a*sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(2*f) + b**(sympy.S(3)/2)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_78():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**5
    F = (3*a - 3*b)*sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(8*f) - (a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*coth(e + f*x)*csch(e + f*x)**3/(4*f) - 3*(a - b)**2*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_79():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**7
    F = -(a - b)*(5*a + b)*sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(16*a*f) + (5*a + b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)*coth(e + f*x)*csch(e + f*x)**3/(24*a*f) - (a + b*cosh(e + f*x)**2 - b)**(sympy.S(5)/2)*coth(e + f*x)*csch(e + f*x)**5/(6*a*f) + (a - b)**2*(5*a + b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(16*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_80():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)**4
    F = b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**5*cosh(e + f*x)/(7*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a - 6*b)*sinh(e + f*x)**3*cosh(e + f*x)/(35*f) + sqrt(a + b*sinh(e + f*x)**2)*(a**2 - 11*a*b + 8*b**2)*sinh(e + f*x)*cosh(e + f*x)/(35*b*f) - sqrt(a + b*sinh(e + f*x)**2)*(a**2 - 11*a*b + 8*b**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(35*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*(a**2 + 4*a*b - 4*b**2)*tanh(e + f*x)/(35*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*(a**2 + 4*a*b - 4*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(35*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_81():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)**2
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)*(3*a - 4*b)*elliptic_f(I*e + I*f*x, b/a)/(15*b*f*sqrt(a + b*sinh(e + f*x)**2)) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)*cosh(e + f*x)/(5*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*sinh(e + f*x)*cosh(e + f*x)/(15*f) - I*sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 - 13*a*b + 8*b**2)*elliptic_e(I*e + I*f*x, b/a)/(15*b*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_82():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)*elliptic_f(I*e + I*f*x, b/a)/(3*f*sqrt(a + b*sinh(e + f*x)**2)) + b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_83():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**2
    F = -a*sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/f + 2*b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (a + b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f - (a + b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_84():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**4
    F = -a*sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)*csch(e + f*x)**2/(3*f) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*tanh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*coth(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - b*(a - 3*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_85():
    f = (a + b*sinh(c + d*x)**2)**(sympy.S(5)/2)
    F = 4*I*a*sqrt(1 + b*sinh(c + d*x)**2/a)*(a - b)*(2*a - b)*elliptic_f(I*c + I*d*x, b/a)/(15*d*sqrt(a + b*sinh(c + d*x)**2)) + b*(a + b*sinh(c + d*x)**2)**(sympy.S(3)/2)*sinh(c + d*x)*cosh(c + d*x)/(5*d) + b*sqrt(a + b*sinh(c + d*x)**2)*(8*a - 4*b)*sinh(c + d*x)*cosh(c + d*x)/(15*d) - I*sqrt(a + b*sinh(c + d*x)**2)*(23*a**2 - 23*a*b + 8*b**2)*elliptic_e(I*c + I*d*x, b/a)/(15*d*sqrt(1 + b*sinh(c + d*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_86():
    f = sqrt(sinh(x)**2 + 1)
    F = sqrt(cosh(x)**2)*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_87():
    f = sqrt(-sinh(x)**2 - 1)
    F = sqrt(-cosh(x)**2)*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_88():
    f = sqrt(1 - sinh(x)**2)
    F = -I*elliptic_e(I*x, -1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_89():
    f = sqrt(sinh(x)**2 - 1)
    F = -I*sqrt(sinh(x)**2 - 1)*elliptic_e(I*x, -1)/sqrt(1 - sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_90():
    f = sqrt(a + b*sinh(x)**2)
    F = -I*sqrt(a + b*sinh(x)**2)*elliptic_e(I*x, b/a)/sqrt(1 + b*sinh(x)**2/a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_91():
    f = (sinh(x)**2 + 1)**(sympy.S(3)/2)
    F = (cosh(x)**2)**(sympy.S(3)/2)*tanh(x)/3 + 2*sqrt(cosh(x)**2)*tanh(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_92():
    f = (-sinh(x)**2 - 1)**(sympy.S(3)/2)
    F = (-cosh(x)**2)**(sympy.S(3)/2)*tanh(x)/3 - 2*sqrt(-cosh(x)**2)*tanh(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_93():
    f = (1 - sinh(x)**2)**(sympy.S(3)/2)
    F = -sqrt(1 - sinh(x)**2)*sinh(x)*cosh(x)/3 - 2*I*elliptic_e(I*x, -1) + 2*I*elliptic_f(I*x, -1)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_94():
    f = (sinh(x)**2 - 1)**(sympy.S(3)/2)
    F = 2*I*sqrt(1 - sinh(x)**2)*elliptic_f(I*x, -1)/(3*sqrt(sinh(x)**2 - 1)) + sqrt(sinh(x)**2 - 1)*sinh(x)*cosh(x)/3 + 2*I*sqrt(sinh(x)**2 - 1)*elliptic_e(I*x, -1)/sqrt(1 - sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_95():
    f = (a + b*sinh(x)**2)**(sympy.S(3)/2)
    F = I*a*sqrt(1 + b*sinh(x)**2/a)*(a - b)*elliptic_f(I*x, b/a)/(3*sqrt(a + b*sinh(x)**2)) + b*sqrt(a + b*sinh(x)**2)*sinh(x)*cosh(x)/3 - 2*I*sqrt(a + b*sinh(x)**2)*(2*a - b)*elliptic_e(I*x, b/a)/(3*sqrt(1 + b*sinh(x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_96():
    f = sinh(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*cosh(e + f*x)**2 - b)*cosh(e + f*x)/(2*b*f) - (a + b)*atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_97():
    f = sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_98():
    f = csch(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = -atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_99():
    f = csch(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*cosh(e + f*x)**2 - b)*coth(e + f*x)*csch(e + f*x)/(2*a*f) + (a + b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_100():
    f = sinh(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*b*f) - sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*tanh(e + f*x)/(3*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_101():
    f = sinh(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(b*f*sqrt(a + b*sinh(e + f*x)**2)) - I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(b*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_102():
    f = 1/sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_103():
    f = csch(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(a*f) - sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/(a*f) - sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_104():
    f = csch(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)*csch(e + f*x)**2/(3*a*f) - b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*tanh(e + f*x)/(3*a**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*coth(e + f*x)/(3*a**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_105():
    f = sinh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*cosh(e + f*x)/(b*f*(a - b)*sqrt(a + b*cosh(e + f*x)**2 - b)) + atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_106():
    f = sinh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = cosh(e + f*x)/(f*(a - b)*sqrt(a + b*cosh(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_107():
    f = csch(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*cosh(e + f*x)/(a*f*(a - b)*sqrt(a + b*cosh(e + f*x)**2 - b)) - atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_108():
    f = csch(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -coth(e + f*x)*csch(e + f*x)/(2*a*f*sqrt(a + b*cosh(e + f*x)**2 - b)) - b*(a - 3*b)*cosh(e + f*x)/(2*a**2*f*(a - b)*sqrt(a + b*cosh(e + f*x)**2 - b)) + (a + 3*b)*atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_109():
    f = sinh(e + f*x)**6/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sinh(e + f*x)**3*cosh(e + f*x)/(b*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(4*a - b)*sinh(e + f*x)*cosh(e + f*x)/(b**2*f*(3*a - 3*b)) - sqrt(a + b*sinh(e + f*x)**2)*(4*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b)) - sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 3*a*b - 2*b**2)*tanh(e + f*x)/(b**3*f*(3*a - 3*b)) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 3*a*b - 2*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(b**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_110():
    f = sinh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sinh(e + f*x)*cosh(e + f*x)/(b*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*tanh(e + f*x)/(b**2*f*(a - b)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_111():
    f = sinh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = sinh(e + f*x)*cosh(e + f*x)/(f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(b*f*sqrt(a + b*sinh(e + f*x)**2)) + I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(b*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_112():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(a*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_113():
    f = csch(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*coth(e + f*x)/(a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) + (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(a**2*f*(a - b)) - (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/(a**2*f*(a - b)) - (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_114():
    f = sinh(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*sinh(e + f*x)**2*cosh(e + f*x)/(b*f*(3*a - 3*b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)) - a*(3*a - 5*b)*cosh(e + f*x)/(3*b**2*f*(a - b)**2*sqrt(a + b*cosh(e + f*x)**2 - b)) + atanh(sqrt(b)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_115():
    f = sinh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sinh(e + f*x)**2*cosh(e + f*x)/(f*(3*a - 3*b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)) - 2*cosh(e + f*x)/(3*f*(a - b)**2*sqrt(a + b*cosh(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_116():
    f = sinh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = cosh(e + f*x)/(f*(3*a - 3*b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)) + 2*cosh(e + f*x)/(3*f*(a - b)**2*sqrt(a + b*cosh(e + f*x)**2 - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_117():
    f = csch(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*cosh(e + f*x)/(3*a*f*(a - b)*(a + b*cosh(e + f*x)**2 - b)**(sympy.S(3)/2)) - b*(5*a - 3*b)*cosh(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*cosh(e + f*x)**2 - b)) - atanh(sqrt(a)*cosh(e + f*x)/sqrt(a + b*cosh(e + f*x)**2 - b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_118():
    f = sinh(e + f*x)**6/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*sinh(e + f*x)**3*cosh(e + f*x)/(b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - 2*a*(2*a - 3*b)*sinh(e + f*x)*cosh(e + f*x)/(3*b**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(4*a - 6*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 13*a*b + 3*b**2)*tanh(e + f*x)/(3*b**3*f*(a - b)**2) - sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 13*a*b + 3*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_119():
    f = sinh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = 2*sqrt(a)*(a - 2*b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*b**(sympy.S(3)/2)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - a*sinh(e + f*x)*cosh(e + f*x)/(b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - (a - 3*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_120():
    f = sinh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sinh(e + f*x)*cosh(e + f*x)/(f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(3*b*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) + (a + b)*sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) + I*(a + b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(3*a*b*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_121():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(3*a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - b*(4*a - 2*b)*sinh(e + f*x)*cosh(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*a**2*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_122():
    f = csch(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*coth(e + f*x)/(3*a*f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - b*(6*a - 4*b)*coth(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - b*sqrt(a + b*sinh(e + f*x)**2)*(6*a - 4*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2) + sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 - 13*a*b + 8*b**2)*tanh(e + f*x)/(3*a**3*f*(a - b)**2) - sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 - 13*a*b + 8*b**2)*coth(e + f*x)/(3*a**3*f*(a - b)**2) - sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 - 13*a*b + 8*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_123():
    f = 1/sqrt(sinh(x)**2 + 1)
    F = cosh(x)*atan(sinh(x))/sqrt(cosh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_124():
    f = 1/sqrt(1 - sinh(x)**2)
    F = -I*elliptic_f(I*x, -1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_125():
    f = 1/sqrt(sinh(x)**2 - 1)
    F = -I*sqrt(1 - sinh(x)**2)*elliptic_f(I*x, -1)/sqrt(sinh(x)**2 - 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_126():
    f = 1/sqrt(-sinh(x)**2 - 1)
    F = cosh(x)*atan(sinh(x))/sqrt(-cosh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_127():
    f = 1/sqrt(a + b*sinh(x)**2)
    F = -I*sqrt(1 + b*sinh(x)**2/a)*elliptic_f(I*x, b/a)/sqrt(a + b*sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_128():
    f = (d*sinh(e + f*x))**m*(a + b*sinh(e + f*x)**2)**p
    F = d*(-sinh(e + f*x)**2)**(sympy.S.Half - m/2)*(d*sinh(e + f*x))**(m - 1)*(a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*appellf1(sympy.S.Half, -p, sympy.S.Half - m/2, sympy.S(3)/2, -b*cosh(e + f*x)**2/(a - b), cosh(e + f*x)**2)/(f*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_129():
    f = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)**5
    F = (a + b*cosh(e + f*x)**2 - b)**(p + 1)*sinh(e + f*x)**2*cosh(e + f*x)/(b*f*(2*p + 5)) - (3*a + 2*b*(p + 2))*(a + b*cosh(e + f*x)**2 - b)**(p + 1)*cosh(e + f*x)/(b**2*f*(2*p + 3)*(2*p + 5)) + (a + b*cosh(e + f*x)**2 - b)**p*(3*a**2 + 4*a*b*(p + 1) + 4*b**2*(p**2 + 3*p + 2))*cosh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*cosh(e + f*x)**2/(a - b))/(b**2*f*(2*p + 3)*(2*p + 5)*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_130():
    f = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)**3
    F = -(a + 2*b*(p + 1))*(a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*cosh(e + f*x)**2/(a - b))/(b*f*(2*p + 3)*(b*cosh(e + f*x)**2/(a - b) + 1)**p) + (a + b*cosh(e + f*x)**2 - b)**(p + 1)*cosh(e + f*x)/(b*f*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_131():
    f = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)
    F = (a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*cosh(e + f*x)**2/(a - b))/(f*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_132():
    f = (a + b*sinh(e + f*x)**2)**p*csch(e + f*x)
    F = -(a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, cosh(e + f*x)**2, -b*cosh(e + f*x)**2/(a - b))/(f*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_133():
    f = (a + b*sinh(e + f*x)**2)**p*csch(e + f*x)**3
    F = (a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, cosh(e + f*x)**2, -b*cosh(e + f*x)**2/(a - b))/(f*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_134():
    f = (a + b*sinh(e + f*x)**2)**p*csch(e + f*x)**5
    F = -(a + b*cosh(e + f*x)**2 - b)**p*cosh(e + f*x)*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, cosh(e + f*x)**2, -b*cosh(e + f*x)**2/(a - b))/(f*(b*cosh(e + f*x)**2/(a - b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_135():
    f = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)**4
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*sinh(e + f*x)**4*tanh(e + f*x)*appellf1(sympy.S(5)/2, sympy.S.Half, -p, sympy.S(7)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(5*f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_136():
    f = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)**2
    F = (a + b*sinh(e + f*x)**2)**p*(sech(e + f*x)**2)**p*tanh(e + f*x)**3*appellf1(sympy.S(3)/2, -p, p + 2, sympy.S(5)/2, (a - b)*tanh(e + f*x)**2/a, tanh(e + f*x)**2)/(3*f*(1 - (a - b)*tanh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_137():
    f = (a + b*sinh(e + f*x)**2)**p*csch(e + f*x)**2
    F = -(a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*appellf1(sympy.S(-1)/2, sympy.S.Half, -p, sympy.S.Half, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)*csch(e + f*x)*sech(e + f*x)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_138():
    f = (a + b*sinh(e + f*x)**2)**p*csch(e + f*x)**4
    F = -(a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*appellf1(sympy.S(-3)/2, sympy.S.Half, -p, sympy.S(-1)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)*csch(e + f*x)**3*sech(e + f*x)/(3*f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_139():
    f = (a + b*sinh(c + d*x)**3)*sinh(c + d*x)**4
    F = 3*a*x/8 + a*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 3*a*sinh(c + d*x)*cosh(c + d*x)/(8*d) + b*cosh(c + d*x)**7/(7*d) - 3*b*cosh(c + d*x)**5/(5*d) + b*cosh(c + d*x)**3/d - b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_140():
    f = (a + b*sinh(c + d*x)**3)*sinh(c + d*x)**3
    F = a*cosh(c + d*x)**3/(3*d) - a*cosh(c + d*x)/d - 5*b*x/16 + b*sinh(c + d*x)**5*cosh(c + d*x)/(6*d) - 5*b*sinh(c + d*x)**3*cosh(c + d*x)/(24*d) + 5*b*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_141():
    f = (a + b*sinh(c + d*x)**3)*sinh(c + d*x)**2
    F = -a*x/2 + a*sinh(c + d*x)*cosh(c + d*x)/(2*d) + b*cosh(c + d*x)**5/(5*d) - 2*b*cosh(c + d*x)**3/(3*d) + b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_142():
    f = (a + b*sinh(c + d*x)**3)*sinh(c + d*x)
    F = a*cosh(c + d*x)/d + 3*b*x/8 + b*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 3*b*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_143():
    f = a + b*sinh(c + d*x)**3
    F = a*x + b*cosh(c + d*x)**3/(3*d) - b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_144():
    f = (a + b*sinh(c + d*x)**3)*csch(c + d*x)
    F = -a*atanh(cosh(c + d*x))/d - b*x/2 + b*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_145():
    f = (a + b*sinh(c + d*x)**3)*csch(c + d*x)**2
    F = -a*coth(c + d*x)/d + b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_146():
    f = (a + b*sinh(c + d*x)**3)*csch(c + d*x)**3
    F = -a*coth(c + d*x)*csch(c + d*x)/(2*d) + a*atanh(cosh(c + d*x))/(2*d) + b*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_147():
    f = (a + b*sinh(c + d*x)**3)*csch(c + d*x)**4
    F = -a*coth(c + d*x)**3/(3*d) + a*coth(c + d*x)/d - b*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_148():
    f = (a + b*sinh(c + d*x)**3)**2*sinh(c + d*x)**3
    F = a**2*cosh(c + d*x)**3/(3*d) - a**2*cosh(c + d*x)/d - 5*a*b*x/8 + a*b*sinh(c + d*x)**5*cosh(c + d*x)/(3*d) - 5*a*b*sinh(c + d*x)**3*cosh(c + d*x)/(12*d) + 5*a*b*sinh(c + d*x)*cosh(c + d*x)/(8*d) + b**2*cosh(c + d*x)**9/(9*d) - 4*b**2*cosh(c + d*x)**7/(7*d) + 6*b**2*cosh(c + d*x)**5/(5*d) - 4*b**2*cosh(c + d*x)**3/(3*d) + b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_149():
    f = (a + b*sinh(c + d*x)**3)**2*sinh(c + d*x)**2
    F = -a**2*x/2 + a**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + 2*a*b*cosh(c + d*x)**5/(5*d) - 4*a*b*cosh(c + d*x)**3/(3*d) + 2*a*b*cosh(c + d*x)/d + 35*b**2*x/128 + b**2*sinh(c + d*x)**7*cosh(c + d*x)/(8*d) - 7*b**2*sinh(c + d*x)**5*cosh(c + d*x)/(48*d) + 35*b**2*sinh(c + d*x)**3*cosh(c + d*x)/(192*d) - 35*b**2*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_150():
    f = (a + b*sinh(c + d*x)**3)**2*sinh(c + d*x)
    F = a**2*cosh(c + d*x)/d + 3*a*b*x/4 + a*b*sinh(c + d*x)**3*cosh(c + d*x)/(2*d) - 3*a*b*sinh(c + d*x)*cosh(c + d*x)/(4*d) + b**2*cosh(c + d*x)**7/(7*d) - 3*b**2*cosh(c + d*x)**5/(5*d) + b**2*cosh(c + d*x)**3/d - b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_151():
    f = (a + b*sinh(c + d*x)**3)**2
    F = a**2*x + 2*a*b*cosh(c + d*x)**3/(3*d) - 2*a*b*cosh(c + d*x)/d - 5*b**2*x/16 + b**2*sinh(c + d*x)**5*cosh(c + d*x)/(6*d) - 5*b**2*sinh(c + d*x)**3*cosh(c + d*x)/(24*d) + 5*b**2*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_152():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)
    F = -a**2*atanh(cosh(c + d*x))/d - a*b*x + a*b*sinh(c + d*x)*cosh(c + d*x)/d + b**2*cosh(c + d*x)**5/(5*d) - 2*b**2*cosh(c + d*x)**3/(3*d) + b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_153():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**2
    F = -a**2*coth(c + d*x)/d + 2*a*b*cosh(c + d*x)/d + 3*b**2*x/8 + b**2*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 3*b**2*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_154():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**3
    F = -a**2*coth(c + d*x)*csch(c + d*x)/(2*d) + a**2*atanh(cosh(c + d*x))/(2*d) + 2*a*b*x + b**2*cosh(c + d*x)**3/(3*d) - b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_155():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) + a**2*coth(c + d*x)/d - 2*a*b*atanh(cosh(c + d*x))/d - b**2*x/2 + b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_156():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**5
    F = -a**2*coth(c + d*x)*csch(c + d*x)**3/(4*d) + 3*a**2*coth(c + d*x)*csch(c + d*x)/(8*d) - 3*a**2*atanh(cosh(c + d*x))/(8*d) - 2*a*b*coth(c + d*x)/d + b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_157():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**6
    F = -a**2*coth(c + d*x)**5/(5*d) + 2*a**2*coth(c + d*x)**3/(3*d) - a**2*coth(c + d*x)/d - a*b*coth(c + d*x)*csch(c + d*x)/d + a*b*atanh(cosh(c + d*x))/d + b**2*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_158():
    f = (a + b*sinh(c + d*x)**3)**2*csch(c + d*x)**7
    F = -a**2*coth(c + d*x)*csch(c + d*x)**5/(6*d) + 5*a**2*coth(c + d*x)*csch(c + d*x)**3/(24*d) - 5*a**2*coth(c + d*x)*csch(c + d*x)/(16*d) + 5*a**2*atanh(cosh(c + d*x))/(16*d) - 2*a*b*coth(c + d*x)**3/(3*d) + 2*a*b*coth(c + d*x)/d - b**2*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_159():
    f = (a + b*sinh(c + d*x)**3)**3*sinh(c + d*x)**2
    F = -a**3*x/2 + a**3*sinh(c + d*x)*cosh(c + d*x)/(2*d) + 3*a**2*b*cosh(c + d*x)**5/(5*d) - 2*a**2*b*cosh(c + d*x)**3/d + 3*a**2*b*cosh(c + d*x)/d + 105*a*b**2*x/128 + 3*a*b**2*sinh(c + d*x)**7*cosh(c + d*x)/(8*d) - 7*a*b**2*sinh(c + d*x)**5*cosh(c + d*x)/(16*d) + 35*a*b**2*sinh(c + d*x)**3*cosh(c + d*x)/(64*d) - 105*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(128*d) + b**3*cosh(c + d*x)**11/(11*d) - 5*b**3*cosh(c + d*x)**9/(9*d) + 10*b**3*cosh(c + d*x)**7/(7*d) - 2*b**3*cosh(c + d*x)**5/d + 5*b**3*cosh(c + d*x)**3/(3*d) - b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_160():
    f = (a + b*sinh(c + d*x)**3)**3*sinh(c + d*x)
    F = a**3*cosh(c + d*x)/d + 9*a**2*b*x/8 + 3*a**2*b*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 9*a**2*b*sinh(c + d*x)*cosh(c + d*x)/(8*d) + 3*a*b**2*cosh(c + d*x)**7/(7*d) - 9*a*b**2*cosh(c + d*x)**5/(5*d) + 3*a*b**2*cosh(c + d*x)**3/d - 3*a*b**2*cosh(c + d*x)/d - 63*b**3*x/256 + b**3*sinh(c + d*x)**9*cosh(c + d*x)/(10*d) - 9*b**3*sinh(c + d*x)**7*cosh(c + d*x)/(80*d) + 21*b**3*sinh(c + d*x)**5*cosh(c + d*x)/(160*d) - 21*b**3*sinh(c + d*x)**3*cosh(c + d*x)/(128*d) + 63*b**3*sinh(c + d*x)*cosh(c + d*x)/(256*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_161():
    f = (a + b*sinh(c + d*x)**3)**3
    F = a**3*x + a**2*b*cosh(c + d*x)**3/d - 3*a**2*b*cosh(c + d*x)/d - 15*a*b**2*x/16 + a*b**2*sinh(c + d*x)**5*cosh(c + d*x)/(2*d) - 5*a*b**2*sinh(c + d*x)**3*cosh(c + d*x)/(8*d) + 15*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(16*d) + b**3*cosh(c + d*x)**9/(9*d) - 4*b**3*cosh(c + d*x)**7/(7*d) + 6*b**3*cosh(c + d*x)**5/(5*d) - 4*b**3*cosh(c + d*x)**3/(3*d) + b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_162():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)
    F = -a**3*atanh(cosh(c + d*x))/d - 3*a**2*b*x/2 + 3*a**2*b*sinh(c + d*x)*cosh(c + d*x)/(2*d) + 3*a*b**2*cosh(c + d*x)**5/(5*d) - 2*a*b**2*cosh(c + d*x)**3/d + 3*a*b**2*cosh(c + d*x)/d + 35*b**3*x/128 + b**3*sinh(c + d*x)**7*cosh(c + d*x)/(8*d) - 7*b**3*sinh(c + d*x)**5*cosh(c + d*x)/(48*d) + 35*b**3*sinh(c + d*x)**3*cosh(c + d*x)/(192*d) - 35*b**3*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_163():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**2
    F = -a**3*coth(c + d*x)/d + 3*a**2*b*cosh(c + d*x)/d + 9*a*b**2*x/8 + 3*a*b**2*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 9*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(8*d) + b**3*cosh(c + d*x)**7/(7*d) - 3*b**3*cosh(c + d*x)**5/(5*d) + b**3*cosh(c + d*x)**3/d - b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_164():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**3
    F = -a**3*coth(c + d*x)*csch(c + d*x)/(2*d) + a**3*atanh(cosh(c + d*x))/(2*d) + 3*a**2*b*x + a*b**2*cosh(c + d*x)**3/d - 3*a*b**2*cosh(c + d*x)/d - 5*b**3*x/16 + b**3*sinh(c + d*x)**5*cosh(c + d*x)/(6*d) - 5*b**3*sinh(c + d*x)**3*cosh(c + d*x)/(24*d) + 5*b**3*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_165():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**4
    F = -a**3*coth(c + d*x)**3/(3*d) + a**3*coth(c + d*x)/d - 3*a**2*b*atanh(cosh(c + d*x))/d - 3*a*b**2*x/2 + 3*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + b**3*cosh(c + d*x)**5/(5*d) - 2*b**3*cosh(c + d*x)**3/(3*d) + b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_166():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**5
    F = -a**3*coth(c + d*x)*csch(c + d*x)**3/(4*d) + 3*a**3*coth(c + d*x)*csch(c + d*x)/(8*d) - 3*a**3*atanh(cosh(c + d*x))/(8*d) - 3*a**2*b*coth(c + d*x)/d + 3*a*b**2*cosh(c + d*x)/d + 3*b**3*x/8 + b**3*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 3*b**3*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_167():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**6
    F = -a**3*coth(c + d*x)**5/(5*d) + 2*a**3*coth(c + d*x)**3/(3*d) - a**3*coth(c + d*x)/d - 3*a**2*b*coth(c + d*x)*csch(c + d*x)/(2*d) + 3*a**2*b*atanh(cosh(c + d*x))/(2*d) + 3*a*b**2*x + b**3*cosh(c + d*x)**3/(3*d) - b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_168():
    f = (a + b*sinh(c + d*x)**3)**3*csch(c + d*x)**7
    F = -a**3*coth(c + d*x)*csch(c + d*x)**5/(6*d) + 5*a**3*coth(c + d*x)*csch(c + d*x)**3/(24*d) - 5*a**3*coth(c + d*x)*csch(c + d*x)/(16*d) + 5*a**3*atanh(cosh(c + d*x))/(16*d) - a**2*b*coth(c + d*x)**3/d + 3*a**2*b*coth(c + d*x)/d - 3*a*b**2*atanh(cosh(c + d*x))/d - b**3*x/2 + b**3*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_169():
    f = sinh(c + d*x)**6/(a + b*sinh(c + d*x)**3)
    F = -2*(-1)**(sympy.S(2)/3)*a**(sympy.S(4)/3)*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**2*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(1)/6)*a**(sympy.S(4)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**2*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*a**(sympy.S(4)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**2*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - a*x/b**2 + cosh(c + d*x)**3/(3*b*d) - cosh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_170():
    f = sinh(c + d*x)**5/(a + b*sinh(c + d*x)**3)
    F = 2*a*atan((-1)**(sympy.S(5)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*I*a*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*a*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - x/(2*b) + sinh(c + d*x)*cosh(c + d*x)/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_171():
    f = sinh(c + d*x)**4/(a + b*sinh(c + d*x)**3)
    F = -2*a**(sympy.S(2)/3)*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(5)/6)*a**(sympy.S(2)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*a**(sympy.S(2)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + cosh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_172():
    f = sinh(c + d*x)**3/(a + b*sinh(c + d*x)**3)
    F = 2*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*b*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(1)/6)*a**(sympy.S(1)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*a**(sympy.S(1)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_173():
    f = sinh(c + d*x)**2/(a + b*sinh(c + d*x)**3)
    F = -2*atan((-1)**(sympy.S(5)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*I*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*b**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_174():
    f = sinh(c + d*x)/(a + b*sinh(c + d*x)**3)
    F = 2*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(5)/6)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_175():
    f = 1/(a + b*sinh(c + d*x)**3)
    F = -2*(-1)**(sympy.S(2)/3)*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(1)/6)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**(sympy.S(2)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_176():
    f = csch(c + d*x)/(a + b*sinh(c + d*x)**3)
    F = 2*b**(sympy.S(1)/3)*atan((-1)**(sympy.S(5)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*I*b**(sympy.S(1)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*b**(sympy.S(1)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) - atanh(cosh(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_177():
    f = csch(c + d*x)**2/(a + b*sinh(c + d*x)**3)
    F = -coth(c + d*x)/(a*d) - 2*b**(sympy.S(2)/3)*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) + 2*(-1)**(sympy.S(5)/6)*b**(sympy.S(2)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*b**(sympy.S(2)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**(sympy.S(4)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_178():
    f = csch(c + d*x)**3/(a + b*sinh(c + d*x)**3)
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d) + atanh(cosh(c + d*x))/(2*a*d) + 2*(-1)**(sympy.S(2)/3)*b*atan((-1)**(sympy.S(1)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3))) - 2*(-1)**(sympy.S(1)/6)*b*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) + 2*b*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**(sympy.S(5)/3)*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_179():
    f = csch(c + d*x)**4/(a + b*sinh(c + d*x)**3)
    F = -coth(c + d*x)**3/(3*a*d) + coth(c + d*x)/(a*d) - 2*b**(sympy.S(4)/3)*atan((-1)**(sympy.S(5)/6)*(I*a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/6)*b**(sympy.S(1)/3))/sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**2*d*sqrt(-(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*I*b**(sympy.S(4)/3)*atanh((-1)**(sympy.S(1)/6)*(a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))/sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3)))/(3*a**2*d*sqrt((-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3) - b**(sympy.S(2)/3))) - 2*b**(sympy.S(4)/3)*atanh((-a**(sympy.S(1)/3)*tanh(c/2 + d*x/2) + b**(sympy.S(1)/3))/sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3)))/(3*a**2*d*sqrt(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))) + b*atanh(cosh(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_180():
    f = 1/(sinh(x)**3 + 1)
    F = -(-1)**(sympy.S(1)/6)*log(-(-1)**(sympy.S(1)/6)*tanh(x/2) + 1 + (-1)**(sympy.S(5)/6))/3 + (-1)**(sympy.S(1)/6)*log((-1)**(sympy.S(1)/3)*tanh(x/2) + 1 + (-1)**(sympy.S(1)/6))/3 - 2*(-1)**(sympy.S(1)/6)*atan(((-1)**(sympy.S(1)/6)*tanh(x/2) + I)/sqrt(1 - (-1)**(sympy.S(1)/3)))/(3*sqrt(1 - (-1)**(sympy.S(1)/3))) - sqrt(2)*atanh(sqrt(2)*(1 - tanh(x/2))/2)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_181():
    f = 1/(1 - sinh(x)**3)
    F = -(-1)**(sympy.S(5)/6)*log((-1)**(sympy.S(2)/3)*tanh(x/2) + 1 + (-1)**(sympy.S(5)/6))/3 + (-1)**(sympy.S(5)/6)*log((-1)**(sympy.S(5)/6)*tanh(x/2) + 1 + (-1)**(sympy.S(1)/6))/3 + sqrt(2)*atanh(sqrt(2)*(tanh(x/2) + 1)/2)/3 - 2*(-1)**(sympy.S(1)/3)*atanh((-(-1)**(sympy.S(1)/3)*tanh(x/2) + 1)/sqrt(1 + (-1)**(sympy.S(2)/3)))/(3*sqrt(1 + (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_182():
    f = (a + b*sinh(c + d*x)**4)*sinh(c + d*x)**4
    F = b*sinh(c + d*x)*cosh(c + d*x)**7/(8*d) - 25*b*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) + x*(3*a/8 + 35*b/128) + (48*a + 163*b)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) - (80*a + 93*b)*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_183():
    f = (a + b*sinh(c + d*x)**4)*sinh(c + d*x)**3
    F = b*cosh(c + d*x)**7/(7*d) - 3*b*cosh(c + d*x)**5/(5*d) - (a + b)*cosh(c + d*x)/d + (a + 3*b)*cosh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_184():
    f = (a + b*sinh(c + d*x)**4)*sinh(c + d*x)**2
    F = b*sinh(c + d*x)*cosh(c + d*x)**5/(6*d) - 13*b*sinh(c + d*x)*cosh(c + d*x)**3/(24*d) + x*(-a/2 - 5*b/16) + (8*a + 11*b)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_185():
    f = (a + b*sinh(c + d*x)**4)*sinh(c + d*x)
    F = b*cosh(c + d*x)**5/(5*d) - 2*b*cosh(c + d*x)**3/(3*d) + (a + b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_186():
    f = a + b*sinh(c + d*x)**4
    F = a*x + 3*b*x/8 + b*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - 3*b*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_187():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)
    F = -a*atanh(cosh(c + d*x))/d + b*cosh(c + d*x)**3/(3*d) - b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_188():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**2
    F = -a*coth(c + d*x)/d - b*x/2 + b*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_189():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**3
    F = -a*coth(c + d*x)*csch(c + d*x)/(2*d) + a*atanh(cosh(c + d*x))/(2*d) + b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_190():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**4
    F = -a*coth(c + d*x)**3/(3*d) + a*coth(c + d*x)/d + b*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_191():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**5
    F = -a*coth(c + d*x)*csch(c + d*x)**3/(4*d) + 3*a*coth(c + d*x)*csch(c + d*x)/(8*d) - (3*a + 8*b)*atanh(cosh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_192():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**6
    F = -a*coth(c + d*x)**5/(5*d) + 2*a*coth(c + d*x)**3/(3*d) - (a + b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_193():
    f = (a + b*sinh(c + d*x)**4)*csch(c + d*x)**7
    F = -a*coth(c + d*x)*csch(c + d*x)**5/(6*d) + 5*a*coth(c + d*x)*csch(c + d*x)**3/(24*d) - (5*a + 8*b)*coth(c + d*x)*csch(c + d*x)/(16*d) + (5*a + 8*b)*atanh(cosh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_194():
    f = (a + b*sinh(c + d*x)**4)**2*sinh(c + d*x)**3
    F = b**2*cosh(c + d*x)**11/(11*d) - 5*b**2*cosh(c + d*x)**9/(9*d) + 2*b*(a + 5*b)*cosh(c + d*x)**7/(7*d) - 2*b*(3*a + 5*b)*cosh(c + d*x)**5/(5*d) - (a + b)**2*cosh(c + d*x)/d + (a + b)*(a + 5*b)*cosh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_195():
    f = (a + b*sinh(c + d*x)**4)**2*sinh(c + d*x)**2
    F = b**2*sinh(c + d*x)*cosh(c + d*x)**9/(10*d) - 41*b**2*sinh(c + d*x)*cosh(c + d*x)**7/(80*d) + b*(160*a + 513*b)*sinh(c + d*x)*cosh(c + d*x)**5/(480*d) - b*(416*a + 447*b)*sinh(c + d*x)*cosh(c + d*x)**3/(384*d) + x*(-a**2/2 - 5*a*b/8 - 63*b**2/256) + (128*a**2 + 352*a*b + 193*b**2)*sinh(c + d*x)*cosh(c + d*x)/(256*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_196():
    f = (a + b*sinh(c + d*x)**4)**2*sinh(c + d*x)
    F = b**2*cosh(c + d*x)**9/(9*d) - 4*b**2*cosh(c + d*x)**7/(7*d) - 4*b*(a + b)*cosh(c + d*x)**3/(3*d) + 2*b*(a + 3*b)*cosh(c + d*x)**5/(5*d) + (a + b)**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_197():
    f = (a + b*sinh(c + d*x)**4)**2
    F = b**2*sinh(c + d*x)*cosh(c + d*x)**7/(8*d) - 25*b**2*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) + b*(96*a + 163*b)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) - b*(160*a + 93*b)*sinh(c + d*x)*cosh(c + d*x)/(128*d) + x*(a**2 + 3*a*b/4 + 35*b**2/128)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_198():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)
    F = -a**2*atanh(cosh(c + d*x))/d + b**2*cosh(c + d*x)**7/(7*d) - 3*b**2*cosh(c + d*x)**5/(5*d) - b*(2*a + b)*cosh(c + d*x)/d + b*(2*a + 3*b)*cosh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_199():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**2
    F = -a**2*coth(c + d*x)/d + b**2*sinh(c + d*x)*cosh(c + d*x)**5/(6*d) - 13*b**2*sinh(c + d*x)*cosh(c + d*x)**3/(24*d) - b*x*(16*a + 5*b)/16 + b*(16*a + 11*b)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_200():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**3
    F = -a**2*coth(c + d*x)*csch(c + d*x)/(2*d) + a**2*atanh(cosh(c + d*x))/(2*d) + b**2*cosh(c + d*x)**5/(5*d) - 2*b**2*cosh(c + d*x)**3/(3*d) + b*(2*a + b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_201():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) + a**2*coth(c + d*x)/d + b**2*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) - 5*b**2*sinh(c + d*x)*cosh(c + d*x)/(8*d) + b*x*(16*a + 3*b)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_202():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**5
    F = -a**2*coth(c + d*x)*csch(c + d*x)**3/(4*d) + 3*a**2*coth(c + d*x)*csch(c + d*x)/(8*d) - a*(3*a + 16*b)*atanh(cosh(c + d*x))/(8*d) + b**2*cosh(c + d*x)**3/(3*d) - b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_203():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**6
    F = -a**2*coth(c + d*x)**5/(5*d) + 2*a**2*coth(c + d*x)**3/(3*d) - a*(a + 2*b)*coth(c + d*x)/d - b**2*x/2 + b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_204():
    f = (a + b*sinh(c + d*x)**4)**2*csch(c + d*x)**7
    F = -a**2*coth(c + d*x)*csch(c + d*x)**5/(6*d) + 5*a**2*coth(c + d*x)*csch(c + d*x)**3/(24*d) - a*(5*a + 16*b)*coth(c + d*x)*csch(c + d*x)/(16*d) + a*(5*a + 16*b)*atanh(cosh(c + d*x))/(16*d) + b**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_205():
    f = (a + b*sinh(c + d*x)**4)**3*sinh(c + d*x)**5
    F = b**3*cosh(c + d*x)**17/(17*d) - 8*b**3*cosh(c + d*x)**15/(15*d) + b**2*(3*a + 28*b)*cosh(c + d*x)**13/(13*d) - 2*b**2*(9*a + 28*b)*cosh(c + d*x)**11/(11*d) - 4*b*(3*a**2 + 15*a*b + 14*b**2)*cosh(c + d*x)**7/(7*d) + b*(3*a**2 + 45*a*b + 70*b**2)*cosh(c + d*x)**9/(9*d) + (a + b)**3*cosh(c + d*x)/d - 2*(a + b)**2*(a + 4*b)*cosh(c + d*x)**3/(3*d) + (a + b)*(a**2 + 17*a*b + 28*b**2)*cosh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_206():
    f = (a + b*sinh(c + d*x)**4)**3*sinh(c + d*x)**3
    F = b**3*cosh(c + d*x)**15/(15*d) - 7*b**3*cosh(c + d*x)**13/(13*d) + 3*b**2*(a + 7*b)*cosh(c + d*x)**11/(11*d) - 5*b**2*(3*a + 7*b)*cosh(c + d*x)**9/(9*d) - 3*b*(a + b)*(3*a + 7*b)*cosh(c + d*x)**5/(5*d) + b*(3*a**2 + 30*a*b + 35*b**2)*cosh(c + d*x)**7/(7*d) - (a + b)**3*cosh(c + d*x)/d + (a + b)**2*(a + 7*b)*cosh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_207():
    f = (a + b*sinh(c + d*x)**4)**3*sinh(c + d*x)
    F = b**3*cosh(c + d*x)**13/(13*d) - 6*b**3*cosh(c + d*x)**11/(11*d) + b**2*(a + 5*b)*cosh(c + d*x)**9/(3*d) - 4*b**2*(3*a + 5*b)*cosh(c + d*x)**7/(7*d) - 2*b*(a + b)**2*cosh(c + d*x)**3/d + 3*b*(a + b)*(a + 5*b)*cosh(c + d*x)**5/(5*d) + (a + b)**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_208():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)
    F = -a**3*atanh(cosh(c + d*x))/d + b**3*cosh(c + d*x)**11/(11*d) - 5*b**3*cosh(c + d*x)**9/(9*d) + b**2*(3*a + 10*b)*cosh(c + d*x)**7/(7*d) - b**2*(9*a + 10*b)*cosh(c + d*x)**5/(5*d) - b*(3*a**2 + 3*a*b + b**2)*cosh(c + d*x)/d + b*(3*a**2 + 9*a*b + 5*b**2)*cosh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_209():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**3
    F = -a**3*coth(c + d*x)*csch(c + d*x)/(2*d) + a**3*atanh(cosh(c + d*x))/(2*d) + b**3*cosh(c + d*x)**9/(9*d) - 4*b**3*cosh(c + d*x)**7/(7*d) + 3*b**2*(a + 2*b)*cosh(c + d*x)**5/(5*d) - 2*b**2*(3*a + 2*b)*cosh(c + d*x)**3/(3*d) + b*(3*a**2 + 3*a*b + b**2)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_210():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**5
    F = -a**3*coth(c + d*x)*csch(c + d*x)**3/(4*d) + 3*a**3*coth(c + d*x)*csch(c + d*x)/(8*d) - 3*a**2*(a + 8*b)*atanh(cosh(c + d*x))/(8*d) + b**3*cosh(c + d*x)**7/(7*d) - 3*b**3*cosh(c + d*x)**5/(5*d) + b**2*(a + b)*cosh(c + d*x)**3/d - b**2*(3*a + b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_211():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**7
    F = -a**3*coth(c + d*x)*csch(c + d*x)**5/(6*d) + 5*a**3*coth(c + d*x)*csch(c + d*x)**3/(24*d) - a**2*(5*a + 24*b)*coth(c + d*x)*csch(c + d*x)/(16*d) + a**2*(5*a + 24*b)*atanh(cosh(c + d*x))/(16*d) + b**3*cosh(c + d*x)**5/(5*d) - 2*b**3*cosh(c + d*x)**3/(3*d) + b**2*(3*a + b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_212():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**9
    F = -a**3*coth(c + d*x)*csch(c + d*x)**7/(8*d) + 7*a**3*coth(c + d*x)*csch(c + d*x)**5/(48*d) - a**2*(35*a + 144*b)*coth(c + d*x)*csch(c + d*x)**3/(192*d) + a**2*(35*a + 144*b)*coth(c + d*x)*csch(c + d*x)/(128*d) - a*(35*a**2 + 144*a*b + 384*b**2)*atanh(cosh(c + d*x))/(128*d) + b**3*cosh(c + d*x)**3/(3*d) - b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_213():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**11
    F = -a**3*coth(c + d*x)*csch(c + d*x)**9/(10*d) + 9*a**3*coth(c + d*x)*csch(c + d*x)**7/(80*d) - a**2*(21*a + 80*b)*coth(c + d*x)*csch(c + d*x)**5/(160*d) + a**2*(21*a + 80*b)*coth(c + d*x)*csch(c + d*x)**3/(128*d) - 3*a*(21*a**2 + 80*a*b + 128*b**2)*coth(c + d*x)*csch(c + d*x)/(256*d) + 3*a*(21*a**2 + 80*a*b + 128*b**2)*atanh(cosh(c + d*x))/(256*d) + b**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_214():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**13
    F = -a**3*coth(c + d*x)*csch(c + d*x)**11/(12*d) + 11*a**3*coth(c + d*x)*csch(c + d*x)**9/(120*d) - 3*a**2*(11*a + 40*b)*coth(c + d*x)*csch(c + d*x)**7/(320*d) + 7*a**2*(11*a + 40*b)*coth(c + d*x)*csch(c + d*x)**5/(640*d) - a*(77*a**2 + 280*a*b + 384*b**2)*coth(c + d*x)*csch(c + d*x)**3/(512*d) + 3*a*(77*a**2 + 280*a*b + 384*b**2)*coth(c + d*x)*csch(c + d*x)/(1024*d) - (231*a**3 + 840*a**2*b + 1152*a*b**2 + 1024*b**3)*atanh(cosh(c + d*x))/(1024*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_215():
    f = (a + b*sinh(c + d*x)**4)**3*sinh(c + d*x)**2
    F = b**3*sinh(c + d*x)*cosh(c + d*x)**13/(14*d) - 85*b**3*sinh(c + d*x)*cosh(c + d*x)**11/(168*d) + b**2*(504*a + 2593*b)*sinh(c + d*x)*cosh(c + d*x)**9/(1680*d) - b**2*(6888*a + 11821*b)*sinh(c + d*x)*cosh(c + d*x)**7/(4480*d) + b*(1920*a**2 + 12312*a*b + 10579*b**2)*sinh(c + d*x)*cosh(c + d*x)**5/(3840*d) - b*(4992*a**2 + 10728*a*b + 5549*b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(3072*d) - x*(1024*a**3 + 1920*a**2*b + 1512*a*b**2 + 429*b**3)/2048 + (1024*a**3 + 4224*a**2*b + 4632*a*b**2 + 1619*b**3)*sinh(c + d*x)*cosh(c + d*x)/(2048*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_216():
    f = (a + b*sinh(c + d*x)**4)**3
    F = b**3*sinh(c + d*x)*cosh(c + d*x)**11/(12*d) - 61*b**3*sinh(c + d*x)*cosh(c + d*x)**9/(120*d) + 3*b**2*(40*a + 139*b)*sinh(c + d*x)*cosh(c + d*x)**7/(320*d) - b**2*(3000*a + 3481*b)*sinh(c + d*x)*cosh(c + d*x)**5/(1920*d) + b*(1152*a**2 + 3912*a*b + 2279*b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(1536*d) - b*(1920*a**2 + 2232*a*b + 793*b**2)*sinh(c + d*x)*cosh(c + d*x)/(1024*d) + x*(1024*a**3 + 1152*a**2*b + 840*a*b**2 + 231*b**3)/1024
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_217():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**2
    F = -a**3*coth(c + d*x)/d + b**3*sinh(c + d*x)*cosh(c + d*x)**9/(10*d) - 41*b**3*sinh(c + d*x)*cosh(c + d*x)**7/(80*d) + b**2*(80*a + 171*b)*sinh(c + d*x)*cosh(c + d*x)**5/(160*d) - b**2*(208*a + 149*b)*sinh(c + d*x)*cosh(c + d*x)**3/(128*d) - 3*b*x*(128*a**2 + 80*a*b + 21*b**2)/256 + b*(384*a**2 + 528*a*b + 193*b**2)*sinh(c + d*x)*cosh(c + d*x)/(256*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_218():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**4
    F = -a**3*coth(c + d*x)**3/(3*d) + a**3*coth(c + d*x)/d + b**3*sinh(c + d*x)*cosh(c + d*x)**7/(8*d) - 25*b**3*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) - 3*b**2*(80*a + 31*b)*sinh(c + d*x)*cosh(c + d*x)/(128*d) + b**2*(144*a + 163*b)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) + b*x*(384*a**2 + 144*a*b + 35*b**2)/128
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_219():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**6
    F = -a**3*coth(c + d*x)**5/(5*d) + 2*a**3*coth(c + d*x)**3/(3*d) - a**2*(a + 3*b)*coth(c + d*x)/d + b**3*sinh(c + d*x)*cosh(c + d*x)**5/(6*d) - 13*b**3*sinh(c + d*x)*cosh(c + d*x)**3/(24*d) - b**2*x*(24*a + 5*b)/16 + b**2*(24*a + 11*b)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_220():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**8
    F = -a**3*coth(c + d*x)**7/(7*d) + 3*a**3*coth(c + d*x)**5/(5*d) - a**2*(a + b)*coth(c + d*x)**3/d + a**2*(a + 3*b)*coth(c + d*x)/d + b**3*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) - 5*b**3*sinh(c + d*x)*cosh(c + d*x)/(8*d) + 3*b**2*x*(8*a + b)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_221():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**10
    F = -a**3*coth(c + d*x)**9/(9*d) + 4*a**3*coth(c + d*x)**7/(7*d) - 3*a**2*(2*a + b)*coth(c + d*x)**5/(5*d) + 2*a**2*(2*a + 3*b)*coth(c + d*x)**3/(3*d) - a*(a**2 + 3*a*b + 3*b**2)*coth(c + d*x)/d - b**3*x/2 + b**3*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_222():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**12
    F = -a**3*coth(c + d*x)**11/(11*d) + 5*a**3*coth(c + d*x)**9/(9*d) - a**2*(10*a + 3*b)*coth(c + d*x)**7/(7*d) + a**2*(10*a + 9*b)*coth(c + d*x)**5/(5*d) + a*(a**2 + 3*a*b + 3*b**2)*coth(c + d*x)/d - a*(5*a**2 + 9*a*b + 3*b**2)*coth(c + d*x)**3/(3*d) + b**3*x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_223():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**14
    F = -a**3*coth(c + d*x)**13/(13*d) + 6*a**3*coth(c + d*x)**11/(11*d) - a**2*(5*a + b)*coth(c + d*x)**9/(3*d) + 4*a**2*(5*a + 3*b)*coth(c + d*x)**7/(7*d) + 2*a*(a + b)**2*coth(c + d*x)**3/d - 3*a*(a + b)*(5*a + b)*coth(c + d*x)**5/(5*d) - (a + b)**3*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_224():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**16
    F = -a**3*coth(c + d*x)**15/(15*d) + 7*a**3*coth(c + d*x)**13/(13*d) - 3*a**2*(7*a + b)*coth(c + d*x)**11/(11*d) + 5*a**2*(7*a + 3*b)*coth(c + d*x)**9/(9*d) + 3*a*(a + b)*(7*a + 3*b)*coth(c + d*x)**5/(5*d) - a*(35*a**2 + 30*a*b + 3*b**2)*coth(c + d*x)**7/(7*d) + (a + b)**3*coth(c + d*x)/d - (a + b)**2*(7*a + b)*coth(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_225():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**18
    F = -a**3*coth(c + d*x)**17/(17*d) + 8*a**3*coth(c + d*x)**15/(15*d) - a**2*(28*a + 3*b)*coth(c + d*x)**13/(13*d) + 2*a**2*(28*a + 9*b)*coth(c + d*x)**11/(11*d) + 4*a*(14*a**2 + 15*a*b + 3*b**2)*coth(c + d*x)**7/(7*d) - a*(70*a**2 + 45*a*b + 3*b**2)*coth(c + d*x)**9/(9*d) - (a + b)**3*coth(c + d*x)/d + 2*(a + b)**2*(4*a + b)*coth(c + d*x)**3/(3*d) - (a + b)*(28*a**2 + 17*a*b + b**2)*coth(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_226():
    f = (a + b*sinh(c + d*x)**4)**3*csch(c + d*x)**20
    F = -a**3*coth(c + d*x)**19/(19*d) + 9*a**3*coth(c + d*x)**17/(17*d) + 21*a**2*(4*a + b)*coth(c + d*x)**13/(13*d) - a**2*(12*a + b)*coth(c + d*x)**15/(5*d) - 3*a*(42*a**2 + 21*a*b + b**2)*coth(c + d*x)**11/(11*d) + a*(42*a**2 + 35*a*b + 5*b**2)*coth(c + d*x)**9/(3*d) + (a + b)**3*coth(c + d*x)/d - (a + b)**2*(3*a + b)*coth(c + d*x)**3/d + (3*a + 3*b)*(12*a**2 + 9*a*b + b**2)*coth(c + d*x)**5/(5*d) - (84*a**3 + 105*a**2*b + 30*a*b**2 + b**3)*coth(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_227():
    f = sinh(c + d*x)**7/(a - b*sinh(c + d*x)**4)
    F = a*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(7)/4)*d*sqrt(sqrt(a) + sqrt(b))) - a*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(7)/4)*d*sqrt(sqrt(a) - sqrt(b))) - cosh(c + d*x)**3/(3*b*d) + cosh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_228():
    f = sinh(c + d*x)**5/(a - b*sinh(c + d*x)**4)
    F = sqrt(a)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(5)/4)*d*sqrt(sqrt(a) + sqrt(b))) + sqrt(a)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(5)/4)*d*sqrt(sqrt(a) - sqrt(b))) - cosh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_229():
    f = sinh(c + d*x)**3/(a - b*sinh(c + d*x)**4)
    F = atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*b**(sympy.S(3)/4)*d*sqrt(sqrt(a) + sqrt(b))) - atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*b**(sympy.S(3)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_230():
    f = sinh(c + d*x)/(a - b*sinh(c + d*x)**4)
    F = atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*sqrt(a)*b**(sympy.S(1)/4)*d*sqrt(sqrt(a) + sqrt(b))) + atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*sqrt(a)*b**(sympy.S(1)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_231():
    f = csch(c + d*x)/(a - b*sinh(c + d*x)**4)
    F = b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cosh(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_232():
    f = csch(c + d*x)**3/(a - b*sinh(c + d*x)**4)
    F = atanh(cosh(c + d*x))/(2*a*d) - 1/(4*a*d*(cosh(c + d*x) + 1)) + 1/(4*a*d*(1 - cosh(c + d*x))) + b**(sympy.S(3)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(3)/2)*d*sqrt(sqrt(a) + sqrt(b))) + b**(sympy.S(3)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(3)/2)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_233():
    f = sinh(c + d*x)**6/(a - b*sinh(c + d*x)**4)
    F = a**(sympy.S(3)/4)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/2)*d*sqrt(sqrt(a) + sqrt(b))) - a**(sympy.S(3)/4)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b**(sympy.S(3)/2)*d*sqrt(sqrt(a) - sqrt(b))) + x/(2*b) + 1/(4*b*d*(tanh(c + d*x) + 1)) - 1/(4*b*d*(1 - tanh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_234():
    f = sinh(c + d*x)**4/(a - b*sinh(c + d*x)**4)
    F = a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b*d*sqrt(sqrt(a) + sqrt(b))) + a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b*d*sqrt(sqrt(a) - sqrt(b))) - x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_235():
    f = sinh(c + d*x)**2/(a - b*sinh(c + d*x)**4)
    F = atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*sqrt(b)*d*sqrt(sqrt(a) + sqrt(b))) - atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*sqrt(b)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_236():
    f = 1/(a - b*sinh(c + d*x)**4)
    F = atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*sqrt(sqrt(a) + sqrt(b))) + atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_237():
    f = csch(c + d*x)**2/(a - b*sinh(c + d*x)**4)
    F = -coth(c + d*x)/(a*d) + sqrt(b)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*d*sqrt(sqrt(a) + sqrt(b))) - sqrt(b)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_238():
    f = csch(c + d*x)**4/(a - b*sinh(c + d*x)**4)
    F = -coth(c + d*x)**3/(3*a*d) + coth(c + d*x)/(a*d) + b*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*d*sqrt(sqrt(a) + sqrt(b))) + b*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(7)/4)*d*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_239():
    f = sinh(c + d*x)**9/(a - b*sinh(c + d*x)**4)**2
    F = -sqrt(a)*(5*sqrt(a) + 6*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*b**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - sqrt(a)*(5*sqrt(a) - 6*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*b**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) + a*(a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(b**2*d*(4*a - 4*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + cosh(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_240():
    f = sinh(c + d*x)**7/(a - b*sinh(c + d*x)**4)**2
    F = -a*(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(b*d*(4*a - 4*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) - (3*sqrt(a) + 4*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*b**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (3*sqrt(a) - 4*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*b**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_241():
    f = sinh(c + d*x)**5/(a - b*sinh(c + d*x)**4)**2
    F = (a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(b*d*(4*a - 4*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) - (sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*sqrt(a)*b**(sympy.S(5)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) - (sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*sqrt(a)*b**(sympy.S(5)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_242():
    f = sinh(c + d*x)**3/(a - b*sinh(c + d*x)**4)**2
    F = -(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(d*(4*a - 4*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*sqrt(a)*b**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*sqrt(a)*b**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_243():
    f = sinh(c + d*x)/(a - b*sinh(c + d*x)**4)**2
    F = (a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(4*a*d*(a - b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + (3*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(3)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (3*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(3)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_244():
    f = csch(c + d*x)/(a - b*sinh(c + d*x)**4)**2
    F = -b*(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(4*a*d*(a - b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**2*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cosh(c + d*x))/(a**2*d) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_245():
    f = sinh(c + d*x)**8/(a - b*sinh(c + d*x)**4)**2
    F = -a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) + sqrt(b))) - a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(2*b**2*d*sqrt(sqrt(a) - sqrt(b))) - a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + a**(sympy.S(1)/4)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) + tanh(c + d*x)**5/(4*b*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - tanh(c + d*x)/(b*d*(4*a - 4*b)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_246():
    f = sinh(c + d*x)**6/(a - b*sinh(c + d*x)**4)**2
    F = tanh(c + d*x)**3*sech(c + d*x)**2/(4*b*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) + tanh(c + d*x)/(b*d*(4*a - 4*b)) - (2*sqrt(a) + 3*sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (2*sqrt(a) - 3*sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(1)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_247():
    f = sinh(c + d*x)**4/(a - b*sinh(c + d*x)**4)**2
    F = tanh(c + d*x)**5/(4*a*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - tanh(c + d*x)/(4*a*d*(a - b)) - atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(3)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_248():
    f = sinh(c + d*x)**2/(a - b*sinh(c + d*x)**4)**2
    F = (a - (a + b)*tanh(c + d*x)**2)*tanh(c + d*x)/(4*a*d*(a - b)*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) + (2*sqrt(a) + sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - (2*sqrt(a) - sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_249():
    f = (a - b*sinh(c + d*x)**4)**(-2)
    F = -b*(1 - 2*tanh(c + d*x)**2)*tanh(c + d*x)/(4*a*d*(a - b)*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) + (4*sqrt(a) + 3*sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + (4*sqrt(a) - 3*sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_250():
    f = csch(c + d*x)**2/(a - b*sinh(c + d*x)**4)**2
    F = b*(a - (a + b)*tanh(c + d*x)**2)*tanh(c + d*x)/(4*a**2*d*(a - b)*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - coth(c + d*x)/(a**2*d) + sqrt(b)*(6*sqrt(a) + 5*sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) - sqrt(b)*(6*sqrt(a) - 5*sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(8*a**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_251():
    f = sinh(c + d*x)**9/(a - b*sinh(c + d*x)**4)**3
    F = a*(a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(b**2*d*(8*a - 8*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) - (9*a**2 - 11*a*b - 10*b**2 - b*(4*a - 10*b)*cosh(c + d*x)**2)*cosh(c + d*x)/(32*b**2*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + (14*sqrt(a)*sqrt(b) + 5*a + 12*b)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*sqrt(a)*b**(sympy.S(9)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-14*sqrt(a)*sqrt(b) + 5*a + 12*b)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*sqrt(a)*b**(sympy.S(9)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_252():
    f = sinh(c + d*x)**7/(a - b*sinh(c + d*x)**4)**3
    F = -a*(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(b*d*(8*a - 8*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) + (5*a - 17*b - (3*a - 9*b)*cosh(c + d*x)**2)*cosh(c + d*x)/(32*b*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) - (3*sqrt(a) + 6*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*sqrt(a)*b**(sympy.S(7)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (3*sqrt(a) - 6*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*sqrt(a)*b**(sympy.S(7)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_253():
    f = sinh(c + d*x)**5/(a - b*sinh(c + d*x)**4)**3
    F = (a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(b*d*(8*a - 8*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) - (a**2 - 11*a*b - 2*b**2 + 2*b*(2*a + b)*cosh(c + d*x)**2)*cosh(c + d*x)/(32*a*b*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) - (10*sqrt(a)*sqrt(b) + 3*a + 4*b)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(5)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (-10*sqrt(a)*sqrt(b) + 3*a + 4*b)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(5)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_254():
    f = sinh(c + d*x)**3/(a - b*sinh(c + d*x)**4)**3
    F = -(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(d*(8*a - 8*b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) - (11*a + b - (5*a + b)*cosh(c + d*x)**2)*cosh(c + d*x)/(32*a*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + (5*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(3)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (5*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(3)/2)*b**(sympy.S(3)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_255():
    f = sinh(c + d*x)/(a - b*sinh(c + d*x)**4)**3
    F = (a - b*cosh(c + d*x)**2 + b)*cosh(c + d*x)/(8*a*d*(a - b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) + (-b*(12*a - 6*b)*cosh(c + d*x)**2 + (a + 2*b)*(7*a - 3*b))*cosh(c + d*x)/(32*a**2*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + (30*sqrt(a)*sqrt(b) + 21*a + 12*b)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(5)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-30*sqrt(a)*sqrt(b) + 21*a + 12*b)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(5)/2)*b**(sympy.S(1)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_256():
    f = csch(c + d*x)/(a - b*sinh(c + d*x)**4)**3
    F = -b*(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(8*a*d*(a - b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)**2) - b*(2 - cosh(c + d*x)**2)*cosh(c + d*x)/(4*a**2*d*(a - b)*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) - b*(11*a + b - (5*a + b)*cosh(c + d*x)**2)*cosh(c + d*x)/(32*a**2*d*(a - b)**2*(a - b*cosh(c + d*x)**4 + 2*b*cosh(c + d*x)**2 - b)) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**3*d*sqrt(sqrt(a) + sqrt(b))) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**3*d*sqrt(sqrt(a) - sqrt(b))) - atanh(cosh(c + d*x))/(a**3*d) + b**(sympy.S(1)/4)*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(8*a**(sympy.S(5)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(3)/2)) + b**(sympy.S(1)/4)*(5*sqrt(a) + 2*sqrt(b))*atanh(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) + sqrt(b)))/(64*a**(sympy.S(5)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - b**(sympy.S(1)/4)*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(8*a**(sympy.S(5)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(3)/2)) - b**(sympy.S(1)/4)*(5*sqrt(a) - 2*sqrt(b))*atan(b**(sympy.S(1)/4)*cosh(c + d*x)/sqrt(sqrt(a) - sqrt(b)))/(64*a**(sympy.S(5)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_257():
    f = sinh(c + d*x)**8/(a - b*sinh(c + d*x)**4)**3
    F = tanh(c + d*x)**9/(8*a*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) - tanh(c + d*x)**5*sech(c + d*x)**2/(32*a*b*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - tanh(c + d*x)**3/(32*a*b*d*(a - b)) - (a + 5*b)*tanh(c + d*x)/(32*a*b*d*(a - b)**2) + (2*sqrt(a) + 5*sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (2*sqrt(a) - 5*sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(3)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_258():
    f = sinh(c + d*x)**6/(a - b*sinh(c + d*x)**4)**3
    F = (a*(a + 3*b) - (a**2 + 6*a*b + b**2)*tanh(c + d*x)**2)*tanh(c + d*x)/(8*d*(a - b)**3*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) + (2*a*(a**2 - a*b - 8*b**2)/(a - b)**3 - (2*a**2 + 15*a*b + 3*b**2)*tanh(c + d*x)**2/(a - b)**2)*tanh(c + d*x)/(32*a*b*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - (10*sqrt(a)*sqrt(b) + 4*a + 3*b)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-10*sqrt(a)*sqrt(b) + 4*a + 3*b)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(5)/4)*b**(sympy.S(3)/2)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_259():
    f = sinh(c + d*x)**4/(a - b*sinh(c + d*x)**4)**3
    F = -b*(3*a + b - (4*a + 4*b)*tanh(c + d*x)**2)*tanh(c + d*x)/(8*d*(a - b)**3*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) - (-(17*a + 3*b)*tanh(c + d*x)**2/(a - b)**2 + (9*a**2 - 24*a*b - b**2)/(a - b)**3)*tanh(c + d*x)/(32*a*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - (6*sqrt(a) + 3*sqrt(b))*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (6*sqrt(a) - 3*sqrt(b))*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(7)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_260():
    f = sinh(c + d*x)**2/(a - b*sinh(c + d*x)**4)**3
    F = b*(a*(a + 3*b) - (a**2 + 6*a*b + b**2)*tanh(c + d*x)**2)*tanh(c + d*x)/(8*a*d*(a - b)**3*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) + (2*a*(5*a**2 - 9*a*b - 4*b**2)/(a - b)**3 - (10*a**2 + 15*a*b - 5*b**2)*tanh(c + d*x)**2/(a - b)**2)*tanh(c + d*x)/(32*a**2*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) + (14*sqrt(a)*sqrt(b) + 12*a + 5*b)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*sqrt(b)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - (-14*sqrt(a)*sqrt(b) + 12*a + 5*b)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*sqrt(b)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_261():
    f = (a - b*sinh(c + d*x)**4)**(-3)
    F = -b**2*(3*a + b - (4*a + 4*b)*tanh(c + d*x)**2)*tanh(c + d*x)/(8*a*d*(a - b)**3*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) - b*(-(33*a - 13*b)*tanh(c + d*x)**2/(a - b)**2 + (17*a**2 - 40*a*b + 7*b**2)/(a - b)**3)*tanh(c + d*x)/(32*a**2*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) + (50*sqrt(a)*sqrt(b) + 32*a + 21*b)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) + (-50*sqrt(a)*sqrt(b) + 32*a + 21*b)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_262():
    f = csch(c + d*x)**2/(a - b*sinh(c + d*x)**4)**3
    F = b**2*(a*(a + 3*b) - (a**2 + 6*a*b + b**2)*tanh(c + d*x)**2)*tanh(c + d*x)/(8*a**2*d*(a - b)**3*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)**2) + b*(2*a**2*(9*a - 17*b)/(a - b)**3 - (18*a**2 + 15*a*b - 13*b**2)*tanh(c + d*x)**2/(a - b)**2)*tanh(c + d*x)/(32*a**3*d*(-2*a*tanh(c + d*x)**2 + a + (a - b)*tanh(c + d*x)**4)) - coth(c + d*x)/(a**3*d) + 3*sqrt(b)*(34*sqrt(a)*sqrt(b) + 20*a + 15*b)*atanh(sqrt(sqrt(a) + sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*d*(sqrt(a) + sqrt(b))**(sympy.S(5)/2)) - 3*sqrt(b)*(-34*sqrt(a)*sqrt(b) + 20*a + 15*b)*atanh(sqrt(sqrt(a) - sqrt(b))*tanh(c + d*x)/a**(sympy.S(1)/4))/(64*a**(sympy.S(13)/4)*d*(sqrt(a) - sqrt(b))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_263():
    f = 1/(1 - sinh(x)**4)
    F = tanh(x)/2 + sqrt(2)*atanh(sqrt(2)*tanh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_264():
    f = 1/(sinh(x)**4 + 1)
    F = sqrt(1 + sqrt(2))*log(sqrt(2)*tanh(x)**2 + sqrt(2 + 2*sqrt(2))*tanh(x) + 1)/8 - sqrt(1 + sqrt(2))*log(2*tanh(x)**2 - 2*sqrt(1 + sqrt(2))*tanh(x) + sqrt(2))/8 - atan((-2*tanh(x) + sqrt(1 + sqrt(2)))/sqrt(-1 + sqrt(2)))/(4*sqrt(1 + sqrt(2))) + atan((2*tanh(x) + sqrt(1 + sqrt(2)))/sqrt(-1 + sqrt(2)))/(4*sqrt(1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_265():
    f = 1/(a + b*sinh(x)**5)
    F = 2*(-1)**(sympy.S(9)/10)*atanh((-1)**(sympy.S(3)/10)*((-1)**(sympy.S(3)/5)*a**(sympy.S(1)/5)*tanh(x/2) + b**(sympy.S(1)/5))/sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(2)/5))) - 2*(-1)**(sympy.S(2)/5)*atan((-1)**(sympy.S(2)/5)*(a**(sympy.S(1)/5)*tanh(x/2) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))/sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(2)/5))) + 2*(-1)**(sympy.S(2)/5)*atan((-(-1)**(sympy.S(2)/5)*a**(sympy.S(1)/5)*tanh(x/2) + b**(sympy.S(1)/5))/sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) - b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(-(-1)**(sympy.S(4)/5)*a**(sympy.S(2)/5) - b**(sympy.S(2)/5))) + 2*(-1)**(sympy.S(1)/5)*atanh(((-1)**(sympy.S(1)/5)*a**(sympy.S(1)/5)*tanh(x/2) + b**(sympy.S(1)/5))/sqrt((-1)**(sympy.S(2)/5)*a**(sympy.S(2)/5) + b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt((-1)**(sympy.S(2)/5)*a**(sympy.S(2)/5) + b**(sympy.S(2)/5))) - 2*atanh((-a**(sympy.S(1)/5)*tanh(x/2) + b**(sympy.S(1)/5))/sqrt(a**(sympy.S(2)/5) + b**(sympy.S(2)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(2)/5) + b**(sympy.S(2)/5)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_266():
    f = 1/(a + b*sinh(x)**6)
    F = atanh(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*tanh(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + atanh(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*tanh(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + atanh(sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*tanh(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_267():
    f = 1/(a + b*sinh(x)**8)
    F = -atanh(sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tanh(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh(sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tanh(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh(sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tanh(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh(sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*tanh(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_268():
    f = 1/(sinh(x)**5 + 1)
    F = -2*(-1)**(sympy.S(3)/5)*atan(((-1)**(sympy.S(3)/5)*tanh(x/2) + 1)/sqrt(-1 + (-1)**(sympy.S(1)/5)))/(5*sqrt(-1 + (-1)**(sympy.S(1)/5))) - 2*(-1)**(sympy.S(2)/5)*atan((-1)**(sympy.S(1)/5)*((-1)**(sympy.S(1)/5)*tanh(x/2) + 1)/sqrt(-(-1)**(sympy.S(2)/5)*(1 + (-1)**(sympy.S(2)/5))))/(5*sqrt(-(-1)**(sympy.S(2)/5)*(1 + (-1)**(sympy.S(2)/5)))) - sqrt(2)*atanh(sqrt(2)*(1 - tanh(x/2))/2)/5 - 2*(-1)**(sympy.S(4)/5)*atanh((-(-1)**(sympy.S(4)/5)*tanh(x/2) + 1)/sqrt(1 - (-1)**(sympy.S(3)/5)))/(5*sqrt(1 - (-1)**(sympy.S(3)/5))) - 2*(-1)**(sympy.S(2)/5)*atanh((-(-1)**(sympy.S(2)/5)*tanh(x/2) + 1)/sqrt(1 + (-1)**(sympy.S(4)/5)))/(5*sqrt(1 + (-1)**(sympy.S(4)/5)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_269():
    f = 1/(sinh(x)**6 + 1)
    F = tanh(x)/3 + atanh(sqrt(1 + (-1)**(sympy.S(1)/3))*tanh(x))/(3*sqrt(1 + (-1)**(sympy.S(1)/3))) + atanh(sqrt(1 - (-1)**(sympy.S(2)/3))*tanh(x))/(3*sqrt(1 - (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_270():
    f = 1/(sinh(x)**8 + 1)
    F = atanh(sqrt(1 - (-1)**(sympy.S(1)/4))*tanh(x))/(4*sqrt(1 - (-1)**(sympy.S(1)/4))) + atanh(sqrt(1 + (-1)**(sympy.S(1)/4))*tanh(x))/(4*sqrt(1 + (-1)**(sympy.S(1)/4))) + atanh(sqrt(1 - (-1)**(sympy.S(3)/4))*tanh(x))/(4*sqrt(1 - (-1)**(sympy.S(3)/4))) + atanh(sqrt(1 + (-1)**(sympy.S(3)/4))*tanh(x))/(4*sqrt(1 + (-1)**(sympy.S(3)/4)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_271():
    f = 1/(1 - sinh(x)**5)
    F = -2*(-1)**(sympy.S(1)/10)*atan(((-1)**(sympy.S(1)/10)*tanh(x/2) + I)/sqrt(1 - (-1)**(sympy.S(1)/5)))/(5*sqrt(1 - (-1)**(sympy.S(1)/5))) + sqrt(2)*atanh(sqrt(2)*(tanh(x/2) + 1)/2)/5 - 2*atanh((-tanh(x/2) + (-1)**(sympy.S(3)/5))/sqrt(1 - (-1)**(sympy.S(1)/5)))/(5*sqrt(1 - (-1)**(sympy.S(1)/5))) + 2*atanh((tanh(x/2) + (-1)**(sympy.S(4)/5))/sqrt(1 - (-1)**(sympy.S(3)/5)))/(5*sqrt(1 - (-1)**(sympy.S(3)/5))) - 2*(-1)**(sympy.S(1)/10)*atanh((-1)**(sympy.S(3)/10)*((-1)**(sympy.S(4)/5)*tanh(x/2) + 1)/sqrt((-1)**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)))/(5*sqrt((-1)**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_272():
    f = 1/(1 - sinh(x)**6)
    F = sqrt(2)*atanh(sqrt(2)*tanh(x))/6 + atanh(sqrt(1 - (-1)**(sympy.S(1)/3))*tanh(x))/(3*sqrt(1 - (-1)**(sympy.S(1)/3))) + atanh(sqrt(1 + (-1)**(sympy.S(2)/3))*tanh(x))/(3*sqrt(1 + (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_273():
    f = 1/(1 - sinh(x)**8)
    F = tanh(x)/4 + sqrt(2)*atanh(sqrt(2)*tanh(x))/8 + atanh(sqrt(1 - I)*tanh(x))/(4*sqrt(1 - I)) + atanh(sqrt(1 + I)*tanh(x))/(4*sqrt(1 + I))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_274():
    f = cosh(x)**5/(a*sinh(x)**2 + a)
    F = sinh(x)**3/(3*a) + sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_275():
    f = cosh(x)**4/(a*sinh(x)**2 + a)
    F = x/(2*a) + sinh(x)*cosh(x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_276():
    f = cosh(x)**3/(a*sinh(x)**2 + a)
    F = sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_277():
    f = cosh(x)**2/(a*sinh(x)**2 + a)
    F = x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_278():
    f = cosh(x)/(a*sinh(x)**2 + a)
    F = atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_279():
    f = sech(x)/(a*sinh(x)**2 + a)
    F = tanh(x)*sech(x)/(2*a) + atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_280():
    f = sech(x)**3/(a*sinh(x)**2 + a)
    F = tanh(x)*sech(x)**3/(4*a) + 3*tanh(x)*sech(x)/(8*a) + 3*atan(sinh(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_281():
    f = (a + b*sinh(c + d*x)**2)*cosh(c + d*x)**4
    F = b*sinh(c + d*x)*cosh(c + d*x)**5/(6*d) + x*(3*a/8 - b/16) + (6*a - b)*sinh(c + d*x)*cosh(c + d*x)**3/(24*d) + (6*a - b)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_282():
    f = (a + b*sinh(c + d*x)**2)*cosh(c + d*x)**3
    F = a*sinh(c + d*x)/d + b*sinh(c + d*x)**5/(5*d) + (a + b)*sinh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_283():
    f = (a + b*sinh(c + d*x)**2)*cosh(c + d*x)**2
    F = b*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + x*(a/2 - b/8) + (4*a - b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_284():
    f = (a + b*sinh(c + d*x)**2)*cosh(c + d*x)
    F = a*sinh(c + d*x)/d + b*sinh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_285():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)
    F = b*sinh(c + d*x)/d + (a - b)*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_286():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)**2
    F = b*x + (a - b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_287():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)**3
    F = (a - b)*tanh(c + d*x)*sech(c + d*x)/(2*d) + (a + b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_288():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)**4
    F = a*tanh(c + d*x)/d - (a - b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_289():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)**5
    F = (a - b)*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + (3*a + b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (3*a + b)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_290():
    f = (a + b*sinh(c + d*x)**2)*sech(c + d*x)**6
    F = a*tanh(c + d*x)/d + (a - b)*tanh(c + d*x)**5/(5*d) - (2*a - b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_291():
    f = (a + b*sinh(c + d*x)**2)**2*cosh(c + d*x)**4
    F = b*(a - (a - b)*tanh(c + d*x)**2)*sinh(c + d*x)*cosh(c + d*x)**7/(8*d) + b*(10*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) + x*(3*a**2/8 - a*b/8 + 3*b**2/128) + (48*a**2 - 16*a*b + 3*b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) + (48*a**2 - 16*a*b + 3*b**2)*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_292():
    f = (a + b*sinh(c + d*x)**2)**2*cosh(c + d*x)**3
    F = a**2*sinh(c + d*x)/d + a*(a + 2*b)*sinh(c + d*x)**3/(3*d) + b**2*sinh(c + d*x)**7/(7*d) + b*(2*a + b)*sinh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_293():
    f = (a + b*sinh(c + d*x)**2)**2*cosh(c + d*x)**2
    F = b*(a - (a - b)*tanh(c + d*x)**2)*sinh(c + d*x)*cosh(c + d*x)**5/(6*d) + b*(8*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)**3/(24*d) + x*(a**2/2 - a*b/4 + b**2/16) + (8*a**2 - 4*a*b + b**2)*sinh(c + d*x)*cosh(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_294():
    f = (a + b*sinh(c + d*x)**2)**2*cosh(c + d*x)
    F = a**2*sinh(c + d*x)/d + 2*a*b*sinh(c + d*x)**3/(3*d) + b**2*sinh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_295():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)
    F = b**2*sinh(c + d*x)**3/(3*d) + b*(2*a - b)*sinh(c + d*x)/d + (a - b)**2*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_296():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**2
    F = b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + b*x*(2*a - 3*b/2) + (a - b)**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_297():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**3
    F = b**2*sinh(c + d*x)/d + (a - b)**2*tanh(c + d*x)*sech(c + d*x)/(2*d) + (a - b)*(a + 3*b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_298():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**4
    F = b**2*x - (a - b)**2*tanh(c + d*x)**3/(3*d) + (a**2 - b**2)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_299():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**5
    F = (a - b)*(a + b*sinh(c + d*x)**2)*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + (3*a**2 - 3*b**2)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (3*a**2 + 2*a*b + 3*b**2)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_300():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**6
    F = a**2*tanh(c + d*x)/d - 2*a*(a - b)*tanh(c + d*x)**3/(3*d) + (a - b)**2*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_301():
    f = (a + b*sinh(c + d*x)**2)**2*sech(c + d*x)**7
    F = (a - b)*(a + b*sinh(c + d*x)**2)*tanh(c + d*x)*sech(c + d*x)**5/(6*d) + (a - b)*(5*a + 3*b)*tanh(c + d*x)*sech(c + d*x)**3/(24*d) + (5*a**2 + 2*a*b + b**2)*tanh(c + d*x)*sech(c + d*x)/(16*d) + (5*a**2 + 2*a*b + b**2)*atan(sinh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_302():
    f = (a + b*sinh(c + d*x)**2)**3*cosh(c + d*x)**4
    F = b*(a - (a - b)*tanh(c + d*x)**2)**2*sinh(c + d*x)*cosh(c + d*x)**9/(10*d) + b*(a*(10*a - b) - (2*a - b)*(5*a - 5*b)*tanh(c + d*x)**2)*sinh(c + d*x)*cosh(c + d*x)**7/(80*d) + b*(44*a**2 - 28*a*b + 5*b**2)*sinh(c + d*x)*cosh(c + d*x)**5/(160*d) + x*(3*a/64 - 3*b/256)*(8*a**2 - 2*a*b + b**2) + (4*a - b)*(8*a**2 - 2*a*b + b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(128*d) + (12*a - 3*b)*(8*a**2 - 2*a*b + b**2)*sinh(c + d*x)*cosh(c + d*x)/(256*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_303():
    f = (a + b*sinh(c + d*x)**2)**3*cosh(c + d*x)**3
    F = a**3*sinh(c + d*x)/d + a**2*(a + 3*b)*sinh(c + d*x)**3/(3*d) + 3*a*b*(a + b)*sinh(c + d*x)**5/(5*d) + b**3*sinh(c + d*x)**9/(9*d) + b**2*(3*a + b)*sinh(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_304():
    f = (a + b*sinh(c + d*x)**2)**3*cosh(c + d*x)**2
    F = b*(a - (a - b)*tanh(c + d*x)**2)**2*sinh(c + d*x)*cosh(c + d*x)**7/(8*d) + b*(a*(8*a - b) - (a - b)*(8*a - 5*b)*tanh(c + d*x)**2)*sinh(c + d*x)*cosh(c + d*x)**5/(48*d) + b*(88*a**2 - 68*a*b + 15*b**2)*sinh(c + d*x)*cosh(c + d*x)**3/(192*d) + x*(a**3/2 - 3*a**2*b/8 + 3*a*b**2/16 - 5*b**3/128) + (64*a**3 - 48*a**2*b + 24*a*b**2 - 5*b**3)*sinh(c + d*x)*cosh(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_305():
    f = (a + b*sinh(c + d*x)**2)**3*cosh(c + d*x)
    F = a**3*sinh(c + d*x)/d + a**2*b*sinh(c + d*x)**3/d + 3*a*b**2*sinh(c + d*x)**5/(5*d) + b**3*sinh(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_306():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)
    F = b**3*sinh(c + d*x)**5/(5*d) + b**2*(3*a - b)*sinh(c + d*x)**3/(3*d) + b*(3*a**2 - 3*a*b + b**2)*sinh(c + d*x)/d + (a - b)**3*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_307():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**2
    F = b**3*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + b**2*(12*a - 9*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) + 3*b*x*(8*a**2 - 12*a*b + 5*b**2)/8 + (a - b)**3*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_308():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**3
    F = b**3*sinh(c + d*x)**3/(3*d) + b**2*(3*a - 2*b)*sinh(c + d*x)/d + (a - b)**3*tanh(c + d*x)*sech(c + d*x)/(2*d) + (a - b)**2*(a + 5*b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_309():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**4
    F = b**3*sinh(c + d*x)*cosh(c + d*x)/(2*d) + b**2*x*(3*a - 5*b/2) - (a - b)**3*tanh(c + d*x)**3/(3*d) + (a - b)**2*(a + 2*b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_310():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**5
    F = b**3*sinh(c + d*x)/d + (a - b)**3*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + 3*(a - b)**2*(a + 3*b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (3*a - 3*b)*(4*b**2 + (a + b)**2)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_311():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**6
    F = b**3*x + (a - b)**3*tanh(c + d*x)**5/(5*d) - (a - b)**2*(2*a + b)*tanh(c + d*x)**3/(3*d) + (a**3 - b**3)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_312():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**7
    F = (a - b)*(a + b*sinh(c + d*x)**2)**2*tanh(c + d*x)*sech(c + d*x)**5/(6*d) + (a - b)*(15*a**2 + 14*a*b + 15*b**2)*tanh(c + d*x)*sech(c + d*x)/(48*d) + (a + b)*(5*a**2 - 2*a*b + 5*b**2)*atan(sinh(c + d*x))/(16*d) + (a + b*sinh(c + d*x)**2)*(5*a**2 - 5*b**2)*tanh(c + d*x)*sech(c + d*x)**3/(24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_313():
    f = (a + b*sinh(c + d*x)**2)**3*sech(c + d*x)**8
    F = a**3*tanh(c + d*x)/d - a**2*(a - b)*tanh(c + d*x)**3/d + 3*a*(a - b)**2*tanh(c + d*x)**5/(5*d) - (a - b)**3*tanh(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_314():
    f = cosh(c + d*x)**7/(a + b*sinh(c + d*x)**2)
    F = sinh(c + d*x)**5/(5*b*d) - (a - 3*b)*sinh(c + d*x)**3/(3*b**2*d) + (a**2 - 3*a*b + 3*b**2)*sinh(c + d*x)/(b**3*d) - (a - b)**3*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_315():
    f = cosh(c + d*x)**6/(a + b*sinh(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)**3/(4*b*d) - (4*a - 7*b)*sinh(c + d*x)*cosh(c + d*x)/(8*b**2*d) + x*(8*a**2 - 20*a*b + 15*b**2)/(8*b**3) - (a - b)**(sympy.S(5)/2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_316():
    f = cosh(c + d*x)**5/(a + b*sinh(c + d*x)**2)
    F = sinh(c + d*x)**3/(3*b*d) - (a - 2*b)*sinh(c + d*x)/(b**2*d) + (a - b)**2*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_317():
    f = cosh(c + d*x)**4/(a + b*sinh(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)/(2*b*d) - x*(2*a - 3*b)/(2*b**2) + (a - b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_318():
    f = cosh(c + d*x)**3/(a + b*sinh(c + d*x)**2)
    F = sinh(c + d*x)/(b*d) - (a - b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_319():
    f = cosh(c + d*x)**2/(a + b*sinh(c + d*x)**2)
    F = x/b - sqrt(a - b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_320():
    f = cosh(c + d*x)/(a + b*sinh(c + d*x)**2)
    F = atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_321():
    f = sech(c + d*x)/(a + b*sinh(c + d*x)**2)
    F = atan(sinh(c + d*x))/(d*(a - b)) - sqrt(b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_322():
    f = sech(c + d*x)**2/(a + b*sinh(c + d*x)**2)
    F = tanh(c + d*x)/(d*(a - b)) - b*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_323():
    f = sech(c + d*x)**3/(a + b*sinh(c + d*x)**2)
    F = (a - 3*b)*atan(sinh(c + d*x))/(2*d*(a - b)**2) + tanh(c + d*x)*sech(c + d*x)/(d*(2*a - 2*b)) + b**(sympy.S(3)/2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_324():
    f = sech(c + d*x)**4/(a + b*sinh(c + d*x)**2)
    F = (a - 2*b)*tanh(c + d*x)/(d*(a - b)**2) - tanh(c + d*x)**3/(d*(3*a - 3*b)) + b**2*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_325():
    f = sech(c + d*x)**5/(a + b*sinh(c + d*x)**2)
    F = tanh(c + d*x)*sech(c + d*x)**3/(d*(4*a - 4*b)) + (3*a - 7*b)*tanh(c + d*x)*sech(c + d*x)/(8*d*(a - b)**2) + (3*a**2 - 10*a*b + 15*b**2)*atan(sinh(c + d*x))/(8*d*(a - b)**3) - b**(sympy.S(5)/2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_326():
    f = sech(c + d*x)**6/(a + b*sinh(c + d*x)**2)
    F = tanh(c + d*x)**5/(d*(5*a - 5*b)) - (2*a - 3*b)*tanh(c + d*x)**3/(3*d*(a - b)**2) + (a**2 - 3*a*b + 3*b**2)*tanh(c + d*x)/(d*(a - b)**3) - b**3*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_327():
    f = cosh(c + d*x)**6/(a + b*sinh(c + d*x)**2)**2
    F = sinh(c + d*x)*cosh(c + d*x)/(2*b*d*(a - (a - b)*tanh(c + d*x)**2)) - x*(4*a - 5*b)/(2*b**3) + (a - b)*(2*a - b)*tanh(c + d*x)/(2*a*b**2*d*(a - (a - b)*tanh(c + d*x)**2)) + (a - b)**(sympy.S(3)/2)*(4*a + b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_328():
    f = cosh(c + d*x)**5/(a + b*sinh(c + d*x)**2)**2
    F = sinh(c + d*x)/(b**2*d) + (a - b)**2*sinh(c + d*x)/(2*a*b**2*d*(a + b*sinh(c + d*x)**2)) - (3*a**2 - 2*a*b - b**2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_329():
    f = cosh(c + d*x)**4/(a + b*sinh(c + d*x)**2)**2
    F = x/b**2 - (a - b)*tanh(c + d*x)/(2*a*b*d*(a - (a - b)*tanh(c + d*x)**2)) - sqrt(a - b)*(2*a + b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_330():
    f = cosh(c + d*x)**3/(a + b*sinh(c + d*x)**2)**2
    F = -(a - b)*sinh(c + d*x)/(2*a*b*d*(a + b*sinh(c + d*x)**2)) + (a + b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_331():
    f = cosh(c + d*x)**2/(a + b*sinh(c + d*x)**2)**2
    F = tanh(c + d*x)/(2*a*d*(a - (a - b)*tanh(c + d*x)**2)) + atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_332():
    f = cosh(c + d*x)/(a + b*sinh(c + d*x)**2)**2
    F = sinh(c + d*x)/(2*a*d*(a + b*sinh(c + d*x)**2)) + atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_333():
    f = sech(c + d*x)/(a + b*sinh(c + d*x)**2)**2
    F = atan(sinh(c + d*x))/(d*(a - b)**2) - b*sinh(c + d*x)/(2*a*d*(a - b)*(a + b*sinh(c + d*x)**2)) - sqrt(b)*(3*a - b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_334():
    f = sech(c + d*x)**2/(a + b*sinh(c + d*x)**2)**2
    F = tanh(c + d*x)/(d*(a - b)**2) + b**2*tanh(c + d*x)/(2*a*d*(a - b)**2*(a - (a - b)*tanh(c + d*x)**2)) - b*(4*a - b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_335():
    f = sech(c + d*x)**3/(a + b*sinh(c + d*x)**2)**2
    F = (a - 5*b)*atan(sinh(c + d*x))/(2*d*(a - b)**3) + tanh(c + d*x)*sech(c + d*x)/(d*(a + b*sinh(c + d*x)**2)*(2*a - 2*b)) + b*(a + b)*sinh(c + d*x)/(2*a*d*(a - b)**2*(a + b*sinh(c + d*x)**2)) + b**(sympy.S(3)/2)*(5*a - b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_336():
    f = sech(c + d*x)**4/(a + b*sinh(c + d*x)**2)**2
    F = (a - 3*b)*tanh(c + d*x)/(d*(a - b)**3) - tanh(c + d*x)**3/(3*d*(a - b)**2) - b**3*tanh(c + d*x)/(2*a*d*(a - b)**3*(a - (a - b)*tanh(c + d*x)**2)) + b**2*(6*a - b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_337():
    f = cosh(c + d*x)**6/(a + b*sinh(c + d*x)**2)**3
    F = x/b**3 - (a - b)*tanh(c + d*x)/(4*a*b*d*(a - (a - b)*tanh(c + d*x)**2)**2) - (a - b)*(4*a + 3*b)*tanh(c + d*x)/(8*a**2*b**2*d*(a - (a - b)*tanh(c + d*x)**2)) - sqrt(a - b)*(8*a**2 + 4*a*b + 3*b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_338():
    f = cosh(c + d*x)**5/(a + b*sinh(c + d*x)**2)**3
    F = (-3/b**2 + 3/a**2)*sinh(c + d*x)/(8*d*(a + b*sinh(c + d*x)**2)) - (a - b)*sinh(c + d*x)*cosh(c + d*x)**2/(4*a*b*d*(a + b*sinh(c + d*x)**2)**2) + (3*a**2 + 2*a*b + 3*b**2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_339():
    f = cosh(c + d*x)**4/(a + b*sinh(c + d*x)**2)**3
    F = tanh(c + d*x)/(4*a*d*(a - (a - b)*tanh(c + d*x)**2)**2) + 3*tanh(c + d*x)/(8*a**2*d*(a - (a - b)*tanh(c + d*x)**2)) + 3*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_340():
    f = cosh(c + d*x)**3/(a + b*sinh(c + d*x)**2)**3
    F = -(a - b)*sinh(c + d*x)/(4*a*b*d*(a + b*sinh(c + d*x)**2)**2) + (a + 3*b)*sinh(c + d*x)/(8*a**2*b*d*(a + b*sinh(c + d*x)**2)) + (a + 3*b)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_341():
    f = cosh(c + d*x)**2/(a + b*sinh(c + d*x)**2)**3
    F = -b*tanh(c + d*x)/(4*a*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)**2) + (4*a - 3*b)*tanh(c + d*x)/(8*a**2*d*(a - b)*(a - (a - b)*tanh(c + d*x)**2)) + (4*a - 3*b)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_342():
    f = cosh(c + d*x)/(a + b*sinh(c + d*x)**2)**3
    F = sinh(c + d*x)/(4*a*d*(a + b*sinh(c + d*x)**2)**2) + 3*sinh(c + d*x)/(8*a**2*d*(a + b*sinh(c + d*x)**2)) + 3*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_343():
    f = sech(c + d*x)/(a + b*sinh(c + d*x)**2)**3
    F = atan(sinh(c + d*x))/(d*(a - b)**3) - b*sinh(c + d*x)/(4*a*d*(a - b)*(a + b*sinh(c + d*x)**2)**2) - b*(7*a - 3*b)*sinh(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*sinh(c + d*x)**2)) - sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_344():
    f = sech(c + d*x)**2/(a + b*sinh(c + d*x)**2)**3
    F = tanh(c + d*x)/(d*(a - b)**3) - b**3*tanh(c + d*x)/(4*a*d*(a - b)**3*(a - (a - b)*tanh(c + d*x)**2)**2) + b**2*(12*a - 3*b)*tanh(c + d*x)/(8*a**2*d*(a - b)**3*(a - (a - b)*tanh(c + d*x)**2)) - 3*b*(8*a**2 - 4*a*b + b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_345():
    f = sech(c + d*x)**3/(a + b*sinh(c + d*x)**2)**3
    F = (a - 7*b)*atan(sinh(c + d*x))/(2*d*(a - b)**4) + tanh(c + d*x)*sech(c + d*x)/(d*(a + b*sinh(c + d*x)**2)**2*(2*a - 2*b)) + b*(2*a + b)*sinh(c + d*x)/(4*a*d*(a - b)**2*(a + b*sinh(c + d*x)**2)**2) + b*(a + 3*b)*(4*a - b)*sinh(c + d*x)/(8*a**2*d*(a - b)**3*(a + b*sinh(c + d*x)**2)) + b**(sympy.S(3)/2)*(35*a**2 - 14*a*b + 3*b**2)*atan(sqrt(b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_346():
    f = sech(c + d*x)**4/(a + b*sinh(c + d*x)**2)**3
    F = (a - 4*b)*tanh(c + d*x)/(d*(a - b)**4) - tanh(c + d*x)**3/(3*d*(a - b)**3) + b**4*tanh(c + d*x)/(4*a*d*(a - b)**4*(a - (a - b)*tanh(c + d*x)**2)**2) - b**3*(16*a - 3*b)*tanh(c + d*x)/(8*a**2*d*(a - b)**4*(a - (a - b)*tanh(c + d*x)**2)) + b**2*(48*a**2 - 16*a*b + 3*b**2)*atanh(sqrt(a - b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a - b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_347():
    f = cosh(x)**2/(1 - sinh(x)**2)
    F = -x + sqrt(2)*atanh(sqrt(2)*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_348():
    f = cosh(x)**3/(1 - sinh(x)**2)
    F = -sinh(x) + 2*atanh(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_349():
    f = cosh(x)**4/(1 - sinh(x)**2)
    F = -5*x/2 - sinh(x)*cosh(x)/2 + 2*sqrt(2)*atanh(sqrt(2)*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_350():
    f = sqrt(a + b*sinh(e + f*x)**2)*cosh(e + f*x)**3
    F = -a*(a - 4*b)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(8*b**(sympy.S(3)/2)*f) - (a - 4*b)*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(8*b*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)/(4*b*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_351():
    f = sqrt(a + b*sinh(e + f*x)**2)*cosh(e + f*x)
    F = a*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*sqrt(b)*f) + sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_352():
    f = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)
    F = sqrt(b)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/f + sqrt(a - b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_353():
    f = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**3
    F = a*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*f*sqrt(a - b)) + sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_354():
    f = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**5
    F = a*(3*a - 4*b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(8*f*(a - b)**(sympy.S(3)/2)) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)*sech(e + f*x)**3/(f*(4*a - 4*b)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*tanh(e + f*x)*sech(e + f*x)/(f*(8*a - 8*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_355():
    f = sqrt(a + b*sinh(e + f*x)**2)*cosh(e + f*x)**4
    F = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)*cosh(e + f*x)/(5*b*f) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - 6*b)*sinh(e + f*x)*cosh(e + f*x)/(15*b*f) - (a - 9*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a**2 - 7*a*b - 3*b**2)*tanh(e + f*x)/(15*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a**2 - 7*a*b - 3*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_356():
    f = sqrt(a + b*sinh(e + f*x)**2)*cosh(e + f*x)**2
    F = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) + 2*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (a + b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(3*b*f) - (a + b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_357():
    f = sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_358():
    f = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**2
    F = sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_359():
    f = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**4
    F = -b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b)) + sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)**2/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_360():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*cosh(e + f*x)**3
    F = -a**2*(a - 6*b)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(16*b**(sympy.S(3)/2)*f) - a*(a - 6*b)*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(16*b*f) - (a - 6*b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)/(24*b*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*sinh(e + f*x)/(6*b*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_361():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*cosh(e + f*x)
    F = 3*a**2*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(8*sqrt(b)*f) + 3*a*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(8*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sinh(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_362():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)
    F = sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*f) + b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(2*f) + (a - b)**(sympy.S(3)/2)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_363():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**3
    F = b**(sympy.S(3)/2)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/f + sqrt(a - b)*(a + 2*b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*f) + (a - b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_364():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**5
    F = 3*a**2*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(8*f*sqrt(a - b)) + 3*a*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)/(8*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)*sech(e + f*x)**3/(4*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_365():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**7
    F = a**2*(5*a - 6*b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(16*f*(a - b)**(sympy.S(3)/2)) + a*sqrt(a + b*sinh(e + f*x)**2)*(5*a - 6*b)*tanh(e + f*x)*sech(e + f*x)/(f*(16*a - 16*b)) + (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*tanh(e + f*x)*sech(e + f*x)**5/(f*(6*a - 6*b)) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(5*a - 6*b)*tanh(e + f*x)*sech(e + f*x)**3/(f*(24*a - 24*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_366():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*cosh(e + f*x)**4
    F = b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)**5/(7*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a - 2*b)*sinh(e + f*x)*cosh(e + f*x)**3/(35*f) + sqrt(a + b*sinh(e + f*x)**2)*(a**2 + 9*a*b - 2*b**2)*sinh(e + f*x)*cosh(e + f*x)/(35*b*f) - sqrt(a + b*sinh(e + f*x)**2)*(a**2 - 18*a*b + b**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(35*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*(a**2 - 6*a*b + b**2)*tanh(e + f*x)/(35*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*(a**2 - 6*a*b + b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(35*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_367():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*cosh(e + f*x)**2
    F = b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)**3/(5*f) + sqrt(a + b*sinh(e + f*x)**2)*(6*a - 2*b)*sinh(e + f*x)*cosh(e + f*x)/(15*f) + sqrt(a + b*sinh(e + f*x)**2)*(9*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 + 7*a*b - 2*b**2)*tanh(e + f*x)/(15*b*f) - sqrt(a + b*sinh(e + f*x)**2)*(3*a**2 + 7*a*b - 2*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(15*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_368():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)*elliptic_f(I*e + I*f*x, b/a)/(3*f*sqrt(a + b*sinh(e + f*x)**2)) + b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_369():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**2
    F = b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f + (a - b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f + (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_370():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**4
    F = -b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (a - b)*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)**2/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 2*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_371():
    f = cosh(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)/(2*b*f) - (a - 2*b)*atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_372():
    f = cosh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_373():
    f = sech(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_374():
    f = sech(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = (a - 2*b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*f*(a - b)**(sympy.S(3)/2)) + sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)/(f*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_375():
    f = cosh(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*b*f) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*tanh(e + f*x)/(3*b**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - (a - 3*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_376():
    f = cosh(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(b*f) - sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_377():
    f = 1/sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_378():
    f = sech(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) - b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_379():
    f = sech(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)**2/(f*(3*a - 3*b)) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 4*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2) - b*(a - 3*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_380():
    f = cosh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(b**(sympy.S(3)/2)*f) - (a - b)*sinh(e + f*x)/(a*b*f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_381():
    f = cosh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = sinh(e + f*x)/(a*f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_382():
    f = sech(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(f*(a - b)**(sympy.S(3)/2)) - b*sinh(e + f*x)/(a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_383():
    f = sech(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = (a - 4*b)*atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(2*f*(a - b)**(sympy.S(5)/2)) + tanh(e + f*x)*sech(e + f*x)/(f*sqrt(a + b*sinh(e + f*x)**2)*(2*a - 2*b)) + b*(a + 2*b)*sinh(e + f*x)/(2*a*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_384():
    f = cosh(e + f*x)**6/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -(a - b)*sinh(e + f*x)*cosh(e + f*x)**3/(a*b*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(4*a - 3*b)*sinh(e + f*x)*cosh(e + f*x)/(3*a*b**2*f) - sqrt(a + b*sinh(e + f*x)**2)*(4*a - 6*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 13*a*b + 3*b**2)*tanh(e + f*x)/(3*a*b**3*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 13*a*b + 3*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*b**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_385():
    f = cosh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -(a - b)*sinh(e + f*x)*cosh(e + f*x)/(a*b*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*tanh(e + f*x)/(a*b**2*f) - sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_386():
    f = cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(sqrt(a)*sqrt(b)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_387():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(a*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_388():
    f = sech(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = tanh(e + f*x)/(f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - 2*b*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2) + sqrt(b)*(a + b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(sqrt(a)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_389():
    f = cosh(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(b**(sympy.S(5)/2)*f) - (a - b)*sinh(e + f*x)*cosh(e + f*x)**2/(3*a*b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - (a - b)*(3*a + 2*b)*sinh(e + f*x)/(3*a**2*b**2*f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_390():
    f = cosh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sinh(e + f*x)*cosh(e + f*x)**2/(3*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + 2*sinh(e + f*x)/(3*a**2*f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_391():
    f = cosh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sinh(e + f*x)/(3*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + 2*sinh(e + f*x)/(3*a**2*f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_392():
    f = sech(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = atan(sqrt(a - b)*sinh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2))/(f*(a - b)**(sympy.S(5)/2)) - b*sinh(e + f*x)/(3*a*f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - b*(5*a - 2*b)*sinh(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_393():
    f = cosh(e + f*x)**6/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -(a - b)*sinh(e + f*x)*cosh(e + f*x)**3/(3*a*b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - (2*a - 2*b)*(2*a + b)*sinh(e + f*x)*cosh(e + f*x)/(3*a**2*b**2*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(4*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*b**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 3*a*b - 2*b**2)*tanh(e + f*x)/(3*a**2*b**3*f) - sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 3*a*b - 2*b**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*b**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_394():
    f = cosh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -(a - b)*sinh(e + f*x)*cosh(e + f*x)/(3*a*b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*b*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (2*a + 2*b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_395():
    f = cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) + (a - 2*b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*a**(sympy.S(3)/2)*sqrt(b)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_396():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(3*a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - b*(4*a - 2*b)*sinh(e + f*x)*cosh(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*a**2*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_397():
    f = sech(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = tanh(e + f*x)/(f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + b*(3*a + b)*sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a - b)**2*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - b*sqrt(a + b*sinh(e + f*x)**2)*(9*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**3) + sqrt(b)*(3*a**2 + 7*a*b - 2*b**2)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*a**(sympy.S(3)/2)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**3*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_398():
    f = (d*cosh(e + f*x))**m*(a + b*sinh(e + f*x)**2)**p
    F = d*(d*cosh(e + f*x))**(m - 1)*(a + b*sinh(e + f*x)**2)**p*(cosh(e + f*x)**2)**(sympy.S.Half - m/2)*sinh(e + f*x)*appellf1(sympy.S.Half, -p, sympy.S.Half - m/2, sympy.S(3)/2, -b*sinh(e + f*x)**2/a, -sinh(e + f*x)**2)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_399():
    f = (a + b*sinh(e + f*x)**2)**p*cosh(e + f*x)**5
    F = (a + b*sinh(e + f*x)**2)**(p + 1)*sinh(e + f*x)*cosh(e + f*x)**2/(b*f*(2*p + 5)) - (a + b*sinh(e + f*x)**2)**(p + 1)*(3*a - b*(2*p + 7))*sinh(e + f*x)/(b**2*f*(2*p + 3)*(2*p + 5)) + (a + b*sinh(e + f*x)**2)**p*(3*a**2 - 2*a*b*(2*p + 5) + b**2*(4*p**2 + 16*p + 15))*sinh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sinh(e + f*x)**2/a)/(b**2*f*(1 + b*sinh(e + f*x)**2/a)**p*(2*p + 3)*(2*p + 5))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_400():
    f = (a + b*sinh(e + f*x)**2)**p*cosh(e + f*x)**3
    F = (a + b*sinh(e + f*x)**2)**(p + 1)*sinh(e + f*x)/(b*f*(2*p + 3)) - (a - b*(2*p + 3))*(a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sinh(e + f*x)**2/a)/(b*f*(1 + b*sinh(e + f*x)**2/a)**p*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_401():
    f = (a + b*sinh(e + f*x)**2)**p*cosh(e + f*x)
    F = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_402():
    f = (a + b*sinh(e + f*x)**2)**p*sech(e + f*x)
    F = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_403():
    f = (a + b*sinh(e + f*x)**2)**p*sech(e + f*x)**3
    F = (a + b*sinh(e + f*x)**2)**p*sinh(e + f*x)*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_404():
    f = (a + b*sinh(e + f*x)**2)**p*cosh(e + f*x)**4
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*tanh(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_405():
    f = (a + b*sinh(e + f*x)**2)**p*cosh(e + f*x)**2
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*tanh(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_406():
    f = (a + b*sinh(e + f*x)**2)**p
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*tanh(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_407():
    f = (a + b*sinh(e + f*x)**2)**p*sech(e + f*x)**2
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*tanh(e + f*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_408():
    f = (a + b*sinh(e + f*x)**2)**p*sech(e + f*x)**4
    F = (a + b*sinh(e + f*x)**2)**p*sqrt(cosh(e + f*x)**2)*tanh(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/2, -p, sympy.S(3)/2, -sinh(e + f*x)**2, -b*sinh(e + f*x)**2/a)/(f*(1 + b*sinh(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_409():
    f = cosh(c + d*x)**5/(a + b*sqrt(sinh(c + d*x)))
    F = -a**3*sinh(c + d*x)**3/(3*b**4*d) - a**3*(a**4 + 2*b**4)*sinh(c + d*x)/(b**8*d) + 2*a**2*sinh(c + d*x)**(sympy.S(7)/2)/(7*b**3*d) + 2*a**2*(a**4 + 2*b**4)*sinh(c + d*x)**(sympy.S(3)/2)/(3*b**7*d) - a*sinh(c + d*x)**4/(4*b**2*d) - a*(a**4 + 2*b**4)*sinh(c + d*x)**2/(2*b**6*d) - 2*a*(a**4 + b**4)**2*log(a + b*sqrt(sinh(c + d*x)))/(b**10*d) + 2*sinh(c + d*x)**(sympy.S(9)/2)/(9*b*d) + (2*a**4 + 4*b**4)*sinh(c + d*x)**(sympy.S(5)/2)/(5*b**5*d) + 2*(a**4 + b**4)**2*sqrt(sinh(c + d*x))/(b**9*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_410():
    f = cosh(c + d*x)**3/(a + b*sqrt(sinh(c + d*x)))
    F = -a**3*sinh(c + d*x)/(b**4*d) + 2*a**2*sinh(c + d*x)**(sympy.S(3)/2)/(3*b**3*d) - a*sinh(c + d*x)**2/(2*b**2*d) - 2*a*(a**4 + b**4)*log(a + b*sqrt(sinh(c + d*x)))/(b**6*d) + 2*sinh(c + d*x)**(sympy.S(5)/2)/(5*b*d) + (2*a**4 + 2*b**4)*sqrt(sinh(c + d*x))/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_411():
    f = cosh(c + d*x)/(a + b*sqrt(sinh(c + d*x)))
    F = -2*a*log(a + b*sqrt(sinh(c + d*x)))/(b**2*d) + 2*sqrt(sinh(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_412():
    f = sech(c + d*x)/(a + b*sqrt(sinh(c + d*x)))
    F = a**3*atan(sinh(c + d*x))/(d*(a**4 + b**4)) - 2*a*b**2*log(a + b*sqrt(sinh(c + d*x)))/(d*(a**4 + b**4)) + a*b**2*log(cosh(c + d*x))/(d*(a**4 + b**4)) - sqrt(2)*b*(a**2 - b**2)*atan(sqrt(2)*sqrt(sinh(c + d*x)) - 1)/(2*d*(a**4 + b**4)) - sqrt(2)*b*(a**2 - b**2)*atan(sqrt(2)*sqrt(sinh(c + d*x)) + 1)/(2*d*(a**4 + b**4)) - sqrt(2)*b*(a**2 + b**2)*log(-sqrt(2)*sqrt(sinh(c + d*x)) + sinh(c + d*x) + 1)/(4*d*(a**4 + b**4)) + sqrt(2)*b*(a**2 + b**2)*log(sqrt(2)*sqrt(sinh(c + d*x)) + sinh(c + d*x) + 1)/(4*d*(a**4 + b**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_413():
    f = cosh(c + d*x)**5/(a + b*sqrt(sinh(c + d*x)))**2
    F = -8*a**3*sinh(c + d*x)**(sympy.S(5)/2)/(5*b**5*d) - 16*a**3*(a**4 + b**4)*sqrt(sinh(c + d*x))/(b**9*d) + a**2*sinh(c + d*x)**3/(b**4*d) + a**2*(7*a**4 + 6*b**4)*sinh(c + d*x)/(b**8*d) - 4*a*sinh(c + d*x)**(sympy.S(7)/2)/(7*b**3*d) - 4*a*(3*a**4 + 2*b**4)*sinh(c + d*x)**(sympy.S(3)/2)/(3*b**7*d) + 2*a*(a**4 + b**4)**2/(b**10*d*(a + b*sqrt(sinh(c + d*x)))) + sinh(c + d*x)**4/(4*b**2*d) + (5*a**4 + 2*b**4)*sinh(c + d*x)**2/(2*b**6*d) + (2*a**4 + 2*b**4)*(9*a**4 + b**4)*log(a + b*sqrt(sinh(c + d*x)))/(b**10*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_414():
    f = cosh(c + d*x)**3/(a + b*sqrt(sinh(c + d*x)))**2
    F = -8*a**3*sqrt(sinh(c + d*x))/(b**5*d) + 3*a**2*sinh(c + d*x)/(b**4*d) - 4*a*sinh(c + d*x)**(sympy.S(3)/2)/(3*b**3*d) + 2*a*(a**4 + b**4)/(b**6*d*(a + b*sqrt(sinh(c + d*x)))) + sinh(c + d*x)**2/(2*b**2*d) + (10*a**4 + 2*b**4)*log(a + b*sqrt(sinh(c + d*x)))/(b**6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_415():
    f = cosh(c + d*x)/(a + b*sqrt(sinh(c + d*x)))**2
    F = 2*a/(b**2*d*(a + b*sqrt(sinh(c + d*x)))) + 2*log(a + b*sqrt(sinh(c + d*x)))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_416():
    f = sech(c + d*x)/(a + b*sqrt(sinh(c + d*x)))**2
    F = a**2*(a**4 - 3*b**4)*atan(sinh(c + d*x))/(d*(a**4 + b**4)**2) + 2*a*b**2/(d*(a + b*sqrt(sinh(c + d*x)))*(a**4 + b**4)) - sqrt(2)*a*b*(a**4 - 2*a**2*b**2 - b**4)*atan(sqrt(2)*sqrt(sinh(c + d*x)) - 1)/(d*(a**4 + b**4)**2) - sqrt(2)*a*b*(a**4 - 2*a**2*b**2 - b**4)*atan(sqrt(2)*sqrt(sinh(c + d*x)) + 1)/(d*(a**4 + b**4)**2) - sqrt(2)*a*b*(a**4 + 2*a**2*b**2 - b**4)*log(-sqrt(2)*sqrt(sinh(c + d*x)) + sinh(c + d*x) + 1)/(2*d*(a**4 + b**4)**2) + sqrt(2)*a*b*(a**4 + 2*a**2*b**2 - b**4)*log(sqrt(2)*sqrt(sinh(c + d*x)) + sinh(c + d*x) + 1)/(2*d*(a**4 + b**4)**2) - 2*b**2*(3*a**4 - b**4)*log(a + b*sqrt(sinh(c + d*x)))/(d*(a**4 + b**4)**2) + b**2*(3*a**4 - b**4)*log(cosh(c + d*x))/(d*(a**4 + b**4)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_417():
    f = cosh(c + d*x)**5/(a + b*sinh(c + d*x)**n)
    F = sinh(c + d*x)**5*hyper((1, 5/n), ((n + 5)/n,), -b*sinh(c + d*x)**n/a)/(5*a*d) + 2*sinh(c + d*x)**3*hyper((1, 3/n), ((n + 3)/n,), -b*sinh(c + d*x)**n/a)/(3*a*d) + sinh(c + d*x)*hyper((1, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_418():
    f = cosh(c + d*x)**3/(a + b*sinh(c + d*x)**n)
    F = sinh(c + d*x)**3*hyper((1, 3/n), ((n + 3)/n,), -b*sinh(c + d*x)**n/a)/(3*a*d) + sinh(c + d*x)*hyper((1, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_419():
    f = cosh(c + d*x)/(a + b*sinh(c + d*x)**n)
    F = sinh(c + d*x)*hyper((1, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_420():
    f = cosh(c + d*x)**5/(a + b*sinh(c + d*x)**n)**2
    F = sinh(c + d*x)**5*hyper((2, 5/n), ((n + 5)/n,), -b*sinh(c + d*x)**n/a)/(5*a**2*d) + 2*sinh(c + d*x)**3*hyper((2, 3/n), ((n + 3)/n,), -b*sinh(c + d*x)**n/a)/(3*a**2*d) + sinh(c + d*x)*hyper((2, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_421():
    f = cosh(c + d*x)**3/(a + b*sinh(c + d*x)**n)**2
    F = sinh(c + d*x)**3*hyper((2, 3/n), ((n + 3)/n,), -b*sinh(c + d*x)**n/a)/(3*a**2*d) + sinh(c + d*x)*hyper((2, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_422():
    f = cosh(c + d*x)/(a + b*sinh(c + d*x)**n)**2
    F = sinh(c + d*x)*hyper((2, 1/n), (1 + 1/n,), -b*sinh(c + d*x)**n/a)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_423():
    f = coth(x)/(1 - sinh(x)**2)
    F = -log(1 - sinh(x)**2)/2 + log(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_424():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)**5
    F = -a**2/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2)) + 2*a/(f*sqrt(a*cosh(e + f*x)**2)) + sqrt(a*cosh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_425():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)**3
    F = a/(f*sqrt(a*cosh(e + f*x)**2)) + sqrt(a*cosh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_426():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)
    F = sqrt(a*cosh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_427():
    f = sqrt(a*sinh(e + f*x)**2 + a)*coth(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/f + sqrt(a*cosh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_428():
    f = sqrt(a*sinh(e + f*x)**2 + a)*coth(e + f*x)**3
    F = -3*sqrt(a)*atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/(2*f) + 3*sqrt(a*cosh(e + f*x)**2)/(2*f) - (a*cosh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**2/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_429():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)**6
    F = -sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)**5/(4*f) - 5*sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)**3/(8*f) + 15*sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/(8*f) - 15*sqrt(a*cosh(e + f*x)**2)*atan(sinh(e + f*x))*sech(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_430():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)**4
    F = -sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)**3/(2*f) + 3*sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/(2*f) - 3*sqrt(a*cosh(e + f*x)**2)*atan(sinh(e + f*x))*sech(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_431():
    f = sqrt(a*sinh(e + f*x)**2 + a)*tanh(e + f*x)**2
    F = sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a*cosh(e + f*x)**2)*atan(sinh(e + f*x))*sech(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_432():
    f = sqrt(a*sinh(e + f*x)**2 + a)*coth(e + f*x)**2
    F = sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)*sech(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_433():
    f = sqrt(a*sinh(e + f*x)**2 + a)*coth(e + f*x)**4
    F = sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)**3*sech(e + f*x)/(3*f) - 2*sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)*sech(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_434():
    f = sqrt(a*sinh(e + f*x)**2 + a)*coth(e + f*x)**6
    F = sqrt(a*cosh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)**5*sech(e + f*x)/(5*f) - sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)**3*sech(e + f*x)/f - 3*sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)*sech(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_435():
    f = tanh(e + f*x)**5/sqrt(a*sinh(e + f*x)**2 + a)
    F = -a**2/(5*f*(a*cosh(e + f*x)**2)**(sympy.S(5)/2)) + 2*a/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2)) - 1/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_436():
    f = tanh(e + f*x)**3/sqrt(a*sinh(e + f*x)**2 + a)
    F = a/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2)) - 1/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_437():
    f = tanh(e + f*x)/sqrt(a*sinh(e + f*x)**2 + a)
    F = -1/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_438():
    f = coth(e + f*x)/sqrt(a*sinh(e + f*x)**2 + a)
    F = -atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_439():
    f = coth(e + f*x)**3/sqrt(a*sinh(e + f*x)**2 + a)
    F = -sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)**2/(2*a*f) - atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_440():
    f = tanh(e + f*x)**4/sqrt(a*sinh(e + f*x)**2 + a)
    F = 3*cosh(e + f*x)*atan(sinh(e + f*x))/(8*f*sqrt(a*cosh(e + f*x)**2)) - tanh(e + f*x)**3/(4*f*sqrt(a*cosh(e + f*x)**2)) - 3*tanh(e + f*x)/(8*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_441():
    f = tanh(e + f*x)**2/sqrt(a*sinh(e + f*x)**2 + a)
    F = cosh(e + f*x)*atan(sinh(e + f*x))/(2*f*sqrt(a*cosh(e + f*x)**2)) - tanh(e + f*x)/(2*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_442():
    f = coth(e + f*x)**2/sqrt(a*sinh(e + f*x)**2 + a)
    F = -coth(e + f*x)/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_443():
    f = coth(e + f*x)**4/sqrt(a*sinh(e + f*x)**2 + a)
    F = -coth(e + f*x)*csch(e + f*x)**2/(3*f*sqrt(a*cosh(e + f*x)**2)) - coth(e + f*x)/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_444():
    f = coth(e + f*x)**6/sqrt(a*sinh(e + f*x)**2 + a)
    F = -coth(e + f*x)*csch(e + f*x)**4/(5*f*sqrt(a*cosh(e + f*x)**2)) - 2*coth(e + f*x)*csch(e + f*x)**2/(3*f*sqrt(a*cosh(e + f*x)**2)) - coth(e + f*x)/(f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_445():
    f = tanh(e + f*x)**5/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -a**2/(7*f*(a*cosh(e + f*x)**2)**(sympy.S(7)/2)) + 2*a/(5*f*(a*cosh(e + f*x)**2)**(sympy.S(5)/2)) - 1/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_446():
    f = tanh(e + f*x)**3/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = a/(5*f*(a*cosh(e + f*x)**2)**(sympy.S(5)/2)) - 1/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_447():
    f = tanh(e + f*x)/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -1/(3*f*(a*cosh(e + f*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_448():
    f = coth(e + f*x)/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = 1/(a*f*sqrt(a*cosh(e + f*x)**2)) - atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_449():
    f = coth(e + f*x)**3/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -sqrt(a*cosh(e + f*x)**2)*csch(e + f*x)**2/(2*a**2*f) + atanh(sqrt(a*cosh(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_450():
    f = tanh(e + f*x)**2/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = cosh(e + f*x)*atan(sinh(e + f*x))/(8*a*f*sqrt(a*cosh(e + f*x)**2)) - tanh(e + f*x)*sech(e + f*x)**2/(4*a*f*sqrt(a*cosh(e + f*x)**2)) + tanh(e + f*x)/(8*a*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_451():
    f = coth(e + f*x)**2/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -cosh(e + f*x)*atan(sinh(e + f*x))/(a*f*sqrt(a*cosh(e + f*x)**2)) - coth(e + f*x)/(a*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_452():
    f = coth(e + f*x)**4/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -coth(e + f*x)*csch(e + f*x)**2/(3*a*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_453():
    f = coth(e + f*x)**6/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -coth(e + f*x)*csch(e + f*x)**4/(5*a*f*sqrt(a*cosh(e + f*x)**2)) - coth(e + f*x)*csch(e + f*x)**2/(3*a*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_454():
    f = coth(e + f*x)**8/(a*sinh(e + f*x)**2 + a)**(sympy.S(3)/2)
    F = -coth(e + f*x)*csch(e + f*x)**6/(7*a*f*sqrt(a*cosh(e + f*x)**2)) - 2*coth(e + f*x)*csch(e + f*x)**4/(5*a*f*sqrt(a*cosh(e + f*x)**2)) - coth(e + f*x)*csch(e + f*x)**2/(3*a*f*sqrt(a*cosh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_455():
    f = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)**5
    F = -(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**4/(f*(4*a - 4*b)) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(8*a - 7*b)*sech(e + f*x)**2/(8*f*(a - b)**2) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 24*a*b + 15*b**2)/(8*f*(a - b)**2) - (8*a**2 - 24*a*b + 15*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_456():
    f = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)**3
    F = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*sech(e + f*x)**2/(f*(2*a - 2*b)) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 3*b)/(f*(2*a - 2*b)) - (2*a - 3*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(2*f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_457():
    f = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)
    F = -sqrt(a - b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/f + sqrt(a + b*sinh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_458():
    f = sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/f + sqrt(a + b*sinh(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_459():
    f = sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)**3
    F = -(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**2/(2*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + b)/(2*a*f) - (2*a + b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_460():
    f = sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)**5
    F = -(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*csch(e + f*x)**4/(4*a*f) - (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(8*a - b)*csch(e + f*x)**2/(8*a**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 + 8*a*b - b**2)/(8*a**2*f) - (8*a**2 + 8*a*b - b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_461():
    f = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)**4
    F = -sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*tanh(e + f*x)/(f*(3*a - 3*b)) - sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)**3/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*tanh(e + f*x)/(f*(3*a - 3*b)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_462():
    f = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)**2
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f - 2*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_463():
    f = sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_464():
    f = sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)**2
    F = 2*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/f - sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/f - 2*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (a + b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_465():
    f = sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)**4
    F = -sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)**3/(3*f) - sqrt(a + b*sinh(e + f*x)**2)*(3*a + b)*coth(e + f*x)/(3*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(7*a + b)*tanh(e + f*x)/(3*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a + 5*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a + b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_466():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)**5
    F = -(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*sech(e + f*x)**4/(f*(4*a - 4*b)) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 - 40*a*b + 35*b**2)/(f*(8*a - 8*b)) + (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*(8*a - 9*b)*sech(e + f*x)**2/(8*f*(a - b)**2) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(8*a**2 - 40*a*b + 35*b**2)/(24*f*(a - b)**2) - (8*a**2 - 40*a*b + 35*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(8*f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_467():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)**3
    F = -sqrt(a - b)*(2*a - 5*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(2*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*sech(e + f*x)**2/(f*(2*a - 2*b)) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(2*a - 5*b)/(f*(6*a - 6*b)) + sqrt(a + b*sinh(e + f*x)**2)*(2*a - 5*b)/(2*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_468():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)
    F = -(a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/f + (a - b)*sqrt(a + b*sinh(e + f*x)**2)/f + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_469():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/f + a*sqrt(a + b*sinh(e + f*x)**2)/f + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_470():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)**3
    F = -sqrt(a)*(2*a + 3*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(2*f) + sqrt(a + b*sinh(e + f*x)**2)*(2*a + 3*b)/(2*f) - (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*csch(e + f*x)**2/(2*a*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(2*a + 3*b)/(6*a*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_471():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)**5
    F = -(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*csch(e + f*x)**4/(4*a*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a**2 + 3*b*(8*a + b))/(8*a*f) - (a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)*(8*a + b)*csch(e + f*x)**2/(8*a**2*f) + (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(8*a**2 + 3*b*(8*a + b))/(24*a**2*f) - (8*a**2 + 3*b*(8*a + b))*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_472():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)**4
    F = (a - 2*b)*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)**2*tanh(e + f*x)/f - (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)**3/(3*f) - sqrt(a + b*sinh(e + f*x)**2)*(3*a - 8*b)*sinh(e + f*x)*cosh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a - 16*b)*tanh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 8*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(8*a - 16*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_473():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)**2
    F = 4*b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*tanh(e + f*x)/f + sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*tanh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_474():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = I*a*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)*elliptic_f(I*e + I*f*x, b/a)/(3*f*sqrt(a + b*sinh(e + f*x)**2)) + b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*f*sqrt(1 + b*sinh(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_475():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)**2
    F = 4*b*sqrt(a + b*sinh(e + f*x)**2)*sinh(e + f*x)*cosh(e + f*x)/(3*f) - (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)/f + sqrt(a + b*sinh(e + f*x)**2)*(7*a + b)*tanh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a + 5*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a + b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_476():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)**4
    F = -(a + b)*sqrt(a + b*sinh(e + f*x)**2)*cosh(e + f*x)**2*coth(e + f*x)/f - (a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*coth(e + f*x)**3/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a + 5*b)*sinh(e + f*x)*cosh(e + f*x)/(3*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a + 8*b)*tanh(e + f*x)/(3*f) - sqrt(a + b*sinh(e + f*x)**2)*(8*a + 8*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + (a + 3*b)*sqrt(a + b*sinh(e + f*x)**2)*(3*a + b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_477():
    f = tanh(e + f*x)**5/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**4/(f*(4*a - 4*b)) + sqrt(a + b*sinh(e + f*x)**2)*(8*a - 5*b)*sech(e + f*x)**2/(8*f*(a - b)**2) - (8*a**2 - 8*a*b + 3*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_478():
    f = tanh(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/(f*(2*a - 2*b)) - (2*a - b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_479():
    f = tanh(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = -atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(f*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_480():
    f = coth(e + f*x)/sqrt(a + b*sinh(e + f*x)**2)
    F = -atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_481():
    f = coth(e + f*x)**3/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**2/(2*a*f) - (2*a - b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_482():
    f = coth(e + f*x)**5/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*csch(e + f*x)**4/(4*a*f) - sqrt(a + b*sinh(e + f*x)**2)*(8*a - 3*b)*csch(e + f*x)**2/(8*a**2*f) - (8*a**2 - 8*a*b + 3*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_483():
    f = tanh(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)*sech(e + f*x)**2/(f*(3*a - 3*b)) + (-4*a + 2*b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_484():
    f = tanh(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_485():
    f = 1/sqrt(a + b*sinh(e + f*x)**2)
    F = -I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(f*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_486():
    f = coth(e + f*x)**2/sqrt(a + b*sinh(e + f*x)**2)
    F = sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(a*f) - sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/(a*f) - sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_487():
    f = coth(e + f*x)**4/sqrt(a + b*sinh(e + f*x)**2)
    F = -sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)*csch(e + f*x)**2/(3*a*f) + (-4*a + 2*b)*sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/(3*a**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(4*a - 2*b)*tanh(e + f*x)/(3*a**2*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(4*a - 2*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_488():
    f = tanh(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -sech(e + f*x)**4/(f*sqrt(a + b*sinh(e + f*x)**2)*(4*a - 4*b)) + (8*a - 3*b)*sech(e + f*x)**2/(8*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) + (8*a**2 + 8*a*b - b**2)/(8*f*(a - b)**3*sqrt(a + b*sinh(e + f*x)**2)) - (8*a**2 + 8*a*b - b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_489():
    f = tanh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = sech(e + f*x)**2/(f*sqrt(a + b*sinh(e + f*x)**2)*(2*a - 2*b)) + (2*a + b)/(2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - (2*a + b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_490():
    f = tanh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = 1/(f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_491():
    f = coth(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = 1/(a*f*sqrt(a + b*sinh(e + f*x)**2)) - atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_492():
    f = coth(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -csch(e + f*x)**2/(2*a*f*sqrt(a + b*sinh(e + f*x)**2)) + (2*a - 3*b)/(2*a**2*f*sqrt(a + b*sinh(e + f*x)**2)) - (2*a - 3*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_493():
    f = coth(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -csch(e + f*x)**4/(4*a*f*sqrt(a + b*sinh(e + f*x)**2)) - (8*a - 5*b)*csch(e + f*x)**2/(8*a**2*f*sqrt(a + b*sinh(e + f*x)**2)) + (8*a**2 - 24*a*b + 15*b**2)/(8*a**3*f*sqrt(a + b*sinh(e + f*x)**2)) - (8*a**2 - 24*a*b + 15*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_494():
    f = tanh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -sqrt(a)*sqrt(b)*(7*a + b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**3*sqrt(a + b*sinh(e + f*x)**2)) - 4*a*tanh(e + f*x)/(3*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) + tanh(e + f*x)*sech(e + f*x)**2/(f*sqrt(a + b*sinh(e + f*x)**2)*(3*a - 3*b)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a + 5*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_495():
    f = tanh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*sqrt(a)*sqrt(b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - tanh(e + f*x)/(f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) + (a + b)*sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_496():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - I*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(I*e + I*f*x, b/a)/(a*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_497():
    f = coth(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = coth(e + f*x)/(a*f*sqrt(a + b*sinh(e + f*x)**2)) + 2*sqrt(a + b*sinh(e + f*x)**2)*tanh(e + f*x)/(a**2*f) - 2*sqrt(a + b*sinh(e + f*x)**2)*coth(e + f*x)/(a**2*f) - 2*sqrt(a + b*sinh(e + f*x)**2)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) + sqrt(a + b*sinh(e + f*x)**2)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(a**2*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_498():
    f = coth(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)
    F = -(a - b)*coth(e + f*x)*csch(e + f*x)**2/(a*b*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*coth(e + f*x)*csch(e + f*x)**2/(3*a**2*b*f) + sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*tanh(e + f*x)/(3*a**3*f) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*coth(e + f*x)/(3*a**3*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_499():
    f = tanh(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -sech(e + f*x)**4/(f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(4*a - 4*b)) + (8*a - b)*sech(e + f*x)**2/(8*f*(a - b)**2*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 + 24*a*b + 3*b**2)/(24*f*(a - b)**3*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 + 24*a*b + 3*b**2)/(8*f*(a - b)**4*sqrt(a + b*sinh(e + f*x)**2)) - (8*a**2 + 24*a*b + 3*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(8*f*(a - b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_500():
    f = tanh(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = sech(e + f*x)**2/(f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(2*a - 2*b)) + (2*a + 3*b)/(6*f*(a - b)**2*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (2*a + 3*b)/(2*f*(a - b)**3*sqrt(a + b*sinh(e + f*x)**2)) - (2*a + 3*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(2*f*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_501():
    f = tanh(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = 1/(f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) + 1/(f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a - b))/(f*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_502():
    f = coth(e + f*x)/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = 1/(3*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + 1/(a**2*f*sqrt(a + b*sinh(e + f*x)**2)) - atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_503():
    f = coth(e + f*x)**3/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -csch(e + f*x)**2/(2*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (2*a - 5*b)/(6*a**2*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (2*a - 5*b)/(2*a**3*f*sqrt(a + b*sinh(e + f*x)**2)) - (2*a - 5*b)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_504():
    f = coth(e + f*x)**5/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -csch(e + f*x)**4/(4*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - (8*a - 7*b)*csch(e + f*x)**2/(8*a**2*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 - 40*a*b + 35*b**2)/(24*a**3*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (8*a**2 - 40*a*b + 35*b**2)/(8*a**4*f*sqrt(a + b*sinh(e + f*x)**2)) - (8*a**2 - 40*a*b + 35*b**2)*atanh(sqrt(a + b*sinh(e + f*x)**2)/sqrt(a))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_505():
    f = tanh(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -8*sqrt(a)*sqrt(b)*(a + b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**4*sqrt(a + b*sinh(e + f*x)**2)) - b*(5*a + 3*b)*sinh(e + f*x)*cosh(e + f*x)/(3*f*(a - b)**3*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + tanh(e + f*x)*sech(e + f*x)**2/(f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 3*b)) - (4*a + 2*b)*tanh(e + f*x)/(3*f*(a - b)**2*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (a + 3*b)*sqrt(a + b*sinh(e + f*x)**2)*(3*a + b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_506():
    f = tanh(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*sinh(e + f*x)*cosh(e + f*x)/(3*f*(a - b)**2*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - tanh(e + f*x)/(f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a + 5*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)**3) - sqrt(b)*(7*a + b)*cosh(e + f*x)*elliptic_e(atan(sqrt(b)*sinh(e + f*x)/sqrt(a)), -a/b + 1)/(3*sqrt(a)*f*sqrt(a*cosh(e + f*x)**2/(a + b*sinh(e + f*x)**2))*(a - b)**3*sqrt(a + b*sinh(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_507():
    f = (a + b*sinh(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*sinh(e + f*x)*cosh(e + f*x)/(3*a*f*(a - b)*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + I*sqrt(1 + b*sinh(e + f*x)**2/a)*elliptic_f(I*e + I*f*x, b/a)/(3*a*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) - b*(4*a - 2*b)*sinh(e + f*x)*cosh(e + f*x)/(3*a**2*f*(a - b)**2*sqrt(a + b*sinh(e + f*x)**2)) - 2*I*sqrt(a + b*sinh(e + f*x)**2)*(2*a - b)*elliptic_e(I*e + I*f*x, b/a)/(3*a**2*f*sqrt(1 + b*sinh(e + f*x)**2/a)*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_508():
    f = coth(e + f*x)**2/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = coth(e + f*x)/(3*a*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) + (3*a - 4*b)*coth(e + f*x)/(3*a**2*f*(a - b)*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*tanh(e + f*x)/(3*a**3*f*(a - b)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*coth(e + f*x)/(3*a**3*f*(a - b)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 4*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b)) - sqrt(a + b*sinh(e + f*x)**2)*(7*a - 8*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**3*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_509():
    f = coth(e + f*x)**4/(a + b*sinh(e + f*x)**2)**(sympy.S(5)/2)
    F = -(a - b)*coth(e + f*x)*csch(e + f*x)**2/(3*a*b*f*(a + b*sinh(e + f*x)**2)**(sympy.S(3)/2)) - (2*a - 6*b)*coth(e + f*x)*csch(e + f*x)**2/(3*a**2*b*f*sqrt(a + b*sinh(e + f*x)**2)) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 8*b)*coth(e + f*x)*csch(e + f*x)**2/(3*a**3*b*f) + sqrt(a + b*sinh(e + f*x)**2)*(8*a - 16*b)*tanh(e + f*x)/(3*a**4*f) - sqrt(a + b*sinh(e + f*x)**2)*(8*a - 16*b)*coth(e + f*x)/(3*a**4*f) + sqrt(a + b*sinh(e + f*x)**2)*(3*a - 8*b)*elliptic_f(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**4*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a)) - sqrt(a + b*sinh(e + f*x)**2)*(8*a - 16*b)*elliptic_e(atan(sinh(e + f*x)), 1 - b/a)*sech(e + f*x)/(3*a**4*f*sqrt((a + b*sinh(e + f*x)**2)*sech(e + f*x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_510():
    f = (d*tanh(e + f*x))**m*(a + b*sinh(e + f*x)**2)**p
    F = (d*tanh(e + f*x))**(m + 1)*(a + b*sinh(e + f*x)**2)**p*(cosh(e + f*x)**2)**(m/2 + sympy.S.Half)*appellf1(m/2 + sympy.S.Half, -p, m/2 + sympy.S.Half, m/2 + sympy.S(3)/2, -b*sinh(e + f*x)**2/a, -sinh(e + f*x)**2)/(d*f*(1 + b*sinh(e + f*x)**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_511():
    f = (a + b*sinh(c + d*x)**2)**p*tanh(c + d*x)**3
    F = (a + b*sinh(c + d*x)**2)**(p + 1)*sech(c + d*x)**2/(d*(2*a - 2*b)) - (a - b*(p + 1))*(a + b*sinh(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sinh(c + d*x)**2)/(a - b))/(2*d*(a - b)**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_512():
    f = (a + b*sinh(c + d*x)**2)**p*tanh(c + d*x)
    F = -(a + b*sinh(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sinh(c + d*x)**2)/(a - b))/(d*(2*a - 2*b)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_513():
    f = (a + b*sinh(c + d*x)**2)**p*coth(c + d*x)
    F = -(a + b*sinh(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sinh(c + d*x)**2/a)/(2*a*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_514():
    f = (a + b*sinh(c + d*x)**2)**p*coth(c + d*x)**3
    F = -(a + b*sinh(c + d*x)**2)**(p + 1)*csch(c + d*x)**2/(2*a*d) - (a + b*p)*(a + b*sinh(c + d*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sinh(c + d*x)**2/a)/(2*a**2*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_515():
    f = (a + b*sinh(c + d*x)**2)**p*tanh(c + d*x)**4
    F = (a + b*sinh(c + d*x)**2)**p*sqrt(cosh(c + d*x)**2)*sinh(c + d*x)**4*tanh(c + d*x)*appellf1(sympy.S(5)/2, sympy.S(5)/2, -p, sympy.S(7)/2, -sinh(c + d*x)**2, -b*sinh(c + d*x)**2/a)/(5*d*(1 + b*sinh(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_516():
    f = (a + b*sinh(c + d*x)**2)**p*tanh(c + d*x)**2
    F = (a + b*sinh(c + d*x)**2)**p*sqrt(cosh(c + d*x)**2)*sinh(c + d*x)**2*tanh(c + d*x)*appellf1(sympy.S(3)/2, sympy.S(3)/2, -p, sympy.S(5)/2, -sinh(c + d*x)**2, -b*sinh(c + d*x)**2/a)/(3*d*(1 + b*sinh(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_517():
    f = (a + b*sinh(c + d*x)**2)**p*coth(c + d*x)**2
    F = -(a + b*sinh(c + d*x)**2)**p*sqrt(cosh(c + d*x)**2)*appellf1(sympy.S(-1)/2, sympy.S(-1)/2, -p, sympy.S.Half, -sinh(c + d*x)**2, -b*sinh(c + d*x)**2/a)*csch(c + d*x)*sech(c + d*x)/(d*(1 + b*sinh(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_518():
    f = (a + b*sinh(c + d*x)**2)**p*coth(c + d*x)**4
    F = -(a + b*sinh(c + d*x)**2)**p*sqrt(cosh(c + d*x)**2)*appellf1(sympy.S(-3)/2, sympy.S(-3)/2, -p, sympy.S(-1)/2, -sinh(c + d*x)**2, -b*sinh(c + d*x)**2/a)*csch(c + d*x)**3*sech(c + d*x)/(3*d*(1 + b*sinh(c + d*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_519():
    f = coth(x)**3/(a + b*sinh(x)**3)
    F = -log(a + b*sinh(x)**3)/(3*a) + log(sinh(x))/a - csch(x)**2/(2*a) - b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*sinh(x))/(3*a**(sympy.S(5)/3)) + b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sinh(x) + b**(sympy.S(2)/3)*sinh(x)**2)/(6*a**(sympy.S(5)/3)) + sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*sinh(x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_520():
    f = coth(x)/sqrt(a + b*sinh(x)**3)
    F = -2*atanh(sqrt(a + b*sinh(x)**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_521():
    f = sqrt(a + b*sinh(x)**3)*coth(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sinh(x)**3)/sqrt(a))/3 + 2*sqrt(a + b*sinh(x)**3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_522():
    f = coth(x)/sqrt(a + b*sinh(x)**n)
    F = -2*atanh(sqrt(a + b*sinh(x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_7_hyper_pow_m_a_plus_b_sinh_pow_n_pow_p_523():
    f = sqrt(a + b*sinh(x)**n)*coth(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sinh(x)**n)/sqrt(a))/n + 2*sqrt(a + b*sinh(x)**n)/n
    assert integrate(f, x) == F

