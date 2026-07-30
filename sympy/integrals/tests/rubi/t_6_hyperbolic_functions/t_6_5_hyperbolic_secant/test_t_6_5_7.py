"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.5 Hyperbolic secant/6.5.7 (d hyper)^m (a+b (c sech)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_1():
    f = (a + b*sech(c + d*x)**2)*sinh(c + d*x)**4
    F = a*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + b*tanh(c + d*x)/d + x*(3*a/8 - 3*b/2) - (5*a - 4*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_2():
    f = (a + b*sech(c + d*x)**2)*sinh(c + d*x)**3
    F = a*cosh(c + d*x)**3/(3*d) + b*sech(c + d*x)/d - (a - b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_3():
    f = (a + b*sech(c + d*x)**2)*sinh(c + d*x)**2
    F = a*sinh(c + d*x)*cosh(c + d*x)/(2*d) - b*tanh(c + d*x)/d + x*(-a/2 + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_4():
    f = (a + b*sech(c + d*x)**2)*sinh(c + d*x)
    F = a*cosh(c + d*x)/d - b*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_5():
    f = (a + b*sech(c + d*x)**2)*csch(c + d*x)
    F = b*sech(c + d*x)/d - (a + b)*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_6():
    f = (a + b*sech(c + d*x)**2)*csch(c + d*x)**2
    F = -b*tanh(c + d*x)/d - (a + b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_7():
    f = (a + b*sech(c + d*x)**2)*csch(c + d*x)**3
    F = -b*sech(c + d*x)/d - (a + b)*coth(c + d*x)*csch(c + d*x)/(2*d) + (a + 3*b)*atanh(cosh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_8():
    f = (a + b*sech(c + d*x)**2)*csch(c + d*x)**4
    F = b*tanh(c + d*x)/d - (a + b)*coth(c + d*x)**3/(3*d) + (a + 2*b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_9():
    f = (a + b*sech(c + d*x)**2)**2*sinh(c + d*x)**4
    F = a**2*sinh(c + d*x)**4*tanh(c + d*x)/(4*d) - a*(a - 8*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) - b**2*tanh(c + d*x)**3/(3*d) + x*(3*a**2/8 - 3*a*b + b**2) - (a**2 - 8*a*b + 4*b**2)*tanh(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_10():
    f = (a + b*sech(c + d*x)**2)**2*sinh(c + d*x)**3
    F = a**2*cosh(c + d*x)**3/(3*d) - a*(a - 2*b)*cosh(c + d*x)/d + b**2*sech(c + d*x)**3/(3*d) + b*(2*a - b)*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_11():
    f = (a + b*sech(c + d*x)**2)**2*sinh(c + d*x)**2
    F = a**2*sinh(c + d*x)**2*tanh(c + d*x)/(2*d) - a*x*(a - 4*b)/2 + a*(a - 4*b)*tanh(c + d*x)/(2*d) + b**2*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_12():
    f = (a + b*sech(c + d*x)**2)**2*sinh(c + d*x)
    F = a**2*cosh(c + d*x)/d - 2*a*b*sech(c + d*x)/d - b**2*sech(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_13():
    f = (a + b*sech(c + d*x)**2)**2*csch(c + d*x)
    F = b**2*sech(c + d*x)**3/(3*d) + b*(2*a + b)*sech(c + d*x)/d - (a + b)**2*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_14():
    f = (a + b*sech(c + d*x)**2)**2*csch(c + d*x)**2
    F = b**2*tanh(c + d*x)**3/(3*d) - 2*b*(a + b)*tanh(c + d*x)/d - (a + b)**2*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_15():
    f = (a + b*sech(c + d*x)**2)**2*csch(c + d*x)**3
    F = b**2*csch(c + d*x)**2*sech(c + d*x)**3/(3*d) - b*(6*a + 5*b)*sech(c + d*x)/(3*d) + (a + b)*(a + 5*b)*atanh(cosh(c + d*x))/(2*d) - (3*a**2 + 6*a*b + 5*b**2)*coth(c + d*x)*csch(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_16():
    f = (a + b*sech(c + d*x)**2)**2*csch(c + d*x)**4
    F = -b**2*tanh(c + d*x)**3/(3*d) + b*(2*a + 3*b)*tanh(c + d*x)/d - (a + b)**2*coth(c + d*x)**3/(3*d) + (a + b)*(a + 3*b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_17():
    f = (a + b*sech(c + d*x)**2)**3*sinh(c + d*x)**4
    F = 3*a*x*(a**2 - 12*a*b + 8*b**2)/8 - 3*a*(a**2 - 12*a*b + 8*b**2)*tanh(c + d*x)/(8*d) - b**2*(15*a - 48*b)*tanh(c + d*x)**5/(40*d) + b*(6*a**2 - 23*a*b - 8*b**2)*tanh(c + d*x)**3/(8*d) - (3*a - 6*b)*(a - b*tanh(c + d*x)**2 + b)**2*sinh(c + d*x)**2*tanh(c + d*x)/(8*d) + (a - b*tanh(c + d*x)**2 + b)**3*sinh(c + d*x)**3*cosh(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_18():
    f = (a + b*sech(c + d*x)**2)**3*sinh(c + d*x)**3
    F = a**3*cosh(c + d*x)**3/(3*d) - a**2*(a - 3*b)*cosh(c + d*x)/d + 3*a*b*(a - b)*sech(c + d*x)/d + b**3*sech(c + d*x)**5/(5*d) + b**2*(3*a - b)*sech(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_19():
    f = (a + b*sech(c + d*x)**2)**3*sinh(c + d*x)
    F = a**3*cosh(c + d*x)/d - 3*a**2*b*sech(c + d*x)/d - a*b**2*sech(c + d*x)**3/d - b**3*sech(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_20():
    f = (a + b*sech(c + d*x)**2)**3*csch(c + d*x)
    F = b**3*sech(c + d*x)**5/(5*d) + b**2*(3*a + b)*sech(c + d*x)**3/(3*d) + b*(3*a**2 + 3*a*b + b**2)*sech(c + d*x)/d - (a + b)**3*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_21():
    f = (a + b*sech(c + d*x)**2)**3*csch(c + d*x)**2
    F = -b**3*tanh(c + d*x)**5/(5*d) + b**2*(a + b)*tanh(c + d*x)**3/d - 3*b*(a + b)**2*tanh(c + d*x)/d - (a + b)**3*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_22():
    f = (a + b*sech(c + d*x)**2)**3*csch(c + d*x)**3
    F = -b**2*(5*a + 7*b)*sech(c + d*x)**5/(10*d) - b*(6*a**2 + 15*a*b + 7*b**2)*sech(c + d*x)**3/(6*d) + (a + b)**2*(a + 7*b)*atanh(cosh(c + d*x))/(2*d) - (a + b)**2*(a + 7*b)*sech(c + d*x)/(2*d) - (a + b)*(a*cosh(c + d*x)**2 + b)**2*csch(c + d*x)**2*sech(c + d*x)**5/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_23():
    f = (a + b*sech(c + d*x)**2)**3*csch(c + d*x)**4
    F = b**3*tanh(c + d*x)**5/(5*d) - b**2*(3*a + 4*b)*tanh(c + d*x)**3/(3*d) + 3*b*(a + b)*(a + 2*b)*tanh(c + d*x)/d - (a + b)**3*coth(c + d*x)**3/(3*d) + (a + b)**2*(a + 4*b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_24():
    f = sinh(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)**3/(4*a*d) - (5*a + 4*b)*sinh(c + d*x)*cosh(c + d*x)/(8*a**2*d) - sqrt(b)*(a + b)**(sympy.S(3)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a**3*d) + x*(3*a**2 + 12*a*b + 8*b**2)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_25():
    f = sinh(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = cosh(c + d*x)**3/(3*a*d) - (a + b)*cosh(c + d*x)/(a**2*d) + sqrt(b)*(a + b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_26():
    f = sinh(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d) + sqrt(b)*sqrt(a + b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a**2*d) - x*(a + 2*b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_27():
    f = sinh(c + d*x)/(a + b*sech(c + d*x)**2)
    F = cosh(c + d*x)/(a*d) - sqrt(b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_28():
    f = csch(c + d*x)/(a + b*sech(c + d*x)**2)
    F = -atanh(cosh(c + d*x))/(d*(a + b)) + sqrt(b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(sqrt(a)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_29():
    f = csch(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = sqrt(b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(d*(a + b)**(sympy.S(3)/2)) - coth(c + d*x)/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_30():
    f = csch(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = -sqrt(a)*sqrt(b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(d*(a + b)**2) + (a - b)*atanh(cosh(c + d*x))/(2*d*(a + b)**2) - coth(c + d*x)*csch(c + d*x)/(d*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_31():
    f = csch(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = -a*sqrt(b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(d*(a + b)**(sympy.S(5)/2)) + a*coth(c + d*x)/(d*(a + b)**2) - coth(c + d*x)**3/(d*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_32():
    f = sinh(c + d*x)**4/(a + b*sech(c + d*x)**2)**2
    F = sinh(c + d*x)*cosh(c + d*x)**3/(4*a*d*(a - b*tanh(c + d*x)**2 + b)) - (5*a + 6*b)*sinh(c + d*x)*cosh(c + d*x)/(8*a**2*d*(a - b*tanh(c + d*x)**2 + b)) - 3*b*(3*a + 4*b)*tanh(c + d*x)/(8*a**3*d*(a - b*tanh(c + d*x)**2 + b)) - 3*sqrt(b)*sqrt(a + b)*(a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**4*d) + x*(3*a**2 + 24*a*b + 24*b**2)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_33():
    f = sinh(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = cosh(c + d*x)**3/(3*a**2*d) - b*(a + b)*cosh(c + d*x)/(2*a**3*d*(a*cosh(c + d*x)**2 + b)) - (a + 2*b)*cosh(c + d*x)/(a**3*d) + sqrt(b)*(3*a + 5*b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(2*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_34():
    f = sinh(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d*(a - b*tanh(c + d*x)**2 + b)) + b*tanh(c + d*x)/(a**2*d*(a - b*tanh(c + d*x)**2 + b)) + sqrt(b)*(3*a + 4*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**3*d*sqrt(a + b)) - x*(a + 4*b)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_35():
    f = sinh(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = -cosh(c + d*x)**3/(2*a*d*(a*cosh(c + d*x)**2 + b)) + 3*cosh(c + d*x)/(2*a**2*d) - 3*sqrt(b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_36():
    f = csch(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = -atanh(cosh(c + d*x))/(d*(a + b)**2) - b*cosh(c + d*x)/(2*a*d*(a + b)*(a*cosh(c + d*x)**2 + b)) + sqrt(b)*(3*a + b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(2*a**(sympy.S(3)/2)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_37():
    f = csch(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = 3*sqrt(b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*d*(a + b)**(sympy.S(5)/2)) + coth(c + d*x)/(d*(2*a + 2*b)*(a - b*tanh(c + d*x)**2 + b)) - 3*coth(c + d*x)/(2*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_38():
    f = csch(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = (a - 3*b)*atanh(cosh(c + d*x))/(2*d*(a + b)**3) - (a - b)*cosh(c + d*x)/(2*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)) - coth(c + d*x)*csch(c + d*x)/(d*(2*a + 2*b)*(a*cosh(c + d*x)**2 + b)) - sqrt(b)*(3*a - b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(2*sqrt(a)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_39():
    f = csch(c + d*x)**4/(a + b*sech(c + d*x)**2)**2
    F = -a*b*tanh(c + d*x)/(2*d*(a + b)**3*(a - b*tanh(c + d*x)**2 + b)) - sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*d*(a + b)**(sympy.S(7)/2)) + (a - b)*coth(c + d*x)/(d*(a + b)**3) - coth(c + d*x)**3/(3*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_40():
    f = sinh(c + d*x)**4/(a + b*sech(c + d*x)**2)**3
    F = sinh(c + d*x)*cosh(c + d*x)**3/(4*a*d*(a - b*tanh(c + d*x)**2 + b)**2) - (5*a + 8*b)*sinh(c + d*x)*cosh(c + d*x)/(8*a**2*d*(a - b*tanh(c + d*x)**2 + b)**2) - b*(7*a + 12*b)*tanh(c + d*x)/(8*a**3*d*(a - b*tanh(c + d*x)**2 + b)**2) - 3*b*(a + 2*b)*tanh(c + d*x)/(2*a**4*d*(a - b*tanh(c + d*x)**2 + b)) - 3*sqrt(b)*(5*a**2 + 20*a*b + 16*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**5*d*sqrt(a + b)) + x*(3*a**2 + 36*a*b + 48*b**2)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_41():
    f = sinh(c + d*x)**3/(a + b*sech(c + d*x)**2)**3
    F = cosh(c + d*x)**3/(3*a**3*d) + b**2*(a + b)*cosh(c + d*x)/(4*a**4*d*(a*cosh(c + d*x)**2 + b)**2) - b*(9*a + 13*b)*cosh(c + d*x)/(8*a**4*d*(a*cosh(c + d*x)**2 + b)) - (a + 3*b)*cosh(c + d*x)/(a**4*d) + 5*sqrt(b)*(3*a + 7*b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(8*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_42():
    f = sinh(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d*(a - b*tanh(c + d*x)**2 + b)**2) + 3*b*tanh(c + d*x)/(4*a**2*d*(a - b*tanh(c + d*x)**2 + b)**2) + b*(11*a + 12*b)*tanh(c + d*x)/(8*a**3*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) + sqrt(b)*(15*a**2 + 40*a*b + 24*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**4*d*(a + b)**(sympy.S(3)/2)) - x*(a + 6*b)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_43():
    f = sinh(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = -cosh(c + d*x)**5/(4*a*d*(a*cosh(c + d*x)**2 + b)**2) - 5*cosh(c + d*x)**3/(8*a**2*d*(a*cosh(c + d*x)**2 + b)) + 15*cosh(c + d*x)/(8*a**3*d) - 15*sqrt(b)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(8*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_44():
    f = csch(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = -atanh(cosh(c + d*x))/(d*(a + b)**3) - b*cosh(c + d*x)**3/(4*a*d*(a + b)*(a*cosh(c + d*x)**2 + b)**2) - b*(7*a + 3*b)*cosh(c + d*x)/(8*a**2*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)) + sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(8*a**(sympy.S(5)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_45():
    f = csch(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = 15*sqrt(b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*d*(a + b)**(sympy.S(7)/2)) + coth(c + d*x)/(d*(4*a + 4*b)*(a - b*tanh(c + d*x)**2 + b)**2) + 5*coth(c + d*x)/(8*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) - 15*coth(c + d*x)/(8*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_46():
    f = csch(c + d*x)**3/(a + b*sech(c + d*x)**2)**3
    F = (a - 5*b)*atanh(cosh(c + d*x))/(2*d*(a + b)**4) - cosh(c + d*x)*coth(c + d*x)**2/(d*(2*a + 2*b)*(a*cosh(c + d*x)**2 + b)**2) + b*(2*a - b)*cosh(c + d*x)/(4*a*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)**2) - (4*a**2 - 9*a*b - b**2)*cosh(c + d*x)/(8*a*d*(a + b)**3*(a*cosh(c + d*x)**2 + b)) - sqrt(b)*(15*a**2 - 10*a*b - b**2)*atan(sqrt(a)*cosh(c + d*x)/sqrt(b))/(8*a**(sympy.S(3)/2)*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_47():
    f = csch(c + d*x)**4/(a + b*sech(c + d*x)**2)**3
    F = -a*b*tanh(c + d*x)/(4*d*(a + b)**3*(a - b*tanh(c + d*x)**2 + b)**2) - sqrt(b)*(15*a - 20*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*d*(a + b)**(sympy.S(9)/2)) - b*(7*a - 4*b)*tanh(c + d*x)/(8*d*(a + b)**4*(a - b*tanh(c + d*x)**2 + b)) + (a - 2*b)*coth(c + d*x)/(d*(a + b)**4) - coth(c + d*x)**3/(3*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_48():
    f = (a + b*sech(c + d*x)**2)*cosh(c + d*x)**4
    F = a*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + x*(3*a + 4*b)/8 + (3*a + 4*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_49():
    f = (a + b*sech(c + d*x)**2)*cosh(c + d*x)**3
    F = a*sinh(c + d*x)**3/(3*d) + (a + b)*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_50():
    f = (a + b*sech(c + d*x)**2)*cosh(c + d*x)**2
    F = a*sinh(c + d*x)*cosh(c + d*x)/(2*d) + x*(a + 2*b)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_51():
    f = (a + b*sech(c + d*x)**2)*cosh(c + d*x)
    F = a*sinh(c + d*x)/d + b*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_52():
    f = (a + b*sech(c + d*x)**2)*sech(c + d*x)
    F = b*tanh(c + d*x)*sech(c + d*x)/(2*d) + (2*a + b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_53():
    f = (a + b*sech(c + d*x)**2)*sech(c + d*x)**3
    F = b*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + (4*a + 3*b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (4*a + 3*b)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_54():
    f = (a + b*sech(c + d*x)**2)**2*cosh(c + d*x)**4
    F = 3*a*(a + 2*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) + a*(a - b*tanh(c + d*x)**2 + b)*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + x*(3*a**2/8 + a*b + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_55():
    f = (a + b*sech(c + d*x)**2)**2*cosh(c + d*x)**3
    F = a**2*sinh(c + d*x)**3/(3*d) + a*(a + 2*b)*sinh(c + d*x)/d + b**2*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_56():
    f = (a + b*sech(c + d*x)**2)**2*cosh(c + d*x)**2
    F = a**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + a*x*(a + 4*b)/2 + b**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_57():
    f = (a + b*sech(c + d*x)**2)**2*cosh(c + d*x)
    F = a**2*sinh(c + d*x)/d + b**2*tanh(c + d*x)*sech(c + d*x)/(2*d) + b*(4*a + b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_58():
    f = (a + b*sech(c + d*x)**2)**2*sech(c + d*x)
    F = 3*b*(2*a + b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + b*(a*sinh(c + d*x)**2 + a + b)*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + (8*a**2 + 8*a*b + 3*b**2)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_59():
    f = (a + b*sech(c + d*x)**2)**2*sech(c + d*x)**2
    F = b**2*tanh(c + d*x)**5/(5*d) - 2*b*(a + b)*tanh(c + d*x)**3/(3*d) + (a + b)**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_60():
    f = (a + b*sech(c + d*x)**2)**2*sech(c + d*x)**3
    F = b*(8*a + 5*b)*tanh(c + d*x)*sech(c + d*x)**3/(24*d) + b*(a*sinh(c + d*x)**2 + a + b)*tanh(c + d*x)*sech(c + d*x)**5/(6*d) + (8*a**2 + 12*a*b + 5*b**2)*tanh(c + d*x)*sech(c + d*x)/(16*d) + (8*a**2 + 12*a*b + 5*b**2)*atan(sinh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_61():
    f = (a + b*sech(c + d*x)**2)**2*sech(c + d*x)**4
    F = -b**2*tanh(c + d*x)**7/(7*d) + b*(2*a + 3*b)*tanh(c + d*x)**5/(5*d) + (a + b)**2*tanh(c + d*x)/d - (a + b)*(a + 3*b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_62():
    f = (a + b*sech(c + d*x)**2)**3*cosh(c + d*x)**4
    F = a**3*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + 3*a**2*(a + 4*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) + 3*a*x*(a**2 + 4*a*b + 8*b**2)/8 + b**3*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_63():
    f = (a + b*sech(c + d*x)**2)**3*cosh(c + d*x)**3
    F = a**3*sinh(c + d*x)**3/(3*d) + a**2*(a + 3*b)*sinh(c + d*x)/d + b**3*tanh(c + d*x)*sech(c + d*x)/(2*d) + b**2*(6*a + b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_64():
    f = (a + b*sech(c + d*x)**2)**3*cosh(c + d*x)**2
    F = a**3*sinh(c + d*x)*cosh(c + d*x)/(2*d) + a**2*x*(a + 6*b)/2 - b**3*tanh(c + d*x)**3/(3*d) + b**2*(3*a + b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_65():
    f = (a + b*sech(c + d*x)**2)**3*cosh(c + d*x)
    F = a**3*sinh(c + d*x)/d + b**3*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + 3*b**2*(4*a + b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + 3*b*(8*a**2 + 4*a*b + b**2)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_66():
    f = (a + b*sech(c + d*x)**2)**3*sech(c + d*x)
    F = 5*b*(2*a + b)*(a*sinh(c + d*x)**2 + a + b)*tanh(c + d*x)*sech(c + d*x)**3/(24*d) + b*(44*a**2 + 44*a*b + 15*b**2)*tanh(c + d*x)*sech(c + d*x)/(48*d) + b*(a*sinh(c + d*x)**2 + a + b)**2*tanh(c + d*x)*sech(c + d*x)**5/(6*d) + (2*a + b)*(8*a**2 + 8*a*b + 5*b**2)*atan(sinh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_67():
    f = (a + b*sech(c + d*x)**2)**3*sech(c + d*x)**2
    F = -b**3*tanh(c + d*x)**7/(7*d) + 3*b**2*(a + b)*tanh(c + d*x)**5/(5*d) - b*(a + b)**2*tanh(c + d*x)**3/d + (a + b)**3*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_68():
    f = (a + b*sech(c + d*x)**2)**3*sech(c + d*x)**3
    F = b*(12*a + 7*b)*(a*sinh(c + d*x)**2 + a + b)*tanh(c + d*x)*sech(c + d*x)**5/(48*d) + b*(72*a**2 + 92*a*b + 35*b**2)*tanh(c + d*x)*sech(c + d*x)**3/(192*d) + b*(a*sinh(c + d*x)**2 + a + b)**2*tanh(c + d*x)*sech(c + d*x)**7/(8*d) + (64*a**3 + 144*a**2*b + 120*a*b**2 + 35*b**3)*tanh(c + d*x)*sech(c + d*x)/(128*d) + (64*a**3 + 144*a**2*b + 120*a*b**2 + 35*b**3)*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_69():
    f = (a + b*sech(c + d*x)**2)**3*sech(c + d*x)**4
    F = b**3*tanh(c + d*x)**9/(9*d) - b**2*(3*a + 4*b)*tanh(c + d*x)**7/(7*d) + 3*b*(a + b)*(a + 2*b)*tanh(c + d*x)**5/(5*d) + (a + b)**3*tanh(c + d*x)/d - (a + b)**2*(a + 4*b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_70():
    f = cosh(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)**3/(4*a*d) + (3*a - 4*b)*sinh(c + d*x)*cosh(c + d*x)/(8*a**2*d) - b**(sympy.S(5)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a**3*d*sqrt(a + b)) + x*(3*a**2 - 4*a*b + 8*b**2)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_71():
    f = cosh(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)**3/(3*a*d) + (a - b)*sinh(c + d*x)/(a**2*d) + b**2*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(a**(sympy.S(5)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_72():
    f = cosh(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d) + b**(sympy.S(3)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a**2*d*sqrt(a + b)) + x*(a - 2*b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_73():
    f = cosh(c + d*x)/(a + b*sech(c + d*x)**2)
    F = sinh(c + d*x)/(a*d) - b*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(a**(sympy.S(3)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_74():
    f = sech(c + d*x)/(a + b*sech(c + d*x)**2)
    F = atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(sqrt(a)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_75():
    f = sech(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(sqrt(b)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_76():
    f = sech(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = -sqrt(a)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(b*d*sqrt(a + b)) + atan(sinh(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_77():
    f = sech(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = -a*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(b**(sympy.S(3)/2)*d*sqrt(a + b)) + tanh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_78():
    f = sech(c + d*x)**5/(a + b*sech(c + d*x)**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(b**2*d*sqrt(a + b)) + tanh(c + d*x)*sech(c + d*x)/(2*b*d) - (2*a - b)*atan(sinh(c + d*x))/(2*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_79():
    f = sech(c + d*x)**6/(a + b*sech(c + d*x)**2)
    F = a**2*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(b**(sympy.S(5)/2)*d*sqrt(a + b)) - tanh(c + d*x)**3/(3*b*d) - (a - b)*tanh(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_80():
    f = cosh(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = sinh(c + d*x)**3/(3*a**2*d) - b**3*sinh(c + d*x)/(2*a**3*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)) + (a - 2*b)*sinh(c + d*x)/(a**3*d) + b**2*(6*a + 5*b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*a**(sympy.S(7)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_81():
    f = cosh(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d*(a - b*tanh(c + d*x)**2 + b)) + b*(a + 2*b)*tanh(c + d*x)/(2*a**2*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) + b**(sympy.S(3)/2)*(5*a + 4*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**3*d*(a + b)**(sympy.S(3)/2)) + x*(a - 4*b)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_82():
    f = cosh(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = b**2*sinh(c + d*x)/(2*a**2*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)) + sinh(c + d*x)/(a**2*d) - b*(4*a + 3*b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_83():
    f = sech(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = -b*sinh(c + d*x)/(2*a*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)) + (2*a + b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_84():
    f = sech(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = tanh(c + d*x)/(d*(2*a + 2*b)*(a - b*tanh(c + d*x)**2 + b)) + atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*sqrt(b)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_85():
    f = sech(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = sinh(c + d*x)/(d*(2*a + 2*b)*(a*sinh(c + d*x)**2 + a + b)) + atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*sqrt(a)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_86():
    f = sech(c + d*x)**4/(a + b*sech(c + d*x)**2)**2
    F = -a*tanh(c + d*x)/(2*b*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) + (a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*b**(sympy.S(3)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_87():
    f = sech(c + d*x)**5/(a + b*sech(c + d*x)**2)**2
    F = -sqrt(a)*(2*a + 3*b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*b**2*d*(a + b)**(sympy.S(3)/2)) - a*sinh(c + d*x)/(2*b*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)) + atan(sinh(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_88():
    f = sech(c + d*x)**6/(a + b*sech(c + d*x)**2)**2
    F = a**2*tanh(c + d*x)/(2*b**2*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) - a*(3*a + 4*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*b**(sympy.S(5)/2)*d*(a + b)**(sympy.S(3)/2)) + tanh(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_89():
    f = sech(c + d*x)**7/(a + b*sech(c + d*x)**2)**2
    F = a**(sympy.S(3)/2)*(4*a + 5*b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(2*b**3*d*(a + b)**(sympy.S(3)/2)) + a*(2*a + b)*sinh(c + d*x)/(2*b**2*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)) + tanh(c + d*x)*sech(c + d*x)/(2*b*d*(a*sinh(c + d*x)**2 + a + b)) - (4*a - b)*atan(sinh(c + d*x))/(2*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_90():
    f = cosh(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = sinh(c + d*x)*cosh(c + d*x)/(2*a*d*(a - b*tanh(c + d*x)**2 + b)**2) + b*(2*a + 3*b)*tanh(c + d*x)/(4*a**2*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) + b*(a + 4*b)*(4*a + 3*b)*tanh(c + d*x)/(8*a**3*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) + b**(sympy.S(3)/2)*(35*a**2 + 56*a*b + 24*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**4*d*(a + b)**(sympy.S(5)/2)) + x*(a - 6*b)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_91():
    f = cosh(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = -b**3*sinh(c + d*x)/(4*a**3*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)**2) + 3*b**2*(4*a + 3*b)*sinh(c + d*x)/(8*a**3*d*(a + b)**2*(a*sinh(c + d*x)**2 + a + b)) + sinh(c + d*x)/(a**3*d) - 3*b*(4*(a + b)**2 + (2*a + b)**2)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(8*a**(sympy.S(7)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_92():
    f = sech(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = -b*sinh(c + d*x)*cosh(c + d*x)**2/(4*a*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)**2) - 3*b*(2*a + b)*sinh(c + d*x)/(8*a**2*d*(a + b)**2*(a*sinh(c + d*x)**2 + a + b)) + (8*a**2 + 8*a*b + 3*b**2)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(8*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_93():
    f = sech(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = tanh(c + d*x)/(d*(4*a + 4*b)*(a - b*tanh(c + d*x)**2 + b)**2) + 3*tanh(c + d*x)/(8*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) + 3*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*sqrt(b)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_94():
    f = sech(c + d*x)**3/(a + b*sech(c + d*x)**2)**3
    F = -b*sinh(c + d*x)/(4*a*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)**2) + (4*a + b)*sinh(c + d*x)/(8*a*d*(a + b)**2*(a*sinh(c + d*x)**2 + a + b)) + (4*a + b)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(8*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_95():
    f = sech(c + d*x)**4/(a + b*sech(c + d*x)**2)**3
    F = -a*tanh(c + d*x)/(4*b*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) + (a + 4*b)*tanh(c + d*x)/(8*b*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) + (a + 4*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*b**(sympy.S(3)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_96():
    f = sech(c + d*x)**5/(a + b*sech(c + d*x)**2)**3
    F = sinh(c + d*x)/(d*(4*a + 4*b)*(a*sinh(c + d*x)**2 + a + b)**2) + 3*sinh(c + d*x)/(8*d*(a + b)**2*(a*sinh(c + d*x)**2 + a + b)) + 3*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(8*sqrt(a)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_97():
    f = sech(c + d*x)**6/(a + b*sech(c + d*x)**2)**3
    F = -a*tanh(c + d*x)*sech(c + d*x)**2/(4*b*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) - 3*a*(a + 2*b)*tanh(c + d*x)/(8*b**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) + (3*a**2 + 8*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*b**(sympy.S(5)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_98():
    f = sech(c + d*x)**7/(a + b*sech(c + d*x)**2)**3
    F = -sqrt(a)*(8*a**2 + 20*a*b + 15*b**2)*atan(sqrt(a)*sinh(c + d*x)/sqrt(a + b))/(8*b**3*d*(a + b)**(sympy.S(5)/2)) - a*sinh(c + d*x)/(4*b*d*(a + b)*(a*sinh(c + d*x)**2 + a + b)**2) - a*(4*a + 7*b)*sinh(c + d*x)/(8*b**2*d*(a + b)**2*(a*sinh(c + d*x)**2 + a + b)) + atan(sinh(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_99():
    f = (a + b*sech(c + d*x)**2)*tanh(c + d*x)**4
    F = a*x - a*tanh(c + d*x)**3/(3*d) - a*tanh(c + d*x)/d + b*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_100():
    f = (a + b*sech(c + d*x)**2)*tanh(c + d*x)**3
    F = a*log(cosh(c + d*x))/d + b*sech(c + d*x)**4/(4*d) + (a - b)*sech(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_101():
    f = (a + b*sech(c + d*x)**2)*tanh(c + d*x)**2
    F = a*x - a*tanh(c + d*x)/d + b*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_102():
    f = (a + b*sech(c + d*x)**2)*tanh(c + d*x)
    F = a*log(cosh(c + d*x))/d - b*sech(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_103():
    f = a + b*sech(c + d*x)**2
    F = a*x + b*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_104():
    f = (a + b*sech(c + d*x)**2)*coth(c + d*x)
    F = -b*log(cosh(c + d*x))/d + (a + b)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_105():
    f = (a + b*sech(c + d*x)**2)*coth(c + d*x)**2
    F = a*x - (a + b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_106():
    f = (a + b*sech(c + d*x)**2)*coth(c + d*x)**3
    F = a*log(sinh(c + d*x))/d - (a + b)*csch(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_107():
    f = (a + b*sech(c + d*x)**2)*coth(c + d*x)**4
    F = a*x - a*coth(c + d*x)/d - (a + b)*coth(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_108():
    f = (a + b*sech(c + d*x)**2)*coth(c + d*x)**5
    F = a*log(sinh(c + d*x))/d - (a + b)*csch(c + d*x)**4/(4*d) - (2*a + b)*csch(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_109():
    f = (a + b*sech(c + d*x)**2)**2*tanh(c + d*x)**4
    F = a**2*x - a**2*tanh(c + d*x)**3/(3*d) - a**2*tanh(c + d*x)/d - b**2*tanh(c + d*x)**7/(7*d) + b*(2*a + b)*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_110():
    f = (a + b*sech(c + d*x)**2)**2*tanh(c + d*x)**3
    F = a**2*log(cosh(c + d*x))/d + a*(a - 2*b)*sech(c + d*x)**2/(2*d) + b**2*sech(c + d*x)**6/(6*d) + b*(2*a - b)*sech(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_111():
    f = (a + b*sech(c + d*x)**2)**2*tanh(c + d*x)**2
    F = a**2*x - a**2*tanh(c + d*x)/d - b**2*tanh(c + d*x)**5/(5*d) + b*(2*a + b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_112():
    f = (a + b*sech(c + d*x)**2)**2*tanh(c + d*x)
    F = a**2*log(cosh(c + d*x))/d - a*b*sech(c + d*x)**2/d - b**2*sech(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_113():
    f = (a + b*sech(c + d*x)**2)**2
    F = a**2*x - b**2*tanh(c + d*x)**3/(3*d) + b*(2*a + b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_114():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)
    F = b**2*sech(c + d*x)**2/(2*d) - b*(2*a + b)*log(cosh(c + d*x))/d + (a + b)**2*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_115():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**2
    F = a**2*x - b**2*tanh(c + d*x)/d - (a + b)**2*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_116():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**3
    F = b**2*log(cosh(c + d*x))/d - (a + b)**2*csch(c + d*x)**2/(2*d) + (a**2 - b**2)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_117():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**4
    F = a**2*x - (a + b)**2*coth(c + d*x)**3/(3*d) - (a**2 - b**2)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_118():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**5
    F = a**2*log(sinh(c + d*x))/d - a*(a + b)*csch(c + d*x)**2/d - (a + b)**2*csch(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_119():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**6
    F = a**2*x - a**2*coth(c + d*x)/d - (a + b)**2*coth(c + d*x)**5/(5*d) - (a**2 - b**2)*coth(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_120():
    f = (a + b*sech(c + d*x)**2)**2*coth(c + d*x)**7
    F = a**2*log(sinh(c + d*x))/d - a*(a + b)*csch(c + d*x)**2/d - (a + b)**2*csch(c + d*x)**4/(4*d) - (a*cosh(c + d*x)**2 + b)**3*csch(c + d*x)**6/(d*(6*a + 6*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_121():
    f = (a + b*sech(c + d*x)**2)**3*tanh(c + d*x)**4
    F = a**3*x - a**3*tanh(c + d*x)**3/(3*d) - a**3*tanh(c + d*x)/d + b**3*tanh(c + d*x)**9/(9*d) - b**2*(3*a + 2*b)*tanh(c + d*x)**7/(7*d) + b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_122():
    f = (a + b*sech(c + d*x)**2)**3*tanh(c + d*x)**3
    F = a**3*log(cosh(c + d*x))/d - 3*a**2*b*sech(c + d*x)**2/(2*d) - 3*a*b**2*sech(c + d*x)**4/(4*d) - b**3*sech(c + d*x)**6/(6*d) + (a*cosh(c + d*x)**2 + b)**4*sech(c + d*x)**8/(8*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_123():
    f = (a + b*sech(c + d*x)**2)**3*tanh(c + d*x)**2
    F = a**3*x - a**3*tanh(c + d*x)/d + b**3*tanh(c + d*x)**7/(7*d) - b**2*(3*a + 2*b)*tanh(c + d*x)**5/(5*d) + b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_124():
    f = (a + b*sech(c + d*x)**2)**3*tanh(c + d*x)
    F = a**3*log(cosh(c + d*x))/d - 3*a**2*b*sech(c + d*x)**2/(2*d) - 3*a*b**2*sech(c + d*x)**4/(4*d) - b**3*sech(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_125():
    f = (a + b*sech(c + d*x)**2)**3
    F = a**3*x + b**3*tanh(c + d*x)**5/(5*d) - b**2*(3*a + 2*b)*tanh(c + d*x)**3/(3*d) + b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_126():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)
    F = b**3*sech(c + d*x)**4/(4*d) + b**2*(3*a + b)*sech(c + d*x)**2/(2*d) - b*(3*a**2 + 3*a*b + b**2)*log(cosh(c + d*x))/d + (a + b)**3*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_127():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**2
    F = a**3*x + b**3*tanh(c + d*x)**3/(3*d) - b**2*(3*a + 2*b)*tanh(c + d*x)/d - (a + b)**3*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_128():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**3
    F = -b**3*sech(c + d*x)**2/(2*d) + b**2*(3*a + 2*b)*log(cosh(c + d*x))/d + (a - 2*b)*(a + b)**2*log(sinh(c + d*x))/d - (a + b)**3*csch(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_129():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**4
    F = a**3*x + b**3*tanh(c + d*x)/d - (a - 2*b)*(a + b)**2*coth(c + d*x)/d - (a + b)**3*coth(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_130():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**5
    F = -b**3*log(cosh(c + d*x))/d - (a + b)**3*csch(c + d*x)**4/(4*d) - (a + b)**2*(2*a - b)*csch(c + d*x)**2/(2*d) + (a**3 + b**3)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_131():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**6
    F = a**3*x - (a - 2*b)*(a + b)**2*coth(c + d*x)**3/(3*d) - (a + b)**3*coth(c + d*x)**5/(5*d) - (a**3 + b**3)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_132():
    f = (a + b*sech(c + d*x)**2)**3*coth(c + d*x)**7
    F = a**3*log(sinh(c + d*x))/d - 3*a**2*(a + b)*csch(c + d*x)**2/(2*d) - 3*a*(a + b)**2*csch(c + d*x)**4/(4*d) - (a + b)**3*csch(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_133():
    f = (a + b*sech(c + d*x)**2)**4
    F = a**4*x - b**4*tanh(c + d*x)**7/(7*d) + b**3*(4*a + 3*b)*tanh(c + d*x)**5/(5*d) - b**2*(6*a**2 + 8*a*b + 3*b**2)*tanh(c + d*x)**3/(3*d) + b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_134():
    f = (a + b*sech(c + d*x)**2)**5
    F = a**5*x + b**5*tanh(c + d*x)**9/(9*d) - b**4*(5*a + 4*b)*tanh(c + d*x)**7/(7*d) + b**3*(10*a**2 + 15*a*b + 6*b**2)*tanh(c + d*x)**5/(5*d) - b**2*(10*a**3 + 20*a**2*b + 15*a*b**2 + 4*b**3)*tanh(c + d*x)**3/(3*d) + b*(5*a**4 + 10*a**3*b + 10*a**2*b**2 + 5*a*b**3 + b**4)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_135():
    f = tanh(c + d*x)**5/(a + b*sech(c + d*x)**2)
    F = -sech(c + d*x)**2/(2*b*d) - (a + 2*b)*log(cosh(c + d*x))/(b**2*d) + (a + b)**2*log(a*cosh(c + d*x)**2 + b)/(2*a*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_136():
    f = tanh(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = tanh(c + d*x)/(b*d) + x/a - (a + b)**(sympy.S(3)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_137():
    f = tanh(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = -log(cosh(c + d*x))/(b*d) + (a + b)*log(a*cosh(c + d*x)**2 + b)/(2*a*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_138():
    f = tanh(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = x/a - sqrt(a + b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_139():
    f = tanh(c + d*x)/(a + b*sech(c + d*x)**2)
    F = log(a*cosh(c + d*x)**2 + b)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_140():
    f = 1/(a + b*sech(c + d*x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a*d*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_141():
    f = coth(c + d*x)/(a + b*sech(c + d*x)**2)
    F = log(sinh(c + d*x))/(d*(a + b)) + b*log(a*cosh(c + d*x)**2 + b)/(2*a*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_142():
    f = coth(c + d*x)**2/(a + b*sech(c + d*x)**2)
    F = -coth(c + d*x)/(d*(a + b)) - b**(sympy.S(3)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a*d*(a + b)**(sympy.S(3)/2)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_143():
    f = coth(c + d*x)**3/(a + b*sech(c + d*x)**2)
    F = -csch(c + d*x)**2/(d*(2*a + 2*b)) + (a + 2*b)*log(sinh(c + d*x))/(d*(a + b)**2) + b**2*log(a*cosh(c + d*x)**2 + b)/(2*a*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_144():
    f = coth(c + d*x)**4/(a + b*sech(c + d*x)**2)
    F = -coth(c + d*x)**3/(d*(3*a + 3*b)) - (a + 2*b)*coth(c + d*x)/(d*(a + b)**2) - b**(sympy.S(5)/2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(a*d*(a + b)**(sympy.S(5)/2)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_145():
    f = tanh(c + d*x)**5/(a + b*sech(c + d*x)**2)**2
    F = (-1/b**2 + a**(-2))*log(a*cosh(c + d*x)**2 + b)/(2*d) + log(cosh(c + d*x))/(b**2*d) + (a + b)**2/(2*a**2*b*d*(a*cosh(c + d*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_146():
    f = tanh(c + d*x)**4/(a + b*sech(c + d*x)**2)**2
    F = -(a + b)*tanh(c + d*x)/(2*a*b*d*(a - b*tanh(c + d*x)**2 + b)) + x/a**2 + (a - 2*b)*sqrt(a + b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**2*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_147():
    f = tanh(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = (a + b)/(2*a**2*d*(a*cosh(c + d*x)**2 + b)) + log(a*cosh(c + d*x)**2 + b)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_148():
    f = tanh(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = -tanh(c + d*x)/(2*a*d*(a - b*tanh(c + d*x)**2 + b)) + x/a**2 - (a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**2*sqrt(b)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_149():
    f = tanh(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = b/(2*a**2*d*(a*cosh(c + d*x)**2 + b)) + log(a*cosh(c + d*x)**2 + b)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_150():
    f = (a + b*sech(c + d*x)**2)**(-2)
    F = -b*tanh(c + d*x)/(2*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) - sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_151():
    f = coth(c + d*x)/(a + b*sech(c + d*x)**2)**2
    F = log(sinh(c + d*x))/(d*(a + b)**2) + b**2/(2*a**2*d*(a + b)*(a*cosh(c + d*x)**2 + b)) + b*(2*a + b)*log(a*cosh(c + d*x)**2 + b)/(2*a**2*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_152():
    f = coth(c + d*x)**2/(a + b*sech(c + d*x)**2)**2
    F = -b*coth(c + d*x)/(2*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) - (2*a - b)*coth(c + d*x)/(2*a*d*(a + b)**2) - b**(sympy.S(3)/2)*(5*a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(5)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_153():
    f = coth(c + d*x)**3/(a + b*sech(c + d*x)**2)**2
    F = -csch(c + d*x)**2/(2*d*(a + b)**2) + (a + 3*b)*log(sinh(c + d*x))/(d*(a + b)**3) + b**3/(2*a**2*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)) + b**2*(3*a + b)*log(a*cosh(c + d*x)**2 + b)/(2*a**2*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_154():
    f = coth(c + d*x)**4/(a + b*sech(c + d*x)**2)**2
    F = -b*coth(c + d*x)**3/(2*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) - (2*a - 3*b)*coth(c + d*x)**3/(6*a*d*(a + b)**2) - (2*a**2 + 6*a*b - b**2)*coth(c + d*x)/(2*a*d*(a + b)**3) - b**(sympy.S(5)/2)*(7*a + 2*b)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(7)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_155():
    f = tanh(c + d*x)**6/(a + b*sech(c + d*x)**2)**3
    F = -(a + b)*tanh(c + d*x)**3/(4*a*b*d*(a - b*tanh(c + d*x)**2 + b)**2) + (a + b)*(3*a - 4*b)*tanh(c + d*x)/(8*a**2*b**2*d*(a - b*tanh(c + d*x)**2 + b)) + x/a**3 - sqrt(a + b)*(3*a**2 - 4*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_156():
    f = tanh(c + d*x)**5/(a + b*sech(c + d*x)**2)**3
    F = -(a + b)**2/(4*a**3*d*(a*cosh(c + d*x)**2 + b)**2) + (a + b)/(a**3*d*(a*cosh(c + d*x)**2 + b)) + log(a*cosh(c + d*x)**2 + b)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_157():
    f = tanh(c + d*x)**4/(a + b*sech(c + d*x)**2)**3
    F = -(a + b)*tanh(c + d*x)/(4*a*b*d*(a - b*tanh(c + d*x)**2 + b)**2) + (a - 4*b)*tanh(c + d*x)/(8*a**2*b*d*(a - b*tanh(c + d*x)**2 + b)) + x/a**3 + (a**2 - 4*a*b - 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*b**(sympy.S(3)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_158():
    f = tanh(c + d*x)**3/(a + b*sech(c + d*x)**2)**3
    F = -b*(a + b)/(4*a**3*d*(a*cosh(c + d*x)**2 + b)**2) + (a + 2*b)/(2*a**3*d*(a*cosh(c + d*x)**2 + b)) + log(a*cosh(c + d*x)**2 + b)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_159():
    f = tanh(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = -tanh(c + d*x)/(4*a*d*(a - b*tanh(c + d*x)**2 + b)**2) - (3*a + 4*b)*tanh(c + d*x)/(8*a**2*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)) + x/a**3 - (3*a**2 + 12*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*sqrt(b)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_160():
    f = tanh(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = -b**2/(4*a**3*d*(a*cosh(c + d*x)**2 + b)**2) + b/(a**3*d*(a*cosh(c + d*x)**2 + b)) + log(a*cosh(c + d*x)**2 + b)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_161():
    f = (a + b*sech(c + d*x)**2)**(-3)
    F = -b*tanh(c + d*x)/(4*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) - b*(7*a + 4*b)*tanh(c + d*x)/(8*a**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) - sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*d*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_162():
    f = coth(c + d*x)/(a + b*sech(c + d*x)**2)**3
    F = log(sinh(c + d*x))/(d*(a + b)**3) - b**3/(4*a**3*d*(a + b)*(a*cosh(c + d*x)**2 + b)**2) + b**2*(3*a + 2*b)/(2*a**3*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)) + b*(3*a**2 + 3*a*b + b**2)*log(a*cosh(c + d*x)**2 + b)/(2*a**3*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_163():
    f = coth(c + d*x)**2/(a + b*sech(c + d*x)**2)**3
    F = -b*coth(c + d*x)/(4*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) - b*(9*a + 4*b)*coth(c + d*x)/(8*a**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) - (8*a**2 - 11*a*b - 4*b**2)*coth(c + d*x)/(8*a**2*d*(a + b)**3) - b**(sympy.S(3)/2)*(35*a**2 + 28*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*d*(a + b)**(sympy.S(7)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_164():
    f = coth(c + d*x)**3/(a + b*sech(c + d*x)**2)**3
    F = -csch(c + d*x)**2/(2*d*(a + b)**3) + (a + 4*b)*log(sinh(c + d*x))/(d*(a + b)**4) - b**4/(4*a**3*d*(a + b)**2*(a*cosh(c + d*x)**2 + b)**2) + b**3*(2*a + b)/(a**3*d*(a + b)**3*(a*cosh(c + d*x)**2 + b)) + b**2*(6*a**2 + 4*a*b + b**2)*log(a*cosh(c + d*x)**2 + b)/(2*a**3*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_165():
    f = coth(c + d*x)**4/(a + b*sech(c + d*x)**2)**3
    F = -b*coth(c + d*x)**3/(4*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**2) - b*(11*a + 4*b)*coth(c + d*x)**3/(8*a**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)) - (8*a**2 - 39*a*b - 12*b**2)*coth(c + d*x)**3/(24*a**2*d*(a + b)**3) - (8*a**3 + 32*a**2*b - 15*a*b**2 - 4*b**3)*coth(c + d*x)/(8*a**2*d*(a + b)**4) - b**(sympy.S(5)/2)*(63*a**2 + 36*a*b + 8*b**2)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(8*a**3*d*(a + b)**(sympy.S(9)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_166():
    f = (a + b*sech(c + d*x)**2)**(-4)
    F = -b*tanh(c + d*x)/(6*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**3) - b*(11*a + 6*b)*tanh(c + d*x)/(24*a**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)**2) - b*(19*a**2 + 22*a*b + 8*b**2)*tanh(c + d*x)/(16*a**3*d*(a + b)**3*(a - b*tanh(c + d*x)**2 + b)) - sqrt(b)*(35*a**3 + 70*a**2*b + 56*a*b**2 + 16*b**3)*atanh(sqrt(b)*tanh(c + d*x)/sqrt(a + b))/(16*a**4*d*(a + b)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_167():
    f = (1 - sech(x)**2)**(sympy.S(3)/2)
    F = -(tanh(x)**2)**(sympy.S(3)/2)*coth(x)/2 + sqrt(tanh(x)**2)*log(cosh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_168():
    f = sqrt(1 - sech(x)**2)
    F = sqrt(tanh(x)**2)*log(cosh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_169():
    f = 1/sqrt(1 - sech(x)**2)
    F = log(sinh(x))*tanh(x)/sqrt(tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_170():
    f = (sech(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-tanh(x)**2)*log(cosh(x))*coth(x) + sqrt(-tanh(x)**2)*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_171():
    f = sqrt(sech(x)**2 - 1)
    F = sqrt(-tanh(x)**2)*log(cosh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_172():
    f = 1/sqrt(sech(x)**2 - 1)
    F = log(sinh(x))*tanh(x)/sqrt(-tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_173():
    f = sqrt(a + b*sech(x)**2)*tanh(x)**5
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b*sech(x)**2) + (a + 2*b)*(a + b*sech(x)**2)**(sympy.S(3)/2)/(3*b**2) - (a + b*sech(x)**2)**(sympy.S(5)/2)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_174():
    f = sqrt(a + b*sech(x)**2)*tanh(x)**4
    F = sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - sqrt(a - b*tanh(x)**2 + b)*tanh(x)**3/4 + (a - 3*b)*sqrt(a - b*tanh(x)**2 + b)*tanh(x)/(8*b) - (a**2 + 6*a*b - 3*b**2)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/(8*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_175():
    f = sqrt(a + b*sech(x)**2)*tanh(x)**3
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b*sech(x)**2) + (a + b*sech(x)**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_176():
    f = sqrt(a + b*sech(x)**2)*tanh(x)**2
    F = sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - sqrt(a - b*tanh(x)**2 + b)*tanh(x)/2 - (a - b)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_177():
    f = sqrt(a + b*sech(x)**2)*tanh(x)
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b*sech(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_178():
    f = sqrt(a + b*sech(x)**2)
    F = sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) + sqrt(b)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_179():
    f = sqrt(a + b*sech(x)**2)*coth(x)
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_180():
    f = sqrt(a + b*sech(x)**2)*coth(x)**2
    F = sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - sqrt(a - b*tanh(x)**2 + b)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_181():
    f = sqrt(a + b*sech(x)**2)*coth(x)**3
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b*sech(x)**2)*coth(x)**2/2 - (2*a + b)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/(2*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_182():
    f = sqrt(a + b*sech(x)**2)*coth(x)**4
    F = sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - (3*a + 2*b)*sqrt(a - b*tanh(x)**2 + b)*coth(x)/(3*a + 3*b) - sqrt(a - b*tanh(x)**2 + b)*coth(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_183():
    f = sqrt(a + b*sech(x)**2)*coth(x)**5
    F = sqrt(a)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - sqrt(a + b*sech(x)**2)*(4*a + 3*b)*coth(x)**2/(8*a + 8*b) - sqrt(a + b*sech(x)**2)*coth(x)**4/4 - (8*a**2 + 12*a*b + 3*b**2)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/(8*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_184():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)*tanh(x)**3
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - a*sqrt(a + b*sech(x)**2) - (a + b*sech(x)**2)**(sympy.S(3)/2)/3 + (a + b*sech(x)**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_185():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)*tanh(x)**2
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) + b*sqrt(a - b*tanh(x)**2 + b)*tanh(x)**3/4 - (5*a + b)*sqrt(a - b*tanh(x)**2 + b)*tanh(x)/8 - (3*a**2 - 6*a*b - b**2)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_186():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)*tanh(x)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) - a*sqrt(a + b*sech(x)**2) - (a + b*sech(x)**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_187():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) + sqrt(b)*(3*a + b)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/2 + b*sqrt(a - b*tanh(x)**2 + b)*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_188():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)*coth(x)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a)) + b*sqrt(a + b*sech(x)**2) - (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_189():
    f = (a + b*sech(x)**2)**(sympy.S(3)/2)*coth(x)**2
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - b**(sympy.S(3)/2)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b)) - (a + b)*sqrt(a - b*tanh(x)**2 + b)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_190():
    f = (a + b*sech(c + d*x)**2)**(sympy.S(5)/2)
    F = a**(sympy.S(5)/2)*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a - b*tanh(c + d*x)**2 + b))/d + sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a - b*tanh(c + d*x)**2 + b))/(8*d) + b*(7*a + 3*b)*sqrt(a - b*tanh(c + d*x)**2 + b)*tanh(c + d*x)/(8*d) + b*(a - b*tanh(c + d*x)**2 + b)**(sympy.S(3)/2)*tanh(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_191():
    f = tanh(x)**5/sqrt(a + b*sech(x)**2)
    F = (a + 2*b)*sqrt(a + b*sech(x)**2)/b**2 - (a + b*sech(x)**2)**(sympy.S(3)/2)/(3*b**2) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_192():
    f = tanh(x)**4/sqrt(a + b*sech(x)**2)
    F = sqrt(a - b*tanh(x)**2 + b)*tanh(x)/(2*b) - (a + 3*b)*atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/(2*b**(sympy.S(3)/2)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_193():
    f = tanh(x)**3/sqrt(a + b*sech(x)**2)
    F = sqrt(a + b*sech(x)**2)/b + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_194():
    f = tanh(x)**2/sqrt(a + b*sech(x)**2)
    F = -atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/sqrt(b) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_195():
    f = tanh(x)/sqrt(a + b*sech(x)**2)
    F = atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_196():
    f = 1/sqrt(a + b*sech(x)**2)
    F = atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_197():
    f = coth(x)/sqrt(a + b*sech(x)**2)
    F = -atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/sqrt(a + b) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_198():
    f = coth(x)**2/sqrt(a + b*sech(x)**2)
    F = -sqrt(a - b*tanh(x)**2 + b)*coth(x)/(a + b) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_199():
    f = coth(x)**3/sqrt(a + b*sech(x)**2)
    F = -sqrt(a + b*sech(x)**2)*coth(x)**2/(2*a + 2*b) - (2*a + 3*b)*atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/(2*(a + b)**(sympy.S(3)/2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_200():
    f = tanh(x)**5/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -sqrt(a + b*sech(x)**2)/b**2 - (a + b)**2/(a*b**2*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_201():
    f = tanh(x)**4/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/b**(sympy.S(3)/2) - (a + b)*tanh(x)/(a*b*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_202():
    f = tanh(x)**3/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -(a + b)/(a*b*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_203():
    f = tanh(x)**2/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -tanh(x)/(a*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_204():
    f = tanh(x)/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -1/(a*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_205():
    f = (a + b*sech(x)**2)**(sympy.S(-3)/2)
    F = -b*tanh(x)/(a*(a + b)*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_206():
    f = coth(x)/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2) - b/(a*(a + b)*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_207():
    f = coth(x)**2/(a + b*sech(x)**2)**(sympy.S(3)/2)
    F = -b*coth(x)/(a*(a + b)*sqrt(a - b*tanh(x)**2 + b)) - (a - b)*sqrt(a - b*tanh(x)**2 + b)*coth(x)/(a*(a + b)**2) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_208():
    f = tanh(x)**6/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -(-1/b**2 + a**(-2))*tanh(x)/sqrt(a - b*tanh(x)**2 + b) - atan(sqrt(b)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/b**(sympy.S(5)/2) - (a + b)*tanh(x)**3/(3*a*b*(a - b*tanh(x)**2 + b)**(sympy.S(3)/2)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_209():
    f = tanh(x)**5/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -(-1/b**2 + a**(-2))/sqrt(a + b*sech(x)**2) - (a + b)**2/(3*a*b**2*(a + b*sech(x)**2)**(sympy.S(3)/2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_210():
    f = tanh(x)**4/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -(a + b)*tanh(x)/(3*a*b*(a - b*tanh(x)**2 + b)**(sympy.S(3)/2)) + (a - 3*b)*tanh(x)/(3*a**2*b*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_211():
    f = tanh(x)**3/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -(a + b)/(3*a*b*(a + b*sech(x)**2)**(sympy.S(3)/2)) - 1/(a**2*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_212():
    f = tanh(x)**2/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -tanh(x)/(3*a*(a - b*tanh(x)**2 + b)**(sympy.S(3)/2)) - (2*a + 3*b)*tanh(x)/(3*a**2*(a + b)*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_213():
    f = tanh(x)/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -1/(3*a*(a + b*sech(x)**2)**(sympy.S(3)/2)) - 1/(a**2*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_214():
    f = (a + b*sech(x)**2)**(sympy.S(-5)/2)
    F = -b*tanh(x)/(3*a*(a + b)*(a - b*tanh(x)**2 + b)**(sympy.S(3)/2)) - b*(5*a + 3*b)*tanh(x)/(3*a**2*(a + b)**2*sqrt(a - b*tanh(x)**2 + b)) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_215():
    f = coth(x)/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b*sech(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2) - b/(3*a*(a + b)*(a + b*sech(x)**2)**(sympy.S(3)/2)) - b*(2*a + b)/(a**2*(a + b)**2*sqrt(a + b*sech(x)**2)) + atanh(sqrt(a + b*sech(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_216():
    f = coth(x)**2/(a + b*sech(x)**2)**(sympy.S(5)/2)
    F = -b*coth(x)/(3*a*(a + b)*(a - b*tanh(x)**2 + b)**(sympy.S(3)/2)) - b*(7*a + 3*b)*coth(x)/(3*a**2*(a + b)**2*sqrt(a - b*tanh(x)**2 + b)) - (a - 3*b)*(3*a + b)*sqrt(a - b*tanh(x)**2 + b)*coth(x)/(3*a**2*(a + b)**3) + atanh(sqrt(a)*tanh(x)/sqrt(a - b*tanh(x)**2 + b))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_7_d_hyper_pow_m_a_plus_b_c_sech_pow_n_pow_p_217():
    f = (a + b*sech(c + d*x)**2)**(sympy.S(-7)/2)
    F = -b*tanh(c + d*x)/(5*a*d*(a + b)*(a - b*tanh(c + d*x)**2 + b)**(sympy.S(5)/2)) - b*(9*a + 5*b)*tanh(c + d*x)/(15*a**2*d*(a + b)**2*(a - b*tanh(c + d*x)**2 + b)**(sympy.S(3)/2)) - b*(33*a**2 + 40*a*b + 15*b**2)*tanh(c + d*x)/(15*a**3*d*(a + b)**3*sqrt(a - b*tanh(c + d*x)**2 + b)) + atanh(sqrt(a)*tanh(c + d*x)/sqrt(a - b*tanh(c + d*x)**2 + b))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F

