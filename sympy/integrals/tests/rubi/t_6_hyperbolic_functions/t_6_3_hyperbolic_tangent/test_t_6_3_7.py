"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.3 Hyperbolic tangent/6.3.7 (d hyper)^m (a+b (c tanh)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_1():
    f = (a + b*tanh(c + d*x)**2)*sinh(c + d*x)**4
    F = -b*tanh(c + d*x)/d + x*(3*a/8 + 15*b/8) + (a + b)*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) - (5*a + 9*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_2():
    f = (a + b*tanh(c + d*x)**2)*sinh(c + d*x)**3
    F = -b*sech(c + d*x)/d + (a + b)*cosh(c + d*x)**3/(3*d) - (a + 2*b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_3():
    f = (a + b*tanh(c + d*x)**2)*sinh(c + d*x)**2
    F = b*tanh(c + d*x)/d + x*(-a/2 - 3*b/2) + (a + b)*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_4():
    f = (a + b*tanh(c + d*x)**2)*sinh(c + d*x)
    F = b*sech(c + d*x)/d + (a + b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_5():
    f = (a + b*tanh(c + d*x)**2)*csch(c + d*x)
    F = -a*atanh(cosh(c + d*x))/d - b*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_6():
    f = (a + b*tanh(c + d*x)**2)*csch(c + d*x)**2
    F = -a*coth(c + d*x)/d + b*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_7():
    f = (a + b*tanh(c + d*x)**2)*csch(c + d*x)**3
    F = -a*coth(c + d*x)*csch(c + d*x)/(2*d) + b*sech(c + d*x)/d + (a - 2*b)*atanh(cosh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_8():
    f = (a + b*tanh(c + d*x)**2)*csch(c + d*x)**4
    F = -a*coth(c + d*x)**3/(3*d) - b*tanh(c + d*x)/d + (a - b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_9():
    f = (a + b*tanh(c + d*x)**2)**2*sinh(c + d*x)**4
    F = -b**2*tanh(c + d*x)**3/(3*d) + x*(3*a**2/8 + 15*a*b/4 + 35*b**2/8) + (a + b)**2*sinh(c + d*x)**4*tanh(c + d*x)/(4*d) - (a + b)*(a + 9*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d) - (a**2 + 10*a*b + 13*b**2)*tanh(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_10():
    f = (a + b*tanh(c + d*x)**2)**2*sinh(c + d*x)**3
    F = b**2*sech(c + d*x)**3/(3*d) - b*(2*a + 3*b)*sech(c + d*x)/d + (a + b)**2*cosh(c + d*x)**3/(3*d) - (a + b)*(a + 3*b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_11():
    f = (a + b*tanh(c + d*x)**2)**2*sinh(c + d*x)**2
    F = b**2*tanh(c + d*x)**3/(3*d) + x*(-a/2 - b/2)*(a + 5*b) + (a + b)**2*sinh(c + d*x)**2*tanh(c + d*x)/(2*d) + (a + b)*(a + 5*b)*tanh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_12():
    f = (a + b*tanh(c + d*x)**2)**2*sinh(c + d*x)
    F = -b**2*sech(c + d*x)**3/(3*d) + 2*b*(a + b)*sech(c + d*x)/d + (a + b)**2*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_13():
    f = (a + b*tanh(c + d*x)**2)**2*csch(c + d*x)
    F = -a**2*atanh(cosh(c + d*x))/d + b**2*sech(c + d*x)**3/(3*d) - b*(2*a + b)*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_14():
    f = (a + b*tanh(c + d*x)**2)**2*csch(c + d*x)**2
    F = -a**2*coth(c + d*x)/d + 2*a*b*tanh(c + d*x)/d + b**2*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_15():
    f = (a + b*tanh(c + d*x)**2)**2*csch(c + d*x)**3
    F = -a**2*csch(c + d*x)**2*sech(c + d*x)/(2*d) + a*(a - 4*b)*atanh(cosh(c + d*x))/(2*d) - a*(a - 4*b)*sech(c + d*x)/(2*d) - b**2*sech(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_16():
    f = (a + b*tanh(c + d*x)**2)**2*csch(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) + a*(a - 2*b)*coth(c + d*x)/d - b**2*tanh(c + d*x)**3/(3*d) - b*(2*a - b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_17():
    f = (a + b*tanh(c + d*x)**2)**3*sinh(c + d*x)**4
    F = -3*b**2*(5*a + 21*b)*tanh(c + d*x)**5/(40*d) - b*(6*a**2 + 35*a*b + 21*b**2)*tanh(c + d*x)**3/(8*d) + x*(3*a/8 + 3*b/8)*(a**2 + 14*a*b + 21*b**2) + (a + b*tanh(c + d*x)**2)**3*sinh(c + d*x)**3*cosh(c + d*x)/(4*d) - (a + b*tanh(c + d*x)**2)**2*(3*a + 9*b)*sinh(c + d*x)**2*tanh(c + d*x)/(8*d) - (3*a + 3*b)*(a**2 + 14*a*b + 21*b**2)*tanh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_18():
    f = (a + b*tanh(c + d*x)**2)**3*sinh(c + d*x)**3
    F = -b**3*sech(c + d*x)**5/(5*d) + b**2*(3*a + 4*b)*sech(c + d*x)**3/(3*d) - 3*b*(a + b)*(a + 2*b)*sech(c + d*x)/d + (a + b)**3*cosh(c + d*x)**3/(3*d) - (a + b)**2*(a + 4*b)*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_19():
    f = (a + b*tanh(c + d*x)**2)**3*sinh(c + d*x)
    F = b**3*sech(c + d*x)**5/(5*d) - b**2*(a + b)*sech(c + d*x)**3/d + 3*b*(a + b)**2*sech(c + d*x)/d + (a + b)**3*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_20():
    f = (a + b*tanh(c + d*x)**2)**3*csch(c + d*x)
    F = -a**3*atanh(cosh(c + d*x))/d - b**3*sech(c + d*x)**5/(5*d) + b**2*(3*a + 2*b)*sech(c + d*x)**3/(3*d) - b*(3*a**2 + 3*a*b + b**2)*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_21():
    f = (a + b*tanh(c + d*x)**2)**3*csch(c + d*x)**2
    F = -a**3*coth(c + d*x)/d + 3*a**2*b*tanh(c + d*x)/d + a*b**2*tanh(c + d*x)**3/d + b**3*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_22():
    f = (a + b*tanh(c + d*x)**2)**3*csch(c + d*x)**3
    F = a**2*(a - 6*b)*atanh(cosh(c + d*x))/(2*d) + b*(33*a - 2*b)*(a - b*sech(c + d*x)**2 + b)*sech(c + d*x)/(30*d) + 7*b*(a - b*sech(c + d*x)**2 + b)**2*sech(c + d*x)/(10*d) + b*(81*a**2 - 28*a*b - 4*b**2)*sech(c + d*x)/(30*d) - (a - b*sech(c + d*x)**2 + b)**3*coth(c + d*x)*csch(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_23():
    f = (a + b*tanh(c + d*x)**2)**3*csch(c + d*x)**4
    F = -a**3*coth(c + d*x)**3/(3*d) + a**2*(a - 3*b)*coth(c + d*x)/d - 3*a*b*(a - b)*tanh(c + d*x)/d - b**3*tanh(c + d*x)**5/(5*d) - b**2*(3*a - b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_24():
    f = sinh(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = a**(sympy.S(3)/2)*sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(d*(a + b)**3) + x*(3*a**2 - 6*a*b - b**2)/(8*(a + b)**3) + sinh(c + d*x)*cosh(c + d*x)**3/(d*(4*a + 4*b)) - (5*a + b)*sinh(c + d*x)*cosh(c + d*x)/(8*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_25():
    f = sinh(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = a*sqrt(b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(d*(a + b)**(sympy.S(5)/2)) - a*cosh(c + d*x)/(d*(a + b)**2) + cosh(c + d*x)**3/(d*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_26():
    f = sinh(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = -sqrt(a)*sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(d*(a + b)**2) - x*(a - b)/(2*(a + b)**2) + sinh(c + d*x)*cosh(c + d*x)/(d*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_27():
    f = sinh(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(d*(a + b)**(sympy.S(3)/2)) + cosh(c + d*x)/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_28():
    f = csch(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = sqrt(b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(a*d*sqrt(a + b)) - atanh(cosh(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_29():
    f = csch(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = -coth(c + d*x)/(a*d) - sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_30():
    f = csch(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d) - sqrt(b)*sqrt(a + b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(a**2*d) + (a + 2*b)*atanh(cosh(c + d*x))/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_31():
    f = csch(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = -coth(c + d*x)**3/(3*a*d) + (a + b)*coth(c + d*x)/(a**2*d) + sqrt(b)*(a + b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_32():
    f = sinh(c + d*x)**4/(a + b*tanh(c + d*x)**2)**2
    F = 3*sqrt(a)*sqrt(b)*(a - b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*d*(a + b)**4) + b*(9*a - 3*b)*tanh(c + d*x)/(8*d*(a + b)**3*(a + b*tanh(c + d*x)**2)) + x*(3*a**2 - 18*a*b + 3*b**2)/(8*(a + b)**4) + sinh(c + d*x)*cosh(c + d*x)**3/(d*(a + b*tanh(c + d*x)**2)*(4*a + 4*b)) - (5*a - b)*sinh(c + d*x)*cosh(c + d*x)/(8*d*(a + b)**2*(a + b*tanh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_33():
    f = sinh(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = a*b*sech(c + d*x)/(2*d*(a + b)**3*(a - b*sech(c + d*x)**2 + b)) + sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(2*d*(a + b)**(sympy.S(7)/2)) - (a - b)*cosh(c + d*x)/(d*(a + b)**3) + cosh(c + d*x)**3/(3*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_34():
    f = sinh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = -b*tanh(c + d*x)/(d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - x*(a - 3*b)/(2*(a + b)**3) + sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)*(2*a + 2*b)) - sqrt(b)*(3*a - b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*sqrt(a)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_35():
    f = sinh(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = -3*sqrt(b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(2*d*(a + b)**(sympy.S(5)/2)) - cosh(c + d*x)/(d*(2*a + 2*b)*(a - b*sech(c + d*x)**2 + b)) + 3*cosh(c + d*x)/(2*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_36():
    f = csch(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = b*sech(c + d*x)/(2*a*d*(a + b)*(a - b*sech(c + d*x)**2 + b)) + sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(2*a**2*d*(a + b)**(sympy.S(3)/2)) - atanh(cosh(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_37():
    f = csch(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = coth(c + d*x)/(2*a*d*(a + b*tanh(c + d*x)**2)) - 3*coth(c + d*x)/(2*a**2*d) - 3*sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_38():
    f = csch(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d*(a - b*sech(c + d*x)**2 + b)) - b*sech(c + d*x)/(a**2*d*(a - b*sech(c + d*x)**2 + b)) - sqrt(b)*(3*a + 4*b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(2*a**3*d*sqrt(a + b)) + (a + 4*b)*atanh(cosh(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_39():
    f = csch(c + d*x)**4/(a + b*tanh(c + d*x)**2)**2
    F = -coth(c + d*x)**3/(3*a**2*d) + b*(a + b)*tanh(c + d*x)/(2*a**3*d*(a + b*tanh(c + d*x)**2)) + (a + 2*b)*coth(c + d*x)/(a**3*d) + sqrt(b)*(3*a + 5*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_40():
    f = sinh(c + d*x)**4/(a + b*tanh(c + d*x)**2)**3
    F = b*(7*a - 5*b)*tanh(c + d*x)/(8*d*(a + b)**3*(a + b*tanh(c + d*x)**2)**2) + b*(3*a - 3*b)*tanh(c + d*x)/(2*d*(a + b)**4*(a + b*tanh(c + d*x)**2)) + x*(3*a**2 - 30*a*b + 15*b**2)/(8*(a + b)**5) + sinh(c + d*x)*cosh(c + d*x)**3/(d*(a + b*tanh(c + d*x)**2)**2*(4*a + 4*b)) - (5*a - 3*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d*(a + b)**2*(a + b*tanh(c + d*x)**2)**2) + 3*sqrt(b)*(5*a**2 - 10*a*b + b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*sqrt(a)*d*(a + b)**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_41():
    f = sinh(c + d*x)**3/(a + b*tanh(c + d*x)**2)**3
    F = a*b*sech(c + d*x)/(4*d*(a + b)**3*(a - b*sech(c + d*x)**2 + b)**2) + sqrt(b)*(15*a - 20*b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(8*d*(a + b)**(sympy.S(9)/2)) + b*(7*a - 4*b)*sech(c + d*x)/(8*d*(a + b)**4*(a - b*sech(c + d*x)**2 + b)) - (a - 2*b)*cosh(c + d*x)/(d*(a + b)**4) + cosh(c + d*x)**3/(3*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_42():
    f = sinh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = -3*b*tanh(c + d*x)/(4*d*(a + b)**2*(a + b*tanh(c + d*x)**2)**2) - x*(a - 5*b)/(2*(a + b)**4) + sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)**2*(2*a + 2*b)) - b*(11*a - b)*tanh(c + d*x)/(8*a*d*(a + b)**3*(a + b*tanh(c + d*x)**2)) - sqrt(b)*(15*a**2 - 10*a*b - b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_43():
    f = sinh(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = -15*sqrt(b)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(8*d*(a + b)**(sympy.S(7)/2)) - cosh(c + d*x)/(d*(4*a + 4*b)*(a - b*sech(c + d*x)**2 + b)**2) - 5*cosh(c + d*x)/(8*d*(a + b)**2*(a - b*sech(c + d*x)**2 + b)) + 15*cosh(c + d*x)/(8*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_44():
    f = csch(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = b*sech(c + d*x)/(4*a*d*(a + b)*(a - b*sech(c + d*x)**2 + b)**2) + b*(7*a + 4*b)*sech(c + d*x)/(8*a**2*d*(a + b)**2*(a - b*sech(c + d*x)**2 + b)) + sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(8*a**3*d*(a + b)**(sympy.S(5)/2)) - atanh(cosh(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_45():
    f = csch(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = coth(c + d*x)/(4*a*d*(a + b*tanh(c + d*x)**2)**2) + 5*coth(c + d*x)/(8*a**2*d*(a + b*tanh(c + d*x)**2)) - 15*coth(c + d*x)/(8*a**3*d) - 15*sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_46():
    f = csch(c + d*x)**3/(a + b*tanh(c + d*x)**2)**3
    F = -coth(c + d*x)*csch(c + d*x)/(2*a*d*(a - b*sech(c + d*x)**2 + b)**2) - 3*b*sech(c + d*x)/(4*a**2*d*(a - b*sech(c + d*x)**2 + b)**2) - b*(11*a + 12*b)*sech(c + d*x)/(8*a**3*d*(a + b)*(a - b*sech(c + d*x)**2 + b)) - sqrt(b)*(15*a**2 + 40*a*b + 24*b**2)*atanh(sqrt(b)*sech(c + d*x)/sqrt(a + b))/(8*a**4*d*(a + b)**(sympy.S(3)/2)) + (a + 6*b)*atanh(cosh(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_47():
    f = csch(c + d*x)**4/(a + b*tanh(c + d*x)**2)**3
    F = b*(a + b)*tanh(c + d*x)/(4*a**3*d*(a + b*tanh(c + d*x)**2)**2) - coth(c + d*x)**3/(3*a**3*d) + b*(7*a + 11*b)*tanh(c + d*x)/(8*a**4*d*(a + b*tanh(c + d*x)**2)) + (a + 3*b)*coth(c + d*x)/(a**4*d) + 5*sqrt(b)*(3*a + 7*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_48():
    f = (a + b*tanh(c + d*x)**3)*sinh(c + d*x)**4
    F = -3*a*tanh(c + d*x)/(8*d) - 3*b*tanh(c + d*x)**2/(2*d) - (a + 8*b*tanh(c + d*x))*sinh(c + d*x)**2*tanh(c + d*x)/(8*d) + (3*a - 24*b)*log(tanh(c + d*x) + 1)/(16*d) - (3*a + 24*b)*log(1 - tanh(c + d*x))/(16*d) + (a*tanh(c + d*x) + b)*sinh(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_49():
    f = (a + b*tanh(c + d*x)**3)*sinh(c + d*x)**3
    F = a*cosh(c + d*x)**3/(3*d) - a*cosh(c + d*x)/d - b*sinh(c + d*x)**3*tanh(c + d*x)**2/(2*d) + 5*b*sinh(c + d*x)**3/(6*d) - 5*b*sinh(c + d*x)/(2*d) + 5*b*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_50():
    f = (a + b*tanh(c + d*x)**3)*sinh(c + d*x)**2
    F = a*tanh(c + d*x)/(2*d) + b*tanh(c + d*x)**2/(2*d) - (a - 4*b)*log(tanh(c + d*x) + 1)/(4*d) + (a + 4*b)*log(1 - tanh(c + d*x))/(4*d) + (a*tanh(c + d*x) + b)*sinh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_51():
    f = (a + b*tanh(c + d*x)**3)*sinh(c + d*x)
    F = a*cosh(c + d*x)/d - b*sinh(c + d*x)*tanh(c + d*x)**2/(2*d) + 3*b*sinh(c + d*x)/(2*d) - 3*b*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_52():
    f = (a + b*tanh(c + d*x)**3)*csch(c + d*x)
    F = -a*atanh(cosh(c + d*x))/d - b*tanh(c + d*x)*sech(c + d*x)/(2*d) + b*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_53():
    f = (a + b*tanh(c + d*x)**3)*csch(c + d*x)**2
    F = -a*coth(c + d*x)/d + b*tanh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_54():
    f = (a + b*tanh(c + d*x)**3)*csch(c + d*x)**3
    F = -a*coth(c + d*x)*csch(c + d*x)/(2*d) + a*atanh(cosh(c + d*x))/(2*d) + b*tanh(c + d*x)*sech(c + d*x)/(2*d) + b*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_55():
    f = (a + b*tanh(c + d*x)**3)*csch(c + d*x)**4
    F = -a*coth(c + d*x)**3/(3*d) + a*coth(c + d*x)/d + b*log(tanh(c + d*x))/d - b*tanh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_56():
    f = (a + b*tanh(c + d*x)**3)**2*sinh(c + d*x)**3
    F = a**2*cosh(c + d*x)**3/(3*d) - a**2*cosh(c + d*x)/d - a*b*sinh(c + d*x)**3*tanh(c + d*x)**2/d + 5*a*b*sinh(c + d*x)**3/(3*d) - 5*a*b*sinh(c + d*x)/d + 5*a*b*atan(sinh(c + d*x))/d + b**2*cosh(c + d*x)**3/(3*d) - 4*b**2*cosh(c + d*x)/d - b**2*sech(c + d*x)**5/(5*d) + 4*b**2*sech(c + d*x)**3/(3*d) - 6*b**2*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_57():
    f = (a + b*tanh(c + d*x)**3)**2*sinh(c + d*x)
    F = a**2*cosh(c + d*x)/d - a*b*sinh(c + d*x)*tanh(c + d*x)**2/d + 3*a*b*sinh(c + d*x)/d - 3*a*b*atan(sinh(c + d*x))/d + b**2*cosh(c + d*x)/d + b**2*sech(c + d*x)**5/(5*d) - b**2*sech(c + d*x)**3/d + 3*b**2*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_58():
    f = (a + b*tanh(c + d*x)**3)**2*csch(c + d*x)
    F = -a**2*atanh(cosh(c + d*x))/d - a*b*tanh(c + d*x)*sech(c + d*x)/d + a*b*atan(sinh(c + d*x))/d - b**2*sech(c + d*x)**5/(5*d) + 2*b**2*sech(c + d*x)**3/(3*d) - b**2*sech(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_59():
    f = (a + b*tanh(c + d*x)**3)**2*csch(c + d*x)**2
    F = -a**2*coth(c + d*x)/d + a*b*tanh(c + d*x)**2/d + b**2*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_60():
    f = (a + b*tanh(c + d*x)**3)**2*csch(c + d*x)**3
    F = -a**2*coth(c + d*x)*csch(c + d*x)/(2*d) + a**2*atanh(cosh(c + d*x))/(2*d) + a*b*tanh(c + d*x)*sech(c + d*x)/d + a*b*atan(sinh(c + d*x))/d + b**2*sech(c + d*x)**5/(5*d) - b**2*sech(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_61():
    f = (a + b*tanh(c + d*x)**3)**2*csch(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) + a**2*coth(c + d*x)/d + 2*a*b*log(tanh(c + d*x))/d - a*b*tanh(c + d*x)**2/d - b**2*tanh(c + d*x)**5/(5*d) + b**2*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_62():
    f = (a + b*tanh(c + d*x)**3)**3*sinh(c + d*x)**3
    F = a**3*cosh(c + d*x)**3/(3*d) - a**3*cosh(c + d*x)/d - 3*a**2*b*sinh(c + d*x)**3*tanh(c + d*x)**2/(2*d) + 5*a**2*b*sinh(c + d*x)**3/(2*d) - 15*a**2*b*sinh(c + d*x)/(2*d) + 15*a**2*b*atan(sinh(c + d*x))/(2*d) + a*b**2*cosh(c + d*x)**3/d - 12*a*b**2*cosh(c + d*x)/d - 3*a*b**2*sech(c + d*x)**5/(5*d) + 4*a*b**2*sech(c + d*x)**3/d - 18*a*b**2*sech(c + d*x)/d - b**3*sinh(c + d*x)**3*tanh(c + d*x)**8/(8*d) - 11*b**3*sinh(c + d*x)**3*tanh(c + d*x)**6/(48*d) - 33*b**3*sinh(c + d*x)**3*tanh(c + d*x)**4/(64*d) - 231*b**3*sinh(c + d*x)**3*tanh(c + d*x)**2/(128*d) + 385*b**3*sinh(c + d*x)**3/(128*d) - 1155*b**3*sinh(c + d*x)/(128*d) + 1155*b**3*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_63():
    f = (a + b*tanh(c + d*x)**3)**3*sinh(c + d*x)
    F = a**3*cosh(c + d*x)/d - 3*a**2*b*sinh(c + d*x)*tanh(c + d*x)**2/(2*d) + 9*a**2*b*sinh(c + d*x)/(2*d) - 9*a**2*b*atan(sinh(c + d*x))/(2*d) + 3*a*b**2*cosh(c + d*x)/d + 3*a*b**2*sech(c + d*x)**5/(5*d) - 3*a*b**2*sech(c + d*x)**3/d + 9*a*b**2*sech(c + d*x)/d - b**3*sinh(c + d*x)*tanh(c + d*x)**8/(8*d) - 3*b**3*sinh(c + d*x)*tanh(c + d*x)**6/(16*d) - 21*b**3*sinh(c + d*x)*tanh(c + d*x)**4/(64*d) - 105*b**3*sinh(c + d*x)*tanh(c + d*x)**2/(128*d) + 315*b**3*sinh(c + d*x)/(128*d) - 315*b**3*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_64():
    f = (a + b*tanh(c + d*x)**3)**3*csch(c + d*x)
    F = -a**3*atanh(cosh(c + d*x))/d - 3*a**2*b*tanh(c + d*x)*sech(c + d*x)/(2*d) + 3*a**2*b*atan(sinh(c + d*x))/(2*d) - 3*a*b**2*sech(c + d*x)**5/(5*d) + 2*a*b**2*sech(c + d*x)**3/d - 3*a*b**2*sech(c + d*x)/d - b**3*tanh(c + d*x)**7*sech(c + d*x)/(8*d) - 7*b**3*tanh(c + d*x)**5*sech(c + d*x)/(48*d) - 35*b**3*tanh(c + d*x)**3*sech(c + d*x)/(192*d) - 35*b**3*tanh(c + d*x)*sech(c + d*x)/(128*d) + 35*b**3*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_65():
    f = (a + b*tanh(c + d*x)**3)**3*csch(c + d*x)**2
    F = -a**3*coth(c + d*x)/d + 3*a**2*b*tanh(c + d*x)**2/(2*d) + 3*a*b**2*tanh(c + d*x)**5/(5*d) + b**3*tanh(c + d*x)**8/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_66():
    f = (a + b*tanh(c + d*x)**3)**3*csch(c + d*x)**3
    F = -a**3*coth(c + d*x)*csch(c + d*x)/(2*d) + a**3*atanh(cosh(c + d*x))/(2*d) + 3*a**2*b*tanh(c + d*x)*sech(c + d*x)/(2*d) + 3*a**2*b*atan(sinh(c + d*x))/(2*d) + 3*a*b**2*sech(c + d*x)**5/(5*d) - a*b**2*sech(c + d*x)**3/d - b**3*tanh(c + d*x)**5*sech(c + d*x)**3/(8*d) - 5*b**3*tanh(c + d*x)**3*sech(c + d*x)**3/(48*d) - 5*b**3*tanh(c + d*x)*sech(c + d*x)**3/(64*d) + 5*b**3*tanh(c + d*x)*sech(c + d*x)/(128*d) + 5*b**3*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_67():
    f = (a + b*tanh(c + d*x)**3)**3*csch(c + d*x)**4
    F = -a**3*coth(c + d*x)**3/(3*d) + a**3*coth(c + d*x)/d + 3*a**2*b*log(tanh(c + d*x))/d - 3*a**2*b*tanh(c + d*x)**2/(2*d) - 3*a*b**2*tanh(c + d*x)**5/(5*d) + a*b**2*tanh(c + d*x)**3/d - b**3*tanh(c + d*x)**8/(8*d) + b**3*tanh(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_68():
    f = sinh(c + d*x)**4/(a + b*tanh(c + d*x)**3)
    F = -sqrt(3)*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(3*a**(sympy.S(4)/3)*b**(sympy.S(2)/3) + a**2 - b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*d*(a**(sympy.S(4)/3) + a**(sympy.S(2)/3)*b**(sympy.S(2)/3) + b**(sympy.S(4)/3))**3) - a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*(2*a**2 + b**2) + a**4 + 7*a**2*b**2 + b**4)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tanh(c + d*x))/(3*d*(a**2 - b**2)**3) + a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*(2*a**2 + b**2) + a**4 + 7*a**2*b**2 + b**4)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tanh(c + d*x) + b**(sympy.S(2)/3)*tanh(c + d*x)**2)/(6*d*(a**2 - b**2)**3) - a**2*b*(a**2 + 2*b**2)*log(a + b*tanh(c + d*x)**3)/(d*(a**2 - b**2)**3) - 3*a*(a - 5*b)*log(1 - tanh(c + d*x))/(16*d*(a + b)**3) + 3*a*(a + 5*b)*log(tanh(c + d*x) + 1)/(16*d*(a - b)**3) - 1/(d*(16*a - 16*b)*(tanh(c + d*x) + 1)**2) + (5*a + b)/(16*d*(a - b)**2*(tanh(c + d*x) + 1)) - (5*a - b)/(16*d*(1 - tanh(c + d*x))*(a + b)**2) + 1/(d*(1 - tanh(c + d*x))**2*(16*a + 16*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_69():
    f = sinh(c + d*x)**3/(a + b*tanh(c + d*x)**3)
    F = sympy.I * sympy.Function('Unintegrable')((((Integer(-1) * sympy.I) * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(3))) * ((Symbol('a') + (Symbol('b') * (sympy.tanh((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_70():
    f = sinh(c + d*x)**2/(a + b*tanh(c + d*x)**3)
    F = sqrt(3)*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(-3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*d*(a**2 - b**2)**2) + a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tanh(c + d*x))/(3*d*(a**2 - b**2)**2) - a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tanh(c + d*x) + b**(sympy.S(2)/3)*tanh(c + d*x)**2)/(6*d*(a**2 - b**2)**2) + b*(2*a**2 + b**2)*log(a + b*tanh(c + d*x)**3)/(3*d*(a**2 - b**2)**2) + (a - 2*b)*log(1 - tanh(c + d*x))/(4*d*(a + b)**2) - 1/(d*(4*a - 4*b)*(tanh(c + d*x) + 1)) - (a + 2*b)*log(tanh(c + d*x) + 1)/(4*d*(a - b)**2) + 1/(d*(1 - tanh(c + d*x))*(4*a + 4*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_71():
    f = sinh(c + d*x)/(a + b*tanh(c + d*x)**3)
    F = (Integer(-1) * sympy.I) * sympy.Function('Unintegrable')(((sympy.I * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * (sympy.tanh((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_72():
    f = csch(c + d*x)/(a + b*tanh(c + d*x)**3)
    F = sympy.I * sympy.Function('Unintegrable')((((Integer(-1) * sympy.I) * sympy.csch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * (sympy.tanh((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_73():
    f = csch(c + d*x)**2/(a + b*tanh(c + d*x)**3)
    F = -coth(c + d*x)/(a*d) + b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(4)/3)*d) - b**(sympy.S(1)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tanh(c + d*x) + b**(sympy.S(2)/3)*tanh(c + d*x)**2)/(6*a**(sympy.S(4)/3)*d) + sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_74():
    f = csch(c + d*x)**3/(a + b*tanh(c + d*x)**3)
    F = (Integer(-1) * sympy.I) * sympy.Function('Unintegrable')(((sympy.I * (sympy.csch((Symbol('c') + (Symbol('d') * x))))**(Integer(3))) * ((Symbol('a') + (Symbol('b') * (sympy.tanh((Symbol('c') + (Symbol('d') * x))))**(Integer(3)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_75():
    f = csch(c + d*x)**4/(a + b*tanh(c + d*x)**3)
    F = -coth(c + d*x)**3/(3*a*d) + coth(c + d*x)/(a*d) + b*log(a + b*tanh(c + d*x)**3)/(3*a**2*d) - b*log(tanh(c + d*x))/(a**2*d) - b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(4)/3)*d) + b**(sympy.S(1)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*tanh(c + d*x) + b**(sympy.S(2)/3)*tanh(c + d*x)**2)/(6*a**(sympy.S(4)/3)*d) - sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*tanh(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_76():
    f = (a + b*tanh(c + d*x)**2)*cosh(c + d*x)**4
    F = x*(3*a/8 - b/8) + (a + b)*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + (3*a - b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_77():
    f = (a + b*tanh(c + d*x)**2)*cosh(c + d*x)**3
    F = a*sinh(c + d*x)/d + (a + b)*sinh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_78():
    f = (a + b*tanh(c + d*x)**2)*cosh(c + d*x)**2
    F = x*(a/2 - b/2) + (a + b)*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_79():
    f = (a + b*tanh(c + d*x)**2)*cosh(c + d*x)
    F = -b*atan(sinh(c + d*x))/d + (a + b)*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_80():
    f = (a + b*tanh(c + d*x)**2)*sech(c + d*x)
    F = -b*tanh(c + d*x)*sech(c + d*x)/(2*d) + (2*a + b)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_81():
    f = (a + b*tanh(c + d*x)**2)*sech(c + d*x)**2
    F = a*tanh(c + d*x)/d + b*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_82():
    f = (a + b*tanh(c + d*x)**2)*sech(c + d*x)**3
    F = -b*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + (4*a + b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (4*a + b)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_83():
    f = (a + b*tanh(c + d*x)**2)*sech(c + d*x)**4
    F = a*tanh(c + d*x)/d - b*tanh(c + d*x)**5/(5*d) - (a - b)*tanh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_84():
    f = (a + b*tanh(c + d*x)**2)**2*cosh(c + d*x)**4
    F = x*(3*a**2/8 - a*b/4 + 3*b**2/8) + (a + b)*(a + b*tanh(c + d*x)**2)*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + (3*a**2 - 3*b**2)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_85():
    f = (a + b*tanh(c + d*x)**2)**2*cosh(c + d*x)**3
    F = b**2*atan(sinh(c + d*x))/d + (a + b)**2*sinh(c + d*x)**3/(3*d) + (a**2 - b**2)*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_86():
    f = (a + b*tanh(c + d*x)**2)**2*cosh(c + d*x)**2
    F = b**2*tanh(c + d*x)/d + x*(a/2 - 3*b/2)*(a + b) + (a + b)**2*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_87():
    f = (a + b*tanh(c + d*x)**2)**2*cosh(c + d*x)
    F = b**2*tanh(c + d*x)*sech(c + d*x)/(2*d) - b*(4*a + 3*b)*atan(sinh(c + d*x))/(2*d) + (a + b)**2*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_88():
    f = (a + b*tanh(c + d*x)**2)**2*sech(c + d*x)
    F = -b*(a + (a + b)*sinh(c + d*x)**2)*tanh(c + d*x)*sech(c + d*x)**3/(4*d) - 3*b*(2*a + b)*tanh(c + d*x)*sech(c + d*x)/(8*d) + (8*a**2 + 8*a*b + 3*b**2)*atan(sinh(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_89():
    f = (a + b*tanh(c + d*x)**2)**2*sech(c + d*x)**2
    F = a**2*tanh(c + d*x)/d + 2*a*b*tanh(c + d*x)**3/(3*d) + b**2*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_90():
    f = (a + b*tanh(c + d*x)**2)**2*sech(c + d*x)**3
    F = -b*(a + (a + b)*sinh(c + d*x)**2)*tanh(c + d*x)*sech(c + d*x)**5/(6*d) - b*(8*a + 3*b)*tanh(c + d*x)*sech(c + d*x)**3/(24*d) + (8*a**2 + 4*a*b + b**2)*tanh(c + d*x)*sech(c + d*x)/(16*d) + (8*a**2 + 4*a*b + b**2)*atan(sinh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_91():
    f = (a + b*tanh(c + d*x)**2)**2*sech(c + d*x)**4
    F = a**2*tanh(c + d*x)/d - a*(a - 2*b)*tanh(c + d*x)**3/(3*d) - b**2*tanh(c + d*x)**7/(7*d) - b*(2*a - b)*tanh(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_92():
    f = (a + b*tanh(c + d*x)**2)**3*cosh(c + d*x)**4
    F = -b**3*tanh(c + d*x)/d + x*(3*a/8 + 3*b/8)*(a**2 - 2*a*b + 5*b**2) + (a + b)**3*sinh(c + d*x)*cosh(c + d*x)**3/(4*d) + (a + b)**2*(3*a - 9*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_93():
    f = (a + b*tanh(c + d*x)**2)**3*cosh(c + d*x)**3
    F = -b**3*tanh(c + d*x)*sech(c + d*x)/(2*d) + b**2*(6*a + 5*b)*atan(sinh(c + d*x))/(2*d) + (a - 2*b)*(a + b)**2*sinh(c + d*x)/d + (a + b)**3*sinh(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_94():
    f = (a + b*tanh(c + d*x)**2)**3*cosh(c + d*x)**2
    F = b**3*tanh(c + d*x)**3/(3*d) + b**2*(3*a + 2*b)*tanh(c + d*x)/d + x*(a/2 - 5*b/2)*(a + b)**2 + (a + b)**3*sinh(c + d*x)*cosh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_95():
    f = (a + b*tanh(c + d*x)**2)**3*cosh(c + d*x)
    F = -b**3*tanh(c + d*x)*sech(c + d*x)**3/(4*d) + 3*b**2*(4*a + 3*b)*tanh(c + d*x)*sech(c + d*x)/(8*d) - 3*b*(4*(a + b)**2 + (2*a + b)**2)*atan(sinh(c + d*x))/(8*d) + (a + b)**3*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_96():
    f = (a + b*tanh(c + d*x)**2)**3*sech(c + d*x)
    F = -b*(a + (a + b)*sinh(c + d*x)**2)**2*tanh(c + d*x)*sech(c + d*x)**5/(6*d) - 5*b*(a + (a + b)*sinh(c + d*x)**2)*(2*a + b)*tanh(c + d*x)*sech(c + d*x)**3/(24*d) - b*(44*a**2 + 44*a*b + 15*b**2)*tanh(c + d*x)*sech(c + d*x)/(48*d) + (2*a + b)*(8*a**2 + 8*a*b + 5*b**2)*atan(sinh(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_97():
    f = (a + b*tanh(c + d*x)**2)**3*sech(c + d*x)**2
    F = a**3*tanh(c + d*x)/d + a**2*b*tanh(c + d*x)**3/d + 3*a*b**2*tanh(c + d*x)**5/(5*d) + b**3*tanh(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_98():
    f = (a + b*tanh(c + d*x)**2)**3*sech(c + d*x)**3
    F = -b*(a + (a + b)*sinh(c + d*x)**2)**2*tanh(c + d*x)*sech(c + d*x)**7/(8*d) - b*(a + (a + b)*sinh(c + d*x)**2)*(12*a + 5*b)*tanh(c + d*x)*sech(c + d*x)**5/(48*d) - b*(72*a**2 + 52*a*b + 15*b**2)*tanh(c + d*x)*sech(c + d*x)**3/(192*d) + (64*a**3 + 48*a**2*b + 24*a*b**2 + 5*b**3)*tanh(c + d*x)*sech(c + d*x)/(128*d) + (64*a**3 + 48*a**2*b + 24*a*b**2 + 5*b**3)*atan(sinh(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_99():
    f = (a + b*tanh(c + d*x)**2)**3*sech(c + d*x)**4
    F = a**3*tanh(c + d*x)/d - a**2*(a - 3*b)*tanh(c + d*x)**3/(3*d) - 3*a*b*(a - b)*tanh(c + d*x)**5/(5*d) - b**3*tanh(c + d*x)**9/(9*d) - b**2*(3*a - b)*tanh(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_100():
    f = cosh(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = x*(3*a**2 + 10*a*b + 15*b**2)/(8*(a + b)**3) + sinh(c + d*x)*cosh(c + d*x)**3/(d*(4*a + 4*b)) + (3*a + 7*b)*sinh(c + d*x)*cosh(c + d*x)/(8*d*(a + b)**2) + b**(sympy.S(5)/2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_101():
    f = cosh(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = sinh(c + d*x)**3/(d*(3*a + 3*b)) + (a + 2*b)*sinh(c + d*x)/(d*(a + b)**2) + b**2*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_102():
    f = cosh(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = x*(a + 3*b)/(2*(a + b)**2) + sinh(c + d*x)*cosh(c + d*x)/(d*(2*a + 2*b)) + b**(sympy.S(3)/2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_103():
    f = cosh(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = sinh(c + d*x)/(d*(a + b)) + b*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_104():
    f = sech(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_105():
    f = sech(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_106():
    f = sech(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = -atan(sinh(c + d*x))/(b*d) + sqrt(a + b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_107():
    f = sech(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = -tanh(c + d*x)/(b*d) + (a + b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_108():
    f = sech(c + d*x)**5/(a + b*tanh(c + d*x)**2)
    F = -tanh(c + d*x)*sech(c + d*x)/(2*b*d) - (2*a + 3*b)*atan(sinh(c + d*x))/(2*b**2*d) + (a + b)**(sympy.S(3)/2)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(sqrt(a)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_109():
    f = sech(c + d*x)**6/(a + b*tanh(c + d*x)**2)
    F = tanh(c + d*x)**3/(3*b*d) - (a + 2*b)*tanh(c + d*x)/(b**2*d) + (a + b)**2*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_110():
    f = cosh(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = sinh(c + d*x)**3/(3*d*(a + b)**2) + (a + 3*b)*sinh(c + d*x)/(d*(a + b)**3) + b**3*sinh(c + d*x)/(2*a*d*(a + b)**3*(a + (a + b)*sinh(c + d*x)**2)) + b**2*(6*a + b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_111():
    f = cosh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = x*(a + 5*b)/(2*(a + b)**3) + sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)*(2*a + 2*b)) - b*(a - b)*tanh(c + d*x)/(2*a*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + b**(sympy.S(3)/2)*(5*a + b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_112():
    f = cosh(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = sinh(c + d*x)/(d*(a + b)**2) + b**2*sinh(c + d*x)/(2*a*d*(a + b)**2*(a + (a + b)*sinh(c + d*x)**2)) + b*(4*a + b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_113():
    f = sech(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = b*sinh(c + d*x)/(2*a*d*(a + b)*(a + (a + b)*sinh(c + d*x)**2)) + (2*a + b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_114():
    f = sech(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = tanh(c + d*x)/(2*a*d*(a + b*tanh(c + d*x)**2)) + atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_115():
    f = sech(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = sinh(c + d*x)/(2*a*d*(a + (a + b)*sinh(c + d*x)**2)) + atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_116():
    f = sech(c + d*x)**4/(a + b*tanh(c + d*x)**2)**2
    F = (a + b)*tanh(c + d*x)/(2*a*b*d*(a + b*tanh(c + d*x)**2)) - (a - b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_117():
    f = sech(c + d*x)**5/(a + b*tanh(c + d*x)**2)**2
    F = atan(sinh(c + d*x))/(b**2*d) + (a + b)*sinh(c + d*x)/(2*a*b*d*(a + (a + b)*sinh(c + d*x)**2)) - sqrt(a + b)*(2*a - b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_118():
    f = sech(c + d*x)**6/(a + b*tanh(c + d*x)**2)**2
    F = tanh(c + d*x)/(b**2*d) + (a + b)**2*tanh(c + d*x)/(2*a*b**2*d*(a + b*tanh(c + d*x)**2)) - (a + b)*(3*a - b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_119():
    f = sech(c + d*x)**7/(a + b*tanh(c + d*x)**2)**2
    F = -tanh(c + d*x)*sech(c + d*x)/(2*b*d*(a + (a + b)*sinh(c + d*x)**2)) + (4*a + 5*b)*atan(sinh(c + d*x))/(2*b**3*d) + (a + b)*(2*a + b)*sinh(c + d*x)/(2*a*b**2*d*(a + (a + b)*sinh(c + d*x)**2)) - (a + b)**(sympy.S(3)/2)*(4*a - b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_120():
    f = cosh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = x*(a + 7*b)/(2*(a + b)**4) + sinh(c + d*x)*cosh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)**2*(2*a + 2*b)) - b*(2*a - b)*tanh(c + d*x)/(4*a*d*(a + b)**2*(a + b*tanh(c + d*x)**2)**2) - b*(a - 3*b)*(4*a + b)*tanh(c + d*x)/(8*a**2*d*(a + b)**3*(a + b*tanh(c + d*x)**2)) + b**(sympy.S(3)/2)*(35*a**2 + 14*a*b + 3*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_121():
    f = cosh(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = sinh(c + d*x)/(d*(a + b)**3) + b**3*sinh(c + d*x)/(4*a*d*(a + b)**3*(a + (a + b)*sinh(c + d*x)**2)**2) + 3*b**2*(4*a + b)*sinh(c + d*x)/(8*a**2*d*(a + b)**3*(a + (a + b)*sinh(c + d*x)**2)) + 3*b*(8*a**2 + 4*a*b + b**2)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_122():
    f = sech(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = b*sinh(c + d*x)*cosh(c + d*x)**2/(4*a*d*(a + b)*(a + (a + b)*sinh(c + d*x)**2)**2) + 3*b*(2*a + b)*sinh(c + d*x)/(8*a**2*d*(a + b)**2*(a + (a + b)*sinh(c + d*x)**2)) + (8*a**2 + 8*a*b + 3*b**2)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_123():
    f = sech(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = tanh(c + d*x)/(4*a*d*(a + b*tanh(c + d*x)**2)**2) + 3*tanh(c + d*x)/(8*a**2*d*(a + b*tanh(c + d*x)**2)) + 3*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_124():
    f = sech(c + d*x)**3/(a + b*tanh(c + d*x)**2)**3
    F = b*sinh(c + d*x)/(4*a*d*(a + b)*(a + (a + b)*sinh(c + d*x)**2)**2) + (4*a + 3*b)*sinh(c + d*x)/(8*a**2*d*(a + b)*(a + (a + b)*sinh(c + d*x)**2)) + (4*a + 3*b)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_125():
    f = sech(c + d*x)**4/(a + b*tanh(c + d*x)**2)**3
    F = (a + b)*tanh(c + d*x)/(4*a*b*d*(a + b*tanh(c + d*x)**2)**2) - (a - 3*b)*tanh(c + d*x)/(8*a**2*b*d*(a + b*tanh(c + d*x)**2)) - (a - 3*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_126():
    f = sech(c + d*x)**5/(a + b*tanh(c + d*x)**2)**3
    F = sinh(c + d*x)/(4*a*d*(a + (a + b)*sinh(c + d*x)**2)**2) + 3*sinh(c + d*x)/(8*a**2*d*(a + (a + b)*sinh(c + d*x)**2)) + 3*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_127():
    f = sech(c + d*x)**6/(a + b*tanh(c + d*x)**2)**3
    F = (-3/b**2 + 3/a**2)*tanh(c + d*x)/(8*d*(a + b*tanh(c + d*x)**2)) + (a + b)*tanh(c + d*x)*sech(c + d*x)**2/(4*a*b*d*(a + b*tanh(c + d*x)**2)**2) + (3*a**2 - 2*a*b + 3*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_128():
    f = sech(c + d*x)**7/(a + b*tanh(c + d*x)**2)**3
    F = -atan(sinh(c + d*x))/(b**3*d) + (a + b)*sinh(c + d*x)/(4*a*b*d*(a + (a + b)*sinh(c + d*x)**2)**2) - (a + b)*(4*a - 3*b)*sinh(c + d*x)/(8*a**2*b**2*d*(a + (a + b)*sinh(c + d*x)**2)) + sqrt(a + b)*(8*a**2 - 4*a*b + 3*b**2)*atan(sqrt(a + b)*sinh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_129():
    f = (a + b*tanh(c + d*x)**2)*tanh(c + d*x)**4
    F = -b*tanh(c + d*x)**5/(5*d) + x*(a + b) - (a + b)*tanh(c + d*x)**3/(3*d) - (a + b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_130():
    f = (a + b*tanh(c + d*x)**2)*tanh(c + d*x)**3
    F = -b*tanh(c + d*x)**4/(4*d) + (a + b)*log(cosh(c + d*x))/d - (a + b)*tanh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_131():
    f = (a + b*tanh(c + d*x)**2)*tanh(c + d*x)**2
    F = -b*tanh(c + d*x)**3/(3*d) + x*(a + b) - (a + b)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_132():
    f = (a + b*tanh(c + d*x)**2)*tanh(c + d*x)
    F = -b*tanh(c + d*x)**2/(2*d) + (a + b)*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_133():
    f = a + b*tanh(c + d*x)**2
    F = a*x + b*x - b*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_134():
    f = (a + b*tanh(c + d*x)**2)*coth(c + d*x)
    F = a*log(sinh(c + d*x))/d + b*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_135():
    f = (a + b*tanh(c + d*x)**2)*coth(c + d*x)**2
    F = -a*coth(c + d*x)/d + x*(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_136():
    f = (a + b*tanh(c + d*x)**2)*coth(c + d*x)**3
    F = -a*coth(c + d*x)**2/(2*d) + (a + b)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_137():
    f = (a + b*tanh(c + d*x)**2)*coth(c + d*x)**4
    F = -a*coth(c + d*x)**3/(3*d) + x*(a + b) - (a + b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_138():
    f = (a + b*tanh(c + d*x)**2)*coth(c + d*x)**5
    F = -a*coth(c + d*x)**4/(4*d) + (a + b)*log(sinh(c + d*x))/d - (a + b)*coth(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_139():
    f = (a + b*tanh(c + d*x)**2)**2*tanh(c + d*x)**4
    F = -b**2*tanh(c + d*x)**7/(7*d) - b*(2*a + b)*tanh(c + d*x)**5/(5*d) + x*(a + b)**2 - (a + b)**2*tanh(c + d*x)**3/(3*d) - (a + b)**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_140():
    f = (a + b*tanh(c + d*x)**2)**2*tanh(c + d*x)**3
    F = -b**2*tanh(c + d*x)**6/(6*d) - b*(2*a + b)*tanh(c + d*x)**4/(4*d) + (a + b)**2*log(cosh(c + d*x))/d - (a + b)**2*tanh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_141():
    f = (a + b*tanh(c + d*x)**2)**2*tanh(c + d*x)**2
    F = -b**2*tanh(c + d*x)**5/(5*d) - b*(2*a + b)*tanh(c + d*x)**3/(3*d) + x*(a + b)**2 - (a + b)**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_142():
    f = (a + b*tanh(c + d*x)**2)**2*tanh(c + d*x)
    F = -b*(a + b)*tanh(c + d*x)**2/(2*d) + (a + b)**2*log(cosh(c + d*x))/d - (a + b*tanh(c + d*x)**2)**2/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_143():
    f = (a + b*tanh(c + d*x)**2)**2
    F = -b**2*tanh(c + d*x)**3/(3*d) - b*(2*a + b)*tanh(c + d*x)/d + x*(a + b)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_144():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)
    F = a**2*log(tanh(c + d*x))/d - b**2*tanh(c + d*x)**2/(2*d) + (a + b)**2*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_145():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**2
    F = -a**2*coth(c + d*x)/d - b**2*tanh(c + d*x)/d + x*(a + b)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_146():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**3
    F = -a**2*coth(c + d*x)**2/(2*d) + a*(a + 2*b)*log(tanh(c + d*x))/d + (a + b)**2*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_147():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**4
    F = -a**2*coth(c + d*x)**3/(3*d) - a*(a + 2*b)*coth(c + d*x)/d + x*(a + b)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_148():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**5
    F = -a**2*coth(c + d*x)**4/(4*d) - a*(a + 2*b)*coth(c + d*x)**2/(2*d) + (a + b)**2*log(cosh(c + d*x))/d + (a + b)**2*log(tanh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_149():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**6
    F = -a**2*coth(c + d*x)**5/(5*d) - a*(a + 2*b)*coth(c + d*x)**3/(3*d) + x*(a + b)**2 - (a + b)**2*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_150():
    f = (a + b*tanh(c + d*x)**2)**2*coth(c + d*x)**7
    F = -a**2*coth(c + d*x)**6/(6*d) - a*(a + 2*b)*coth(c + d*x)**4/(4*d) + (a + b)**2*log(cosh(c + d*x))/d + (a + b)**2*log(tanh(c + d*x))/d - (a + b)**2*coth(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_151():
    f = (a + b*tanh(c + d*x)**2)**3*tanh(c + d*x)**4
    F = -b**3*tanh(c + d*x)**9/(9*d) - b**2*(3*a + b)*tanh(c + d*x)**7/(7*d) - b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)**5/(5*d) + x*(a + b)**3 - (a + b)**3*tanh(c + d*x)**3/(3*d) - (a + b)**3*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_152():
    f = (a + b*tanh(c + d*x)**2)**3*tanh(c + d*x)**3
    F = -b**3*tanh(c + d*x)**8/(8*d) - b**2*(3*a + b)*tanh(c + d*x)**6/(6*d) - b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)**4/(4*d) + (a + b)**3*log(cosh(c + d*x))/d - (a + b)**3*tanh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_153():
    f = (a + b*tanh(c + d*x)**2)**3*tanh(c + d*x)**2
    F = -b**3*tanh(c + d*x)**7/(7*d) - b**2*(3*a + b)*tanh(c + d*x)**5/(5*d) - b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)**3/(3*d) + x*(a + b)**3 - (a + b)**3*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_154():
    f = (a + b*tanh(c + d*x)**2)**3*tanh(c + d*x)
    F = -b*(a + b)**2*tanh(c + d*x)**2/(2*d) + (a + b)**3*log(cosh(c + d*x))/d - (a + b)*(a + b*tanh(c + d*x)**2)**2/(4*d) - (a + b*tanh(c + d*x)**2)**3/(6*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_155():
    f = (a + b*tanh(c + d*x)**2)**3
    F = -b**3*tanh(c + d*x)**5/(5*d) - b**2*(3*a + b)*tanh(c + d*x)**3/(3*d) - b*(3*a**2 + 3*a*b + b**2)*tanh(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_156():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)
    F = a**3*log(tanh(c + d*x))/d - b**3*tanh(c + d*x)**4/(4*d) - b**2*(3*a + b)*tanh(c + d*x)**2/(2*d) + (a + b)**3*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_157():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**2
    F = -a**3*coth(c + d*x)/d - b**3*tanh(c + d*x)**3/(3*d) - b**2*(3*a + b)*tanh(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_158():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**3
    F = -a**3*coth(c + d*x)**2/(2*d) + a**2*(a + 3*b)*log(tanh(c + d*x))/d - b**3*tanh(c + d*x)**2/(2*d) + (a + b)**3*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_159():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**4
    F = -a**3*coth(c + d*x)**3/(3*d) - a**2*(a + 3*b)*coth(c + d*x)/d - b**3*tanh(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_160():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**5
    F = -a**3*coth(c + d*x)**4/(4*d) - a**2*(a + 3*b)*coth(c + d*x)**2/(2*d) + a*(a**2 + 3*a*b + 3*b**2)*log(tanh(c + d*x))/d + (a + b)**3*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_161():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**6
    F = -a**3*coth(c + d*x)**5/(5*d) - a**2*(a + 3*b)*coth(c + d*x)**3/(3*d) - a*(a**2 + 3*a*b + 3*b**2)*coth(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_162():
    f = (a + b*tanh(c + d*x)**2)**3*coth(c + d*x)**7
    F = -a**3*coth(c + d*x)**6/(6*d) - a**2*(a + 3*b)*coth(c + d*x)**4/(4*d) - a*(a**2 + 3*a*b + 3*b**2)*coth(c + d*x)**2/(2*d) + (a + b)**3*log(cosh(c + d*x))/d + (a + b)**3*log(tanh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_163():
    f = (a + b*tanh(c + d*x)**2)**4
    F = -b**4*tanh(c + d*x)**7/(7*d) - b**3*(4*a + b)*tanh(c + d*x)**5/(5*d) - b**2*(6*a**2 + 4*a*b + b**2)*tanh(c + d*x)**3/(3*d) - b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*tanh(c + d*x)/d + x*(a + b)**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_164():
    f = (a + b*tanh(c + d*x)**2)**5
    F = -b**5*tanh(c + d*x)**9/(9*d) - b**4*(5*a + b)*tanh(c + d*x)**7/(7*d) - b**3*(10*a**2 + 5*a*b + b**2)*tanh(c + d*x)**5/(5*d) - b**2*(10*a**3 + 10*a**2*b + 5*a*b**2 + b**3)*tanh(c + d*x)**3/(3*d) - b*(5*a**4 + 10*a**3*b + 10*a**2*b**2 + 5*a*b**3 + b**4)*tanh(c + d*x)/d + x*(a + b)**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_165():
    f = tanh(c + d*x)**5/(a + b*tanh(c + d*x)**2)
    F = a**2*log(a + b*tanh(c + d*x)**2)/(2*b**2*d*(a + b)) + log(cosh(c + d*x))/(d*(a + b)) - tanh(c + d*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_166():
    f = tanh(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = a**(sympy.S(3)/2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(b**(sympy.S(3)/2)*d*(a + b)) + x/(a + b) - tanh(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_167():
    f = tanh(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = -a*log(a + b*tanh(c + d*x)**2)/(2*b*d*(a + b)) + log(cosh(c + d*x))/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_168():
    f = tanh(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = -sqrt(a)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(b)*d*(a + b)) + x/(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_169():
    f = tanh(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = log(a + b*tanh(c + d*x)**2)/(d*(2*a + 2*b)) + log(cosh(c + d*x))/(d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_170():
    f = 1/(a + b*tanh(c + d*x)**2)
    F = x/(a + b) + sqrt(b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(sqrt(a)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_171():
    f = coth(c + d*x)/(a + b*tanh(c + d*x)**2)
    F = log(cosh(c + d*x))/(d*(a + b)) - b*log(a + b*tanh(c + d*x)**2)/(2*a*d*(a + b)) + log(tanh(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_172():
    f = coth(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = x/(a + b) - coth(c + d*x)/(a*d) - b**(sympy.S(3)/2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(3)/2)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_173():
    f = coth(c + d*x)**3/(a + b*tanh(c + d*x)**2)
    F = log(cosh(c + d*x))/(d*(a + b)) - coth(c + d*x)**2/(2*a*d) + b**2*log(a + b*tanh(c + d*x)**2)/(2*a**2*d*(a + b)) + (a - b)*log(tanh(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_174():
    f = coth(c + d*x)**4/(a + b*tanh(c + d*x)**2)
    F = x/(a + b) - coth(c + d*x)**3/(3*a*d) - (a - b)*coth(c + d*x)/(a**2*d) + b**(sympy.S(5)/2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(a**(sympy.S(5)/2)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_175():
    f = tanh(c + d*x)**5/(a + b*tanh(c + d*x)**2)**2
    F = -a**2/(2*b**2*d*(a + b)*(a + b*tanh(c + d*x)**2)) - a*(a + 2*b)*log(a + b*tanh(c + d*x)**2)/(2*b**2*d*(a + b)**2) + log(cosh(c + d*x))/(d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_176():
    f = tanh(c + d*x)**4/(a + b*tanh(c + d*x)**2)**2
    F = -sqrt(a)*(a + 3*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*b**(sympy.S(3)/2)*d*(a + b)**2) + a*tanh(c + d*x)/(2*b*d*(a + b)*(a + b*tanh(c + d*x)**2)) + x/(a + b)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_177():
    f = tanh(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = a/(2*b*d*(a + b)*(a + b*tanh(c + d*x)**2)) + log(a + b*tanh(c + d*x)**2)/(2*d*(a + b)**2) + log(cosh(c + d*x))/(d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_178():
    f = tanh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = x/(a + b)**2 - tanh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)*(2*a + 2*b)) - (a - b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*sqrt(a)*sqrt(b)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_179():
    f = tanh(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = -1/(d*(a + b*tanh(c + d*x)**2)*(2*a + 2*b)) + log(a + b*tanh(c + d*x)**2)/(2*d*(a + b)**2) + log(cosh(c + d*x))/(d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_180():
    f = (a + b*tanh(c + d*x)**2)**(-2)
    F = x/(a + b)**2 + b*tanh(c + d*x)/(2*a*d*(a + b)*(a + b*tanh(c + d*x)**2)) + sqrt(b)*(3*a + b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_181():
    f = coth(c + d*x)/(a + b*tanh(c + d*x)**2)**2
    F = log(cosh(c + d*x))/(d*(a + b)**2) + b/(2*a*d*(a + b)*(a + b*tanh(c + d*x)**2)) - b*(2*a + b)*log(a + b*tanh(c + d*x)**2)/(2*a**2*d*(a + b)**2) + log(tanh(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_182():
    f = coth(c + d*x)**2/(a + b*tanh(c + d*x)**2)**2
    F = x/(a + b)**2 + b*coth(c + d*x)/(2*a*d*(a + b)*(a + b*tanh(c + d*x)**2)) - (2*a + 3*b)*coth(c + d*x)/(2*a**2*d*(a + b)) - b**(sympy.S(3)/2)*(5*a + 3*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(5)/2)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_183():
    f = coth(c + d*x)**3/(a + b*tanh(c + d*x)**2)**2
    F = log(cosh(c + d*x))/(d*(a + b)**2) - b**2/(2*a**2*d*(a + b)*(a + b*tanh(c + d*x)**2)) - coth(c + d*x)**2/(2*a**2*d) + b**2*(3*a + 2*b)*log(a + b*tanh(c + d*x)**2)/(2*a**3*d*(a + b)**2) + (a - 2*b)*log(tanh(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_184():
    f = coth(c + d*x)**4/(a + b*tanh(c + d*x)**2)**2
    F = x/(a + b)**2 + b*coth(c + d*x)**3/(2*a*d*(a + b)*(a + b*tanh(c + d*x)**2)) - (2*a + 5*b)*coth(c + d*x)**3/(6*a**2*d*(a + b)) - (2*a**2 - 2*a*b - 5*b**2)*coth(c + d*x)/(2*a**3*d*(a + b)) + b**(sympy.S(5)/2)*(7*a + 5*b)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(2*a**(sympy.S(7)/2)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_185():
    f = tanh(c + d*x)**6/(a + b*tanh(c + d*x)**2)**3
    F = -sqrt(a)*(3*a**2 + 10*a*b + 15*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*b**(sympy.S(5)/2)*d*(a + b)**3) + a*tanh(c + d*x)**3/(4*b*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + a*(3*a + 7*b)*tanh(c + d*x)/(8*b**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + x/(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_186():
    f = tanh(c + d*x)**5/(a + b*tanh(c + d*x)**2)**3
    F = -a**2/(4*b**2*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + a*(a + 2*b)/(2*b**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + log(a + b*tanh(c + d*x)**2)/(2*d*(a + b)**3) + log(cosh(c + d*x))/(d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_187():
    f = tanh(c + d*x)**4/(a + b*tanh(c + d*x)**2)**3
    F = a*tanh(c + d*x)/(4*b*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + x/(a + b)**3 - (a + 5*b)*tanh(c + d*x)/(8*b*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - (a**2 + 6*a*b - 3*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*sqrt(a)*b**(sympy.S(3)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_188():
    f = tanh(c + d*x)**3/(a + b*tanh(c + d*x)**2)**3
    F = a/(4*b*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) - 1/(2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + log(a + b*tanh(c + d*x)**2)/(2*d*(a + b)**3) + log(cosh(c + d*x))/(d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_189():
    f = tanh(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = x/(a + b)**3 - tanh(c + d*x)/(d*(a + b*tanh(c + d*x)**2)**2*(4*a + 4*b)) - (3*a - b)*tanh(c + d*x)/(8*a*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - (3*a**2 - 6*a*b - b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(3)/2)*sqrt(b)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_190():
    f = tanh(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = -1/(d*(a + b*tanh(c + d*x)**2)**2*(4*a + 4*b)) - 1/(2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + log(a + b*tanh(c + d*x)**2)/(2*d*(a + b)**3) + log(cosh(c + d*x))/(d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_191():
    f = (a + b*tanh(c + d*x)**2)**(-3)
    F = x/(a + b)**3 + b*tanh(c + d*x)/(4*a*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + b*(7*a + 3*b)*tanh(c + d*x)/(8*a**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) + sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_192():
    f = coth(c + d*x)/(a + b*tanh(c + d*x)**2)**3
    F = log(cosh(c + d*x))/(d*(a + b)**3) + b/(4*a*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + b*(2*a + b)/(2*a**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - b*(3*a**2 + 3*a*b + b**2)*log(a + b*tanh(c + d*x)**2)/(2*a**3*d*(a + b)**3) + log(tanh(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_193():
    f = coth(c + d*x)**2/(a + b*tanh(c + d*x)**2)**3
    F = x/(a + b)**3 + b*coth(c + d*x)/(4*a*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + b*(9*a + 5*b)*coth(c + d*x)/(8*a**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - (8*a**2 + 27*a*b + 15*b**2)*coth(c + d*x)/(8*a**3*d*(a + b)**2) - b**(sympy.S(3)/2)*(35*a**2 + 42*a*b + 15*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(7)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_194():
    f = coth(c + d*x)**3/(a + b*tanh(c + d*x)**2)**3
    F = log(cosh(c + d*x))/(d*(a + b)**3) - b**2/(4*a**2*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) - b**2*(3*a + 2*b)/(2*a**3*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - coth(c + d*x)**2/(2*a**3*d) + b**2*(6*a**2 + 8*a*b + 3*b**2)*log(a + b*tanh(c + d*x)**2)/(2*a**4*d*(a + b)**3) + (a - 3*b)*log(tanh(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_195():
    f = coth(c + d*x)**4/(a + b*tanh(c + d*x)**2)**3
    F = x/(a + b)**3 + b*coth(c + d*x)**3/(4*a*d*(a + b)*(a + b*tanh(c + d*x)**2)**2) + b*(11*a + 7*b)*coth(c + d*x)**3/(8*a**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)) - (8*a**2 + 55*a*b + 35*b**2)*coth(c + d*x)**3/(24*a**3*d*(a + b)**2) - (8*a**3 - 8*a**2*b - 55*a*b**2 - 35*b**3)*coth(c + d*x)/(8*a**4*d*(a + b)**2) + b**(sympy.S(5)/2)*(63*a**2 + 90*a*b + 35*b**2)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(8*a**(sympy.S(9)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_196():
    f = (a + b*tanh(c + d*x)**2)**(-4)
    F = x/(a + b)**4 + b*tanh(c + d*x)/(6*a*d*(a + b)*(a + b*tanh(c + d*x)**2)**3) + b*(11*a + 5*b)*tanh(c + d*x)/(24*a**2*d*(a + b)**2*(a + b*tanh(c + d*x)**2)**2) + b*(19*a**2 + 16*a*b + 5*b**2)*tanh(c + d*x)/(16*a**3*d*(a + b)**3*(a + b*tanh(c + d*x)**2)) + sqrt(b)*(35*a**3 + 35*a**2*b + 21*a*b**2 + 5*b**3)*atan(sqrt(b)*tanh(c + d*x)/sqrt(a))/(16*a**(sympy.S(7)/2)*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_197():
    f = sqrt(1 - tanh(x)**2)
    F = asin(tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_198():
    f = sqrt(tanh(x)**2 - 1)
    F = -atanh(tanh(x)/sqrt(-sech(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_199():
    f = (1 - tanh(x)**2)**(sympy.S(3)/2)
    F = sqrt(sech(x)**2)*tanh(x)/2 + asin(tanh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_200():
    f = (tanh(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-sech(x)**2)*tanh(x)/2 + atanh(tanh(x)/sqrt(-sech(x)**2))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_201():
    f = 1/sqrt(1 - tanh(x)**2)
    F = tanh(x)/sqrt(sech(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_202():
    f = 1/sqrt(tanh(x)**2 - 1)
    F = tanh(x)/sqrt(-sech(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_203():
    f = sqrt(a + b*tanh(x)**2)*tanh(x)**5
    F = sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - sqrt(a + b*tanh(x)**2) + (a - b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)/(3*b**2) - (a + b*tanh(x)**2)**(sympy.S(5)/2)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_204():
    f = sqrt(a + b*tanh(x)**2)*tanh(x)**4
    F = sqrt(a + b)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2)) - sqrt(a + b*tanh(x)**2)*tanh(x)**3/4 - (a + 4*b)*sqrt(a + b*tanh(x)**2)*tanh(x)/(8*b) + (a**2 - 4*a*b - 8*b**2)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(8*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_205():
    f = sqrt(a + b*tanh(x)**2)*tanh(x)**3
    F = sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - sqrt(a + b*tanh(x)**2) - (a + b*tanh(x)**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_206():
    f = sqrt(a + b*tanh(x)**2)*tanh(x)**2
    F = sqrt(a + b)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2)) - sqrt(a + b*tanh(x)**2)*tanh(x)/2 - (a + 2*b)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_207():
    f = sqrt(a + b*tanh(x)**2)*tanh(x)
    F = sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - sqrt(a + b*tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_208():
    f = sqrt(a + b*tanh(x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2)) + sqrt(a + b)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_209():
    f = sqrt(a + b*tanh(x)**2)*coth(x)
    F = -sqrt(a)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a)) + sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_210():
    f = sqrt(a + b*tanh(x)**2)*coth(x)**2
    F = sqrt(a + b)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2)) - sqrt(a + b*tanh(x)**2)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_211():
    f = sqrt(a + b*tanh(x)**2)*coth(x)**3
    F = sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - sqrt(a + b*tanh(x)**2)*coth(x)**2/2 - (2*a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_212():
    f = sqrt(a + b*tanh(x)**2)*coth(x)**4
    F = sqrt(a + b)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2)) - sqrt(a + b*tanh(x)**2)*coth(x)**3/3 - sqrt(a + b*tanh(x)**2)*(3*a + b)*coth(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_213():
    f = sqrt(a + b*tanh(x)**2)*coth(x)**5
    F = sqrt(a + b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - sqrt(a + b*tanh(x)**2)*coth(x)**4/4 - sqrt(a + b*tanh(x)**2)*(4*a + b)*coth(x)**2/(8*a) - (8*a**2 + 4*a*b - b**2)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_214():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)*tanh(x)**3
    F = (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - (a + b)*sqrt(a + b*tanh(x)**2) - (a + b*tanh(x)**2)**(sympy.S(3)/2)/3 - (a + b*tanh(x)**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_215():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)*tanh(x)**2
    F = -b*sqrt(a + b*tanh(x)**2)*tanh(x)**3/4 - (5*a/8 + b/2)*sqrt(a + b*tanh(x)**2)*tanh(x) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2)) - (3*a**2 + 12*a*b + 8*b**2)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_216():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)*tanh(x)
    F = (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b)) - (a + b)*sqrt(a + b*tanh(x)**2) - (a + b*tanh(x)**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_217():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = -sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/2 - b*sqrt(a + b*tanh(x)**2)*tanh(x)/2 + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_218():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)*coth(x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a)) - b*sqrt(a + b*tanh(x)**2) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_219():
    f = (a + b*tanh(x)**2)**(sympy.S(3)/2)*coth(x)**2
    F = -a*sqrt(a + b*tanh(x)**2)*coth(x) - b**(sympy.S(3)/2)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2)) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_220():
    f = sqrt(tanh(x)**2 + 1)
    F = -asinh(tanh(x)) + sqrt(2)*atanh(sqrt(2)*tanh(x)/sqrt(tanh(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_221():
    f = sqrt(-tanh(x)**2 - 1)
    F = atan(tanh(x)/sqrt(-tanh(x)**2 - 1)) - sqrt(2)*atan(sqrt(2)*tanh(x)/sqrt(-tanh(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_222():
    f = (tanh(x)**2 + 1)**(sympy.S(3)/2)
    F = -sqrt(tanh(x)**2 + 1)*tanh(x)/2 - 5*asinh(tanh(x))/2 + 2*sqrt(2)*atanh(sqrt(2)*tanh(x)/sqrt(tanh(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_223():
    f = (-tanh(x)**2 - 1)**(sympy.S(3)/2)
    F = sqrt(-tanh(x)**2 - 1)*tanh(x)/2 - 5*atan(tanh(x)/sqrt(-tanh(x)**2 - 1))/2 + 2*sqrt(2)*atan(sqrt(2)*tanh(x)/sqrt(-tanh(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_224():
    f = tanh(x)**5/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/sqrt(a + b) + (a - b)*sqrt(a + b*tanh(x)**2)/b**2 - (a + b*tanh(x)**2)**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_225():
    f = tanh(x)**4/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/sqrt(a + b) - sqrt(a + b*tanh(x)**2)*tanh(x)/(2*b) + (a - 2*b)*atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_226():
    f = tanh(x)**3/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/sqrt(a + b) - sqrt(a + b*tanh(x)**2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_227():
    f = tanh(x)**2/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/sqrt(a + b) - atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_228():
    f = tanh(x)/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/sqrt(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_229():
    f = 1/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/sqrt(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_230():
    f = coth(x)/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/sqrt(a + b) - atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_231():
    f = coth(x)**2/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/sqrt(a + b) - sqrt(a + b*tanh(x)**2)*coth(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_232():
    f = coth(x)**3/sqrt(a + b*tanh(x)**2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/sqrt(a + b) - sqrt(a + b*tanh(x)**2)*coth(x)**2/(2*a) - (2*a - b)*atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_233():
    f = tanh(x)**5/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = -a**2/(b**2*(a + b)*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2) - sqrt(a + b*tanh(x)**2)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_234():
    f = tanh(x)**4/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = a*tanh(x)/(b*(a + b)*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(3)/2) - atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_235():
    f = tanh(x)**3/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = a/(b*(a + b)*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_236():
    f = tanh(x)**2/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = -tanh(x)/((a + b)*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_237():
    f = tanh(x)/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = -1/((a + b)*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_238():
    f = (a + b*tanh(x)**2)**(sympy.S(-3)/2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(3)/2) + b*tanh(x)/(a*(a + b)*sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_239():
    f = coth(x)/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2) + b/(a*(a + b)*sqrt(a + b*tanh(x)**2)) - atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_240():
    f = coth(x)**2/(a + b*tanh(x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(3)/2) + b*coth(x)/(a*(a + b)*sqrt(a + b*tanh(x)**2)) - (a + 2*b)*sqrt(a + b*tanh(x)**2)*coth(x)/(a**2*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_241():
    f = tanh(x)**6/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = a*tanh(x)**3/(3*b*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + a*(a + 2*b)*tanh(x)/(b**2*(a + b)**2*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(5)/2) - atanh(sqrt(b)*tanh(x)/sqrt(a + b*tanh(x)**2))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_242():
    f = tanh(x)**5/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = -a**2/(3*b**2*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + a*(a + 2*b)/(b**2*(a + b)**2*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_243():
    f = tanh(x)**4/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = a*tanh(x)/(3*b*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(5)/2) - (a + 4*b)*tanh(x)/(3*b*(a + b)**2*sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_244():
    f = tanh(x)**3/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = a/(3*b*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) - 1/((a + b)**2*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_245():
    f = tanh(x)**2/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = -tanh(x)/((a + b*tanh(x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) + atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(5)/2) - (2*a - b)*tanh(x)/(3*a*(a + b)**2*sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_246():
    f = tanh(x)/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = -1/((a + b*tanh(x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) - 1/((a + b)**2*sqrt(a + b*tanh(x)**2)) + atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_247():
    f = (a + b*tanh(x)**2)**(sympy.S(-5)/2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(5)/2) + b*tanh(x)/(3*a*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + b*(5*a + 2*b)*tanh(x)/(3*a**2*(a + b)**2*sqrt(a + b*tanh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_248():
    f = coth(x)/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b*tanh(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2) + b/(3*a*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + b*(2*a + b)/(a**2*(a + b)**2*sqrt(a + b*tanh(x)**2)) - atanh(sqrt(a + b*tanh(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_249():
    f = coth(x)**2/(a + b*tanh(x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b)*tanh(x)/sqrt(a + b*tanh(x)**2))/(a + b)**(sympy.S(5)/2) + b*coth(x)/(3*a*(a + b)*(a + b*tanh(x)**2)**(sympy.S(3)/2)) + b*(7*a + 4*b)*coth(x)/(3*a**2*(a + b)**2*sqrt(a + b*tanh(x)**2)) - (a + 4*b)*sqrt(a + b*tanh(x)**2)*(3*a + 2*b)*coth(x)/(3*a**3*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_250():
    f = 1/sqrt(tanh(x)**2 + 1)
    F = sqrt(2)*atanh(sqrt(2)*tanh(x)/sqrt(tanh(x)**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_251():
    f = 1/sqrt(-tanh(x)**2 - 1)
    F = sqrt(2)*atan(sqrt(2)*tanh(x)/sqrt(-tanh(x)**2 - 1))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_252():
    f = 1/(tanh(x)**3 + 1)
    F = x/2 - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(x))/3)/9 - 1/(6*tanh(x) + 6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_253():
    f = (a + b*tanh(x)**4)**(sympy.S(3)/2)*tanh(x)
    F = -sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*tanh(x)**2/sqrt(a + b*tanh(x)**4))/4 + (a + b)**(sympy.S(3)/2)*atanh((a + b*tanh(x)**2)/(sqrt(a + b)*sqrt(a + b*tanh(x)**4)))/2 - (a + b*tanh(x)**4)**(sympy.S(3)/2)/6 - sqrt(a + b*tanh(x)**4)*(a/2 + b*tanh(x)**2/4 + b/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_254():
    f = sqrt(a + b*tanh(x)**4)*tanh(x)
    F = -sqrt(b)*atanh(sqrt(b)*tanh(x)**2/sqrt(a + b*tanh(x)**4))/2 + sqrt(a + b)*atanh((a + b*tanh(x)**2)/(sqrt(a + b)*sqrt(a + b*tanh(x)**4)))/2 - sqrt(a + b*tanh(x)**4)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_255():
    f = tanh(x)/sqrt(a + b*tanh(x)**4)
    F = atanh((a + b*tanh(x)**2)/(sqrt(a + b)*sqrt(a + b*tanh(x)**4)))/(2*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_256():
    f = tanh(x)/(a + b*tanh(x)**4)**(sympy.S(3)/2)
    F = atanh((a + b*tanh(x)**2)/(sqrt(a + b)*sqrt(a + b*tanh(x)**4)))/(2*(a + b)**(sympy.S(3)/2)) - (a - b*tanh(x)**2)/(2*a*(a + b)*sqrt(a + b*tanh(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_7_d_hyper_pow_m_a_plus_b_c_tanh_pow_n_pow_p_257():
    f = tanh(x)/(a + b*tanh(x)**4)**(sympy.S(5)/2)
    F = atanh((a + b*tanh(x)**2)/(sqrt(a + b)*sqrt(a + b*tanh(x)**4)))/(2*(a + b)**(sympy.S(5)/2)) - (a - b*tanh(x)**2)/(6*a*(a + b)*(a + b*tanh(x)**4)**(sympy.S(3)/2)) - (3*a**2 - b*(5*a + 2*b)*tanh(x)**2)/(6*a**2*(a + b)**2*sqrt(a + b*tanh(x)**4))
    assert integrate(f, x) == F

