"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.6 Hyperbolic cosecant/6.6.7 (d hyper)^m (a+b (c csch)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_1():
    f = (a + b*csch(c + d*x)**2)**4
    F = a**4*x - b**4*coth(c + d*x)**7/(7*d) - b**3*(4*a - 3*b)*coth(c + d*x)**5/(5*d) - b**2*(6*a**2 - 8*a*b + 3*b**2)*coth(c + d*x)**3/(3*d) - b*(2*a - b)*(2*a**2 - 2*a*b + b**2)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_2():
    f = (a + b*csch(c + d*x)**2)**3
    F = a**3*x - b**3*coth(c + d*x)**5/(5*d) - b**2*(3*a - 2*b)*coth(c + d*x)**3/(3*d) - b*(3*a**2 - 3*a*b + b**2)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_3():
    f = (a + b*csch(c + d*x)**2)**2
    F = a**2*x - b**2*coth(c + d*x)**3/(3*d) - b*(2*a - b)*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_4():
    f = a + b*csch(c + d*x)**2
    F = a*x - b*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_5():
    f = 1/(a + b*csch(c + d*x)**2)
    F = -sqrt(b)*atan(sqrt(a - b)*tanh(c + d*x)/sqrt(b))/(a*d*sqrt(a - b)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_6():
    f = (a + b*csch(c + d*x)**2)**(-2)
    F = b*coth(c + d*x)/(2*a*d*(a - b)*(a + b*coth(c + d*x)**2 - b)) - sqrt(b)*(3*a - 2*b)*atan(sqrt(a - b)*tanh(c + d*x)/sqrt(b))/(2*a**2*d*(a - b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_7():
    f = (a + b*csch(c + d*x)**2)**(-3)
    F = b*coth(c + d*x)/(4*a*d*(a - b)*(a + b*coth(c + d*x)**2 - b)**2) + b*(7*a - 4*b)*coth(c + d*x)/(8*a**2*d*(a - b)**2*(a + b*coth(c + d*x)**2 - b)) - sqrt(b)*(15*a**2 - 20*a*b + 8*b**2)*atan(sqrt(a - b)*tanh(c + d*x)/sqrt(b))/(8*a**3*d*(a - b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_8():
    f = (a + b*csch(c + d*x)**2)**(-4)
    F = b*coth(c + d*x)/(6*a*d*(a - b)*(a + b*coth(c + d*x)**2 - b)**3) + b*(11*a - 6*b)*coth(c + d*x)/(24*a**2*d*(a - b)**2*(a + b*coth(c + d*x)**2 - b)**2) + b*(19*a**2 - 22*a*b + 8*b**2)*coth(c + d*x)/(16*a**3*d*(a - b)**3*(a + b*coth(c + d*x)**2 - b)) - sqrt(b)*(35*a**3 - 70*a**2*b + 56*a*b**2 - 16*b**3)*atan(sqrt(a - b)*tanh(c + d*x)/sqrt(b))/(16*a**4*d*(a - b)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_9():
    f = (a + b*csch(c + d*x)**2)**(sympy.S(5)/2)
    F = a**(sympy.S(5)/2)*atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/d - sqrt(b)*(15*a**2 - 10*a*b + 3*b**2)*atanh(sqrt(b)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/(8*d) - b*(7*a - 3*b)*sqrt(a + b*coth(c + d*x)**2 - b)*coth(c + d*x)/(8*d) - b*(a + b*coth(c + d*x)**2 - b)**(sympy.S(3)/2)*coth(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_10():
    f = (a + b*csch(c + d*x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/d - sqrt(b)*(3*a - b)*atanh(sqrt(b)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/(2*d) - b*sqrt(a + b*coth(c + d*x)**2 - b)*coth(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_11():
    f = sqrt(a + b*csch(c + d*x)**2)
    F = sqrt(a)*atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/d - sqrt(b)*atanh(sqrt(b)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_12():
    f = (a + b*csch(c + d*x)**2)**(sympy.S(-3)/2)
    F = b*coth(c + d*x)/(a*d*(a - b)*sqrt(a + b*coth(c + d*x)**2 - b)) + atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_13():
    f = (a + b*csch(c + d*x)**2)**(sympy.S(-5)/2)
    F = b*coth(c + d*x)/(3*a*d*(a - b)*(a + b*coth(c + d*x)**2 - b)**(sympy.S(3)/2)) + b*(5*a - 3*b)*coth(c + d*x)/(3*a**2*d*(a - b)**2*sqrt(a + b*coth(c + d*x)**2 - b)) + atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_14():
    f = (a + b*csch(c + d*x)**2)**(sympy.S(-7)/2)
    F = b*coth(c + d*x)/(5*a*d*(a - b)*(a + b*coth(c + d*x)**2 - b)**(sympy.S(5)/2)) + b*(9*a - 5*b)*coth(c + d*x)/(15*a**2*d*(a - b)**2*(a + b*coth(c + d*x)**2 - b)**(sympy.S(3)/2)) + b*(33*a**2 - 40*a*b + 15*b**2)*coth(c + d*x)/(15*a**3*d*(a - b)**3*sqrt(a + b*coth(c + d*x)**2 - b)) + atanh(sqrt(a)*coth(c + d*x)/sqrt(a + b*coth(c + d*x)**2 - b))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_15():
    f = (csch(x)**2 + 1)**(sympy.S(3)/2)
    F = -(coth(x)**2)**(sympy.S(3)/2)*tanh(x)/2 + sqrt(coth(x)**2)*log(sinh(x))*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_16():
    f = sqrt(csch(x)**2 + 1)
    F = sqrt(coth(x)**2)*log(sinh(x))*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_17():
    f = 1/sqrt(csch(x)**2 + 1)
    F = log(cosh(x))*coth(x)/sqrt(coth(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_18():
    f = (1 - csch(x)**2)**(sympy.S(3)/2)
    F = sqrt(2 - coth(x)**2)*coth(x)/2 + 2*asin(sqrt(2)*coth(x)/2) + atanh(coth(x)/sqrt(2 - coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_19():
    f = sqrt(1 - csch(x)**2)
    F = asin(sqrt(2)*coth(x)/2) + atanh(coth(x)/sqrt(2 - coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_20():
    f = 1/sqrt(1 - csch(x)**2)
    F = atanh(coth(x)/sqrt(2 - coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_21():
    f = (csch(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(coth(x)**2 - 2)*coth(x)/2 + atan(coth(x)/sqrt(coth(x)**2 - 2)) + 2*atanh(coth(x)/sqrt(coth(x)**2 - 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_22():
    f = sqrt(csch(x)**2 - 1)
    F = -atan(coth(x)/sqrt(coth(x)**2 - 2)) - atanh(coth(x)/sqrt(coth(x)**2 - 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_23():
    f = 1/sqrt(csch(x)**2 - 1)
    F = atan(coth(x)/sqrt(coth(x)**2 - 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_24():
    f = (-csch(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-coth(x)**2)*log(sinh(x))*tanh(x) + sqrt(-coth(x)**2)*coth(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_25():
    f = sqrt(-csch(x)**2 - 1)
    F = sqrt(-coth(x)**2)*log(sinh(x))*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_7_d_hyper_pow_m_a_plus_b_c_csch_pow_n_pow_p_26():
    f = 1/sqrt(-csch(x)**2 - 1)
    F = log(cosh(x))*coth(x)/sqrt(-coth(x)**2)
    assert integrate(f, x) == F

