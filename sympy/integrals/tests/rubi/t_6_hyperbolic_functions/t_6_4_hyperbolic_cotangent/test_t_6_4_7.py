"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.4 Hyperbolic cotangent/6.4.7 (d hyper)^m (a+b (c coth)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d = symbols('a b c d')

def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_1():
    f = (a + b*coth(c + d*x)**2)**5
    F = -b**5*coth(c + d*x)**9/(9*d) - b**4*(5*a + b)*coth(c + d*x)**7/(7*d) - b**3*(10*a**2 + 5*a*b + b**2)*coth(c + d*x)**5/(5*d) - b**2*(10*a**3 + 10*a**2*b + 5*a*b**2 + b**3)*coth(c + d*x)**3/(3*d) - b*(5*a**4 + 10*a**3*b + 10*a**2*b**2 + 5*a*b**3 + b**4)*coth(c + d*x)/d + x*(a + b)**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_2():
    f = (a + b*coth(c + d*x)**2)**4
    F = -b**4*coth(c + d*x)**7/(7*d) - b**3*(4*a + b)*coth(c + d*x)**5/(5*d) - b**2*(6*a**2 + 4*a*b + b**2)*coth(c + d*x)**3/(3*d) - b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*coth(c + d*x)/d + x*(a + b)**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_3():
    f = (a + b*coth(c + d*x)**2)**3
    F = -b**3*coth(c + d*x)**5/(5*d) - b**2*(3*a + b)*coth(c + d*x)**3/(3*d) - b*(3*a**2 + 3*a*b + b**2)*coth(c + d*x)/d + x*(a + b)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_4():
    f = (a + b*coth(c + d*x)**2)**2
    F = -b**2*coth(c + d*x)**3/(3*d) - b*(2*a + b)*coth(c + d*x)/d + x*(a + b)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_5():
    f = 1/(a + b*coth(c + d*x)**2)
    F = x/(a + b) - sqrt(b)*atan(sqrt(a)*tanh(c + d*x)/sqrt(b))/(sqrt(a)*d*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_6():
    f = (a + b*coth(c + d*x)**2)**(-2)
    F = x/(a + b)**2 + b*coth(c + d*x)/(2*a*d*(a + b)*(a + b*coth(c + d*x)**2)) - sqrt(b)*(3*a + b)*atan(sqrt(a)*tanh(c + d*x)/sqrt(b))/(2*a**(sympy.S(3)/2)*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_7():
    f = (a + b*coth(c + d*x)**2)**(-3)
    F = x/(a + b)**3 + b*coth(c + d*x)/(4*a*d*(a + b)*(a + b*coth(c + d*x)**2)**2) + b*(7*a + 3*b)*coth(c + d*x)/(8*a**2*d*(a + b)**2*(a + b*coth(c + d*x)**2)) - sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atan(sqrt(a)*tanh(c + d*x)/sqrt(b))/(8*a**(sympy.S(5)/2)*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_8():
    f = (a + b*coth(c + d*x)**2)**(-4)
    F = x/(a + b)**4 + b*coth(c + d*x)/(6*a*d*(a + b)*(a + b*coth(c + d*x)**2)**3) + b*(11*a + 5*b)*coth(c + d*x)/(24*a**2*d*(a + b)**2*(a + b*coth(c + d*x)**2)**2) + b*(19*a**2 + 16*a*b + 5*b**2)*coth(c + d*x)/(16*a**3*d*(a + b)**3*(a + b*coth(c + d*x)**2)) - sqrt(b)*(35*a**3 + 35*a**2*b + 21*a*b**2 + 5*b**3)*atan(sqrt(a)*tanh(c + d*x)/sqrt(b))/(16*a**(sympy.S(7)/2)*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_9():
    f = 1/(1 - 2*coth(x)**2)
    F = -x + sqrt(2)*atanh(sqrt(2)*tanh(x)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_10():
    f = sqrt(1 - coth(x)**2)
    F = asin(coth(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_11():
    f = sqrt(coth(x)**2 - 1)
    F = -atanh(coth(x)/sqrt(csch(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_12():
    f = (1 - coth(x)**2)**(sympy.S(3)/2)
    F = sqrt(-csch(x)**2)*coth(x)/2 + asin(coth(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_13():
    f = (coth(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(csch(x)**2)*coth(x)/2 + atanh(coth(x)/sqrt(csch(x)**2))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_14():
    f = 1/sqrt(1 - coth(x)**2)
    F = coth(x)/sqrt(-csch(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_15():
    f = 1/sqrt(coth(x)**2 - 1)
    F = coth(x)/sqrt(csch(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_16():
    f = sqrt(a + b*coth(x)**2)*coth(x)**3
    F = sqrt(a + b)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b)) - sqrt(a + b*coth(x)**2) - (a + b*coth(x)**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_17():
    f = sqrt(a + b*coth(x)**2)*coth(x)**2
    F = sqrt(a + b)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2)) - sqrt(a + b*coth(x)**2)*coth(x)/2 - (a + 2*b)*atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_18():
    f = sqrt(a + b*coth(x)**2)*coth(x)
    F = sqrt(a + b)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b)) - sqrt(a + b*coth(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_19():
    f = sqrt(a + b*coth(x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2)) + sqrt(a + b)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_20():
    f = sqrt(a + b*coth(x)**2)*tanh(x)
    F = -sqrt(a)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a)) + sqrt(a + b)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_21():
    f = sqrt(a + b*coth(x)**2)*tanh(x)**2
    F = sqrt(a + b)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2)) - sqrt(a + b*coth(x)**2)*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_22():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)*coth(x)**3
    F = (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b)) - (a + b)*sqrt(a + b*coth(x)**2) - (a + b*coth(x)**2)**(sympy.S(3)/2)/3 - (a + b*coth(x)**2)**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_23():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)*coth(x)**2
    F = -b*sqrt(a + b*coth(x)**2)*coth(x)**3/4 - (5*a/8 + b/2)*sqrt(a + b*coth(x)**2)*coth(x) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2)) - (3*a**2 + 12*a*b + 8*b**2)*atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_24():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)*coth(x)
    F = (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b)) - (a + b)*sqrt(a + b*coth(x)**2) - (a + b*coth(x)**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_25():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)
    F = -sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2))/2 - b*sqrt(a + b*coth(x)**2)*coth(x)/2 + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_26():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)*tanh(x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a)) - b*sqrt(a + b*coth(x)**2) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_27():
    f = (a + b*coth(x)**2)**(sympy.S(3)/2)*tanh(x)**2
    F = -a*sqrt(a + b*coth(x)**2)*tanh(x) - b**(sympy.S(3)/2)*atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2)) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_28():
    f = sqrt(coth(x)**2 + 1)
    F = -asinh(coth(x)) + sqrt(2)*atanh(sqrt(2)*coth(x)/sqrt(coth(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_29():
    f = sqrt(-coth(x)**2 - 1)
    F = atan(coth(x)/sqrt(-coth(x)**2 - 1)) - sqrt(2)*atan(sqrt(2)*coth(x)/sqrt(-coth(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_30():
    f = (coth(x)**2 + 1)**(sympy.S(3)/2)
    F = -sqrt(coth(x)**2 + 1)*coth(x)/2 - 5*asinh(coth(x))/2 + 2*sqrt(2)*atanh(sqrt(2)*coth(x)/sqrt(coth(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_31():
    f = (-coth(x)**2 - 1)**(sympy.S(3)/2)
    F = sqrt(-coth(x)**2 - 1)*coth(x)/2 - 5*atan(coth(x)/sqrt(-coth(x)**2 - 1))/2 + 2*sqrt(2)*atan(sqrt(2)*coth(x)/sqrt(-coth(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_32():
    f = coth(x)**3/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/sqrt(a + b) - sqrt(a + b*coth(x)**2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_33():
    f = coth(x)**2/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/sqrt(a + b) - atanh(sqrt(b)*coth(x)/sqrt(a + b*coth(x)**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_34():
    f = coth(x)/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/sqrt(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_35():
    f = 1/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/sqrt(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_36():
    f = tanh(x)/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/sqrt(a + b) - atanh(sqrt(a + b*coth(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_37():
    f = tanh(x)**2/sqrt(a + b*coth(x)**2)
    F = atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/sqrt(a + b) - sqrt(a + b*coth(x)**2)*tanh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_38():
    f = coth(x)**3/(a + b*coth(x)**2)**(sympy.S(3)/2)
    F = a/(b*(a + b)*sqrt(a + b*coth(x)**2)) + atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_39():
    f = coth(x)**2/(a + b*coth(x)**2)**(sympy.S(3)/2)
    F = -coth(x)/((a + b)*sqrt(a + b*coth(x)**2)) + atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_40():
    f = coth(x)/(a + b*coth(x)**2)**(sympy.S(3)/2)
    F = -1/((a + b)*sqrt(a + b*coth(x)**2)) + atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_41():
    f = tanh(x)/(a + b*coth(x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(3)/2) + b/(a*(a + b)*sqrt(a + b*coth(x)**2)) - atanh(sqrt(a + b*coth(x)**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_42():
    f = tanh(x)**2/(a + b*coth(x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/(a + b)**(sympy.S(3)/2) + b*tanh(x)/(a*(a + b)*sqrt(a + b*coth(x)**2)) - (a + 2*b)*sqrt(a + b*coth(x)**2)*tanh(x)/(a**2*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_43():
    f = coth(x)**3/(a + b*coth(x)**2)**(sympy.S(5)/2)
    F = a/(3*b*(a + b)*(a + b*coth(x)**2)**(sympy.S(3)/2)) - 1/((a + b)**2*sqrt(a + b*coth(x)**2)) + atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_44():
    f = coth(x)**2/(a + b*coth(x)**2)**(sympy.S(5)/2)
    F = -coth(x)/((a + b*coth(x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) + atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/(a + b)**(sympy.S(5)/2) - (2*a - b)*coth(x)/(3*a*(a + b)**2*sqrt(a + b*coth(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_45():
    f = coth(x)/(a + b*coth(x)**2)**(sympy.S(5)/2)
    F = -1/((a + b*coth(x)**2)**(sympy.S(3)/2)*(3*a + 3*b)) - 1/((a + b)**2*sqrt(a + b*coth(x)**2)) + atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_46():
    f = tanh(x)/(a + b*coth(x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b*coth(x)**2)/sqrt(a + b))/(a + b)**(sympy.S(5)/2) + b/(3*a*(a + b)*(a + b*coth(x)**2)**(sympy.S(3)/2)) + b*(2*a + b)/(a**2*(a + b)**2*sqrt(a + b*coth(x)**2)) - atanh(sqrt(a + b*coth(x)**2)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_47():
    f = tanh(x)**2/(a + b*coth(x)**2)**(sympy.S(5)/2)
    F = atanh(sqrt(a + b)*coth(x)/sqrt(a + b*coth(x)**2))/(a + b)**(sympy.S(5)/2) + b*tanh(x)/(3*a*(a + b)*(a + b*coth(x)**2)**(sympy.S(3)/2)) + b*(7*a + 4*b)*tanh(x)/(3*a**2*(a + b)**2*sqrt(a + b*coth(x)**2)) - (a + 4*b)*sqrt(a + b*coth(x)**2)*(3*a + 2*b)*tanh(x)/(3*a**3*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_48():
    f = 1/sqrt(coth(x)**2 + 1)
    F = sqrt(2)*atanh(sqrt(2)*coth(x)/sqrt(coth(x)**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_49():
    f = 1/sqrt(-coth(x)**2 - 1)
    F = sqrt(2)*atan(sqrt(2)*coth(x)/sqrt(-coth(x)**2 - 1))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_50():
    f = 1/(coth(x)**3 + 1)
    F = x/2 - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*coth(x))/3)/9 - 1/(6*coth(x) + 6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_51():
    f = sqrt(a + b*coth(x)**4)*coth(x)
    F = -sqrt(b)*atanh(sqrt(b)*coth(x)**2/sqrt(a + b*coth(x)**4))/2 + sqrt(a + b)*atanh((a + b*coth(x)**2)/(sqrt(a + b)*sqrt(a + b*coth(x)**4)))/2 - sqrt(a + b*coth(x)**4)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_52():
    f = coth(x)/sqrt(a + b*coth(x)**4)
    F = atanh((a + b*coth(x)**2)/(sqrt(a + b)*sqrt(a + b*coth(x)**4)))/(2*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_7_d_hyper_pow_m_a_plus_b_c_coth_pow_n_pow_p_53():
    f = coth(x)/(a + b*coth(x)**4)**(sympy.S(3)/2)
    F = atanh((a + b*coth(x)**2)/(sqrt(a + b)*sqrt(a + b*coth(x)**4)))/(2*(a + b)**(sympy.S(3)/2)) - (a - b*coth(x)**2)/(2*a*(a + b)*sqrt(a + b*coth(x)**4))
    assert integrate(f, x) == F

