"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.2 Hyperbolic cosine/6.2.7 hyper^m (a+b cosh^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, n = symbols('a b n')

def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_1():
    f = sinh(x)**4/(-a*cosh(x)**2 + a)
    F = x/(2*a) - sinh(x)*cosh(x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_2():
    f = sinh(x)**3/(-a*cosh(x)**2 + a)
    F = -cosh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_3():
    f = sinh(x)**2/(-a*cosh(x)**2 + a)
    F = -x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_4():
    f = csch(x)**2/(-a*cosh(x)**2 + a)
    F = coth(x)**3/(3*a) - coth(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_5():
    f = csch(x)**4/(-a*cosh(x)**2 + a)
    F = coth(x)**5/(5*a) - 2*coth(x)**3/(3*a) + coth(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_6():
    f = sinh(x)**7/(a + b*cosh(x)**2)
    F = cosh(x)**5/(5*b) - (a + 3*b)*cosh(x)**3/(3*b**2) + (a**2 + 3*a*b + 3*b**2)*cosh(x)/b**3 - (a + b)**3*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_7():
    f = sinh(x)**5/(a + b*cosh(x)**2)
    F = cosh(x)**3/(3*b) - (a + 2*b)*cosh(x)/b**2 + (a + b)**2*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_8():
    f = sinh(x)**3/(a + b*cosh(x)**2)
    F = cosh(x)/b - (a + b)*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_9():
    f = sinh(x)/(a + b*cosh(x)**2)
    F = atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_10():
    f = csch(x)/(a + b*cosh(x)**2)
    F = -atanh(cosh(x))/(a + b) - sqrt(b)*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_11():
    f = csch(x)**3/(a + b*cosh(x)**2)
    F = -coth(x)*csch(x)/(2*a + 2*b) + (a + 3*b)*atanh(cosh(x))/(2*(a + b)**2) + b**(sympy.S(3)/2)*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_12():
    f = csch(x)**5/(a + b*cosh(x)**2)
    F = -coth(x)*csch(x)**3/(4*a + 4*b) + (3*a + 7*b)*coth(x)*csch(x)/(8*(a + b)**2) - (3*a**2 + 10*a*b + 15*b**2)*atanh(cosh(x))/(8*(a + b)**3) - b**(sympy.S(5)/2)*atan(sqrt(b)*cosh(x)/sqrt(a))/(sqrt(a)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_13():
    f = sinh(x)**6/(a + b*cosh(x)**2)
    F = sinh(x)**3*cosh(x)/(4*b) - (4*a + 7*b)*sinh(x)*cosh(x)/(8*b**2) + x*(8*a**2 + 20*a*b + 15*b**2)/(8*b**3) - (a + b)**(sympy.S(5)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_14():
    f = sinh(x)**4/(a + b*cosh(x)**2)
    F = sinh(x)*cosh(x)/(2*b) - x*(2*a + 3*b)/(2*b**2) + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_15():
    f = sinh(x)**2/(a + b*cosh(x)**2)
    F = x/b - sqrt(a + b)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_16():
    f = 1/(a + b*cosh(x)**2)
    F = atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_17():
    f = csch(x)**4/(a + b*cosh(x)**2)
    F = -coth(x)**3/(3*a + 3*b) + (a + 2*b)*coth(x)/(a + b)**2 + b**2*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_18():
    f = csch(x)**6/(a + b*cosh(x)**2)
    F = -coth(x)**5/(5*a + 5*b) + (2*a + 3*b)*coth(x)**3/(3*(a + b)**2) - (a**2 + 3*a*b + 3*b**2)*coth(x)/(a + b)**3 - b**3*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_19():
    f = sinh(x)/(4 - 3*cosh(x)**3)
    F = -6**(sympy.S(2)/3)*log(-3**(sympy.S(1)/3)*cosh(x) + 2**(sympy.S(2)/3))/36 + 6**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*cosh(x)**2 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*cosh(x) + 2*2**(sympy.S(1)/3))/72 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/6)*atan(sqrt(3)*(6**(sympy.S(1)/3)*cosh(x) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_20():
    f = cosh(x)**7/(a + b*cosh(x)**2)
    F = -a**3*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(b**(sympy.S(7)/2)*sqrt(a + b)) + sinh(x)**5/(5*b) - (a - 2*b)*sinh(x)**3/(3*b**2) + (a**2 - a*b + b**2)*sinh(x)/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_21():
    f = cosh(x)**6/(a + b*cosh(x)**2)
    F = -a**(sympy.S(5)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(b**3*sqrt(a + b)) + sinh(x)*cosh(x)**3/(4*b) - (4*a - 3*b)*sinh(x)*cosh(x)/(8*b**2) + x*(8*a**2 - 4*a*b + 3*b**2)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_22():
    f = cosh(x)**5/(a + b*cosh(x)**2)
    F = a**2*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(b**(sympy.S(5)/2)*sqrt(a + b)) + sinh(x)**3/(3*b) - (a - b)*sinh(x)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_23():
    f = cosh(x)**4/(a + b*cosh(x)**2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(b**2*sqrt(a + b)) + sinh(x)*cosh(x)/(2*b) - x*(2*a - b)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_24():
    f = cosh(x)**3/(a + b*cosh(x)**2)
    F = -a*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(b**(sympy.S(3)/2)*sqrt(a + b)) + sinh(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_25():
    f = cosh(x)**2/(a + b*cosh(x)**2)
    F = -sqrt(a)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(b*sqrt(a + b)) + x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_26():
    f = cosh(x)/(a + b*cosh(x)**2)
    F = atan(sqrt(b)*sinh(x)/sqrt(a + b))/(sqrt(b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_27():
    f = 1/(a + b*cosh(x)**2)
    F = atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(sqrt(a)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_28():
    f = sech(x)/(a + b*cosh(x)**2)
    F = -sqrt(b)*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(a*sqrt(a + b)) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_29():
    f = sech(x)**2/(a + b*cosh(x)**2)
    F = tanh(x)/a - b*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(a**(sympy.S(3)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_30():
    f = sech(x)**3/(a + b*cosh(x)**2)
    F = tanh(x)*sech(x)/(2*a) + b**(sympy.S(3)/2)*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(a**2*sqrt(a + b)) + (a - 2*b)*atan(sinh(x))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_31():
    f = sech(x)**4/(a + b*cosh(x)**2)
    F = -tanh(x)**3/(3*a) + (a - b)*tanh(x)/a**2 + b**2*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(a**(sympy.S(5)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_32():
    f = sech(x)**5/(a + b*cosh(x)**2)
    F = tanh(x)*sech(x)**3/(4*a) + (3*a - 4*b)*tanh(x)*sech(x)/(8*a**2) - b**(sympy.S(5)/2)*atan(sqrt(b)*sinh(x)/sqrt(a + b))/(a**3*sqrt(a + b)) + (3*a**2 - 4*a*b + 8*b**2)*atan(sinh(x))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_33():
    f = (a + b*cosh(x)**2)**(-2)
    F = -b*sinh(x)*cosh(x)/(2*a*(a + b)*(a + b*cosh(x)**2)) + (2*a + b)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(2*a**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_34():
    f = (a + b*cosh(x)**2)**(-3)
    F = -b*sinh(x)*cosh(x)/(4*a*(a + b)*(a + b*cosh(x)**2)**2) - 3*b*(2*a + b)*sinh(x)*cosh(x)/(8*a**2*(a + b)**2*(a + b*cosh(x)**2)) + (8*a**2 + 8*a*b + 3*b**2)*atanh(sqrt(a)*tanh(x)/sqrt(a + b))/(8*a**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_35():
    f = 1/(cosh(x)**2 + 1)
    F = sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_36():
    f = (cosh(x)**2 + 1)**(-2)
    F = 3*sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/8 - sinh(x)*cosh(x)/(4*cosh(x)**2 + 4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_37():
    f = (cosh(x)**2 + 1)**(-3)
    F = 19*sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/64 - 9*sinh(x)*cosh(x)/(32*cosh(x)**2 + 32) - sinh(x)*cosh(x)/(8*(cosh(x)**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_38():
    f = 1/(1 - cosh(x)**2)
    F = coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_39():
    f = (1 - cosh(x)**2)**(-2)
    F = -coth(x)**3/3 + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_40():
    f = (1 - cosh(x)**2)**(-3)
    F = coth(x)**5/5 - 2*coth(x)**3/3 + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_41():
    f = sqrt(a + b*cosh(x)**2)
    F = -I*sqrt(a + b*cosh(x)**2)*elliptic_e(I*x + pi/2, -b/a)/sqrt(1 + b*cosh(x)**2/a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_42():
    f = sqrt(cosh(x)**2 + 1)
    F = -I*elliptic_e(I*x + pi/2, -1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_43():
    f = sqrt(1 - cosh(x)**2)
    F = sqrt(-sinh(x)**2)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_44():
    f = sqrt(cosh(x)**2 - 1)
    F = sqrt(sinh(x)**2)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_45():
    f = sqrt(-cosh(x)**2 - 1)
    F = -I*sqrt(-cosh(x)**2 - 1)*elliptic_e(I*x + pi/2, -1)/sqrt(cosh(x)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_46():
    f = (a + b*cosh(x)**2)**(sympy.S(3)/2)
    F = I*a*sqrt(1 + b*cosh(x)**2/a)*(a + b)*elliptic_f(I*x + pi/2, -b/a)/(3*sqrt(a + b*cosh(x)**2)) + b*sqrt(a + b*cosh(x)**2)*sinh(x)*cosh(x)/3 - 2*I*sqrt(a + b*cosh(x)**2)*(2*a + b)*elliptic_e(I*x + pi/2, -b/a)/(3*sqrt(1 + b*cosh(x)**2/a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_47():
    f = (cosh(x)**2 + 1)**(sympy.S(3)/2)
    F = sqrt(cosh(x)**2 + 1)*sinh(x)*cosh(x)/3 - 2*I*elliptic_e(I*x + pi/2, -1) + 2*I*elliptic_f(I*x + pi/2, -1)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_48():
    f = (1 - cosh(x)**2)**(sympy.S(3)/2)
    F = (-sinh(x)**2)**(sympy.S(3)/2)*coth(x)/3 + 2*sqrt(-sinh(x)**2)*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_49():
    f = (cosh(x)**2 - 1)**(sympy.S(3)/2)
    F = (sinh(x)**2)**(sympy.S(3)/2)*coth(x)/3 - 2*sqrt(sinh(x)**2)*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_50():
    f = (-cosh(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-cosh(x)**2 - 1)*sinh(x)*cosh(x)/3 + 2*I*sqrt(-cosh(x)**2 - 1)*elliptic_e(I*x + pi/2, -1)/sqrt(cosh(x)**2 + 1) + 2*I*sqrt(cosh(x)**2 + 1)*elliptic_f(I*x + pi/2, -1)/(3*sqrt(-cosh(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_51():
    f = 1/sqrt(a + b*cosh(x)**2)
    F = -I*sqrt(1 + b*cosh(x)**2/a)*elliptic_f(I*x + pi/2, -b/a)/sqrt(a + b*cosh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_52():
    f = 1/sqrt(cosh(x)**2 + 1)
    F = -I*elliptic_f(I*x + pi/2, -1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_53():
    f = 1/sqrt(1 - cosh(x)**2)
    F = -sinh(x)*atanh(cosh(x))/sqrt(-sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_54():
    f = 1/sqrt(cosh(x)**2 - 1)
    F = -sinh(x)*atanh(cosh(x))/sqrt(sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_55():
    f = 1/sqrt(-cosh(x)**2 - 1)
    F = -I*sqrt(cosh(x)**2 + 1)*elliptic_f(I*x + pi/2, -1)/sqrt(-cosh(x)**2 - 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_56():
    f = 1/(a + b*cosh(x)**3)
    F = 2*atanh(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + 2*atanh(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + 2*atanh(sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_57():
    f = 1/(a - b*cosh(x)**3)
    F = 2*atanh(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + 2*atanh(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + 2*atanh(sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3))*tanh(x/2)/sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_58():
    f = 1/(cosh(x)**3 + 1)
    F = -2*(-3)**(sympy.S(3)/4)*atan((-3)**(sympy.S(1)/4)*tanh(x/2))/(3*(3 + 3*(-1)**(sympy.S(2)/3))) - 2*(-3)**(sympy.S(3)/4)*atanh((-3)**(sympy.S(1)/4)*tanh(x/2))/(3*(3 - 3*(-1)**(sympy.S(1)/3))) + sinh(x)/(3*cosh(x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_59():
    f = 1/(1 - cosh(x)**3)
    F = -2*(-3)**(sympy.S(1)/4)*atan((-3)**(sympy.S(3)/4)*tanh(x/2)/3)/(3*(1 - (-1)**(sympy.S(2)/3))) - 2*(-3)**(sympy.S(1)/4)*atanh((-3)**(sympy.S(3)/4)*tanh(x/2)/3)/(3*(1 + (-1)**(sympy.S(1)/3))) - sinh(x)/(3 - 3*cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_60():
    f = 1/(a - b*cosh(x)**4)
    F = atanh(a**(sympy.S(1)/4)*tanh(x)/sqrt(sqrt(a) + sqrt(b)))/(2*a**(sympy.S(3)/4)*sqrt(sqrt(a) + sqrt(b))) + atanh(a**(sympy.S(1)/4)*tanh(x)/sqrt(sqrt(a) - sqrt(b)))/(2*a**(sympy.S(3)/4)*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_61():
    f = 1/(cosh(x)**4 + 1)
    F = sqrt(1 + sqrt(2))*log(sqrt(2)*coth(x)**2 + sqrt(2 + 2*sqrt(2))*coth(x) + 1)/8 - sqrt(1 + sqrt(2))*log(2*coth(x)**2 - 2*sqrt(1 + sqrt(2))*coth(x) + sqrt(2))/8 - atan((-2*coth(x) + sqrt(1 + sqrt(2)))/sqrt(-1 + sqrt(2)))/(4*sqrt(1 + sqrt(2))) + atan((2*coth(x) + sqrt(1 + sqrt(2)))/sqrt(-1 + sqrt(2)))/(4*sqrt(1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_62():
    f = 1/(1 - cosh(x)**4)
    F = coth(x)/2 + sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_63():
    f = 1/(a + b*cosh(x)**5)
    F = 2*atanh(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_64():
    f = 1/(a + b*cosh(x)**6)
    F = atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_65():
    f = 1/(a + b*cosh(x)**8)
    F = -atanh((-a)**(sympy.S(1)/8)*tanh(x)/sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*(-a)**(sympy.S(7)/8)*sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh((-a)**(sympy.S(1)/8)*tanh(x)/sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*(-a)**(sympy.S(7)/8)*sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh((-a)**(sympy.S(1)/8)*tanh(x)/sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*(-a)**(sympy.S(7)/8)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atanh((-a)**(sympy.S(1)/8)*tanh(x)/sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*(-a)**(sympy.S(7)/8)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_66():
    f = 1/(a - b*cosh(x)**5)
    F = 2*atanh(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))) + 2*atanh(sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5))*tanh(x/2)/sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_67():
    f = 1/(a - b*cosh(x)**6)
    F = atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) + atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) + atanh(a**(sympy.S(1)/6)*tanh(x)/sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_68():
    f = 1/(a - b*cosh(x)**8)
    F = atanh(a**(sympy.S(1)/8)*tanh(x)/sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4)))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4))) + atanh(a**(sympy.S(1)/8)*tanh(x)/sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4)))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4))) + atanh(a**(sympy.S(1)/8)*tanh(x)/sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4)))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4))) + atanh(a**(sympy.S(1)/8)*tanh(x)/sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4)))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_69():
    f = 1/(cosh(x)**5 + 1)
    F = -2*atan(tanh(x/2)/sqrt(-(1 - (-1)**(sympy.S(1)/5))/(1 + (-1)**(sympy.S(1)/5))))/(5*sqrt(-1 + (-1)**(sympy.S(2)/5))) - 2*sqrt(-(1 + (-1)**(sympy.S(3)/5))/(1 - (-1)**(sympy.S(3)/5)))*atan(sqrt(-(1 + (-1)**(sympy.S(3)/5))/(1 - (-1)**(sympy.S(3)/5)))*tanh(x/2))/(5 + 5*(-1)**(sympy.S(3)/5)) + 2*atanh(sqrt((1 - (-1)**(sympy.S(2)/5))/(1 + (-1)**(sympy.S(2)/5)))*tanh(x/2))/(5*sqrt(1 - (-1)**(sympy.S(4)/5))) + 2*atanh(sqrt((1 - (-1)**(sympy.S(4)/5))/(1 + (-1)**(sympy.S(4)/5)))*tanh(x/2))/(5*sqrt(1 + (-1)**(sympy.S(3)/5))) + sinh(x)/(5*cosh(x) + 5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_70():
    f = 1/(cosh(x)**6 + 1)
    F = sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/6 + atanh(tanh(x)/sqrt(1 - (-1)**(sympy.S(1)/3)))/(3*sqrt(1 - (-1)**(sympy.S(1)/3))) + atanh(tanh(x)/sqrt(1 + (-1)**(sympy.S(2)/3)))/(3*sqrt(1 + (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_71():
    f = 1/(cosh(x)**8 + 1)
    F = atanh(tanh(x)/sqrt(1 - (-1)**(sympy.S(1)/4)))/(4*sqrt(1 - (-1)**(sympy.S(1)/4))) + atanh(tanh(x)/sqrt(1 + (-1)**(sympy.S(1)/4)))/(4*sqrt(1 + (-1)**(sympy.S(1)/4))) + atanh(tanh(x)/sqrt(1 - (-1)**(sympy.S(3)/4)))/(4*sqrt(1 - (-1)**(sympy.S(3)/4))) + atanh(tanh(x)/sqrt(1 + (-1)**(sympy.S(3)/4)))/(4*sqrt(1 + (-1)**(sympy.S(3)/4)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_72():
    f = 1/(1 - cosh(x)**5)
    F = -2*atan(tanh(x/2)/sqrt(-(1 - (-1)**(sympy.S(2)/5))/(1 + (-1)**(sympy.S(2)/5))))/(5*sqrt(-1 + (-1)**(sympy.S(4)/5))) + 2*atan(sqrt(-(1 + (-1)**(sympy.S(4)/5))/(1 - (-1)**(sympy.S(4)/5)))*tanh(x/2))/(5*sqrt(-1 - (-1)**(sympy.S(3)/5))) + 2*atanh(sqrt((1 - (-1)**(sympy.S(1)/5))/(1 + (-1)**(sympy.S(1)/5)))*tanh(x/2))/(5*sqrt(1 - (-1)**(sympy.S(2)/5))) + 2*atanh(sqrt((1 - (-1)**(sympy.S(3)/5))/(1 + (-1)**(sympy.S(3)/5)))*tanh(x/2))/(5*sqrt(1 + (-1)**(sympy.S(1)/5))) - sinh(x)/(5 - 5*cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_73():
    f = 1/(1 - cosh(x)**6)
    F = coth(x)/3 + atanh(tanh(x)/sqrt(1 + (-1)**(sympy.S(1)/3)))/(3*sqrt(1 + (-1)**(sympy.S(1)/3))) + atanh(tanh(x)/sqrt(1 - (-1)**(sympy.S(2)/3)))/(3*sqrt(1 - (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_74():
    f = 1/(1 - cosh(x)**8)
    F = coth(x)/4 + sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/8 + atanh(tanh(x)/sqrt(1 - I))/(4*sqrt(1 - I)) + atanh(tanh(x)/sqrt(1 + I))/(4*sqrt(1 + I))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_75():
    f = tanh(x)/(cosh(x)**2 + 1)
    F = -log(cosh(x)**2 + 1)/2 + log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_76():
    f = sqrt(a + b*cosh(x)**2)*tanh(x)
    F = -sqrt(a)*atanh(sqrt(a + b*cosh(x)**2)/sqrt(a)) + sqrt(a + b*cosh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_77():
    f = tanh(x)/sqrt(a + b*cosh(x)**2)
    F = -atanh(sqrt(a + b*cosh(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_78():
    f = tanh(x)/sqrt(cosh(x)**2 + 1)
    F = -atanh(sqrt(cosh(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_79():
    f = tanh(x)/sqrt(1 - cosh(x)**2)
    F = -atanh(sqrt(-sinh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_80():
    f = tanh(x)**3/(a + b*cosh(x)**3)
    F = -log(a + b*cosh(x)**3)/(3*a) + log(cosh(x))/a + sech(x)**2/(2*a) + b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*cosh(x))/(3*a**(sympy.S(5)/3)) - b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cosh(x) + b**(sympy.S(2)/3)*cosh(x)**2)/(6*a**(sympy.S(5)/3)) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*cosh(x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_81():
    f = tanh(x)/sqrt(a + b*cosh(x)**3)
    F = -2*atanh(sqrt(a + b*cosh(x)**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_82():
    f = sqrt(a + b*cosh(x)**3)*tanh(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*cosh(x)**3)/sqrt(a))/3 + 2*sqrt(a + b*cosh(x)**3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_83():
    f = tanh(x)/sqrt(a + b*cosh(x)**n)
    F = -2*atanh(sqrt(a + b*cosh(x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_7_hyper_pow_m_a_plus_b_cosh_pow_n_pow_p_84():
    f = sqrt(a + b*cosh(x)**n)*tanh(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*cosh(x)**n)/sqrt(a))/n + 2*sqrt(a + b*cosh(x)**n)/n
    assert integrate(f, x) == F

