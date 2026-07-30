"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.1 Linear/1.1.1.6 P(x) (a+b x)^m (c+d x)^n (e+f x)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f = symbols('A B C a b c d e f')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_1():
    f = (e + f*x)**3*sqrt(-d*x + 1)*sqrt(d*x + 1)*(A + B*x + C*x**2)
    F = -C*(e + f*x)**4*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(7*d**2*f) + (e + f*x)**3*(-7*B*f + 3*C*e)*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(42*d**2*f) + x*sqrt(-d**2*x**2 + 1)*(8*A*d**4*e**3 + 6*A*d**2*e*f**2 + 6*B*d**2*e**2*f + B*f**3 + 2*C*d**2*e**3 + 3*C*e*f**2)/(16*d**4) - (e + f*x)**2*(-C*(3*d**2*e**2 - 8*f**2) + 7*d**2*f*(2*A*f + B*e))*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(70*d**4*f) + (8*A*d**4*e**3 + 6*A*d**2*e*f**2 + 6*B*d**2*e**2*f + B*f**3 + 2*C*d**2*e**3 + 3*C*e*f**2)*asin(d*x)/(16*d**5) + (-d**2*x**2 + 1)**(sympy.S(3)/2)*(8*C*(3*d**4*e**4 - 30*d**2*e**2*f**2 - 8*f**4) + 3*d**2*f*x*(-98*A*d**2*e*f**2 - 14*B*d**2*e**2*f - 35*B*f**3 + 6*C*d**2*e**3 - 41*C*e*f**2) - 56*d**2*f*(2*A*f*(6*d**2*e**2 + f**2) + B*(d**2*e**3 + 6*e*f**2)))/(840*d**6*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_2():
    f = (e + f*x)**2*sqrt(-d*x + 1)*sqrt(d*x + 1)*(A + B*x + C*x**2)
    F = -C*(e + f*x)**3*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(6*d**2*f) + (e + f*x)**2*(-2*B*f + C*e)*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(10*d**2*f) + x*(C*(2*d**2*e**2 + f**2) + 2*d**2*(A*(4*d**2*e**2 + f**2) + 2*B*e*f))*sqrt(-d**2*x**2 + 1)/(16*d**4) + (-d**2*x**2 + 1)**(sympy.S(3)/2)*(8*C*(d**2*e**3 - 4*e*f**2) - 3*f*x*(-2*d**2*e*(-2*B*f + C*e) + f**2*(10*A*d**2 + 5*C)) - 16*f*(5*A*d**2*e*f + B*(d**2*e**2 + f**2)))/(120*d**4*f) + (C*(2*d**2*e**2 + f**2) + 2*d**2*(A*(4*d**2*e**2 + f**2) + 2*B*e*f))*asin(d*x)/(16*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_3():
    f = (e + f*x)*sqrt(-d*x + 1)*sqrt(d*x + 1)*(A + B*x + C*x**2)
    F = -C*(e + f*x)**2*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(5*d**2*f) + x*sqrt(-d**2*x**2 + 1)*(4*A*d**2*e + B*f + C*e)/(8*d**2) + (4*A*d**2*e + B*f + C*e)*asin(d*x)/(8*d**3) - (-d**2*x**2 + 1)**(sympy.S(3)/2)*(-4*C*(3*d**2*e**2 - 2*f**2) - 3*d**2*f*x*(-5*B*f + 3*C*e) + 20*d**2*f*(A*f + B*e))/(60*d**4*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_4():
    f = sqrt(-d*x + 1)*sqrt(d*x + 1)*(A + B*x + C*x**2)
    F = -B*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(3*d**2) - C*x*(-d**2*x**2 + 1)**(sympy.S(3)/2)/(4*d**2) + x*(4*A*d**2 + C)*sqrt(-d**2*x**2 + 1)/(8*d**2) + (4*A*d**2 + C)*asin(d*x)/(8*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_5():
    f = (A + B*x + C*x**2)/((e + f*x)*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -C*sqrt(-d**2*x**2 + 1)/(d**2*f) + (A*f**2 - B*e*f + C*e**2)*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(f**2*sqrt(d**2*e**2 - f**2)) - (-B*f + C*e)*asin(d*x)/(d*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_6():
    f = (A + B*x + C*x**2)/((e + f*x)**2*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = C*asin(d*x)/(d*f**2) + sqrt(-d**2*x**2 + 1)*(A*f**2 - B*e*f + C*e**2)/(f*(e + f*x)*(d**2*e**2 - f**2)) - (-A*d**2*e*f**2 + B*f**3 + C*d**2*e**3 - 2*C*e*f**2)*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(f**2*(d**2*e**2 - f**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_7():
    f = (A + B*x + C*x**2)/((e + f*x)**3*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = (C*(d**2*e**2 + 2*f**2) - d**2*(-A*(2*d**2*e**2 + f**2) + 3*B*e*f))*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(2*(d**2*e**2 - f**2)**(sympy.S(5)/2)) - sqrt(-d**2*x**2 + 1)*(-3*A*d**2*e*f**2 + B*d**2*e**2*f + 2*B*f**3 + C*d**2*e**3 - 4*C*e*f**2)/(2*f*(e + f*x)*(d**2*e**2 - f**2)**2) + sqrt(-d**2*x**2 + 1)*(A*f**2 - B*e*f + C*e**2)/(2*f*(e + f*x)**2*(d**2*e**2 - f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_8():
    f = (e + f*x)**3*(A + B*x + C*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -C*(e + f*x)**4*sqrt(-d**2*x**2 + 1)/(5*d**2*f) + (e + f*x)**3*(-5*B*f + C*e)*sqrt(-d**2*x**2 + 1)/(20*d**2*f) - (e + f*x)**2*sqrt(-d**2*x**2 + 1)*(-3*d**2*e*(-5*B*f + C*e) + f**2*(20*A*d**2 + 16*C))/(60*d**4*f) + (8*A*d**4*e**3 + 12*A*d**2*e*f**2 + 12*B*d**2*e**2*f + 3*B*f**3 + 4*C*d**2*e**3 + 9*C*e*f**2)*asin(d*x)/(8*d**5) + sqrt(-d**2*x**2 + 1)*(4*C*(3*d**4*e**4 - 52*d**2*e**2*f**2 - 16*f**4) + d**2*f*x*(-100*A*d**2*e*f**2 - 30*B*d**2*e**2*f - 45*B*f**3 + 6*C*d**2*e**3 - 71*C*e*f**2) - 20*d**2*f*(4*A*f*(4*d**2*e**2 + f**2) + 3*B*(d**2*e**3 + 4*e*f**2)))/(120*d**6*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_9():
    f = (e + f*x)**2*(A + B*x + C*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -C*(e + f*x)**3*sqrt(-d**2*x**2 + 1)/(4*d**2*f) + (e + f*x)**2*(-4*B*f + C*e)*sqrt(-d**2*x**2 + 1)/(12*d**2*f) + sqrt(-d**2*x**2 + 1)*(4*C*(d**2*e**3 - 8*e*f**2) - f*x*(-2*d**2*e*(-4*B*f + C*e) + f**2*(12*A*d**2 + 9*C)) - 16*f*(3*A*d**2*e*f + B*(d**2*e**2 + f**2)))/(24*d**4*f) + (C*(4*d**2*e**2 + 3*f**2) + 4*d**2*(A*(2*d**2*e**2 + f**2) + 2*B*e*f))*asin(d*x)/(8*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_10():
    f = (e + f*x)*(A + B*x + C*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -C*(e + f*x)**2*sqrt(-d**2*x**2 + 1)/(3*d**2*f) + (2*A*d**2*e + B*f + C*e)*asin(d*x)/(2*d**3) - sqrt(-d**2*x**2 + 1)*(-2*C*(d**2*e**2 - 2*f**2) - d**2*f*x*(-3*B*f + C*e) + 6*d**2*f*(A*f + B*e))/(6*d**4*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_11():
    f = (A + B*x + C*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -B*sqrt(-d**2*x**2 + 1)/d**2 - C*x*sqrt(-d**2*x**2 + 1)/(2*d**2) + (2*A*d**2 + C)*asin(d*x)/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_12():
    f = (A + B*x + C*x**2)/((e + f*x)*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -C*sqrt(-d**2*x**2 + 1)/(d**2*f) + (A*f**2 - B*e*f + C*e**2)*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(f**2*sqrt(d**2*e**2 - f**2)) - (-B*f + C*e)*asin(d*x)/(d*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_13():
    f = (A + B*x + C*x**2)/((e + f*x)**2*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = C*asin(d*x)/(d*f**2) + sqrt(-d**2*x**2 + 1)*(A*f**2 - B*e*f + C*e**2)/(f*(e + f*x)*(d**2*e**2 - f**2)) - (-A*d**2*e*f**2 + B*f**3 + C*d**2*e**3 - 2*C*e*f**2)*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(f**2*(d**2*e**2 - f**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_14():
    f = (A + B*x + C*x**2)/((e + f*x)**3*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = (C*(d**2*e**2 + 2*f**2) - d**2*(-A*(2*d**2*e**2 + f**2) + 3*B*e*f))*atan((d**2*e*x + f)/(sqrt(d**2*e**2 - f**2)*sqrt(-d**2*x**2 + 1)))/(2*(d**2*e**2 - f**2)**(sympy.S(5)/2)) - sqrt(-d**2*x**2 + 1)*(-3*A*d**2*e*f**2 + B*d**2*e**2*f + 2*B*f**3 + C*d**2*e**3 - 4*C*e*f**2)/(2*f*(e + f*x)*(d**2*e**2 - f**2)**2) + sqrt(-d**2*x**2 + 1)*(A*f**2 - B*e*f + C*e**2)/(2*f*(e + f*x)**2*(d**2*e**2 - f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_15():
    f = x*(a + b*x + c*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = b*asin(d*x)/(2*d**3) - c*x**2*sqrt(-d**2*x**2 + 1)/(3*d**2) - sqrt(-d**2*x**2 + 1)*(6*a*d**2 + 3*b*d**2*x + 4*c)/(6*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_16():
    f = (a + b*x + c*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -b*sqrt(-d**2*x**2 + 1)/d**2 - c*x*sqrt(-d**2*x**2 + 1)/(2*d**2) + (2*a*d**2 + c)*asin(d*x)/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_17():
    f = (a + b*x + c*x**2)/(x*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -a*atanh(sqrt(-d**2*x**2 + 1)) + b*asin(d*x)/d - c*sqrt(-d**2*x**2 + 1)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_18():
    f = (a + b*x + c*x**2)/(x**2*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -a*sqrt(-d**2*x**2 + 1)/x - b*atanh(sqrt(-d**2*x**2 + 1)) + c*asin(d*x)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_19():
    f = (a + b*x + c*x**2)/(x**3*sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -a*sqrt(-d**2*x**2 + 1)/(2*x**2) - b*sqrt(-d**2*x**2 + 1)/x - (a*d**2/2 + c)*atanh(sqrt(-d**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_20():
    f = sqrt(a + b*x)*(e + f*x)**3*sqrt(a*c - b*c*x)*(A + B*x + C*x**2)
    F = -C*sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**4*sqrt(a*c - b*c*x)/(7*b**2*f) + a**2*sqrt(c)*sqrt(a + b*x)*(A*(6*a**2*b**2*e*f**2 + 8*b**4*e**3) + a**2*(a**2*f**2*(B*f + 3*C*e) + 2*b**2*e**2*(3*B*f + C*e)))*sqrt(a*c - b*c*x)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(16*b**5*sqrt(a**2*c - b**2*c*x**2)) + sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**3*(-7*B*f + 3*C*e)*sqrt(a*c - b*c*x)/(42*b**2*f) + x*sqrt(a + b*x)*(A*(6*a**2*b**2*e*f**2 + 8*b**4*e**3) + a**2*(a**2*f**2*(B*f + 3*C*e) + 2*b**2*e**2*(3*B*f + C*e)))*sqrt(a*c - b*c*x)/(16*b**4) - sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**2*sqrt(a*c - b*c*x)*(8*C*a**2*f**2 - b**2*(3*C*e**2 - 7*f*(2*A*f + B*e)))/(70*b**4*f) - sqrt(a + b*x)*(a**2 - b**2*x**2)*sqrt(a*c - b*c*x)*(64*C*a**4*f**4 + 16*a**2*b**2*f**2*(15*C*e**2 + 7*f*(A*f + 3*B*e)) - 8*b**4*e**2*(3*C*e**2 - 7*f*(12*A*f + B*e)) + 3*b**2*f*x*(a**2*f**2*(35*B*f + 41*C*e) - 2*b**2*e*(3*C*e**2 - 7*f*(7*A*f + B*e))))/(840*b**6*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_21():
    f = sqrt(a + b*x)*(e + f*x)**2*sqrt(a*c - b*c*x)*(A + B*x + C*x**2)
    F = -C*sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**3*sqrt(a*c - b*c*x)/(6*b**2*f) + a**2*sqrt(c)*sqrt(a + b*x)*(2*A*(a**2*b**2*f**2 + 4*b**4*e**2) + a**2*(C*a**2*f**2 + 2*b**2*e*(2*B*f + C*e)))*sqrt(a*c - b*c*x)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(16*b**5*sqrt(a**2*c - b**2*c*x**2)) + sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**2*(-2*B*f + C*e)*sqrt(a*c - b*c*x)/(10*b**2*f) + x*sqrt(a + b*x)*(2*A*(a**2*b**2*f**2 + 4*b**4*e**2) + a**2*(C*a**2*f**2 + 2*b**2*e*(2*B*f + C*e)))*sqrt(a*c - b*c*x)/(16*b**4) - sqrt(a + b*x)*(a**2 - b**2*x**2)*sqrt(a*c - b*c*x)*(16*a**2*f**2*(B*f + 2*C*e) - 8*b**2*e*(C*e**2 - 2*f*(5*A*f + B*e)) + 3*f*x*(5*C*a**2*f**2 - b**2*(2*C*e**2 - 2*f*(5*A*f + 2*B*e))))/(120*b**4*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_22():
    f = sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x)*(A + B*x + C*x**2)
    F = -C*sqrt(a + b*x)*(a**2 - b**2*x**2)*(e + f*x)**2*sqrt(a*c - b*c*x)/(5*b**2*f) + a**2*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(4*A*b**2*e + a**2*(B*f + C*e))*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(8*b**3*sqrt(a**2*c - b**2*c*x**2)) + x*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(4*A*b**2*e + a**2*(B*f + C*e))/(8*b**2) - sqrt(a + b*x)*(a**2 - b**2*x**2)*sqrt(a*c - b*c*x)*(8*C*a**2*f**2 - 3*b**2*f*x*(-5*B*f + 3*C*e) - 4*b**2*(3*C*e**2 - 5*f*(A*f + B*e)))/(60*b**4*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_23():
    f = sqrt(a + b*x)*sqrt(a*c - b*c*x)*(A + B*x + C*x**2)
    F = -B*sqrt(a + b*x)*(a**2 - b**2*x**2)*sqrt(a*c - b*c*x)/(3*b**2) - C*x*sqrt(a + b*x)*(a**2 - b**2*x**2)*sqrt(a*c - b*c*x)/(4*b**2) + a**2*sqrt(c)*sqrt(a + b*x)*(4*A*b**2 + C*a**2)*sqrt(a*c - b*c*x)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(8*b**3*sqrt(a**2*c - b**2*c*x**2)) + x*(A/2 + C*a**2/(8*b**2))*sqrt(a + b*x)*sqrt(a*c - b*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_24():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x))
    F = -C*(a**2 - b**2*x**2)/(b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + sqrt(a**2*c - b**2*c*x**2)*(A*f**2 - B*e*f + C*e**2)*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)*sqrt(-a**2*f**2 + b**2*e**2)) - (-B*f + C*e)*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(b*sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_25():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)**2*sqrt(a*c - b*c*x))
    F = C*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(b*sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + f*(A + e*(-B*f + C*e)/f**2)*(a**2 - b**2*x**2)/(sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)) + sqrt(a**2*c - b**2*c*x**2)*(a**2*f**2*(-B*f + 2*C*e) - b**2*(-A*e*f**2 + C*e**3))*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_26():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)**3*sqrt(a*c - b*c*x))
    F = f*(A + e*(-B*f + C*e)/f**2)*(a**2 - b**2*x**2)/(sqrt(a + b*x)*(e + f*x)**2*sqrt(a*c - b*c*x)*(-2*a**2*f**2 + 2*b**2*e**2)) + (a**2 - b**2*x**2)*(2*a**2*f**2*(-B*f + 2*C*e) - b**2*e*(C*e**2 + f*(-3*A*f + B*e)))/(2*f*sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**2) + (A*(a**2*b**2*f**2 + 2*b**4*e**2) + a**2*(2*C*a**2*f**2 + b**2*e*(-3*B*f + C*e)))*sqrt(a**2*c - b**2*c*x**2)*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(2*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_27():
    f = (e + f*x)**3*(A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(a*c - b*c*x))
    F = -C*(a**2 - b**2*x**2)*(e + f*x)**4/(5*b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + (a**2 - b**2*x**2)*(e + f*x)**3*(-5*B*f + C*e)/(20*b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) - (a**2 - b**2*x**2)*(e + f*x)**2*(16*C*a**2*f**2 - b**2*(3*C*e**2 - 5*f*(4*A*f + 3*B*e)))/(60*b**4*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + (4*A*(3*a**2*b**2*e*f**2 + 2*b**4*e**3) + a**2*(3*a**2*f**2*(B*f + 3*C*e) + 4*b**2*e**2*(3*B*f + C*e)))*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(8*b**5*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x)) - (a**2 - b**2*x**2)*(64*C*a**4*f**4 + 16*a**2*b**2*f**2*(13*C*e**2 + 5*f*(A*f + 3*B*e)) - 4*b**4*e**2*(3*C*e**2 - 5*f*(16*A*f + 3*B*e)) + b**2*f*x*(a**2*f**2*(45*B*f + 71*C*e) - 2*b**2*e*(3*C*e**2 - 5*f*(10*A*f + 3*B*e))))/(120*b**6*f*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_28():
    f = (e + f*x)**2*(A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(a*c - b*c*x))
    F = -C*(a**2 - b**2*x**2)*(e + f*x)**3/(4*b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + (a**2 - b**2*x**2)*(e + f*x)**2*(-4*B*f + C*e)/(12*b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) - (a**2 - b**2*x**2)*(16*a**2*f**2*(B*f + 2*C*e) - 4*b**2*e*(C*e**2 - 4*f*(3*A*f + B*e)) + f*x*(9*C*a**2*f**2 - b**2*(2*C*e**2 - 4*f*(3*A*f + 2*B*e))))/(24*b**4*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + (4*A*(a**2*b**2*f**2 + 2*b**4*e**2) + a**2*(3*C*a**2*f**2 + 4*b**2*e*(2*B*f + C*e)))*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(8*b**5*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_29():
    f = (e + f*x)*(A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(a*c - b*c*x))
    F = -C*(a**2 - b**2*x**2)*(e + f*x)**2/(3*b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + sqrt(a**2*c - b**2*c*x**2)*(2*A*b**2*e + a**2*(B*f + C*e))*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(2*b**3*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x)) - (a**2 - b**2*x**2)*(4*C*a**2*f**2 - b**2*f*x*(-3*B*f + C*e) - 2*b**2*(C*e**2 - 3*f*(A*f + B*e)))/(6*b**4*f*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_30():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(a*c - b*c*x))
    F = -B*(a**2 - b**2*x**2)/(b**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)) - C*x*(a**2 - b**2*x**2)/(2*b**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + (2*A*b**2 + C*a**2)*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(2*b**3*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_31():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x))
    F = -C*(a**2 - b**2*x**2)/(b**2*f*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + sqrt(a**2*c - b**2*c*x**2)*(A*f**2 - B*e*f + C*e**2)*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)*sqrt(-a**2*f**2 + b**2*e**2)) - (-B*f + C*e)*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(b*sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_32():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)**2*sqrt(a*c - b*c*x))
    F = C*sqrt(a**2*c - b**2*c*x**2)*atan(b*sqrt(c)*x/sqrt(a**2*c - b**2*c*x**2))/(b*sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)) + f*(A + e*(-B*f + C*e)/f**2)*(a**2 - b**2*x**2)/(sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)) + sqrt(a**2*c - b**2*c*x**2)*(a**2*f**2*(-B*f + 2*C*e) - b**2*(-A*e*f**2 + C*e**3))*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(sqrt(c)*f**2*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_33():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*(e + f*x)**3*sqrt(a*c - b*c*x))
    F = f*(A + e*(-B*f + C*e)/f**2)*(a**2 - b**2*x**2)/(sqrt(a + b*x)*(e + f*x)**2*sqrt(a*c - b*c*x)*(-2*a**2*f**2 + 2*b**2*e**2)) + (a**2 - b**2*x**2)*(2*a**2*f**2*(-B*f + 2*C*e) - b**2*e*(C*e**2 + f*(-3*A*f + B*e)))/(2*f*sqrt(a + b*x)*(e + f*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**2) + (A*(a**2*b**2*f**2 + 2*b**4*e**2) + a**2*(2*C*a**2*f**2 + b**2*e*(-3*B*f + C*e)))*sqrt(a**2*c - b**2*c*x**2)*atan(sqrt(c)*(a**2*f + b**2*e*x)/(sqrt(a**2*c - b**2*c*x**2)*sqrt(-a**2*f**2 + b**2*e**2)))/(2*sqrt(c)*sqrt(a + b*x)*sqrt(a*c - b*c*x)*(-a**2*f**2 + b**2*e**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_34():
    f = (a + b*x)**2*sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)
    F = C*(a + b*x)**3*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)/(6*b*d*f) + (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(8*a**2*d**2*f**2*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - 8*a*b*d*f*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))) + b**2*(C*(21*c**4*f**4 + 28*c**3*d*e*f**3 + 30*c**2*d**2*e**2*f**2 + 28*c*d**3*e**3*f + 21*d**4*e**4) + 4*d*f*(2*A*d*f*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) - B*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3))))/(256*d**5*f**4) + sqrt(c + d*x)*sqrt(e + f*x)*(-c*f + d*e)*(8*a**2*d**2*f**2*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - 8*a*b*d*f*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))) + b**2*(C*(21*c**4*f**4 + 28*c**3*d*e*f**3 + 30*c**2*d**2*e**2*f**2 + 28*c*d**3*e**3*f + 21*d**4*e**4) + 4*d*f*(2*A*d*f*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) - B*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3))))/(512*d**5*f**5) - (-c*f + d*e)**2*(8*a**2*d**2*f**2*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - 8*a*b*d*f*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))) + b**2*(C*(21*c**4*f**4 + 28*c**3*d*e*f**3 + 30*c**2*d**2*e**2*f**2 + 28*c*d**3*e**3*f + 21*d**4*e**4) + 4*d*f*(2*A*d*f*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) - B*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(512*d**(sympy.S(11)/2)*f**(sympy.S(11)/2)) - (a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(2*C*a*d*f - b*(4*B*d*f - 3*C*(c*f + d*e)))/(20*b*d**2*f**2) - (c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(64*C*a**3*d**3*f**3 - 8*a**2*b*d**2*f**2*(16*B*d*f - 7*C*(c*f + d*e)) - 8*a*b**2*d*f*(C*(35*c**2*f**2 + 38*c*d*e*f + 35*d**2*e**2) + 10*d*f*(8*A*d*f - 5*B*(c*f + d*e))) + b**3*(7*C*(15*c**3*f**3 + 17*c**2*d*e*f**2 + 17*c*d**2*e**2*f + 15*d**3*e**3) + 4*d*f*(50*A*d*f*(c*f + d*e) - B*(35*c**2*f**2 + 38*c*d*e*f + 35*d**2*e**2))) + 6*b*d*f*x*(10*b*d*f*(-4*A*b*d*f + C*a*c*f + C*a*d*e + 2*C*b*c*e) + (4*a*d*f - 7*b*(c*f + d*e))*(2*C*a*d*f - b*(4*B*d*f - 3*C*(c*f + d*e)))))/(960*b*d**4*f**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_35():
    f = (a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)
    F = C*(a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)/(5*b*d*f) + (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(2*a*d*f*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - b*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))))/(64*d**4*f**3) + sqrt(c + d*x)*sqrt(e + f*x)*(-c*f + d*e)*(2*a*d*f*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - b*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))))/(128*d**4*f**4) - (-c*f + d*e)**2*(2*a*d*f*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e))) - b*(C*(7*c**3*f**3 + 9*c**2*d*e*f**2 + 9*c*d**2*e**2*f + 7*d**3*e**3) + 2*d*f*(8*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(128*d**(sympy.S(9)/2)*f**(sympy.S(9)/2)) - (c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(48*C*a**2*d**2*f**2 - 10*a*b*d*f*(8*B*d*f - 5*C*(c*f + d*e)) - b**2*(C*(35*c**2*f**2 + 38*c*d*e*f + 35*d**2*e**2) + 10*d*f*(8*A*d*f - 5*B*(c*f + d*e))) + 6*b*d*f*x*(6*C*a*d*f - b*(10*B*d*f - 7*C*(c*f + d*e))))/(240*b*d**3*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_36():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)
    F = C*(c + d*x)**(sympy.S(5)/2)*(e + f*x)**(sympy.S(3)/2)/(4*d**2*f) - (c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(-8*B*d*f + 11*C*c*f + 5*C*d*e)/(24*d**2*f**2) + (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e)))/(32*d**3*f**2) + sqrt(c + d*x)*sqrt(e + f*x)*(C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e)))*(-c*f + d*e)/(64*d**3*f**3) - (C*(5*c**2*f**2 + 6*c*d*e*f + 5*d**2*e**2) + 8*d*f*(2*A*d*f - B*(c*f + d*e)))*(-c*f + d*e)**2*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(64*d**(sympy.S(7)/2)*f**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_37():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)
    F = C*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)/(3*b*d*f) - sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(2*C*a*d*f + b*(-2*B*d*f + C*c*f + C*d*e))/(4*b**2*d*f**2) + sqrt(c + d*x)*sqrt(e + f*x)*(4*b*d*f*(2*A*b*d*f - C*a*(c*f + d*e)) + (2*C*a*d*f + b*(-2*B*d*f + C*c*f + C*d*e))*(4*a*d*f - b*c*f + b*d*e))/(8*b**3*d**2*f**2) - (2*A*b**2 - 2*a*(B*b - C*a))*sqrt(-a*d + b*c)*sqrt(-a*f + b*e)*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/b**4 - (16*C*a**3*d**3*f**3 - 8*a**2*b*d**2*f**2*(2*B*d*f + C*c*f + C*d*e) - 2*a*b**2*d*f*(C*(-c*f + d*e)**2 - 4*d*f*(2*A*d*f + B*c*f + B*d*e)) - b**3*(C*(-c*f + d*e)**2*(c*f + d*e) - 2*d*f*(-4*A*d*f*(c*f + d*e) + B*(-c*f + d*e)**2)))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(8*b**4*d**(sympy.S(5)/2)*f**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_38():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**2
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(A*b**2 - a*(B*b - C*a))/(b*(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(3*C*a**2*d*f - a*b*(2*B*d*f + C*c*f + C*d*e) + b**2*(2*A*d*f + C*c*e))/(2*b**2*f*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*sqrt(e + f*x)*(12*C*a**2*d*f**2 - a*b*f*(8*B*d*f + C*c*f + 7*C*d*e) + b**2*(-C*e*(-c*f + d*e) + 4*d*f*(A*f + B*e)))/(4*b**3*d*f*(-a*f + b*e)) + (6*C*a**3*d*f - a**2*b*(4*B*d*f + 5*C*(c*f + d*e)) + a*b**2*(2*A*d*f + 3*B*c*f + 3*B*d*e + 4*C*c*e) - b**3*(A*c*f + A*d*e + 2*B*c*e))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(b**4*sqrt(-a*d + b*c)*sqrt(-a*f + b*e)) + (24*C*a**2*d**2*f**2 - 8*a*b*d*f*(2*B*d*f + C*c*f + C*d*e) - b**2*(C*(-c*f + d*e)**2 - 4*d*f*(2*A*d*f + B*c*f + B*d*e)))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(4*b**4*d**(sympy.S(3)/2)*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_39():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**3
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(A*b**2 - a*(B*b - C*a))/(2*b*(a + b*x)**2*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(6*C*a**3*d*f - a**2*b*(2*B*d*f + 7*C*(c*f + d*e)) + a*b**2*(-2*A*d*f + 3*B*c*f + 3*B*d*e + 8*C*c*e) - b**3*(-A*c*f - A*d*e + 4*B*c*e))/(4*b**2*(a + b*x)*(-a*d + b*c)*(-a*f + b*e)**2) - sqrt(c + d*x)*sqrt(e + f*x)*(12*C*a**3*d*f**2 - a**2*b*f*(4*B*d*f + 11*C*c*f + 17*C*d*e) + a*b**2*(B*f*(3*c*f + 5*d*e) + 4*C*e*(4*c*f + d*e)) - b**3*(A*d*e*f + c*(-A*f**2 + 4*B*e*f + 4*C*e**2)))/(4*b**3*(-a*d + b*c)*(-a*f + b*e)**2) - (24*C*a**4*d**2*f**2 - 8*a**3*b*d*f*(B*d*f + 5*C*(c*f + d*e)) + 3*a**2*b**2*(4*B*d*f*(c*f + d*e) + C*(5*c**2*f**2 + 22*c*d*e*f + 5*d**2*e**2)) - 3*a*b**3*(B*d**2*e**2 + c**2*f*(B*f + 8*C*e) + 2*c*d*e*(3*B*f + 4*C*e)) - b**4*(A*d**2*e**2 - c**2*(-A*f**2 + 4*B*e*f + 8*C*e**2) - 2*c*d*e*(A*f + 2*B*e)))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(4*b**4*(-a*d + b*c)**(sympy.S(3)/2)*(-a*f + b*e)**(sympy.S(3)/2)) - (6*C*a*d*f - b*(2*B*d*f + C*c*f + C*d*e))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(b**4*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_40():
    f = (a + b*x)**2*sqrt(c + d*x)*(A + B*x + C*x**2)/sqrt(e + f*x)
    F = C*(a + b*x)**3*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(5*b*d*f) - sqrt(c + d*x)*sqrt(e + f*x)*(16*a**2*d**2*f**2*(-C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(-4*A*d*f + B*c*f + 3*B*d*e)) + 4*a*b*d*f*(C*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3) + 8*d*f*(2*A*d*f*(c*f + 3*d*e) - B*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2))) - b**2*(C*(7*c**4*f**4 + 12*c**3*d*e*f**3 + 18*c**2*d**2*e**2*f**2 + 28*c*d**3*e**3*f + 63*d**4*e**4) + 2*d*f*(8*A*d*f*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) - B*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3))))/(128*d**4*f**5) + (-c*f + d*e)*(16*a**2*d**2*f**2*(-C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(-4*A*d*f + B*c*f + 3*B*d*e)) + 4*a*b*d*f*(C*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3) + 8*d*f*(2*A*d*f*(c*f + 3*d*e) - B*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2))) - b**2*(C*(7*c**4*f**4 + 12*c**3*d*e*f**3 + 18*c**2*d**2*e**2*f**2 + 28*c*d**3*e**3*f + 63*d**4*e**4) + 2*d*f*(8*A*d*f*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) - B*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(128*d**(sympy.S(9)/2)*f**(sympy.S(11)/2)) - (a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(4*C*a*d*f + b*(-10*B*d*f + 7*C*c*f + 9*C*d*e))/(40*b*d**2*f**2) - (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(96*C*a**3*d**3*f**3 + 8*a**2*b*d**2*f**2*(-30*B*d*f + 9*C*c*f + 23*C*d*e) + 20*a*b**2*d*f*(-C*(15*c**2*f**2 + 22*c*d*e*f + 35*d**2*e**2) + 8*d*f*(-6*A*d*f + 3*B*c*f + 5*B*d*e)) + b**3*(C*(105*c**3*f**3 + 145*c**2*d*e*f**2 + 203*c*d**2*e**2*f + 315*d**3*e**3) + 10*d*f*(8*A*d*f*(3*c*f + 5*d*e) - B*(15*c**2*f**2 + 22*c*d*e*f + 35*d**2*e**2))) + 4*b*d*f*x*(8*b*d*f*(-10*A*b*d*f + C*a*c*f + 3*C*a*d*e + 6*C*b*c*e) - (4*C*a*d*f + b*(-10*B*d*f + 7*C*c*f + 9*C*d*e))*(-4*a*d*f + 5*b*c*f + 7*b*d*e)))/(960*b*d**4*f**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_41():
    f = (a + b*x)*sqrt(c + d*x)*(A + B*x + C*x**2)/sqrt(e + f*x)
    F = C*(a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(4*b*d*f) - sqrt(c + d*x)*sqrt(e + f*x)*(8*a*d*f*(-C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(-4*A*d*f + B*c*f + 3*B*d*e)) + b*(C*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3) + 8*d*f*(2*A*d*f*(c*f + 3*d*e) - B*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2))))/(64*d**3*f**4) + (-c*f + d*e)*(8*a*d*f*(-C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(-4*A*d*f + B*c*f + 3*B*d*e)) + b*(C*(5*c**3*f**3 + 9*c**2*d*e*f**2 + 15*c*d**2*e**2*f + 35*d**3*e**3) + 8*d*f*(2*A*d*f*(c*f + 3*d*e) - B*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(64*d**(sympy.S(7)/2)*f**(sympy.S(9)/2)) - (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(24*C*a**2*d**2*f**2 + 8*a*b*d*f*(-6*B*d*f + 3*C*c*f + 5*C*d*e) + b**2*(-C*(15*c**2*f**2 + 22*c*d*e*f + 35*d**2*e**2) + 8*d*f*(-6*A*d*f + 3*B*c*f + 5*B*d*e)) + 4*b*d*f*x*(4*C*a*d*f + b*(-8*B*d*f + 5*C*c*f + 7*C*d*e)))/(96*b*d**3*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_42():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/sqrt(e + f*x)
    F = C*(c + d*x)**(sympy.S(5)/2)*sqrt(e + f*x)/(3*d**2*f) - (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(-6*B*d*f + 7*C*c*f + 5*C*d*e)/(12*d**2*f**2) + sqrt(c + d*x)*sqrt(e + f*x)*(C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(4*A*d*f - B*(c*f + 3*d*e)))/(8*d**2*f**3) - (C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(4*A*d*f - B*(c*f + 3*d*e)))*(-c*f + d*e)*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(8*d**(sympy.S(5)/2)*f**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_43():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)*sqrt(e + f*x))
    F = C*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(2*b*d*f) - sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a*d*f + b*(-4*B*d*f + C*c*f + 3*C*d*e))/(4*b**2*d*f**2) - (2*A*b**2 - 2*a*(B*b - C*a))*sqrt(-a*d + b*c)*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(b**3*sqrt(-a*f + b*e)) + (2*b*d*f*(4*A*b*d*f - C*a*(c*f + 3*d*e)) + (4*C*a*d*f + b*(-4*B*d*f + C*c*f + 3*C*d*e))*(2*a*d*f - b*c*f + b*d*e))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(4*b**3*d**(sympy.S(3)/2)*f**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_44():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**2*sqrt(e + f*x))
    F = -(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(b*(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*sqrt(e + f*x)*(2*C*a**2*d*f - a*b*(B*d*f + C*c*f + C*d*e) + b**2*(A*d*f + C*c*e))/(b**2*f*(-a*d + b*c)*(-a*f + b*e)) + (4*C*a**3*d*f - a**2*b*(2*B*d*f + 3*C*c*f + 5*C*d*e) + a*b**2*(B*c*f + 3*B*d*e + 4*C*c*e) - b**3*(-A*c*f + A*d*e + 2*B*c*e))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(b**3*sqrt(-a*d + b*c)*(-a*f + b*e)**(sympy.S(3)/2)) - (4*C*a*d*f + b*(-2*B*d*f - C*c*f + C*d*e))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(b**3*sqrt(d)*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_45():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**3*sqrt(e + f*x))
    F = 2*C*sqrt(d)*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(b**3*sqrt(f)) - (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(2*b*(a + b*x)**2*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**3*d*f - C*a**2*b*(5*c*f + 7*d*e) + a*b**2*(-4*A*d*f + B*c*f + 3*B*d*e + 8*C*c*e) - b**3*(-3*A*c*f - A*d*e + 4*B*c*e))/(4*b**2*(a + b*x)*(-a*d + b*c)*(-a*f + b*e)**2) - (8*C*a**4*d**2*f**2 - 4*C*a**3*b*d*f*(3*c*f + 5*d*e) + 3*C*a**2*b**2*(c**2*f**2 + 10*c*d*e*f + 5*d**2*e**2) - a*b**3*(c**2*f*(-B*f + 8*C*e) + 2*c*d*(2*A*f**2 - B*e*f + 12*C*e**2) + d**2*e*(-4*A*f + 3*B*e)) - b**4*(A*d**2*e**2 - c**2*(3*A*f**2 - 4*B*e*f + 8*C*e**2) - 2*c*d*e*(-A*f + 2*B*e)))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(4*b**3*(-a*d + b*c)**(sympy.S(3)/2)*(-a*f + b*e)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_46():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**4*sqrt(e + f*x))
    F = -(-c*f + d*e)*(-a**2*(-C*(c**2*f**2 + 2*c*d*e*f + 5*d**2*e**2) + 2*d*f*(-4*A*d*f + B*c*f + 3*B*d*e)) + a*b*(-c**2*f*(-B*f + 4*C*e) - 2*c*d*(6*A*f**2 - 7*B*e*f + 6*C*e**2) + d**2*e*(-4*A*f + B*e)) + b**2*(A*d**2*e**2 + c**2*(5*A*f**2 - 6*B*e*f + 8*C*e**2) - 2*c*d*e*(-A*f + B*e)))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(8*(-a*d + b*c)**(sympy.S(5)/2)*(-a*f + b*e)**(sympy.S(7)/2)) - (c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(3*b*(a + b*x)**3*(-a*d + b*c)*(-a*f + b*e)) - sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a**4*d**2*f**2 - 2*a**3*b*d*f*(-2*B*d*f + 7*C*c*f + 13*C*d*e) - a**2*b**2*(-C*(3*c**2*f**2 + 44*c*d*e*f + 33*d**2*e**2) + 4*d*f*(-2*A*d*f + B*c*f + 4*B*d*e)) - a*b**3*(3*c**2*f*(-B*f + 4*C*e) + 2*c*d*(13*A*f**2 - 14*B*e*f + 30*C*e**2) + d**2*e*(-10*A*f + 3*B*e)) - b**4*(3*A*d**2*e**2 - 3*c**2*(5*A*f**2 - 6*B*e*f + 8*C*e**2) - 2*c*d*e*(-2*A*f + 3*B*e)))/(24*b**2*(a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)**3) + sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**3*d*f - a**2*b*(-2*B*d*f + 7*C*c*f + 9*C*d*e) + a*b**2*(-8*A*d*f + B*c*f + 3*B*d*e + 12*C*c*e) - b**3*(-5*A*c*f - 3*A*d*e + 6*B*c*e))/(12*b**2*(a + b*x)**2*(-a*d + b*c)*(-a*f + b*e)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_47():
    f = (a + b*x)**2*(A + B*x + C*x**2)/(sqrt(c + d*x)*sqrt(e + f*x))
    F = C*(a + b*x)**3*sqrt(c + d*x)*sqrt(e + f*x)/(4*b*d*f) + (16*a**2*d**2*f**2*(C*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) + 4*d*f*(2*A*d*f - B*(c*f + d*e))) - 16*a*b*d*f*(C*(5*c**3*f**3 + 3*c**2*d*e*f**2 + 3*c*d**2*e**2*f + 5*d**3*e**3) + 2*d*f*(4*A*d*f*(c*f + d*e) - B*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2))) + b**2*(C*(35*c**4*f**4 + 20*c**3*d*e*f**3 + 18*c**2*d**2*e**2*f**2 + 20*c*d**3*e**3*f + 35*d**4*e**4) + 8*d*f*(2*A*d*f*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) - B*(5*c**3*f**3 + 3*c**2*d*e*f**2 + 3*c*d**2*e**2*f + 5*d**3*e**3))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(64*d**(sympy.S(9)/2)*f**(sympy.S(9)/2)) - (a + b*x)**2*sqrt(c + d*x)*sqrt(e + f*x)*(2*C*a*d*f - b*(8*B*d*f - 7*C*(c*f + d*e)))/(24*b*d**2*f**2) - sqrt(c + d*x)*sqrt(e + f*x)*(32*C*a**3*d**3*f**3 - 8*a**2*b*d**2*f**2*(16*B*d*f - 11*C*(c*f + d*e)) - 16*a*b**2*d*f*(C*(15*c**2*f**2 + 14*c*d*e*f + 15*d**2*e**2) + 6*d*f*(4*A*d*f - 3*B*(c*f + d*e))) + b**3*(5*C*(21*c**3*f**3 + 19*c**2*d*e*f**2 + 19*c*d**2*e**2*f + 21*d**3*e**3) + 8*d*f*(18*A*d*f*(c*f + d*e) - B*(15*c**2*f**2 + 14*c*d*e*f + 15*d**2*e**2))) + 2*b*d*f*x*(6*b*d*f*(-8*A*b*d*f + C*a*c*f + C*a*d*e + 6*C*b*c*e) + (4*a*d*f - 5*b*(c*f + d*e))*(2*C*a*d*f - b*(8*B*d*f - 7*C*(c*f + d*e)))))/(192*b*d**4*f**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_48():
    f = (a + b*x)*(A + B*x + C*x**2)/(sqrt(c + d*x)*sqrt(e + f*x))
    F = C*(a + b*x)**2*sqrt(c + d*x)*sqrt(e + f*x)/(3*b*d*f) + (2*a*d*f*(C*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) + 4*d*f*(2*A*d*f - B*(c*f + d*e))) - b*(C*(5*c**3*f**3 + 3*c**2*d*e*f**2 + 3*c*d**2*e**2*f + 5*d**3*e**3) + 2*d*f*(4*A*d*f*(c*f + d*e) - B*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2))))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(8*d**(sympy.S(7)/2)*f**(sympy.S(7)/2)) - sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a**2*d**2*f**2 - 6*a*b*d*f*(4*B*d*f - 3*C*(c*f + d*e)) - b**2*(C*(15*c**2*f**2 + 14*c*d*e*f + 15*d**2*e**2) + 6*d*f*(4*A*d*f - 3*B*(c*f + d*e))) + 2*b*d*f*x*(2*C*a*d*f - b*(6*B*d*f - 5*C*(c*f + d*e))))/(24*b*d**3*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_49():
    f = (A + B*x + C*x**2)/(sqrt(c + d*x)*sqrt(e + f*x))
    F = C*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(2*d**2*f) - sqrt(c + d*x)*sqrt(e + f*x)*(-4*B*d*f + 5*C*c*f + 3*C*d*e)/(4*d**2*f**2) + (C*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) + 4*d*f*(2*A*d*f - B*(c*f + d*e)))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(4*d**(sympy.S(5)/2)*f**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_50():
    f = (A + B*x + C*x**2)/((a + b*x)*sqrt(c + d*x)*sqrt(e + f*x))
    F = C*sqrt(c + d*x)*sqrt(e + f*x)/(b*d*f) - (2*A*b**2 - 2*a*(B*b - C*a))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(b**2*sqrt(-a*d + b*c)*sqrt(-a*f + b*e)) - (2*C*a*d*f + b*(-2*B*d*f + C*c*f + C*d*e))*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(b**2*d**(sympy.S(3)/2)*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_51():
    f = (A + B*x + C*x**2)/((a + b*x)**2*sqrt(c + d*x)*sqrt(e + f*x))
    F = 2*C*atanh(sqrt(f)*sqrt(c + d*x)/(sqrt(d)*sqrt(e + f*x)))/(b**2*sqrt(d)*sqrt(f)) - sqrt(c + d*x)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(b*(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) + (2*C*a**3*d*f - 3*C*a**2*b*(c*f + d*e) + a*b**2*(-2*A*d*f + B*c*f + B*d*e + 4*C*c*e) - b**3*(-A*c*f - A*d*e + 2*B*c*e))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(b**2*(-a*d + b*c)**(sympy.S(3)/2)*(-a*f + b*e)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_52():
    f = (A + B*x + C*x**2)/((a + b*x)**3*sqrt(c + d*x)*sqrt(e + f*x))
    F = -(a**2*(C*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) + 4*d*f*(2*A*d*f - B*(c*f + d*e))) + a*b*(-c**2*f*(-B*f + 8*C*e) - 2*c*d*(4*A*f**2 - 7*B*e*f + 4*C*e**2) + d**2*e*(-8*A*f + B*e)) + b**2*(3*A*d**2*e**2 + c**2*(3*A*f**2 - 4*B*e*f + 8*C*e**2) - 2*c*d*e*(-A*f + 2*B*e)))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(4*(-a*d + b*c)**(sympy.S(5)/2)*(-a*f + b*e)**(sympy.S(5)/2)) + sqrt(c + d*x)*sqrt(e + f*x)*(2*C*a**3*d*f + a**2*b*(2*B*d*f - 5*C*(c*f + d*e)) + a*b**2*(-6*A*d*f + B*c*f + B*d*e + 8*C*c*e) - b**3*(-3*A*(c*f + d*e) + 4*B*c*e))/(4*b*(a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)**2) - sqrt(c + d*x)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(2*b*(a + b*x)**2*(-a*d + b*c)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_53():
    f = (A + B*x + C*x**2)/((a + b*x)**4*sqrt(c + d*x)*sqrt(e + f*x))
    F = (-2*a**3*d*f*(C*(3*c**2*f**2 + 2*c*d*e*f + 3*d**2*e**2) + 4*d*f*(2*A*d*f - B*(c*f + d*e))) + a**2*b*(C*(c**3*f**3 + 23*c**2*d*e*f**2 + 23*c*d**2*e**2*f + d**3*e**3) + 4*d*f*(6*A*d*f*(c*f + d*e) - B*(c**2*f**2 + 10*c*d*e*f + d**2*e**2))) + a*b**2*(-c**3*f**2*(-B*f + 4*C*e) - c**2*d*f*(18*A*f**2 - 23*B*e*f + 40*C*e**2) - c*d**2*e*(12*A*f**2 - 23*B*e*f + 4*C*e**2) + d**3*e**2*(-18*A*f + B*e)) + b**3*(5*A*d**3*e**3 + c**3*f*(5*A*f**2 - 6*B*e*f + 8*C*e**2) + c**2*d*e*(3*A*f**2 - 4*B*e*f + 8*C*e**2) - 3*c*d**2*e**2*(-A*f + 2*B*e)))*atanh(sqrt(c + d*x)*sqrt(-a*f + b*e)/(sqrt(e + f*x)*sqrt(-a*d + b*c)))/(8*(-a*d + b*c)**(sympy.S(7)/2)*(-a*f + b*e)**(sympy.S(7)/2)) + sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**4*d**2*f**2 + 8*a**3*b*d*f*(B*d*f - 2*C*(c*f + d*e)) - a**2*b**2*(C*(3*c**2*f**2 - 34*c*d*e*f + 3*d**2*e**2) + 2*d*f*(22*A*d*f - 5*B*(c*f + d*e))) - a*b**3*(-3*c**2*f*(-B*f + 4*C*e) - 2*c*d*(22*A*f**2 - 29*B*e*f + 6*C*e**2) + d**2*e*(-44*A*f + 3*B*e)) - b**4*(15*A*d**2*e**2 + 3*c**2*(5*A*f**2 - 6*B*e*f + 8*C*e**2) - 2*c*d*e*(-7*A*f + 9*B*e)))/(24*b*(a + b*x)*(-a*d + b*c)**3*(-a*f + b*e)**3) + sqrt(c + d*x)*sqrt(e + f*x)*(2*C*a**3*d*f + a**2*b*(4*B*d*f - 7*C*(c*f + d*e)) + a*b**2*(-10*A*d*f + B*c*f + B*d*e + 12*C*c*e) - b**3*(-5*A*(c*f + d*e) + 6*B*c*e))/(12*b*(a + b*x)**2*(-a*d + b*c)**2*(-a*f + b*e)**2) - sqrt(c + d*x)*sqrt(e + f*x)*(A*b**2 - a*(B*b - C*a))/(3*b*(a + b*x)**3*(-a*d + b*c)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_54():
    f = sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)
    F = 2*C*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)/(9*b*d*f) - sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(4*C*a*d*f - 2*b*(3*B*d*f - 2*C*(c*f + d*e)))/(21*b*d**2*f**2) - sqrt(a + b*x)*sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(14*b*d*f*(-3*A*b*d*f + C*a*c*f + C*a*d*e + C*b*c*e) + 2*(a*d*f - 4*b*(c*f + d*e))*(2*C*a*d*f - b*(3*B*d*f - 2*C*(c*f + d*e))))/(105*b**2*d**2*f**3) + sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(16*C*a**3*d**3*f**3 + 6*a**2*b*d**2*f**2*(-4*B*d*f - C*c*f + C*d*e) - 6*a*b**2*d*f**2*(B*d*(-2*c*f + d*e) + f*(-7*A*d**2 + C*c**2)) - 2*b**3*(C*(-8*c**3*f**3 - 3*c**2*d*e*f**2 + 16*d**3*e**3) + 3*d*f*(7*A*d*f*(-c*f + 2*d*e) - B*(-4*c**2*f**2 - c*d*e*f + 8*d**2*e**2))))/(315*b**3*d**3*f**3) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(-c*f + d*e)*(8*C*a**3*d**3*f**3 + 3*a**2*b*d**2*f**2*(-4*B*d*f - C*c*f + C*d*e) - 3*a*b**2*d*f**2*(B*d*(-2*c*f + d*e) + f*(-7*A*d**2 + C*c**2)) - b**3*(C*(-8*c**3*f**3 - 3*c**2*d*e*f**2 + 16*d**3*e**3) + 3*d*f*(7*A*d*f*(-c*f + 2*d*e) - B*(-4*c**2*f**2 - c*d*e*f + 8*d**2*e**2))))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(315*b**4*d**(sympy.S(7)/2)*f**4*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(16*C*a**4*d**4*f**4 - 8*a**3*b*d**3*f**3*(3*B*d*f + C*c*f + C*d*e) + 3*a**2*b**2*d**2*f**2*(-2*C*(c**2*f**2 - c*d*e*f + d**2*e**2) + d*f*(14*A*d*f + 5*B*c*f + 5*B*d*e)) - a*b**3*d*f*(C*(8*c**3*f**3 - 6*c**2*d*e*f**2 - 6*c*d**2*e**2*f + 8*d**3*e**3) + 3*d*f*(14*A*d*f*(c*f + d*e) - B*(5*c**2*f**2 - 6*c*d*e*f + 5*d**2*e**2))) + b**4*(2*C*(8*c**4*f**4 - 4*c**3*d*e*f**3 - 3*c**2*d**2*e**2*f**2 - 4*c*d**3*e**3*f + 8*d**4*e**4) + 3*d*f*(14*A*d*f*(c**2*f**2 - c*d*e*f + d**2*e**2) - B*(8*c**3*f**3 - 5*c**2*d*e*f**2 - 5*c*d**2*e**2*f + 8*d**3*e**3))))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(315*b**4*d**(sympy.S(7)/2)*f**4*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_55():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/sqrt(a + b*x)
    F = 2*C*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)/(7*b*d*f) - sqrt(a + b*x)*sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(12*C*a*d*f - 2*b*(7*B*d*f - 4*C*(c*f + d*e)))/(35*b**2*d*f**2) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(10*b*d*f*(3*C*a*(c*f + d*e) + b*(-7*A*d*f + C*c*e)) - 2*(6*C*a*d*f - b*(7*B*d*f - 4*C*(c*f + d*e)))*(4*a*d*f - b*c*f + 2*b*d*e))/(105*b**3*d**2*f**2) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(-c*f + d*e)*(24*C*a**2*d**2*f**2 + a*b*d*f*(-28*B*d*f - 5*C*c*f + 13*C*d*e) - b**2*(-C*(-4*c**2*f**2 - c*d*e*f + 8*d**2*e**2) + 7*d*f*(-5*A*d*f - B*c*f + 2*B*d*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**4*d**(sympy.S(5)/2)*f**3*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(3*b*d*f*(5*b*c*f*(3*C*a*(c*f + d*e) + b*(-7*A*d*f + C*c*e)) - (6*C*a*d*f - b*(7*B*d*f - 4*C*(c*f + d*e)))*(3*a*c*f + a*d*e + b*c*e)) + (b*d*e - 2*f*(a*d + b*c))*(5*b*d*f*(3*C*a*(c*f + d*e) + b*(-7*A*d*f + C*c*e)) - (6*C*a*d*f - b*(7*B*d*f - 4*C*(c*f + d*e)))*(4*a*d*f - b*c*f + 2*b*d*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**4*d**(sympy.S(5)/2)*f**3*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_56():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**(sympy.S(3)/2)
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(b*sqrt(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(a + b*x)*sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(12*C*a**2*d*f - 2*a*b*(5*B*d*f + C*c*f + C*d*e) + 2*b**2*(5*A*d*f + C*c*e))/(5*b**2*f*(-a*d + b*c)*(-a*f + b*e)) + sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(48*C*a**2*d*f**2 - 2*a*b*f*(20*B*d*f + C*c*f + 7*C*d*e) + 2*b**2*(-C*e*(-c*f + 2*d*e) + 5*d*f*(3*A*f + B*e)))/(15*b**3*d*f*(-a*f + b*e)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-c*f + d*e)*(24*C*a**2*d*f**2 - a*b*f*(20*B*d*f + C*c*f + 7*C*d*e) + b**2*(-C*e*(-c*f + 2*d*e) + 5*d*f*(3*A*f + B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**4*d**(sympy.S(3)/2)*f**2*sqrt(c + d*x)*sqrt(e + f*x)) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(48*C*a**2*d**2*f**2 - 8*a*b*d*f*(5*B*d*f + C*c*f + C*d*e) + b**2*(-2*C*(c**2*f**2 - c*d*e*f + d**2*e**2) + 5*d*f*(6*A*d*f + B*c*f + B*d*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**4*d**(sympy.S(3)/2)*f**2*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_57():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**(sympy.S(5)/2)
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*f + b*e)) - sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(2*B*b - 4*C*a)/(b**2*sqrt(a + b*x)*(-a*f + b*e)) + sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(16*C*a**2*d*f - 2*a*b*(4*B*d*f + 7*C*c*f + C*d*e) + 2*b**2*(A*d*f + 3*B*c*f + C*c*e))/(3*b**3*(-a*d + b*c)*(-a*f + b*e)) + sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(-2*c*f + 2*d*e)*(8*C*a**2*d*f - a*b*(4*B*d*f + 7*C*c*f + C*d*e) + b**2*(A*d*f + 3*B*c*f + C*c*e))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**4*sqrt(d)*f*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(a*d - b*c)) + sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(32*C*a**3*d**2*f**2 - 16*a**2*b*d*f*(B*d*f + 2*C*(c*f + d*e)) + 2*a*b**2*(C*(c**2*f**2 + 16*c*d*e*f + d**2*e**2) + d*f*(2*A*d*f + 7*B*c*f + 7*B*d*e)) - 2*b**3*(A*d**2*e*f + C*c**2*e*f + c*d*(A*f**2 + 6*B*e*f + C*e**2)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**4*sqrt(d)*f*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*sqrt(a*d - b*c)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_58():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**(sympy.S(7)/2)
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(12*C*a**3*d*f - 2*a**2*b*(B*d*f + 8*C*(c*f + d*e)) + 2*a*b**2*(-4*A*d*f + 3*B*c*f + 3*B*d*e + 10*C*c*e) - 2*b**3*(-2*A*(c*f + d*e) + 5*B*c*e))/(15*b**2*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*f + b*e)**2) + sqrt(c + d*x)*sqrt(e + f*x)*(48*C*a**3*d**2*f - 2*a**2*b*d*(4*B*d*f + 41*C*c*f + 23*C*d*e) + 2*a*b**2*(15*C*c**2*f + c*(6*B*d*f + 40*C*d*e) + d**2*(-A*f + 3*B*e)) - 2*b**3*(-2*A*d**2*e + 15*C*c**2*e + c*d*(A*f + 5*B*e)))/(15*b**3*sqrt(a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(48*C*a**4*d**2*f**2 - 8*a**3*b*d*f*(B*d*f + 11*C*(c*f + d*e)) + a**2*b**2*(2*C*(19*c**2*f**2 + 81*c*d*e*f + 19*d**2*e**2) - d*f*(2*A*d*f - 13*B*(c*f + d*e))) - a*b**3*(c**2*f*(3*B*f + 70*C*e) + 2*c*d*(-A*f**2 + 11*B*e*f + 35*C*e**2) + d**2*e*(-2*A*f + 3*B*e)) - b**4*(2*A*d**2*e**2 - c**2*(-2*A*f**2 + 5*B*e*f + 30*C*e**2) - c*d*e*(2*A*f + 5*B*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**4*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e)**2) + sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(-2*c*f + 2*d*e)*(24*C*a**3*d**2*f - a**2*b*d*(4*B*d*f + 41*C*c*f + 23*C*d*e) + a*b**2*(15*C*c**2*f + c*(6*B*d*f + 40*C*d*e) + d**2*(-A*f + 3*B*e)) - b**3*(-2*A*d**2*e + 15*C*c**2*e + c*d*(A*f + 5*B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**4*sqrt(d)*sqrt(c + d*x)*sqrt(e + f*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_59():
    f = sqrt(c + d*x)*sqrt(e + f*x)*(A + B*x + C*x**2)/(a + b*x)**(sympy.S(9)/2)
    F = -(c + d*x)**(sympy.S(3)/2)*(e + f*x)**(sympy.S(3)/2)*(2*A*b**2 - 2*a*(B*b - C*a))/(7*b*(a + b*x)**(sympy.S(7)/2)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(c + d*x)*(e + f*x)**(sympy.S(3)/2)*(12*C*a**3*d*f + 2*a**2*b*(B*d*f - 10*C*(c*f + d*e)) + 2*a*b**2*(-8*A*d*f + 3*B*c*f + 3*B*d*e + 14*C*c*e) - 2*b**3*(-4*A*(c*f + d*e) + 7*B*c*e))/(35*b**2*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*f + b*e)**2) + sqrt(c + d*x)*sqrt(e + f*x)*(96*C*a**5*d**3*f**3 + 16*a**4*b*d**2*f**2*(B*d*f - 16*C*(c*f + d*e)) + 2*a**3*b**2*d*f*(C*(103*c**2*f**2 + 344*c*d*e*f + 103*d**2*e**2) + d*f*(6*A*d*f - 19*B*(c*f + d*e))) - 6*a**2*b**3*(C*(5*c**3*f**3 + 94*c**2*d*e*f**2 + 94*c*d**2*e**2*f + 5*d**3*e**3) + d*f*(3*A*d*f*(c*f + d*e) - B*(3*c**2*f**2 + 16*c*d*e*f + 3*d**2*e**2))) - 2*a*b**4*(-6*c**3*f**2*(-B*f + 7*C*e) - c**2*d*f*(238*C*e**2 - 19*f*(-A*f + B*e)) - c*d**2*e*(42*C*e**2 - f*(20*A*f + 19*B*e)) + d**3*e**2*(-19*A*f + 6*B*e)) - 2*b**5*(8*A*d**3*e**3 + c**3*f*(8*A*f**2 - 14*B*e*f + 35*C*e**2) + c**2*d*e*(-5*A*f**2 + 14*B*e*f + 35*C*e**2) - c*d**2*e**2*(5*A*f + 14*B*e)))/(105*b**3*sqrt(a + b*x)*(-a*d + b*c)**3*(-a*f + b*e)**3) - sqrt(c + d*x)*sqrt(e + f*x)*(48*C*a**4*d**2*f**2 - 2*a**3*b*d*f*(-4*B*d*f + 43*C*c*f + 61*C*d*e) - 6*a**2*b**2*(-C*(5*c**2*f**2 + 37*c*d*e*f + 15*d**2*e**2) + d*f*(-A*d*f + 2*B*c*f + 3*B*d*e)) - 6*a*b**3*(2*c**2*f*(-B*f + 7*C*e) + c*d*(5*A*f**2 - 5*B*e*f + 28*C*e**2) + d**2*e*(-3*A*f + B*e)) - 2*b**4*(4*A*d**2*e**2 - c**2*(8*A*f**2 - 14*B*e*f + 35*C*e**2) - c*d*e*(-A*f + 7*B*e)))/(105*b**3*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**2*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(-c*f + d*e)*(24*C*a**4*d**2*f**2 - a**3*b*d*f*(-4*B*d*f + 61*C*c*f + 43*C*d*e) - 3*a**2*b**2*(-C*(15*c**2*f**2 + 37*c*d*e*f + 5*d**2*e**2) + d*f*(-A*d*f + 3*B*c*f + 2*B*d*e)) + 3*a*b**3*(-c**2*f*(B*f + 28*C*e) - c*d*(-3*A*f**2 - 5*B*e*f + 14*C*e**2) + d**2*e*(-5*A*f + 2*B*e)) + b**4*(8*A*d**2*e**2 + c**2*(-4*A*f**2 + 7*B*e*f + 35*C*e**2) - c*d*e*(A*f + 14*B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**4*sqrt(c + d*x)*sqrt(e + f*x)*(a*d - b*c)**(sympy.S(5)/2)*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(48*C*a**5*d**3*f**3 + 8*a**4*b*d**2*f**2*(B*d*f - 16*C*(c*f + d*e)) + a**3*b**2*d*f*(C*(103*c**2*f**2 + 344*c*d*e*f + 103*d**2*e**2) + d*f*(6*A*d*f - 19*B*(c*f + d*e))) - 3*a**2*b**3*(C*(5*c**3*f**3 + 94*c**2*d*e*f**2 + 94*c*d**2*e**2*f + 5*d**3*e**3) + d*f*(3*A*d*f*(c*f + d*e) - B*(3*c**2*f**2 + 16*c*d*e*f + 3*d**2*e**2))) - a*b**4*(-6*c**3*f**2*(-B*f + 7*C*e) - c**2*d*f*(238*C*e**2 - 19*f*(-A*f + B*e)) - c*d**2*e*(42*C*e**2 - f*(20*A*f + 19*B*e)) + d**3*e**2*(-19*A*f + 6*B*e)) - b**5*(8*A*d**3*e**3 + c**3*f*(8*A*f**2 - 14*B*e*f + 35*C*e**2) + c**2*d*e*(-5*A*f**2 + 14*B*e*f + 35*C*e**2) - c*d**2*e**2*(5*A*f + 14*B*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**4*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(a*d - b*c)**(sympy.S(5)/2)*(-a*f + b*e)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_60():
    f = (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*(A + B*x + C*x**2)/sqrt(e + f*x)
    F = 2*C*(a + b*x)**(sympy.S(5)/2)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(9*b*d*f) - (a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(8*C*a*d*f + 2*b*(-9*B*d*f + 6*C*c*f + 8*C*d*e))/(63*b*d**2*f**2) - sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(14*b*d*f*(-9*A*b*d*f + C*a*c*f + 3*C*a*d*e + 5*C*b*c*e) - 2*(4*C*a*d*f + b*(-9*B*d*f + 6*C*c*f + 8*C*d*e))*(-3*a*d*f + 4*b*c*f + 6*b*d*e))/(315*b*d**3*f**3) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(10*b*d*f*(7*a*d*f*(-9*A*b*d*f + C*a*c*f + 3*C*a*d*e + 5*C*b*c*e) - (4*C*a*d*f + b*(-9*B*d*f + 6*C*c*f + 8*C*d*e))*(a*c*f + 3*a*d*e + 3*b*c*e)) + 2*(a*d*f - 2*b*(c*f + 2*d*e))*(7*b*d*f*(-9*A*b*d*f + C*a*c*f + 3*C*a*d*e + 5*C*b*c*e) - (4*C*a*d*f + b*(-9*B*d*f + 6*C*c*f + 8*C*d*e))*(-3*a*d*f + 4*b*c*f + 6*b*d*e)))/(945*b**2*d**3*f**4) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(-c*f + d*e)*(4*C*a**3*d**3*f**3 + 3*a**2*b*d**2*f**2*(-3*B*d*f - C*c*f + 3*C*d*e) - 3*a*b**2*d*f*(-5*C*(c**2*f**2 + 2*c*d*e*f + 8*d**2*e**2) + 3*d*f*(-21*A*d*f + 3*B*c*f + 16*B*d*e)) - b**3*(C*(8*c**3*f**3 + 15*c**2*d*e*f**2 + 24*c*d**2*e**2*f + 128*d**3*e**3) + 3*d*f*(7*A*d*f*(c*f + 8*d*e) - 4*B*(c**2*f**2 + 2*c*d*e*f + 12*d**2*e**2))))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(315*b**3*d**(sympy.S(7)/2)*f**5*sqrt(c + d*x)*sqrt(e + f*x)) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(8*C*a**4*d**4*f**4 + a**3*b*d**3*f**3*(-18*B*d*f - 7*C*c*f + 11*C*d*e) - 3*a**2*b**2*d**2*f**2*(-C*(-3*c**2*f**2 - 5*c*d*e*f + 9*d**2*e**2) + 3*d*f*(-7*A*d*f - 3*B*c*f + 4*B*d*e)) - a*b**3*d*f*(2*C*(-16*c**3*f**3 - 18*c**2*d*e*f**2 - 33*c*d**2*e**2*f + 92*d**3*e**3) + 3*d*f*(7*A*d*f*(-7*c*f + 13*d*e) - B*(-19*c**2*f**2 - 29*c*d*e*f + 72*d**2*e**2))) + b**4*(C*(-16*c**4*f**4 - 16*c**3*d*e*f**3 - 21*c**2*d**2*e**2*f**2 - 40*c*d**3*e**3*f + 128*d**4*e**4) + 3*d*f*(7*A*d*f*(-2*c**2*f**2 - 3*c*d*e*f + 8*d**2*e**2) - B*(-8*c**3*f**3 - 9*c**2*d*e*f**2 - 16*c*d**2*e**2*f + 48*d**3*e**3))))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(315*b**3*d**(sympy.S(7)/2)*f**5*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_61():
    f = sqrt(a + b*x)*sqrt(c + d*x)*(A + B*x + C*x**2)/sqrt(e + f*x)
    F = 2*C*(a + b*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(7*b*d*f) - sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(8*C*a*d*f + 2*b*(-7*B*d*f + 4*C*c*f + 6*C*d*e))/(35*b*d**2*f**2) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(10*b*d*f*(-7*A*b*d*f + C*a*c*f + 3*C*a*d*e + 3*C*b*c*e) + 2*(a*d*f - 2*b*(c*f + 2*d*e))*(4*C*a*d*f + b*(-7*B*d*f + 4*C*c*f + 6*C*d*e)))/(105*b**2*d**2*f**3) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(-c*f + d*e)*(4*C*a**2*d**2*f**2 + a*b*d*f*(-7*B*d*f - 2*C*c*f + 8*C*d*e) - b**2*(-4*C*(c**2*f**2 + 2*c*d*e*f + 12*d**2*e**2) + 7*d*f*(-10*A*d*f + B*c*f + 8*B*d*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**3*d**(sympy.S(5)/2)*f**4*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(3*b*d*f*(5*a*d*f*(-7*A*b*d*f + C*a*c*f + 3*C*a*d*e + 3*C*b*c*e) - (4*C*a*d*f + b*(-7*B*d*f + 4*C*c*f + 6*C*d*e))*(a*c*f + 3*a*d*e + b*c*e)) + (b*c*f - 2*d*(a*f + b*e))*(5*b*d*f*(-7*A*b*d*f + C*a*c*f + 3*C*a*d*e + 3*C*b*c*e) + (a*d*f - 2*b*(c*f + 2*d*e))*(4*C*a*d*f + b*(-7*B*d*f + 4*C*c*f + 6*C*d*e))))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**3*d**(sympy.S(5)/2)*f**4*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_62():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(e + f*x))
    F = 2*C*sqrt(a + b*x)*(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(5*b*d*f) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a*d*f + 2*b*(-5*B*d*f + 2*C*c*f + 4*C*d*e))/(15*b**2*d*f**2) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-c*f + d*e)*(4*C*a**2*d*f**2 + a*b*f*(-5*B*d*f - C*c*f + 3*C*d*e) - b**2*(-C*e*(c*f + 8*d*e) + 5*d*f*(-3*A*f + 2*B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**3*d**(sympy.S(3)/2)*f**3*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(3*b*d*f*(-5*A*b*d*f + C*a*c*f + 3*C*a*d*e + C*b*c*e) - (4*C*a*d*f + b*(-5*B*d*f + 2*C*c*f + 4*C*d*e))*(2*a*d*f - b*c*f + 2*b*d*e))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**3*d**(sympy.S(3)/2)*f**3*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_63():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**(sympy.S(3)/2)*sqrt(e + f*x))
    F = -(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(b*sqrt(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) + sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a**2*d*f - 2*a*b*(3*B*d*f + C*c*f + C*d*e) + 2*b**2*(3*A*d*f + C*c*e))/(3*b**2*f*(-a*d + b*c)*(-a*f + b*e)) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-c*f + d*e)*(-3*B*b*f + 4*C*a*f + 2*C*b*e)*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**3*sqrt(d)*f**2*sqrt(c + d*x)*sqrt(e + f*x)) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(8*C*a**2*d*f**2 - a*b*f*(6*B*d*f + C*c*f + 3*C*d*e) + b**2*(-C*e*(-c*f + 2*d*e) + 3*d*f*(A*f + B*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**3*sqrt(d)*f**2*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_64():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**(sympy.S(5)/2)*sqrt(e + f*x))
    F = -(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*f + b*e)) - sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a**2*f - 2*a*b*(B*f + 6*C*e) + 2*b**2*(-2*A*f + 3*B*e))/(3*b**2*sqrt(a + b*x)*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(8*C*a**3*d*f**2 - a**2*b*f*(2*B*d*f + 7*C*c*f + 13*C*d*e) + a*b**2*(3*C*e*(4*c*f + d*e) + f*(-A*d*f + B*c*f + 4*B*d*e)) - b**3*(A*d*e*f + c*(-2*A*f**2 + 3*B*e*f + 3*C*e**2)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**3*f*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*sqrt(a*d - b*c)*(-a*f + b*e)**2) + sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(-2*c*f + 2*d*e)*(4*C*a**2*d*f - a*b*(B*d*f + 3*C*(c*f + d*e)) + b**2*(A*d*f + 3*C*c*e))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**3*sqrt(d)*f*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(a*d - b*c)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_65():
    f = sqrt(c + d*x)*(A + B*x + C*x**2)/((a + b*x)**(sympy.S(7)/2)*sqrt(e + f*x))
    F = -(c + d*x)**(sympy.S(3)/2)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*f + b*e)) - sqrt(c + d*x)*sqrt(e + f*x)*(16*C*a**4*d**2*f**2 - 2*a**3*b*d*f*(-2*B*d*f + 13*C*c*f + 23*C*d*e) - 2*a**2*b**2*(-C*(3*c**2*f**2 + 37*c*d*e*f + 23*d**2*e**2) + d*f*(-3*A*d*f + 2*B*c*f + 7*B*d*e)) - 2*a*b**3*(2*c**2*f*(-B*f + 5*C*e) + c*d*(40*C*e**2 - 13*f*(-A*f + B*e)) + d**2*e*(-7*A*f + 3*B*e)) - 2*b**4*(2*A*d**2*e**2 - c**2*(8*A*f**2 - 10*B*e*f + 15*C*e**2) - c*d*e*(-3*A*f + 5*B*e)))/(15*b**2*sqrt(a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)**3) + sqrt(c + d*x)*sqrt(e + f*x)*(8*C*a**3*d*f - 2*a**2*b*(-B*d*f + 6*C*c*f + 8*C*d*e) + 2*a*b**2*(-6*A*d*f + B*c*f + 3*B*d*e + 10*C*c*e) - 2*b**3*(-4*A*c*f - 2*A*d*e + 5*B*c*e))/(15*b**2*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(-c*f + d*e)*(4*C*a**3*d*f - a**2*b*(-B*d*f + 6*C*c*f + 8*C*d*e) + a*b**2*(-6*A*d*f + B*c*f + 3*B*d*e + 10*C*c*e) - b**3*(-4*A*c*f - 2*A*d*e + 5*B*c*e))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**3*sqrt(c + d*x)*sqrt(e + f*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(8*C*a**4*d**2*f**2 - a**3*b*d*f*(-2*B*d*f + 13*C*c*f + 23*C*d*e) - a**2*b**2*(-C*(3*c**2*f**2 + 37*c*d*e*f + 23*d**2*e**2) + d*f*(-3*A*d*f + 2*B*c*f + 7*B*d*e)) - a*b**3*(2*c**2*f*(-B*f + 5*C*e) + c*d*(40*C*e**2 - 13*f*(-A*f + B*e)) + d**2*e*(-7*A*f + 3*B*e)) - b**4*(2*A*d**2*e**2 - c**2*(8*A*f**2 - 10*B*e*f + 15*C*e**2) - c*d*e*(-3*A*f + 5*B*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**3*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_66():
    f = (a + b*x)**(sympy.S(3)/2)*(A + B*x + C*x**2)/(sqrt(c + d*x)*sqrt(e + f*x))
    F = 2*C*(a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*sqrt(e + f*x)/(7*b*d*f) - (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a*d*f - 2*b*(7*B*d*f - 6*C*(c*f + d*e)))/(35*b*d**2*f**2) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(10*b*d*f*(-7*A*b*d*f + C*a*c*f + C*a*d*e + 5*C*b*c*e) + 2*(3*a*d*f - 4*b*(c*f + d*e))*(2*C*a*d*f - b*(7*B*d*f - 6*C*(c*f + d*e))))/(105*b*d**3*f**3) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(3*C*a**2*d**2*f**2*(-c*f + d*e) - 3*a*b*d*f*(-C*(11*c**2*f**2 + 8*c*d*e*f + 16*d**2*e**2) + 7*d*f*(-5*A*d*f + 2*B*c*f + 3*B*d*e)) - b**2*(C*(24*c**3*f**3 + 17*c**2*d*e*f**2 + 16*c*d**2*e**2*f + 48*d**3*e**3) + 7*d*f*(5*A*d*f*(c*f + 2*d*e) - B*(4*c**2*f**2 + 3*c*d*e*f + 8*d**2*e**2))))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**2*d**(sympy.S(7)/2)*f**4*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(3*b*d*f*(5*a*d*f*(-7*A*b*d*f + C*a*c*f + C*a*d*e + 5*C*b*c*e) - (2*C*a*d*f - b*(7*B*d*f - 6*C*(c*f + d*e)))*(a*c*f + a*d*e + 3*b*c*e)) + (a*d*f - 2*b*(c*f + d*e))*(5*b*d*f*(-7*A*b*d*f + C*a*c*f + C*a*d*e + 5*C*b*c*e) + (3*a*d*f - 4*b*(c*f + d*e))*(2*C*a*d*f - b*(7*B*d*f - 6*C*(c*f + d*e)))))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(105*b**2*d**(sympy.S(7)/2)*f**4*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_67():
    f = sqrt(a + b*x)*(A + B*x + C*x**2)/(sqrt(c + d*x)*sqrt(e + f*x))
    F = 2*C*(a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*sqrt(e + f*x)/(5*b*d*f) - sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a*d*f - 2*b*(5*B*d*f - 4*C*(c*f + d*e)))/(15*b*d**2*f**2) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(-a*f + b*e)*(C*a*d*f*(-c*f + d*e) - b*(-C*(4*c**2*f**2 + 3*c*d*e*f + 8*d**2*e**2) + 5*d*f*(-3*A*d*f + B*c*f + 2*B*d*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**2*d**(sympy.S(5)/2)*f**3*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(3*b*d*f*(-5*A*b*d*f + C*a*c*f + C*a*d*e + 3*C*b*c*e) + (a*d*f - 2*b*(c*f + d*e))*(2*C*a*d*f - b*(5*B*d*f - 4*C*(c*f + d*e))))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**2*d**(sympy.S(5)/2)*f**3*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_68():
    f = (A + B*x + C*x**2)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x))
    F = 2*C*sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)/(3*b*d*f) + 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(a*d - b*c)*(C*a*f*(-c*f + d*e) - b*(-C*e*(c*f + 2*d*e) + 3*d*f*(-A*f + B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**2*d**(sympy.S(3)/2)*f**2*sqrt(c + d*x)*sqrt(e + f*x)) - 2*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*sqrt(a*d - b*c)*(2*C*a*d*f - b*(3*B*d*f - 2*C*(c*f + d*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**2*d**(sympy.S(3)/2)*f**2*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_69():
    f = (A + B*x + C*x**2)/((a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)*sqrt(e + f*x))
    F = -sqrt(c + d*x)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(b*sqrt(a + b*x)*(-a*d + b*c)*(-a*f + b*e)) - sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(2*C*a*(-c*f + d*e) - 2*b*(A*d*f - B*c*f + C*c*e))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(b**2*sqrt(d)*f*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(a*d - b*c)) - sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(4*C*a**2*d*f - 2*a*b*(B*d*f + C*c*f + C*d*e) + 2*b**2*(A*d*f + C*c*e))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(b**2*sqrt(d)*f*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*sqrt(a*d - b*c)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_70():
    f = (A + B*x + C*x**2)/((a + b*x)**(sympy.S(5)/2)*sqrt(c + d*x)*sqrt(e + f*x))
    F = sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**3*d*f + 2*a**2*b*(B*d*f - 4*C*(c*f + d*e)) + 2*a*b**2*(-4*A*d*f + B*c*f + B*d*e + 6*C*c*e) - 2*b**3*(-2*A*(c*f + d*e) + 3*B*c*e))/(3*b*sqrt(a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)**2) - sqrt(c + d*x)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(3*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)*(-a*f + b*e)) - 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(2*C*a**3*d*f + a**2*b*(B*d*f - 4*C*(c*f + d*e)) + a*b**2*(-4*A*d*f + B*c*f + B*d*e + 6*C*c*e) - b**3*(-2*A*(c*f + d*e) + 3*B*c*e))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**2*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e)**2) - sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(2*C*a**2*d*(-c*f + d*e) + 2*a*b*(-B*d*(2*c*f + d*e) + f*(3*A*d**2 + 3*C*c**2)) - 2*b**2*(A*c*d*f + 2*A*d**2*e - 3*B*c*d*e + 3*C*c**2*e))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(3*b**2*sqrt(d)*sqrt(c + d*x)*sqrt(e + f*x)*(a*d - b*c)**(sympy.S(3)/2)*(-a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_6_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_e_plus_f_x_pow_p_71():
    f = (A + B*x + C*x**2)/((a + b*x)**(sympy.S(7)/2)*sqrt(c + d*x)*sqrt(e + f*x))
    F = sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**4*d**2*f**2 + 2*a**3*b*d*f*(3*B*d*f - 7*C*(c*f + d*e)) - 2*a**2*b**2*(C*(3*c**2*f**2 - 13*c*d*e*f + 3*d**2*e**2) + d*f*(23*A*d*f - 7*B*(c*f + d*e))) - 2*a*b**3*(-2*c**2*f*(-B*f + 5*C*e) - c*d*(23*A*f**2 - 33*B*e*f + 10*C*e**2) + d**2*e*(-23*A*f + 2*B*e)) - 2*b**4*(8*A*d**2*e**2 + c**2*(8*A*f**2 - 10*B*e*f + 15*C*e**2) - c*d*e*(-7*A*f + 10*B*e)))/(15*b*sqrt(a + b*x)*(-a*d + b*c)**3*(-a*f + b*e)**3) + sqrt(c + d*x)*sqrt(e + f*x)*(4*C*a**3*d*f + 6*a**2*b*(B*d*f - 2*C*(c*f + d*e)) + 2*a*b**2*(-8*A*d*f + B*c*f + B*d*e + 10*C*c*e) - 2*b**3*(-4*A*(c*f + d*e) + 5*B*c*e))/(15*b*(a + b*x)**(sympy.S(3)/2)*(-a*d + b*c)**2*(-a*f + b*e)**2) - sqrt(c + d*x)*sqrt(e + f*x)*(2*A*b**2 - 2*a*(B*b - C*a))/(5*b*(a + b*x)**(sympy.S(5)/2)*(-a*d + b*c)*(-a*f + b*e)) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*(C*a**3*d*f*(-c*f + d*e) - 3*a**2*b*(-C*(3*c**2*f**2 + c*d*e*f + d**2*e**2) + d*f*(-5*A*d*f + 3*B*c*f + 2*B*d*e)) + a*b**2*(-c**2*f*(-B*f + 20*C*e) - c*d*(11*A*f**2 - 27*B*e*f + 10*C*e**2) + d**2*e*(-19*A*f + 2*B*e)) + b**3*(8*A*d**2*e**2 + c**2*(4*A*f**2 - 5*B*e*f + 15*C*e**2) - c*d*e*(-3*A*f + 10*B*e)))*elliptic_f(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**2*sqrt(c + d*x)*sqrt(e + f*x)*(a*d - b*c)**(sympy.S(5)/2)*(-a*f + b*e)**2) + 2*sqrt(d)*sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*(2*C*a**4*d**2*f**2 + a**3*b*d*f*(3*B*d*f - 7*C*(c*f + d*e)) - a**2*b**2*(C*(3*c**2*f**2 - 13*c*d*e*f + 3*d**2*e**2) + d*f*(23*A*d*f - 7*B*(c*f + d*e))) - a*b**3*(-2*c**2*f*(-B*f + 5*C*e) - c*d*(23*A*f**2 - 33*B*e*f + 10*C*e**2) + d**2*e*(-23*A*f + 2*B*e)) - b**4*(8*A*d**2*e**2 + c**2*(8*A*f**2 - 10*B*e*f + 15*C*e**2) - c*d*e*(-7*A*f + 10*B*e)))*elliptic_e(asin(sqrt(d)*sqrt(a + b*x)/sqrt(a*d - b*c)), f*(-a*d + b*c)/(d*(-a*f + b*e)))/(15*b**2*sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)*(a*d - b*c)**(sympy.S(5)/2)*(-a*f + b*e)**3)
    assert integrate(f, x) == F

