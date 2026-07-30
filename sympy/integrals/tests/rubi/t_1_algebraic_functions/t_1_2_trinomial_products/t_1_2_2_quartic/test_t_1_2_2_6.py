"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.6 P(x) (d x)^m (a+b x^2+c x^4)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, g, m, p = symbols('A B C a b c d e f g m p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1():
    f = x**2*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)
    F = A*a*x**3/3 + B*a*x**4/4 + B*b*x**6/6 + B*c*x**8/8 + C*c*x**9/9 + x**7*(A*c/7 + C*b/7) + x**5*(A*b/5 + C*a/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_2():
    f = x*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)
    F = A*a*x**2/2 + B*a*x**3/3 + B*b*x**5/5 + B*c*x**7/7 + C*c*x**8/8 + x**6*(A*c/6 + C*b/6) + x**4*(A*b/4 + C*a/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_3():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)
    F = A*a*x + B*a*x**2/2 + B*b*x**4/4 + B*c*x**6/6 + C*c*x**7/7 + x**5*(A*c/5 + C*b/5) + x**3*(A*b/3 + C*a/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_4():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x
    F = A*a*log(x) + B*a*x + B*b*x**3/3 + B*c*x**5/5 + C*c*x**6/6 + x**4*(A*c/4 + C*b/4) + x**2*(A*b/2 + C*a/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_5():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**2
    F = -A*a/x + B*a*log(x) + B*b*x**2/2 + B*c*x**4/4 + C*c*x**5/5 + x**3*(A*c/3 + C*b/3) + x*(A*b + C*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_6():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**3
    F = -A*a/(2*x**2) - B*a/x + B*b*x + B*c*x**3/3 + C*c*x**4/4 + x**2*(A*c/2 + C*b/2) + (A*b + C*a)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_7():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**4
    F = -A*a/(3*x**3) - B*a/(2*x**2) + B*b*log(x) + B*c*x**2/2 + C*c*x**3/3 + x*(A*c + C*b) - (A*b + C*a)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_8():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**5
    F = -A*a/(4*x**4) - B*a/(3*x**3) - B*b/x + B*c*x + C*c*x**2/2 + (A*c + C*b)*log(x) - (A*b + C*a)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_9():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**6
    F = -A*a/(5*x**5) - B*a/(4*x**4) - B*b/(2*x**2) + B*c*log(x) + C*c*x - (A*c + C*b)/x - (A*b + C*a)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_10():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)/x**7
    F = -A*a/(6*x**6) - B*a/(5*x**5) - B*b/(3*x**3) - B*c/x + C*c*log(x) - (A*c + C*b)/(2*x**2) - (A*b + C*a)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_11():
    f = x**2*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2
    F = A*a**2*x**3/3 + B*a**2*x**4/4 + B*a*b*x**6/3 + B*b*c*x**10/5 + B*c**2*x**12/12 + B*x**8*(2*a*c + b**2)/8 + C*c**2*x**13/13 + a*x**5*(2*A*b + C*a)/5 + c*x**11*(A*c + 2*C*b)/11 + x**9*(2*A*b*c/9 + C*(2*a*c + b**2)/9) + x**7*(A*(2*a*c + b**2)/7 + 2*C*a*b/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_12():
    f = x*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2
    F = A*a**2*x**2/2 + B*a**2*x**3/3 + 2*B*a*b*x**5/5 + 2*B*b*c*x**9/9 + B*c**2*x**11/11 + B*x**7*(2*a*c + b**2)/7 + C*c**2*x**12/12 + a*x**4*(2*A*b + C*a)/4 + c*x**10*(A*c + 2*C*b)/10 + x**8*(A*b*c/4 + C*(2*a*c + b**2)/8) + x**6*(A*(2*a*c + b**2)/6 + C*a*b/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_13():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2
    F = A*a**2*x + B*a**2*x**2/2 + B*a*b*x**4/2 + B*b*c*x**8/4 + B*c**2*x**10/10 + B*x**6*(2*a*c + b**2)/6 + C*c**2*x**11/11 + a*x**3*(2*A*b + C*a)/3 + c*x**9*(A*c + 2*C*b)/9 + x**7*(2*A*b*c/7 + C*(2*a*c + b**2)/7) + x**5*(A*(2*a*c + b**2)/5 + 2*C*a*b/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_14():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x
    F = A*a**2*log(x) + B*a**2*x + 2*B*a*b*x**3/3 + 2*B*b*c*x**7/7 + B*c**2*x**9/9 + B*x**5*(2*a*c + b**2)/5 + C*c**2*x**10/10 + a*x**2*(2*A*b + C*a)/2 + c*x**8*(A*c + 2*C*b)/8 + x**6*(A*b*c/3 + C*(2*a*c + b**2)/6) + x**4*(A*(2*a*c + b**2)/4 + C*a*b/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_15():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**2
    F = -A*a**2/x + B*a**2*log(x) + B*a*b*x**2 + B*b*c*x**6/3 + B*c**2*x**8/8 + B*x**4*(2*a*c + b**2)/4 + C*c**2*x**9/9 + a*x*(2*A*b + C*a) + c*x**7*(A*c + 2*C*b)/7 + x**5*(2*A*b*c/5 + C*(2*a*c + b**2)/5) + x**3*(A*(2*a*c + b**2)/3 + 2*C*a*b/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_16():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**3
    F = -A*a**2/(2*x**2) - B*a**2/x + 2*B*a*b*x + 2*B*b*c*x**5/5 + B*c**2*x**7/7 + B*x**3*(2*a*c + b**2)/3 + C*c**2*x**8/8 + a*(2*A*b + C*a)*log(x) + c*x**6*(A*c + 2*C*b)/6 + x**4*(A*b*c/2 + C*(2*a*c + b**2)/4) + x**2*(A*(2*a*c + b**2)/2 + C*a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_17():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**4
    F = -A*a**2/(3*x**3) - B*a**2/(2*x**2) + 2*B*a*b*log(x) + B*b*c*x**4/2 + B*c**2*x**6/6 + B*x**2*(2*a*c + b**2)/2 + C*c**2*x**7/7 - a*(2*A*b + C*a)/x + c*x**5*(A*c + 2*C*b)/5 + x**3*(2*A*b*c/3 + C*(2*a*c + b**2)/3) + x*(A*(2*a*c + b**2) + 2*C*a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_18():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**5
    F = -A*a**2/(4*x**4) - B*a**2/(3*x**3) - 2*B*a*b/x + 2*B*b*c*x**3/3 + B*c**2*x**5/5 + B*x*(2*a*c + b**2) + C*c**2*x**6/6 - a*(2*A*b + C*a)/(2*x**2) + c*x**4*(A*c + 2*C*b)/4 + x**2*(A*b*c + C*(2*a*c + b**2)/2) + (A*(2*a*c + b**2) + 2*C*a*b)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_19():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**6
    F = -A*a**2/(5*x**5) - B*a**2/(4*x**4) - B*a*b/x**2 + B*b*c*x**2 + B*c**2*x**4/4 + B*(2*a*c + b**2)*log(x) + C*c**2*x**5/5 - a*(2*A*b + C*a)/(3*x**3) + c*x**3*(A*c + 2*C*b)/3 + x*(2*A*b*c + C*(2*a*c + b**2)) - (A*(2*a*c + b**2) + 2*C*a*b)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_20():
    f = (A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2/x**7
    F = -A*a**2/(6*x**6) - B*a**2/(5*x**5) - 2*B*a*b/(3*x**3) + 2*B*b*c*x + B*c**2*x**3/3 - B*(2*a*c + b**2)/x + C*c**2*x**4/4 - a*(2*A*b + C*a)/(4*x**4) + c*x**2*(A*c + 2*C*b)/2 + (2*A*b*c + C*(2*a*c + b**2))*log(x) - (A*(2*a*c + b**2) + 2*C*a*b)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_21():
    f = x**4*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = -B*b*log(a + b*x**2 + c*x**4)/(4*c**2) + B*x**2/(2*c) - B*(-2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2)) + C*x**3/(3*c) + x*(A*c - C*b)/c**2 - sqrt(2)*(A*b*c + C*a*c - C*b**2 + (A*c*(-2*a*c + b**2) - C*b*(-3*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*(A*b*c + C*a*c - C*b**2 - (A*c*(-2*a*c + b**2) - C*b*(-3*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_22():
    f = x**3*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = B*x/c - sqrt(2)*B*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*B*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + C*x**2/(2*c) + (A*c - C*b)*log(a + b*x**2 + c*x**4)/(4*c**2) + (A*b*c + 2*C*a*c - C*b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_23():
    f = x**2*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + B*log(a + b*x**2 + c*x**4)/(4*c) + C*x/c + sqrt(2)*(A*c - C*b + (A*b*c + 2*C*a*c - C*b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(A*c - C*b - (A*b*c - C*(-2*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_24():
    f = x*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*B*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*B*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + C*log(a + b*x**2 + c*x**4)/(4*c) - (2*A*c - C*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_25():
    f = (A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = -B*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2) + sqrt(2)*(C - (2*A*c - C*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(C + (2*A*c - C*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_26():
    f = (A + B*x + C*x**2)/(x*(a + b*x**2 + c*x**4))
    F = A*log(x)/a - A*log(a + b*x**2 + c*x**4)/(4*a) - sqrt(2)*B*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*B*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + (A*b - 2*C*a)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_27():
    f = (A + B*x + C*x**2)/(x**2*(a + b*x**2 + c*x**4))
    F = -A/(a*x) + B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2)) + B*log(x)/a - B*log(a + b*x**2 + c*x**4)/(4*a) - sqrt(2)*sqrt(c)*(A - (A*b - 2*C*a)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(A + (A*b - 2*C*a)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_28():
    f = (A + B*x + C*x**2)/(x**3*(a + b*x**2 + c*x**4))
    F = -A/(2*a*x**2) - sqrt(2)*B*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*B*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2))) - B/(a*x) - (A*b - C*a)*log(x)/a**2 + (A*b - C*a)*log(a + b*x**2 + c*x**4)/(4*a**2) - (A*(-2*a*c + b**2) - C*a*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_29():
    f = x**4*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = 2*B*a*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*x**2*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x**3*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + x*(2*A*c - C*b)/(2*c*(-4*a*c + b**2)) + sqrt(2)*(A*b*c + C*(-6*a*c + b**2) + (A*c*(4*a*c + b**2) + C*b*(-8*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(A*b*c + C*(-6*a*c + b**2) - (A*c*(4*a*c + b**2) + C*b*(-8*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_30():
    f = x**3*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = B*x*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*B*(b - (4*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*B*(4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (A*b - 2*C*a)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (a*(2*A*c - C*b) + x**2*(A*b*c + 2*C*a*c - C*b**2))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_31():
    f = x**2*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_32():
    f = x*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = -sqrt(2)*B*sqrt(c)*(2*b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*B*sqrt(c)*(2*b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - B*x*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + (2*A*c - C*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_33():
    f = (A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = 2*B*c*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - B*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(A*b - 2*C*a - (-12*A*a*c + A*b**2 + 4*C*a*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(A*b - 2*C*a + (A*(-12*a*c + b**2) + 4*C*a*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-2*A*a*c + A*b**2 - C*a*b + c*x**2*(A*b - 2*C*a))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_34():
    f = (A + B*x + C*x**2)/(x*(a + b*x**2 + c*x**4)**2)
    F = A*log(x)/a**2 - A*log(a + b*x**2 + c*x**4)/(4*a**2) - sqrt(2)*B*sqrt(c)*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*B*sqrt(c)*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + B*x*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (A*(-2*a*c + b**2) - C*a*b + c*x**2*(A*b - 2*C*a))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (A*(-6*a*b*c + b**3) + 4*C*a**2*c)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_35():
    f = (A + B*x + C*x**2)/(x**2*(a + b*x**2 + c*x**4)**2)
    F = B*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + B*b*(-6*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + B*log(x)/a**2 - B*log(a + b*x**2 + c*x**4)/(4*a**2) + (A*(-2*a*c + b**2) - C*a*b + c*x**2*(A*b - 2*C*a))/(2*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*sqrt(c)*(-10*A*a*c + 3*A*b**2 - C*a*b - (A*(-16*a*b*c + 3*b**3) - C*a*(-12*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(A*(-16*a*b*c - 10*a*c*sqrt(-4*a*c + b**2) + 3*b**3 + 3*b**2*sqrt(-4*a*c + b**2)) - C*a*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*A*a*c + 3*A*b**2 - C*a*b)/(2*a**2*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_36():
    f = (A + B*x + C*x**2)/(x**3*(a + b*x**2 + c*x**4)**2)
    F = B*(-2*a*c + b**2 + b*c*x**2)/(2*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*B*sqrt(c)*(-16*a*b*c + 3*b**3 - (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*B*sqrt(c)*(-16*a*b*c + 3*b**3 + (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - B*(-10*a*c + 3*b**2)/(2*a**2*x*(-4*a*c + b**2)) + (A*(-2*a*c + b**2) - C*a*b + c*x**2*(A*b - 2*C*a))/(2*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-6*A*a*c + 2*A*b**2 - C*a*b)/(2*a**2*x**2*(-4*a*c + b**2)) - (2*A*b - C*a)*log(x)/a**3 + (2*A*b - C*a)*log(a + b*x**2 + c*x**4)/(4*a**3) - (2*A*(6*a**2*c**2 - 6*a*b**2*c + b**4) - C*a*b*(-6*a*c + b**2))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_37():
    f = (d*x)**m*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**3
    F = A*a**3*(d*x)**(m + 1)/(d*(m + 1)) + B*a**3*(d*x)**(m + 2)/(d**2*(m + 2)) + 3*B*a**2*b*(d*x)**(m + 4)/(d**4*(m + 4)) + 3*B*a*(d*x)**(m + 6)*(a*c + b**2)/(d**6*(m + 6)) + 3*B*b*c**2*(d*x)**(m + 12)/(d**12*(m + 12)) + B*b*(d*x)**(m + 8)*(6*a*c + b**2)/(d**8*(m + 8)) + B*c**3*(d*x)**(m + 14)/(d**14*(m + 14)) + 3*B*c*(d*x)**(m + 10)*(a*c + b**2)/(d**10*(m + 10)) + C*c**3*(d*x)**(m + 15)/(d**15*(m + 15)) + a**2*(d*x)**(m + 3)*(3*A*b + C*a)/(d**3*(m + 3)) + 3*a*(d*x)**(m + 5)*(A*(a*c + b**2) + C*a*b)/(d**5*(m + 5)) + c**2*(d*x)**(m + 13)*(A*c + 3*C*b)/(d**13*(m + 13)) + 3*c*(d*x)**(m + 11)*(A*b*c + C*(a*c + b**2))/(d**11*(m + 11)) + (d*x)**(m + 7)*(A*(6*a*b*c + b**3) + 3*C*a*(a*c + b**2))/(d**7*(m + 7)) + (d*x)**(m + 9)*(3*A*c*(a*c + b**2) + C*b*(6*a*c + b**2))/(d**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_38():
    f = (d*x)**m*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)**2
    F = A*a**2*(d*x)**(m + 1)/(d*(m + 1)) + B*a**2*(d*x)**(m + 2)/(d**2*(m + 2)) + 2*B*a*b*(d*x)**(m + 4)/(d**4*(m + 4)) + 2*B*b*c*(d*x)**(m + 8)/(d**8*(m + 8)) + B*c**2*(d*x)**(m + 10)/(d**10*(m + 10)) + B*(d*x)**(m + 6)*(2*a*c + b**2)/(d**6*(m + 6)) + C*c**2*(d*x)**(m + 11)/(d**11*(m + 11)) + a*(d*x)**(m + 3)*(2*A*b + C*a)/(d**3*(m + 3)) + c*(d*x)**(m + 9)*(A*c + 2*C*b)/(d**9*(m + 9)) + (d*x)**(m + 5)*(A*(2*a*c + b**2) + 2*C*a*b)/(d**5*(m + 5)) + (d*x)**(m + 7)*(2*A*b*c + C*(2*a*c + b**2))/(d**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_39():
    f = (d*x)**m*(A + B*x + C*x**2)*(a + b*x**2 + c*x**4)
    F = A*a*(d*x)**(m + 1)/(d*(m + 1)) + B*a*(d*x)**(m + 2)/(d**2*(m + 2)) + B*b*(d*x)**(m + 4)/(d**4*(m + 4)) + B*c*(d*x)**(m + 6)/(d**6*(m + 6)) + C*c*(d*x)**(m + 7)/(d**7*(m + 7)) + (d*x)**(m + 3)*(A*b + C*a)/(d**3*(m + 3)) + (d*x)**(m + 5)*(A*c + C*b)/(d**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_40():
    f = (d*x)**m*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)
    F = -2*B*c*(d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d**2*(b + sqrt(-4*a*c + b**2))*(m + 2)*sqrt(-4*a*c + b**2)) + 2*B*c*(d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(d**2*(b - sqrt(-4*a*c + b**2))*(m + 2)*sqrt(-4*a*c + b**2)) + (d*x)**(m + 1)*(C - (2*A*c - C*b)/sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(d*(b + sqrt(-4*a*c + b**2))*(m + 1)) + (d*x)**(m + 1)*(C + (2*A*c - C*b)/sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(d*(b - sqrt(-4*a*c + b**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_41():
    f = (d*x)**m*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = B*c*(d*x)**(m + 2)*(4*a*c*(2 - m) + b*m*(b - sqrt(-4*a*c + b**2)))*hyper((1, m/2 + 1), (m/2 + 2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(2*a*d**2*(b + sqrt(-4*a*c + b**2))*(m + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) - B*c*(d*x)**(m + 2)*(4*a*c*(2 - m) + b*m*(b + sqrt(-4*a*c + b**2)))*hyper((1, m/2 + 1), (m/2 + 2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(2*a*d**2*(b - sqrt(-4*a*c + b**2))*(m + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) + B*(d*x)**(m + 2)*(-2*a*c + b**2 + b*c*x**2)/(2*a*d**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - c*(d*x)**(m + 1)*(A*(-4*a*c*(3 - m) + b**2*(1 - m) - b*(1 - m)*sqrt(-4*a*c + b**2)) + 2*C*a*(2*b + (1 - m)*sqrt(-4*a*c + b**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(2*a*d*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + c*(d*x)**(m + 1)*(A*(-4*a*c*(3 - m) + b**2*(1 - m) + b*(1 - m)*sqrt(-4*a*c + b**2)) + 2*C*a*(2*b - (1 - m)*sqrt(-4*a*c + b**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(2*a*d*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + (d*x)**(m + 1)*(A*(-2*a*c + b**2) - C*a*b + c*x**2*(A*b - 2*C*a))/(2*a*d*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_42():
    f = x**2*(A + B*x + C*x**2)/(a + b*x**2 + c*x**4)**2
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_43():
    f = x*(A*x + B*x**2 + C*x**3)/(a + b*x**2 + c*x**4)**2
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_44():
    f = (A*x**2 + B*x**3 + C*x**4)/(a + b*x**2 + c*x**4)**2
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_45():
    f = (A*x**3 + B*x**4 + C*x**5)/(x*(a + b*x**2 + c*x**4)**2)
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_46():
    f = (A*x**4 + B*x**5 + C*x**6)/(x**2*(a + b*x**2 + c*x**4)**2)
    F = -B*b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + B*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(A*b - 2*C*a + x**2*(2*A*c - C*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*A*c - C*b + (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*A*c - C*b - (4*A*b*c - C*(4*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_47():
    f = x**7*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**8/(8*c) + x**6*(-b*f + c*e)/(6*c**2) + x**4*(b**2*f + c**2*d - c*(a*f + b*e))/(4*c**3) + x**2*(-a*c**2*e - b**3*f + b**2*c*e - b*c*(-2*a*f + c*d))/(2*c**4) - (-2*a*b*c**2*e + a*c**2*(-a*f + c*d) - b**4*f + b**3*c*e - b**2*c*(-3*a*f + c*d))*log(a + b*x**2 + c*x**4)/(4*c**5) - (2*a**2*c**3*e - 4*a*b**2*c**2*e + a*b*c**2*(-5*a*f + 3*c*d) - b**5*f + b**4*c*e - b**3*c*(-5*a*f + c*d))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**5*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_48():
    f = x**5*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**6/(6*c) + x**4*(-b*f + c*e)/(4*c**2) + x**2*(b**2*f + c**2*d - c*(a*f + b*e))/(2*c**3) + (-a*c**2*e - b**3*f + b**2*c*e - b*c*(-2*a*f + c*d))*log(a + b*x**2 + c*x**4)/(4*c**4) + (-3*a*b*c**2*e + 2*a*c**2*(-a*f + c*d) - b**4*f + b**3*c*e - b**2*c*(-4*a*f + c*d))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_49():
    f = x**3*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**4/(4*c) + x**2*(-b*f + c*e)/(2*c**2) + (b**2*f + c**2*d - c*(a*f + b*e))*log(a + b*x**2 + c*x**4)/(4*c**3) - (-2*a*c**2*e - b**3*f + b**2*c*e - b*c*(-3*a*f + c*d))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_50():
    f = x*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**2/(2*c) + (-b*f + c*e)*log(a + b*x**2 + c*x**4)/(4*c**2) - (-2*a*c*f + b**2*f - b*c*e + 2*c**2*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_51():
    f = (d + e*x**2 + f*x**4)/(x*(a + b*x**2 + c*x**4))
    F = d*log(x)/a - (-a*f + c*d)*log(a + b*x**2 + c*x**4)/(4*a*c) + (a*b*f - 2*a*c*e + b*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_52():
    f = (d + e*x**2 + f*x**4)/(x**3*(a + b*x**2 + c*x**4))
    F = -d/(2*a*x**2) - (-a*e + b*d)*log(x)/a**2 + (-a*e + b*d)*log(a + b*x**2 + c*x**4)/(4*a**2) - (-a*b*e - 2*a*(-a*f + c*d) + b**2*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_53():
    f = (d + e*x**2 + f*x**4)/(x**5*(a + b*x**2 + c*x**4))
    F = -d/(4*a*x**4) + (-a*e + b*d)/(2*a**2*x**2) + (-a*b*e - a*(-a*f + c*d) + b**2*d)*log(x)/a**3 - (-a*b*e - a*(-a*f + c*d) + b**2*d)*log(a + b*x**2 + c*x**4)/(4*a**3) + (2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_54():
    f = (d + e*x**2 + f*x**4)/(x**7*(a + b*x**2 + c*x**4))
    F = -d/(6*a*x**6) + (-a*e + b*d)/(4*a**2*x**4) - (-a*b*e - a*(-a*f + c*d) + b**2*d)/(2*a**3*x**2) - (a**2*c*e - a*b**2*e - a*b*(-a*f + 2*c*d) + b**3*d)*log(x)/a**4 + (a**2*c*e - a*b**2*e - a*b*(-a*f + 2*c*d) + b**3*d)*log(a + b*x**2 + c*x**4)/(4*a**4) - (3*a**2*b*c*e + 2*a**2*c*(-a*f + c*d) - a*b**3*e - a*b**2*(-a*f + 4*c*d) + b**4*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_55():
    f = x**4*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**5/(5*c) + x**3*(-b*f + c*e)/(3*c**2) + x*(b**2*f + c**2*d - c*(a*f + b*e))/c**3 + sqrt(2)*(-a*c**2*e - b**3*f + b**2*c*e - b*c*(-2*a*f + c*d) + (-3*a*b*c**2*e + 2*a*c**2*(-a*f + c*d) - b**4*f + b**3*c*e - b**2*c*(-4*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-a*c**2*e - b**3*f + b**2*c*e - b*c*(-2*a*f + c*d) - (-3*a*b*c**2*e + 2*a*c**2*(-a*f + c*d) - b**4*f + b**3*c*e - b**2*c*(-4*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_56():
    f = x**2*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x**3/(3*c) + x*(-b*f + c*e)/c**2 + sqrt(2)*(-a*c*f + b**2*f - b*c*e + c**2*d - (-2*a*c**2*e - b**3*f + b**2*c*e - b*c*(-3*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-a*c*f + b**2*f - b*c*e + c**2*d + (-2*a*c**2*e - b**3*f + b**2*c*e - b*c*(-3*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_57():
    f = (d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)
    F = f*x/c + sqrt(2)*(-b*f + c*e - (-2*a*c*f + b**2*f - b*c*e + 2*c**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-b*f + c*e + (b**2*f + 2*c**2*d - c*(2*a*f + b*e))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_58():
    f = (d + e*x**2 + f*x**4)/(x**2*(a + b*x**2 + c*x**4))
    F = -d/(a*x) - sqrt(2)*(-a*f + c*d - (a*b*f - 2*a*c*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*(-a*f + c*d + (a*b*f - 2*a*c*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_59():
    f = (d + e*x**2 + f*x**4)/(x**4*(a + b*x**2 + c*x**4))
    F = -d/(3*a*x**3) - sqrt(2)*sqrt(c)*(-a*(-2*a*f + 2*c*d - e*sqrt(-4*a*c + b**2)) + b**2*d - b*(a*e + d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-a*e + b*d + (-a*b*e - 2*a*(-a*f + c*d) + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b - sqrt(-4*a*c + b**2))) + (-a*e + b*d)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_60():
    f = (d + e*x**2 + f*x**4)/(x**6*(a + b*x**2 + c*x**4))
    F = -d/(5*a*x**5) + (-a*e + b*d)/(3*a**2*x**3) - sqrt(2)*sqrt(c)*(-a*b*e - a*(-a*f + c*d) + b**2*d - (2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**3*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(-a*b*e - a*(-a*f + c*d) + b**2*d + (2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**3*sqrt(b - sqrt(-4*a*c + b**2))) - (-a*b*e - a*(-a*f + c*d) + b**2*d)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_61():
    f = x**7*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = x**6*(2*a*c*e - b*(a*f + c*d) - x**2*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + x**4*(3*b**2*f + 4*c**2*d - 2*c*(4*a*f + b*e))/(4*c**2*(-4*a*c + b**2)) + x**2*(-6*a*c**2*e - 3*b**3*f + 2*b**2*c*e - b*c*(-11*a*f + c*d))/(2*c**3*(-4*a*c + b**2)) + (3*b**2*f + c**2*d - 2*c*(a*f + b*e))*log(a + b*x**2 + c*x**4)/(4*c**4) - (12*a**2*c**3*e - 12*a*b**2*c**2*e + 6*a*b*c**2*(-5*a*f + c*d) - 3*b**5*f + 2*b**4*c*e - b**3*c*(-20*a*f + c*d))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**4*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_62():
    f = x**5*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = x**4*(2*a*c*e - b*(a*f + c*d) - x**2*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + x**2*(2*b**2*f + 2*c**2*d - c*(6*a*f + b*e))/(2*c**2*(-4*a*c + b**2)) + (-2*b*f + c*e)*log(a + b*x**2 + c*x**4)/(4*c**3) - (12*a**2*c**2*f - 2*a*c*(6*b**2*f - 3*b*c*e + 2*c**2*d) - b**3*(-2*b*f + c*e))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_63():
    f = x**3*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = x**2*(2*a*c*e - b*(a*f + c*d) - x**2*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + f*log(a + b*x**2 + c*x**4)/(4*c**2) + (4*a*c**2*e + b**3*f - 2*b*c*(3*a*f + c*d))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_64():
    f = x*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = (2*a*f - b*e + 2*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a*c*e - b*(a*f + c*d) - x**2*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_65():
    f = (d + e*x**2 + f*x**4)/(x*(a + b*x**2 + c*x**4)**2)
    F = (-a*b*e - 2*a*(-a*f + c*d) + b**2*d + x**2*(a*b*f - 2*a*c*e + b*c*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + d*log(x)/a**2 - d*log(a + b*x**2 + c*x**4)/(4*a**2) + (4*a**2*c*e - 2*a*b*(a*f + 3*c*d) + b**3*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_66():
    f = (d + e*x**2 + f*x**4)/(x**3*(a + b*x**2 + c*x**4)**2)
    F = -d/(2*a**2*x**2) - (2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d + c*x**2*(-a*b*e - 2*a*(-a*f + c*d) + b**2*d))/(2*a**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-a*e + 2*b*d)*log(x)/a**3 + (-a*e + 2*b*d)*log(a + b*x**2 + c*x**4)/(4*a**3) - (6*a**2*b*c*e + 4*a**2*c*(-a*f + 3*c*d) - a*b**3*e - 12*a*b**2*c*d + 2*b**4*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_67():
    f = (d + e*x**2 + f*x**4)/(x**5*(a + b*x**2 + c*x**4)**2)
    F = -d/(4*a**2*x**4) + (3*a**2*b*c*e + 2*a**2*c*(-a*f + c*d) - a*b**3*e - a*b**2*(-a*f + 4*c*d) + b**4*d + c*x**2*(2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d))/(2*a**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (-a*e + 2*b*d)/(2*a**3*x**2) + (-2*a*b*e - a*(-a*f + 2*c*d) + 3*b**2*d)*log(x)/a**4 - (-2*a*b*e - a*(-a*f + 2*c*d) + 3*b**2*d)*log(a + b*x**2 + c*x**4)/(4*a**4) + (-12*a**3*c**2*e + 12*a**2*b**2*c*e + 6*a**2*b*c*(-a*f + 5*c*d) - 2*a*b**4*e - a*b**3*(-a*f + 20*c*d) + 3*b**5*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**4*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_68():
    f = x**6*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = f*x**3/(3*c**2) + x*(-2*b*f + c*e)/c**3 + x*(a*(-2*a*c**2*e - b**3*f + b**2*c*e - b*c*(-3*a*f + c*d)) + x**2*(-3*a*b*c**2*e + 2*a*c**2*(-a*f + c*d) - b**4*f + b**3*c*e - b**2*c*(-4*a*f + c*d)))/(2*c**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(-13*a*b*c**2*e + 2*a*c**2*(-7*a*f + 3*c*d) - 5*b**4*f + 3*b**3*c*e - b**2*c*(-24*a*f + c*d) + (20*a**2*c**3*e - 19*a*b**2*c**2*e + 4*a*b*c**2*(-13*a*f + 2*c*d) - 5*b**5*f + 3*b**4*c*e - b**3*c*(-34*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(7)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(-13*a*b*c**2*e + 2*a*c**2*(-7*a*f + 3*c*d) - 5*b**4*f + 3*b**3*c*e - b**2*c*(-24*a*f + c*d) - (20*a**2*c**3*e - 19*a*b**2*c**2*e + 4*a*b*c**2*(-13*a*f + 2*c*d) - 5*b**5*f + 3*b**4*c*e - b**3*c*(-34*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(7)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_69():
    f = x**4*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = f*x/c**2 + x*(a*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d) - x**2*(-2*a*c**2*e - b**3*f + b**2*c*e - b*c*(-3*a*f + c*d)))/(2*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(-6*a*c**2*e - 3*b**3*f + b**2*c*e + b*c*(13*a*f + c*d) + (-8*a*b*c**2*e + 4*a*c**2*(-5*a*f + c*d) - 3*b**4*f + b**3*c*e + b**2*c*(19*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(-6*a*c**2*e - 3*b**3*f + b**2*c*e + b*c*(13*a*f + c*d) - (-8*a*b*c**2*e + 4*a*c**2*(-5*a*f + c*d) - 3*b**4*f + b**3*c*e + b**2*c*(19*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_70():
    f = x**2*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = -x*(a*b*f - 2*a*c*e + b*c*d + x**2*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(6*a*f - b**2*f/c - b*e + 2*c*d - (4*a*c**2*e + b**3*f + b**2*c*e - 4*b*c*(2*a*f + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(6*a*f - b**2*f/c - b*e + 2*c*d + (4*a*c**2*e + b**3*f + b**2*c*e - 4*b*c*(2*a*f + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_71():
    f = (d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2
    F = x*(-a*b*e - 2*a*(-a*f + c*d) + b**2*d + x**2*(a*b*f - 2*a*c*e + b*c*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b*f - 2*a*c*e + b*c*d - (4*a*b*c*e - 4*a*c*(a*f + 3*c*d) + b**2*(-a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b*f - 2*a*c*e + b*c*d + (4*a*b*c*e - 4*a*c*(a*f + 3*c*d) + b**2*(-a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_72():
    f = (d + e*x**2 + f*x**4)/(x**2*(a + b*x**2 + c*x**4)**2)
    F = -sqrt(2)*sqrt(c)*(-a*b*e - 2*a*(-a*f + 5*c*d) + 3*b**2*d - (12*a**2*c*e - a*b**2*e - 4*a*b*(a*f + 4*c*d) + 3*b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(-a*b*e - 2*a*(-a*f + 5*c*d) + 3*b**2*d + (12*a**2*c*e - a*b**2*e - 4*a*b*(a*f + 4*c*d) + 3*b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - d/(a**2*x) - x*(a*(a*(b*f + 2*c*e) - b*(b*e + 3*c*d) + b**3*d/a) + c*x**2*(-a*b*e - 2*a*(-a*f + c*d) + b**2*d))/(2*a**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_73():
    f = (d + e*x**2 + f*x**4)/(x**4*(a + b*x**2 + c*x**4)**2)
    F = -d/(3*a**2*x**3) - sqrt(2)*sqrt(c)*(2*a**2*c*(-6*a*f + 14*c*d - 5*e*sqrt(-4*a*c + b**2)) - a*b**2*(-a*f + 29*c*d - 3*e*sqrt(-4*a*c + b**2)) + a*b*(16*a*c*e - a*f*sqrt(-4*a*c + b**2) + 19*c*d*sqrt(-4*a*c + b**2)) + 5*b**4*d - b**3*(3*a*e + 5*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(2*a**2*c*(-6*a*f + 14*c*d + 5*e*sqrt(-4*a*c + b**2)) - a*b**2*(-a*f + 29*c*d + 3*e*sqrt(-4*a*c + b**2)) - a*b*(-16*a*c*e - a*f*sqrt(-4*a*c + b**2) + 19*c*d*sqrt(-4*a*c + b**2)) + 5*b**4*d + b**3*(-3*a*e + 5*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + x*(a**2*(-2*a*c*f + b**2*f + 3*b*c*e + 2*c**2*d - b**2*(b*e + 4*c*d)/a + b**4*d/a**2) + c*x**2*(2*a**2*c*e - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d))/(2*a**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (-a*e + 2*b*d)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_74():
    f = x**9*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**8/8 - 9*x**6/2 + 49*x**4/2 - 293*x**2/2 + (415*x**2 + 414)/(2*x**4 + 6*x**2 + 4) + 2*log(x**2 + 1) + 392*log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_75():
    f = x**7*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**6/6 - 27*x**4/4 + 49*x**2 - (207*x**2 + 206)/(2*x**4 + 6*x**2 + 4) - 5*log(x**2 + 1)/2 - 144*log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_76():
    f = x**5*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**4/4 - 27*x**2/2 + (103*x**2 + 102)/(2*x**4 + 6*x**2 + 4) + 3*log(x**2 + 1) + 46*log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_77():
    f = x**3*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**2/2 - (51*x**2 + 50)/(2*x**4 + 6*x**2 + 4) - 7*log(x**2 + 1)/2 - 10*log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_78():
    f = x*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = (25*x**2 + 24)/(2*x**4 + 6*x**2 + 4) + 4*log(x**2 + 1) - 3*log(x**2 + 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_79():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x*(x**4 + 3*x**2 + 2)**2)
    F = -(12*x**2 + 11)/(2*x**4 + 6*x**2 + 4) + log(x) - 9*log(x**2 + 1)/2 + 4*log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_80():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**3*(x**4 + 3*x**2 + 2)**2)
    F = (11*x**2 + 9)/(4*x**4 + 12*x**2 + 8) - 11*log(x)/4 + 5*log(x**2 + 1) - 29*log(x**2 + 2)/8 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_81():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**5*(x**4 + 3*x**2 + 2)**2)
    F = -(9*x**2 + 5)/(8*x**4 + 24*x**2 + 16) + 23*log(x)/4 - 11*log(x**2 + 1)/2 + 21*log(x**2 + 2)/8 + 11/(8*x**2) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_82():
    f = x**8*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**7/7 - 27*x**5/5 + 98*x**3/3 - x*(207*x**2 + 206)/(2*x**4 + 6*x**2 + 4) - 293*x + 9*atan(x)/2 + 340*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_83():
    f = x**6*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = x**5 - 9*x**3 + x*(103*x**2 + 102)/(2*x**4 + 6*x**2 + 4) + 98*x - 11*atan(x)/2 - 118*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_84():
    f = x**4*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = 5*x**3/3 - x*(51*x**2 + 50)/(2*x**4 + 6*x**2 + 4) - 27*x + 13*atan(x)/2 + 33*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_85():
    f = x**2*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = x*(25*x**2 + 24)/(2*x**4 + 6*x**2 + 4) + 5*x - 15*atan(x)/2 - 7*sqrt(2)*atan(sqrt(2)*x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_86():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**2
    F = -x*(12*x**2 + 11)/(2*x**4 + 6*x**2 + 4) + 17*atan(x)/2 - 19*sqrt(2)*atan(sqrt(2)*x/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_87():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**2*(x**4 + 3*x**2 + 2)**2)
    F = x*(11*x**2 + 9)/(4*x**4 + 12*x**2 + 8) - 19*atan(x)/2 + 45*sqrt(2)*atan(sqrt(2)*x/2)/8 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_88():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4*(x**4 + 3*x**2 + 2)**2)
    F = -x*(9*x**2 + 5)/(8*x**4 + 24*x**2 + 16) + 21*atan(x)/2 - 71*sqrt(2)*atan(sqrt(2)*x/2)/16 + 11/(4*x) - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_89():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**6*(x**4 + 3*x**2 + 2)**2)
    F = -x*(3 - 5*x**2)/(16*x**4 + 48*x**2 + 32) - 23*atan(x)/2 + 97*sqrt(2)*atan(sqrt(2)*x/2)/32 - 23/(4*x) + 11/(12*x**3) - 1/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_90():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**8*(x**4 + 3*x**2 + 2)**2)
    F = x*(3*x**2 + 19)/(32*x**4 + 96*x**2 + 64) + 25*atan(x)/2 - 123*sqrt(2)*atan(sqrt(2)*x/2)/64 + 137/(16*x) - 23/(12*x**3) + 11/(20*x**5) - 1/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_91():
    f = x**10*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = x**5 - 14*x**3 + x*(415*x**2 + 414)/(4*(x**4 + 3*x**2 + 2)**2) + x*(1669*x**2 + 824)/(8*x**4 + 24*x**2 + 16) + 214*x + 477*atan(x)/8 - 351*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_92():
    f = x**8*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = 5*x**3/3 + x*(24 - 409*x**2)/(8*x**4 + 24*x**2 + 16) - x*(207*x**2 + 206)/(4*(x**4 + 3*x**2 + 2)**2) - 42*x - 449*atan(x)/8 + 219*sqrt(2)*atan(sqrt(2)*x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_93():
    f = x**6*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = -x*(15*x**2 + 244)/(8*x**4 + 24*x**2 + 16) + x*(103*x**2 + 102)/(4*(x**4 + 3*x**2 + 2)**2) + 5*x + 413*atan(x)/8 - 191*sqrt(2)*atan(sqrt(2)*x/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_94():
    f = x**4*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = -x*(51*x**2 + 50)/(4*(x**4 + 3*x**2 + 2)**2) + x*(125*x**2 + 254)/(8*x**4 + 24*x**2 + 16) - 369*atan(x)/8 + 267*sqrt(2)*atan(sqrt(2)*x/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_95():
    f = x**2*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = x*(25*x**2 + 24)/(4*(x**4 + 3*x**2 + 2)**2) - x*(130*x**2 + 211)/(8*x**4 + 24*x**2 + 16) + 317*atan(x)/8 - 447*sqrt(2)*atan(sqrt(2)*x/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_96():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 3*x**2 + 2)**3
    F = -x*(12*x**2 + 11)/(4*(x**4 + 3*x**2 + 2)**2) + x*(217*x**2 + 335)/(16*x**4 + 48*x**2 + 32) - 257*atan(x)/8 + 731*sqrt(2)*atan(sqrt(2)*x/2)/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_97():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**2*(x**4 + 3*x**2 + 2)**3)
    F = x*(11*x**2 + 9)/(8*(x**4 + 3*x**2 + 2)**2) - x*(347*x**2 + 547)/(32*x**4 + 96*x**2 + 64) + 189*atan(x)/8 - 1119*sqrt(2)*atan(sqrt(2)*x/2)/64 - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_98():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4*(x**4 + 3*x**2 + 2)**3)
    F = -x*(9*x**2 + 5)/(16*(x**4 + 3*x**2 + 2)**2) + x*(571*x**2 + 951)/(64*x**4 + 192*x**2 + 128) - 113*atan(x)/8 + 1611*sqrt(2)*atan(sqrt(2)*x/2)/128 + 17/(8*x) - 1/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_99():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**6*(x**4 + 3*x**2 + 2)**3)
    F = -x*(3 - 5*x**2)/(32*(x**4 + 3*x**2 + 2)**2) - x*(999*x**2 + 1771)/(128*x**4 + 384*x**2 + 256) + 29*atan(x)/8 - 2207*sqrt(2)*atan(sqrt(2)*x/2)/256 - 93/(16*x) + 17/(24*x**3) - 1/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_100():
    f = x**9*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**8/8 - 17*x**6/6 + 19*x**4/4 + 19*x**2 - (175*x**2 + 375)/(8*x**4 + 16*x**2 + 24) - 183*log(x**4 + 2*x**2 + 3)/4 + 201*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_101():
    f = x**7*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**6/6 - 17*x**4/4 + 19*x**2/2 + (125*x**2 + 75)/(8*x**4 + 16*x**2 + 24) + 19*log(x**4 + 2*x**2 + 3)/2 - 455*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_102():
    f = x**5*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**4/4 - 17*x**2/2 + (75 - 25*x**2)/(8*x**4 + 16*x**2 + 24) + 19*log(x**4 + 2*x**2 + 3)/4 + 203*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_103():
    f = x**3*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**2/2 - (25*x**2 + 75)/(8*x**4 + 16*x**2 + 24) - 17*log(x**4 + 2*x**2 + 3)/4 - 17*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_104():
    f = x*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = (25*x**2 + 25)/(8*x**4 + 16*x**2 + 24) + 5*log(x**4 + 2*x**2 + 3)/4 - 23*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_105():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x*(x**4 + 2*x**2 + 3)**2)
    F = (25 - 25*x**2)/(24*x**4 + 48*x**2 + 72) + 4*log(x)/9 - log(x**4 + 2*x**2 + 3)/9 + 89*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_106():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**3*(x**4 + 2*x**2 + 3)**2)
    F = -(25*x**2 + 125)/(72*x**4 + 144*x**2 + 216) - 13*log(x)/27 + 13*log(x**4 + 2*x**2 + 3)/108 - 71*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/432 - 2/(9*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_107():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**5*(x**4 + 2*x**2 + 3)**2)
    F = (125*x**2 + 175)/(216*x**4 + 432*x**2 + 648) + 13*log(x)/27 - 13*log(x**4 + 2*x**2 + 3)/108 + 125*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/432 + 13/(54*x**2) - 1/(9*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_108():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**7*(x**4 + 2*x**2 + 3)**2)
    F = (25 - 175*x**2)/(648*x**4 + 1296*x**2 + 1944) + 61*log(x)/243 - 61*log(x**4 + 2*x**2 + 3)/972 - 1237*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/3888 - 13/(54*x**2) + 13/(108*x**4) - 2/(27*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_109():
    f = x**8*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**7/7 - 17*x**5/5 + 19*x**3/3 + 25*x*(5*x**2 + 3)/(8*x**4 + 16*x**2 + 24) + 38*x - sqrt(sympy.S(-262771)/2 + 618291*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + sqrt(sympy.S(-262771)/2 + 618291*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + sqrt(sympy.S(262771)/2 + 618291*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16 - sqrt(sympy.S(262771)/2 + 618291*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_110():
    f = x**6*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = x**5 - 17*x**3/3 + 25*x*(3 - x**2)/(8*x**4 + 16*x**2 + 24) + 19*x + 3*sqrt(sympy.S(26007)/2 + 15033*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 - 3*sqrt(sympy.S(26007)/2 + 15033*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + 3*sqrt(sympy.S(-26007)/2 + 15033*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16 - 3*sqrt(sympy.S(-26007)/2 + 15033*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_111():
    f = x**4*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 5*x**3/3 - 25*x*(x**2 + 3)/(8*x**4 + 16*x**2 + 24) - 17*x - sqrt(sympy.S(-14395)/2 + 26499*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + sqrt(sympy.S(-14395)/2 + 26499*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 - sqrt(sympy.S(14395)/2 + 26499*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16 + sqrt(sympy.S(14395)/2 + 26499*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_112():
    f = x**2*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 25*x*(x**2 + 1)/(8*x**4 + 16*x**2 + 24) + 5*x - sqrt(sympy.S(-19291)/6 + 12899*sqrt(3)/6)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + sqrt(sympy.S(-19291)/6 + 12899*sqrt(3)/6)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/32 + sqrt(sympy.S(19291)/6 + 12899*sqrt(3)/6)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16 - sqrt(sympy.S(19291)/6 + 12899*sqrt(3)/6)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_113():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**2
    F = 25*x*(1 - x**2)/(24*x**4 + 48*x**2 + 72) + sqrt(sympy.S(11567)/6 + 4299*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/96 - sqrt(sympy.S(11567)/6 + 4299*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/96 - sqrt(sympy.S(-11567)/6 + 4299*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/48 + sqrt(sympy.S(-11567)/6 + 4299*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/48
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_114():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**2*(x**4 + 2*x**2 + 3)**2)
    F = -25*x*(x**2 + 5)/(72*x**4 + 144*x**2 + 216) - sqrt(sympy.S(965)/6 + 233*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/96 + sqrt(sympy.S(965)/6 + 233*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/96 + sqrt(sympy.S(-965)/6 + 233*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/48 - sqrt(sympy.S(-965)/6 + 233*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/48 - 4/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_115():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4*(x**4 + 2*x**2 + 3)**2)
    F = 25*x*(5*x**2 + 7)/(216*x**4 + 432*x**2 + 648) + sqrt(sympy.S(-6073)/6 + 18891*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/864 - sqrt(sympy.S(-6073)/6 + 18891*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/864 - sqrt(sympy.S(6073)/6 + 18891*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/432 + sqrt(sympy.S(6073)/6 + 18891*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/432 + 13/(27*x) - 4/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_116():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**6*(x**4 + 2*x**2 + 3)**2)
    F = 25*x*(1 - 7*x**2)/(648*x**4 + 1296*x**2 + 1944) - sqrt(sympy.S(1139381)/6 + 229473*sqrt(3)/2)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/2592 + sqrt(sympy.S(1139381)/6 + 229473*sqrt(3)/2)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/2592 + sqrt(sympy.S(-1139381)/6 + 229473*sqrt(3)/2)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/1296 - sqrt(sympy.S(-1139381)/6 + 229473*sqrt(3)/2)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/1296 - 13/(27*x) + 13/(81*x**3) - 4/(45*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_117():
    f = x**10*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = x**5 - 9*x**3 - 25*x*(7*x**2 + 15)/(16*(x**4 + 2*x**2 + 3)**2) + x*(252*x**2 + 3305)/(64*x**4 + 128*x**2 + 192) + 58*x + 3*sqrt(8595619 + 7678611*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - 3*sqrt(8595619 + 7678611*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 + 3*sqrt(-8595619 + 7678611*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256 - 3*sqrt(-8595619 + 7678611*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_118():
    f = x**8*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = 5*x**3/3 + 25*x*(5*x**2 + 3)/(16*(x**4 + 2*x**2 + 3)**2) - x*(835*x**2 + 1468)/(64*x**4 + 128*x**2 + 192) - 27*x - 21*sqrt(-34271 + 22721*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 + 21*sqrt(-34271 + 22721*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - 21*sqrt(34271 + 22721*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256 + 21*sqrt(34271 + 22721*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_119():
    f = x**6*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = 25*x*(3 - x**2)/(16*(x**4 + 2*x**2 + 3)**2) + 7*x*(58*x**2 + 11)/(64*x**4 + 128*x**2 + 192) + 5*x - sqrt(-827621 + 1176531*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 + sqrt(-827621 + 1176531*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 + sqrt(827621 + 1176531*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256 - sqrt(827621 + 1176531*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_120():
    f = x**4*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = x*(238 - 59*x**2)/(64*x**4 + 128*x**2 + 192) - 25*x*(x**2 + 3)/(16*(x**4 + 2*x**2 + 3)**2) + sqrt(146505 + 98481*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - sqrt(146505 + 98481*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - sqrt(-146505 + 98481*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256 + sqrt(-146505 + 98481*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_121():
    f = x**2*(5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = 25*x*(x**2 + 1)/(16*(x**4 + 2*x**2 + 3)**2) - x*(88*x**2 + 353)/(192*x**4 + 384*x**2 + 576) - 11*sqrt(sympy.S(1825)/3 + 363*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/1536 + 11*sqrt(sympy.S(1825)/3 + 363*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/1536 - 11*sqrt(sympy.S(-1825)/3 + 363*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/768 + 11*sqrt(sympy.S(-1825)/3 + 363*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/768
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_122():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4 + 2*x**2 + 3)**3
    F = 25*x*(1 - x**2)/(48*(x**4 + 2*x**2 + 3)**2) + x*(51*x**2 + 64)/(192*x**4 + 384*x**2 + 576) + sqrt(sympy.S(1291)/3 + 1019*sqrt(3)/3)*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - sqrt(sympy.S(1291)/3 + 1019*sqrt(3)/3)*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/512 - sqrt(sympy.S(-1291)/3 + 1019*sqrt(3)/3)*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256 + sqrt(sympy.S(-1291)/3 + 1019*sqrt(3)/3)*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_123():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**2*(x**4 + 2*x**2 + 3)**3)
    F = -25*x*(x**2 + 5)/(144*(x**4 + 2*x**2 + 3)**2) - x*(242*x**2 + 325)/(1728*x**4 + 3456*x**2 + 5184) - sqrt(sympy.S(-59711)/3 + 18387*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/4608 + sqrt(sympy.S(-59711)/3 + 18387*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/4608 + sqrt(sympy.S(59711)/3 + 18387*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/2304 - sqrt(sympy.S(59711)/3 + 18387*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/2304 - 4/(27*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_124():
    f = (5*x**6 + 3*x**4 + x**2 + 4)/(x**4*(x**4 + 2*x**2 + 3)**3)
    F = 25*x*(5*x**2 + 7)/(432*(x**4 + 2*x**2 + 3)**2) + x*(1025*x**2 + 1474)/(5184*x**4 + 10368*x**2 + 15552) + sqrt(sympy.S(-10004741)/3 + 3746817*sqrt(3))*log(x**2 - x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/41472 - sqrt(sympy.S(-10004741)/3 + 3746817*sqrt(3))*log(x**2 + x*sqrt(-2 + 2*sqrt(3)) + sqrt(3))/41472 - sqrt(sympy.S(10004741)/3 + 3746817*sqrt(3))*atan((-2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/20736 + sqrt(sympy.S(10004741)/3 + 3746817*sqrt(3))*atan((2*x + sqrt(-2 + 2*sqrt(3)))/sqrt(2 + 2*sqrt(3)))/20736 + 7/(27*x) - 4/(81*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_125():
    f = x*(d + e*x**2 + f*x**4 + g*x**6)/(a + b*x**2 + c*x**4)
    F = g*x**4/(4*c) + x**2*(-b*g + c*f)/(2*c**2) + (b**2*g + c**2*e - c*(a*g + b*f))*log(a + b*x**2 + c*x**4)/(4*c**3) - (-b**3*g + b*c*(3*a*g + b*f) + 2*c**3*d - c**2*(2*a*f + b*e))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_126():
    f = x**4*(d + e*x**2 + f*x**4 + g*x**6)/(a + b*x**2 + c*x**4)**2
    F = g*x**3/(3*c**2) + x*(-2*b*g + c*f)/c**3 + x*(a*(-b**3*g + b*c*(3*a*g + b*f) + 2*c**3*d - c**2*(2*a*f + b*e)) + x**2*(2*a*c**2*(-a*g + c*e) - b**4*g + b**3*c*f - b**2*c*(-4*a*g + c*e) + b*c**2*(-3*a*f + c*d)))/(2*c**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(2*a*c**2*(-7*a*g + 3*c*e) - 5*b**4*g + 3*b**3*c*f - b**2*c*(-24*a*g + c*e) - b*c**2*(13*a*f + c*d) + (4*a*b*c**2*(-13*a*g + 2*c*e) - 4*a*c**3*(-5*a*f + c*d) - 5*b**5*g + 3*b**4*c*f - b**3*c*(-34*a*g + c*e) - b**2*c**2*(19*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(7)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(2*a*c**2*(-7*a*g + 3*c*e) - 5*b**4*g + 3*b**3*c*f - b**2*c*(-24*a*g + c*e) - b*c**2*(13*a*f + c*d) - (4*a*b*c**2*(-13*a*g + 2*c*e) - 4*a*c**3*(-5*a*f + c*d) - 5*b**5*g + 3*b**4*c*f - b**3*c*(-34*a*g + c*e) - b**2*c**2*(19*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(7)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_127():
    f = x**2*(d + e*x**2 + f*x**4 + g*x**6)/(a + b*x**2 + c*x**4)**2
    F = g*x/c**2 - x*(-a*b**2*g - 2*a*c*(-a*g + c*e) + b*c*(a*f + c*d) + x**2*(-b**3*g + b*c*(3*a*g + b*f) + 2*c**3*d - c**2*(2*a*f + b*e)))/(2*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(3*b**3*g - b*c*(13*a*g + b*f) + 2*c**3*d - c**2*(-6*a*f + b*e) - (4*a*c**2*(-5*a*g + c*e) - 3*b**4*g + b**3*c*f + b**2*c*(19*a*g + c*e) - 4*b*c**2*(2*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(3*b**3*g - b*c*(13*a*g + b*f) + 2*c**3*d - c**2*(-6*a*f + b*e) + (4*a*c**2*(-5*a*g + c*e) - 3*b**4*g + b**3*c*f + b**2*c*(19*a*g + c*e) - 4*b*c**2*(2*a*f + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_128():
    f = (d + e*x**2 + f*x**4 + g*x**6)/(a + b*x**2 + c*x**4)**2
    F = x*(c*(-a*b*(a*g + c*e)/c - 2*a*(-a*f + c*d) + b**2*d) + x**2*(-a*b**2*g - 2*a*c*(-a*g + c*e) + b*c*(a*f + c*d)))/(2*a*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b**2*g/c - 2*a*(3*a*g + c*e) + b*(a*f + c*d) - (-a*b**3*g + 4*a*b*c*(2*a*g + c*e) - 4*a*c**2*(a*f + 3*c*d) + b**2*c*(-a*f + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b**2*g/c - 2*a*(3*a*g + c*e) + b*(a*f + c*d) + (-a*b**3*g + 4*a*b*c*(2*a*g + c*e) - 4*a*c**2*(a*f + 3*c*d) + b**2*c*(-a*f + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_129():
    f = (d + e*x**2 + f*x**4 + g*x**6)/(x**2*(a + b*x**2 + c*x**4)**2)
    F = -d/(a**2*x) - x*(a*(-2*a**2*g + a*(b*f + 2*c*e) - b*(b*e + 3*c*d) + b**3*d/a) + x**2*(-a*b*(a*g + c*e) - 2*a*c*(-a*f + c*d) + b**2*c*d))/(2*a**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*(-a*b*(a*g + c*e) - 2*a*c*(-a*f + 5*c*d) + 3*b**2*c*d - (4*a**2*c*(a*g + 3*c*e) - a*b**2*(-a*g + c*e) - 4*a*b*c*(a*f + 4*c*d) + 3*b**3*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(-a*b*(a*g + c*e) - 2*a*c*(-a*f + 5*c*d) + 3*b**2*c*d + (4*a**2*c*(a*g + 3*c*e) - a*b**2*(-a*g + c*e) - 4*a*b*c*(a*f + 4*c*d) + 3*b**3*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_130():
    f = (d + e*x**2 + f*x**4 + g*x**6)/(x**4*(a + b*x**2 + c*x**4)**2)
    F = -d/(3*a**2*x**3) + sqrt(2)*sqrt(c)*(2*a**2*(-a*g + 5*c*e) - 3*a*b**2*e - a*b*(-a*f + 19*c*d) + 5*b**3*d - (4*a**2*b*(a*g + 4*c*e) + 4*a**2*c*(-3*a*f + 7*c*d) - 3*a*b**3*e - a*b**2*(-a*f + 29*c*d) + 5*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(2*a**2*(-a*g + 5*c*e) - 3*a*b**2*e - a*b*(-a*f + 19*c*d) + 5*b**3*d + (4*a**2*b*(a*g + 4*c*e) + 4*a**2*c*(-3*a*f + 7*c*d) - 3*a*b**3*e - a*b**2*(-a*f + 29*c*d) + 5*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(a**2*(-a*(b*g + 2*c*f) + b**2*f + 3*b*c*e + 2*c**2*d - b**2*(b*e + 4*c*d)/a + b**4*d/a**2) + c*x**2*(2*a**2*(-a*g + c*e) - a*b**2*e - a*b*(-a*f + 3*c*d) + b**3*d))/(2*a**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (-a*e + 2*b*d)/(a**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_131():
    f = x**2*(a + b*x**2 + c*x**4)**p*(3*a + b*x**2*(2*p + 5) + c*x**4*(4*p + 7))
    F = x**3*(a + b*x**2 + c*x**4)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_132():
    f = (a + b*x**2 + c*x**4)/(x**4*sqrt(d - e*x)*sqrt(d + e*x))
    F = -a*(d**2 - e**2*x**2)/(3*d**2*x**3*sqrt(d - e*x)*sqrt(d + e*x)) + c*sqrt(d**2 - e**2*x**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(e*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(2*a*e**2 + 3*b*d**2)/(3*d**4*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_133():
    f = (a + b*x**2 + c*x**4)/(x**6*sqrt(d - e*x)*sqrt(d + e*x))
    F = -a*(d**2 - e**2*x**2)/(5*d**2*x**5*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(4*a*e**2 + 5*b*d**2)/(15*d**4*x**3*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(8*a*e**4 + 10*b*d**2*e**2 + 15*c*d**4)/(15*d**6*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_134():
    f = (a + b*x**2 + c*x**4)/(x**8*sqrt(d - e*x)*sqrt(d + e*x))
    F = -a*(d**2 - e**2*x**2)/(7*d**2*x**7*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(6*a*e**2 + 7*b*d**2)/(35*d**4*x**5*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(24*a*e**4 + 28*b*d**2*e**2 + 35*c*d**4)/(105*d**6*x**3*sqrt(d - e*x)*sqrt(d + e*x)) - 2*e**2*(d**2 - e**2*x**2)*(24*a*e**4 + 28*b*d**2*e**2 + 35*c*d**4)/(105*d**8*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_6_P_x_d_x_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_135():
    f = (a + b*x**2 + c*x**4)/(x**10*sqrt(d - e*x)*sqrt(d + e*x))
    F = -a*(d**2 - e**2*x**2)/(9*d**2*x**9*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(8*a*e**2 + 9*b*d**2)/(63*d**4*x**7*sqrt(d - e*x)*sqrt(d + e*x)) - (d**2 - e**2*x**2)*(16*a*e**4 + 18*b*d**2*e**2 + 21*c*d**4)/(105*d**6*x**5*sqrt(d - e*x)*sqrt(d + e*x)) - 4*e**2*(d**2 - e**2*x**2)*(16*a*e**4 + 18*b*d**2*e**2 + 21*c*d**4)/(315*d**8*x**3*sqrt(d - e*x)*sqrt(d + e*x)) - 8*e**4*(d**2 - e**2*x**2)*(16*a*e**4 + 18*b*d**2*e**2 + 21*c*d**4)/(315*d**10*x*sqrt(d - e*x)*sqrt(d + e*x))
    assert integrate(f, x) == F

