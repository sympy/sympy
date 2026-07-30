"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.4 Improper/1.1.4.2 (c x)^m (a x^j+b x^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, j, m, n, p, q = symbols('a b c j m n p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_1():
    f = x**2*(a*x + b*x**3)
    F = a*x**4/4 + b*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_2():
    f = x*(a*x + b*x**3)
    F = a*x**3/3 + b*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_3():
    f = a*x + b*x**3
    F = a*x**2/2 + b*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_4():
    f = (a*x + b*x**3)/x
    F = a*x + b*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_5():
    f = (a*x + b*x**3)/x**2
    F = a*log(x) + b*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_6():
    f = x**2*(a*x + b*x**3)**2
    F = a**2*x**5/5 + 2*a*b*x**7/7 + b**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_7():
    f = x*(a*x + b*x**3)**2
    F = a**2*x**4/4 + a*b*x**6/3 + b**2*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_8():
    f = (a*x + b*x**3)**2
    F = a**2*x**3/3 + 2*a*b*x**5/5 + b**2*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_9():
    f = (a*x + b*x**3)**2/x
    F = (a + b*x**2)**3/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_10():
    f = (a*x + b*x**3)**2/x**2
    F = a**2*x + 2*a*b*x**3/3 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_11():
    f = (3*x**3 - 4*x)**6
    F = 729*x**19/19 - 5832*x**17/17 + 1296*x**15 - 34560*x**13/13 + 34560*x**11/11 - 2048*x**9 + 4096*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_12():
    f = x**4/(a*x + b*x**3)
    F = -a*log(a + b*x**2)/(2*b**2) + x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_13():
    f = x**3/(a*x + b*x**3)
    F = -sqrt(a)*atan(sqrt(b)*x/sqrt(a))/b**(sympy.S(3)/2) + x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_14():
    f = x**2/(a*x + b*x**3)
    F = log(a + b*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_15():
    f = x/(a*x + b*x**3)
    F = atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_16():
    f = 1/(a*x + b*x**3)
    F = log(x)/a - log(a + b*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_17():
    f = 1/(x*(a*x + b*x**3))
    F = -1/(a*x) - sqrt(b)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_18():
    f = 1/(x**2*(a*x + b*x**3))
    F = -1/(2*a*x**2) - b*log(x)/a**2 + b*log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_19():
    f = 1/(x**3*(a*x + b*x**3))
    F = -1/(3*a*x**3) + b/(a**2*x) + b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_20():
    f = 1/(x**4*(a*x + b*x**3))
    F = -1/(4*a*x**4) + b/(2*a**2*x**2) + b**2*log(x)/a**3 - b**2*log(a + b*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_21():
    f = x**2/(a*x + b*x**3)**2
    F = x/(2*a*(a + b*x**2)) + atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_22():
    f = x/(a*x + b*x**3)**2
    F = 1/(2*a*(a + b*x**2)) + log(x)/a**2 - log(a + b*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_23():
    f = (a*x + b*x**3)**(-2)
    F = 1/(2*a*x*(a + b*x**2)) - 3/(2*a**2*x) - 3*sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_24():
    f = 1/(x*(a*x + b*x**3)**2)
    F = -b/(2*a**2*(a + b*x**2)) - 1/(2*a**2*x**2) - 2*b*log(x)/a**3 + b*log(a + b*x**2)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_25():
    f = 1/(x**2*(a*x + b*x**3)**2)
    F = 1/(2*a*x**3*(a + b*x**2)) - 5/(6*a**2*x**3) + 5*b/(2*a**3*x) + 5*b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_26():
    f = x**5/(-x**3 + x)
    F = -x**3/3 - x + atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_27():
    f = x**4/(-x**3 + x)
    F = -x**2/2 - log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_28():
    f = x**3/(-x**3 + x)
    F = -x + atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_29():
    f = x**2/(-x**3 + x)
    F = -log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_30():
    f = x/(-x**3 + x)
    F = atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_31():
    f = 1/(-x**3 + x)
    F = log(x) - log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_32():
    f = 1/(x*(-x**3 + x))
    F = atanh(x) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_33():
    f = 1/(x**2*(-x**3 + x))
    F = log(x) - log(1 - x**2)/2 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_34():
    f = 1/(x**3*(-x**3 + x))
    F = atanh(x) - 1/x - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_35():
    f = 1/(x**4*(-x**3 + x))
    F = log(x) - log(1 - x**2)/2 - 1/(2*x**2) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_36():
    f = 1/(b*x**3 + x)
    F = log(x) - log(b*x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_37():
    f = 1/(b*x**3 - x)
    F = -log(x) + log(-b*x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_38():
    f = x**3*sqrt(a*x + b*x**3)
    F = 10*a**(sympy.S(11)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a*x + b*x**3)) - 20*a**2*sqrt(a*x + b*x**3)/(231*b**2) + 4*a*x**2*sqrt(a*x + b*x**3)/(77*b) + 2*x**4*sqrt(a*x + b*x**3)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_39():
    f = x**2*sqrt(a*x + b*x**3)
    F = 4*a**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 2*a**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 4*a**2*x*(a + b*x**2)/(15*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 4*a*x*sqrt(a*x + b*x**3)/(45*b) + 2*x**3*sqrt(a*x + b*x**3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_40():
    f = x*sqrt(a*x + b*x**3)
    F = -2*a**(sympy.S(7)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(5)/4)*sqrt(a*x + b*x**3)) + 4*a*sqrt(a*x + b*x**3)/(21*b) + 2*x**2*sqrt(a*x + b*x**3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_41():
    f = sqrt(a*x + b*x**3)
    F = -4*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 2*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 4*a*x*(a + b*x**2)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 2*x*sqrt(a*x + b*x**3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_42():
    f = sqrt(a*x + b*x**3)/x
    F = 2*a**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*sqrt(a*x + b*x**3)) + 2*sqrt(a*x + b*x**3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_43():
    f = sqrt(a*x + b*x**3)/x**2
    F = -4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a*x + b*x**3) + 2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a*x + b*x**3) + 4*sqrt(b)*x*(a + b*x**2)/((sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - 2*sqrt(a*x + b*x**3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_44():
    f = sqrt(a*x + b*x**3)/x**3
    F = -2*sqrt(a*x + b*x**3)/(3*x**2) + 2*b**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(3*a**(sympy.S(1)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_45():
    f = sqrt(a*x + b*x**3)/x**4
    F = -2*sqrt(a*x + b*x**3)/(5*x**3) + 4*b**(sympy.S(3)/2)*x*(a + b*x**2)/(5*a*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - 4*b*sqrt(a*x + b*x**3)/(5*a*x) - 4*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 2*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_46():
    f = x**2*(a*x + b*x**3)**(sympy.S(3)/2)
    F = 4*a**(sympy.S(15)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(231*b**(sympy.S(9)/4)*sqrt(a*x + b*x**3)) - 8*a**3*sqrt(a*x + b*x**3)/(231*b**2) + 8*a**2*x**2*sqrt(a*x + b*x**3)/(385*b) + 4*a*x**4*sqrt(a*x + b*x**3)/55 + 2*x**3*(a*x + b*x**3)**(sympy.S(3)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_47():
    f = x*(a*x + b*x**3)**(sympy.S(3)/2)
    F = 8*a**(sympy.S(13)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 4*a**(sympy.S(13)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 8*a**3*x*(a + b*x**2)/(65*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 8*a**2*x*sqrt(a*x + b*x**3)/(195*b) + 4*a*x**3*sqrt(a*x + b*x**3)/39 + 2*x**2*(a*x + b*x**3)**(sympy.S(3)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_48():
    f = (a*x + b*x**3)**(sympy.S(3)/2)
    F = -4*a**(sympy.S(11)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(77*b**(sympy.S(5)/4)*sqrt(a*x + b*x**3)) + 8*a**2*sqrt(a*x + b*x**3)/(77*b) + 12*a*x**2*sqrt(a*x + b*x**3)/77 + 2*x*(a*x + b*x**3)**(sympy.S(3)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_49():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x
    F = -8*a**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 4*a**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 8*a**2*x*(a + b*x**2)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 4*a*x*sqrt(a*x + b*x**3)/15 + 2*(a*x + b*x**3)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_50():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**2
    F = 4*a**(sympy.S(7)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(7*b**(sympy.S(1)/4)*sqrt(a*x + b*x**3)) + 4*a*sqrt(a*x + b*x**3)/7 + 2*(a*x + b*x**3)**(sympy.S(3)/2)/(7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_51():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**3
    F = -24*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a*x + b*x**3)) + 12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a*x + b*x**3)) + 24*a*sqrt(b)*x*(a + b*x**2)/((5*sqrt(a) + 5*sqrt(b)*x)*sqrt(a*x + b*x**3)) + 12*b*x*sqrt(a*x + b*x**3)/5 - 2*(a*x + b*x**3)**(sympy.S(3)/2)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_52():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**4
    F = 4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(3*sqrt(a*x + b*x**3)) + 4*b*sqrt(a*x + b*x**3)/3 - 2*(a*x + b*x**3)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_53():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**5
    F = -24*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a*x + b*x**3)) + 12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a*x + b*x**3)) + 24*b**(sympy.S(3)/2)*x*(a + b*x**2)/((5*sqrt(a) + 5*sqrt(b)*x)*sqrt(a*x + b*x**3)) - 12*b*sqrt(a*x + b*x**3)/(5*x) - 2*(a*x + b*x**3)**(sympy.S(3)/2)/(5*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_54():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**6
    F = -4*b*sqrt(a*x + b*x**3)/(7*x**2) - 2*(a*x + b*x**3)**(sympy.S(3)/2)/(7*x**5) + 4*b**(sympy.S(7)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(7*a**(sympy.S(1)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_55():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**7
    F = -4*b*sqrt(a*x + b*x**3)/(15*x**3) - 2*(a*x + b*x**3)**(sympy.S(3)/2)/(9*x**6) + 8*b**(sympy.S(5)/2)*x*(a + b*x**2)/(15*a*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - 8*b**2*sqrt(a*x + b*x**3)/(15*a*x) - 8*b**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 4*b**(sympy.S(9)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_56():
    f = (a*x + b*x**3)**(sympy.S(3)/2)/x**8
    F = -12*b*sqrt(a*x + b*x**3)/(77*x**4) - 2*(a*x + b*x**3)**(sympy.S(3)/2)/(11*x**7) - 8*b**2*sqrt(a*x + b*x**3)/(77*a*x**2) - 4*b**(sympy.S(11)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(77*a**(sympy.S(5)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_57():
    f = x**4/sqrt(a*x + b*x**3)
    F = 5*a**(sympy.S(7)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(21*b**(sympy.S(9)/4)*sqrt(a*x + b*x**3)) - 10*a*sqrt(a*x + b*x**3)/(21*b**2) + 2*x**2*sqrt(a*x + b*x**3)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_58():
    f = x**3/sqrt(a*x + b*x**3)
    F = 6*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 3*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 6*a*x*(a + b*x**2)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 2*x*sqrt(a*x + b*x**3)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_59():
    f = x**2/sqrt(a*x + b*x**3)
    F = -a**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(5)/4)*sqrt(a*x + b*x**3)) + 2*sqrt(a*x + b*x**3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_60():
    f = x/sqrt(a*x + b*x**3)
    F = -2*a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + 2*x*(a + b*x**2)/(sqrt(b)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_61():
    f = 1/sqrt(a*x + b*x**3)
    F = sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_62():
    f = 1/(x*sqrt(a*x + b*x**3))
    F = 2*sqrt(b)*x*(a + b*x**2)/(a*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - 2*sqrt(a*x + b*x**3)/(a*x) - 2*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) + b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_63():
    f = 1/(x**2*sqrt(a*x + b*x**3))
    F = -2*sqrt(a*x + b*x**3)/(3*a*x**2) - b**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(3*a**(sympy.S(5)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_64():
    f = 1/(x**3*sqrt(a*x + b*x**3))
    F = -2*sqrt(a*x + b*x**3)/(5*a*x**3) - 6*b**(sympy.S(3)/2)*x*(a + b*x**2)/(5*a**2*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 6*b*sqrt(a*x + b*x**3)/(5*a**2*x) + 6*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - 3*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(7)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_65():
    f = x**7/(a*x + b*x**3)**(sympy.S(3)/2)
    F = 15*a**(sympy.S(7)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(14*b**(sympy.S(13)/4)*sqrt(a*x + b*x**3)) - 15*a*sqrt(a*x + b*x**3)/(7*b**3) - x**5/(b*sqrt(a*x + b*x**3)) + 9*x**2*sqrt(a*x + b*x**3)/(7*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_66():
    f = x**6/(a*x + b*x**3)**(sympy.S(3)/2)
    F = 21*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(11)/4)*sqrt(a*x + b*x**3)) - 21*a**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(10*b**(sympy.S(11)/4)*sqrt(a*x + b*x**3)) - 21*a*x*(a + b*x**2)/(5*b**(sympy.S(5)/2)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - x**4/(b*sqrt(a*x + b*x**3)) + 7*x*sqrt(a*x + b*x**3)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_67():
    f = x**5/(a*x + b*x**3)**(sympy.S(3)/2)
    F = -5*a**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(9)/4)*sqrt(a*x + b*x**3)) - x**3/(b*sqrt(a*x + b*x**3)) + 5*sqrt(a*x + b*x**3)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_68():
    f = x**4/(a*x + b*x**3)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) + 3*a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) - x**2/(b*sqrt(a*x + b*x**3)) + 3*x*(a + b*x**2)/(b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_69():
    f = x**3/(a*x + b*x**3)**(sympy.S(3)/2)
    F = -x/(b*sqrt(a*x + b*x**3)) + sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_70():
    f = x**2/(a*x + b*x**3)**(sympy.S(3)/2)
    F = x**2/(a*sqrt(a*x + b*x**3)) - x*(a + b*x**2)/(a*sqrt(b)*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3)) - sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_71():
    f = x/(a*x + b*x**3)**(sympy.S(3)/2)
    F = x/(a*sqrt(a*x + b*x**3)) + sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_72():
    f = (a*x + b*x**3)**(sympy.S(-3)/2)
    F = 1/(a*sqrt(a*x + b*x**3)) + 3*sqrt(b)*x*(a + b*x**2)/(a**2*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) - 3*sqrt(a*x + b*x**3)/(a**2*x) - 3*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(7)/4)*sqrt(a*x + b*x**3)) + 3*b**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(7)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_73():
    f = 1/(x*(a*x + b*x**3)**(sympy.S(3)/2))
    F = 1/(a*x*sqrt(a*x + b*x**3)) - 5*sqrt(a*x + b*x**3)/(3*a**2*x**2) - 5*b**(sympy.S(3)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(6*a**(sympy.S(9)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_74():
    f = 1/(x**2*(a*x + b*x**3)**(sympy.S(3)/2))
    F = 1/(a*x**2*sqrt(a*x + b*x**3)) - 7*sqrt(a*x + b*x**3)/(5*a**2*x**3) - 21*b**(sympy.S(3)/2)*x*(a + b*x**2)/(5*a**3*(sqrt(a) + sqrt(b)*x)*sqrt(a*x + b*x**3)) + 21*b*sqrt(a*x + b*x**3)/(5*a**3*x) + 21*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_e(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(11)/4)*sqrt(a*x + b*x**3)) - 21*b**(sympy.S(5)/4)*sqrt(x)*sqrt((a + b*x**2)/(sqrt(a) + sqrt(b)*x)**2)*(sqrt(a) + sqrt(b)*x)*elliptic_f(2*atan(b**(sympy.S(1)/4)*sqrt(x)/a**(sympy.S(1)/4)), sympy.S.Half)/(10*a**(sympy.S(11)/4)*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_75():
    f = x**(sympy.S(29)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -9*a*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x + b*x**3))/(2*b**(sympy.S(11)/2)) - x**(sympy.S(25)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - 9*x**(sympy.S(19)/2)/(35*b**2*(a*x + b*x**3)**(sympy.S(5)/2)) - 3*x**(sympy.S(13)/2)/(5*b**3*(a*x + b*x**3)**(sympy.S(3)/2)) - 3*x**(sympy.S(7)/2)/(b**4*sqrt(a*x + b*x**3)) + 9*sqrt(x)*sqrt(a*x + b*x**3)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_76():
    f = x**(sympy.S(27)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(23)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - 8*x**(sympy.S(17)/2)/(35*b**2*(a*x + b*x**3)**(sympy.S(5)/2)) - 16*x**(sympy.S(11)/2)/(35*b**3*(a*x + b*x**3)**(sympy.S(3)/2)) - 64*x**(sympy.S(5)/2)/(35*b**4*sqrt(a*x + b*x**3)) + 128*sqrt(a*x + b*x**3)/(35*b**5*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_77():
    f = x**(sympy.S(25)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(21)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - x**(sympy.S(15)/2)/(5*b**2*(a*x + b*x**3)**(sympy.S(5)/2)) - x**(sympy.S(9)/2)/(3*b**3*(a*x + b*x**3)**(sympy.S(3)/2)) - x**(sympy.S(3)/2)/(b**4*sqrt(a*x + b*x**3)) + atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x + b*x**3))/b**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_78():
    f = x**(sympy.S(23)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(19)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - 6*x**(sympy.S(13)/2)/(35*b**2*(a*x + b*x**3)**(sympy.S(5)/2)) - 8*x**(sympy.S(7)/2)/(35*b**3*(a*x + b*x**3)**(sympy.S(3)/2)) - 16*sqrt(x)/(35*b**4*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_79():
    f = x**(sympy.S(21)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(21)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_80():
    f = x**(sympy.S(19)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(15)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - 4*x**(sympy.S(9)/2)/(35*b**2*(a*x + b*x**3)**(sympy.S(5)/2)) - 8*x**(sympy.S(3)/2)/(105*b**3*(a*x + b*x**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_81():
    f = x**(sympy.S(17)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(17)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 2*x**(sympy.S(15)/2)/(35*a**2*(a*x + b*x**3)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_82():
    f = x**(sympy.S(15)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(11)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2)) - 2*x**(sympy.S(5)/2)/(35*b**2*(a*x + b*x**3)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_83():
    f = x**(sympy.S(13)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(13)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 4*x**(sympy.S(11)/2)/(35*a**2*(a*x + b*x**3)**(sympy.S(5)/2)) + 8*x**(sympy.S(9)/2)/(105*a**3*(a*x + b*x**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_84():
    f = x**(sympy.S(11)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = -x**(sympy.S(7)/2)/(7*b*(a*x + b*x**3)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_85():
    f = x**(sympy.S(9)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(9)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 6*x**(sympy.S(7)/2)/(35*a**2*(a*x + b*x**3)**(sympy.S(5)/2)) + 8*x**(sympy.S(5)/2)/(35*a**3*(a*x + b*x**3)**(sympy.S(3)/2)) + 16*x**(sympy.S(3)/2)/(35*a**4*sqrt(a*x + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_86():
    f = x**(sympy.S(7)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(7)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + x**(sympy.S(5)/2)/(5*a**2*(a*x + b*x**3)**(sympy.S(5)/2)) + x**(sympy.S(3)/2)/(3*a**3*(a*x + b*x**3)**(sympy.S(3)/2)) + sqrt(x)/(a**4*sqrt(a*x + b*x**3)) - atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**3))/a**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_87():
    f = x**(sympy.S(5)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(5)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 8*x**(sympy.S(3)/2)/(35*a**2*(a*x + b*x**3)**(sympy.S(5)/2)) + 16*sqrt(x)/(35*a**3*(a*x + b*x**3)**(sympy.S(3)/2)) + 64/(35*a**4*sqrt(x)*sqrt(a*x + b*x**3)) - 128*sqrt(a*x + b*x**3)/(35*a**5*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_88():
    f = x**(sympy.S(3)/2)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = x**(sympy.S(3)/2)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 9*sqrt(x)/(35*a**2*(a*x + b*x**3)**(sympy.S(5)/2)) + 3/(5*a**3*sqrt(x)*(a*x + b*x**3)**(sympy.S(3)/2)) + 3/(a**4*x**(sympy.S(3)/2)*sqrt(a*x + b*x**3)) - 9*sqrt(a*x + b*x**3)/(2*a**5*x**(sympy.S(5)/2)) + 9*b*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**3))/(2*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_89():
    f = sqrt(x)/(a*x + b*x**3)**(sympy.S(9)/2)
    F = sqrt(x)/(7*a*(a*x + b*x**3)**(sympy.S(7)/2)) + 2/(7*a**2*sqrt(x)*(a*x + b*x**3)**(sympy.S(5)/2)) + 16/(21*a**3*x**(sympy.S(3)/2)*(a*x + b*x**3)**(sympy.S(3)/2)) + 32/(7*a**4*x**(sympy.S(5)/2)*sqrt(a*x + b*x**3)) - 128*sqrt(a*x + b*x**3)/(21*a**5*x**(sympy.S(7)/2)) + 256*b*sqrt(a*x + b*x**3)/(21*a**6*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_90():
    f = 1/(sqrt(x)*(a*x + b*x**3)**(sympy.S(9)/2))
    F = 1/(7*a*sqrt(x)*(a*x + b*x**3)**(sympy.S(7)/2)) + 11/(35*a**2*x**(sympy.S(3)/2)*(a*x + b*x**3)**(sympy.S(5)/2)) + 33/(35*a**3*x**(sympy.S(5)/2)*(a*x + b*x**3)**(sympy.S(3)/2)) + 33/(5*a**4*x**(sympy.S(7)/2)*sqrt(a*x + b*x**3)) - 33*sqrt(a*x + b*x**3)/(4*a**5*x**(sympy.S(9)/2)) + 99*b*sqrt(a*x + b*x**3)/(8*a**6*x**(sympy.S(5)/2)) - 99*b**2*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**3))/(8*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_91():
    f = 1/(x**(sympy.S(3)/2)*(a*x + b*x**3)**(sympy.S(9)/2))
    F = 1/(7*a*x**(sympy.S(3)/2)*(a*x + b*x**3)**(sympy.S(7)/2)) + 12/(35*a**2*x**(sympy.S(5)/2)*(a*x + b*x**3)**(sympy.S(5)/2)) + 8/(7*a**3*x**(sympy.S(7)/2)*(a*x + b*x**3)**(sympy.S(3)/2)) + 64/(7*a**4*x**(sympy.S(9)/2)*sqrt(a*x + b*x**3)) - 384*sqrt(a*x + b*x**3)/(35*a**5*x**(sympy.S(11)/2)) + 512*b*sqrt(a*x + b*x**3)/(35*a**6*x**(sympy.S(7)/2)) - 1024*b**2*sqrt(a*x + b*x**3)/(35*a**7*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_92():
    f = x**4/sqrt(a*x + b*x**4)
    F = -a*atanh(sqrt(b)*x**2/sqrt(a*x + b*x**4))/(3*b**(sympy.S(3)/2)) + x*sqrt(a*x + b*x**4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_93():
    f = x/sqrt(a*x + b*x**4)
    F = 2*atanh(sqrt(b)*x**2/sqrt(a*x + b*x**4))/(3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_94():
    f = 1/(x**2*sqrt(a*x + b*x**4))
    F = -2*sqrt(a*x + b*x**4)/(3*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_95():
    f = 1/(x**5*sqrt(a*x + b*x**4))
    F = -2*sqrt(a*x + b*x**4)/(9*a*x**5) + 4*b*sqrt(a*x + b*x**4)/(9*a**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_96():
    f = 1/(x**8*sqrt(a*x + b*x**4))
    F = -2*sqrt(a*x + b*x**4)/(15*a*x**8) + 8*b*sqrt(a*x + b*x**4)/(45*a**2*x**5) - 16*b**2*sqrt(a*x + b*x**4)/(45*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_97():
    f = x**3/sqrt(a*x + b*x**4)
    F = -3**(sympy.S(3)/4)*a**(sympy.S(2)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(12*b*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) + sqrt(a*x + b*x**4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_98():
    f = 1/sqrt(a*x + b*x**4)
    F = 3**(sympy.S(3)/4)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(1)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_99():
    f = 1/(x**3*sqrt(a*x + b*x**4))
    F = -2*sqrt(a*x + b*x**4)/(5*a*x**3) - 2*3**(sympy.S(3)/4)*b*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(15*a**(sympy.S(4)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_100():
    f = x**5/sqrt(a*x + b*x**4)
    F = 5*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) + 3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(5 - 5*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(48*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) - a*x*(5 + 5*sqrt(3))*(a + b*x**3)/(8*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x + b*x**4)) + x**2*sqrt(a*x + b*x**4)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_101():
    f = x**2/sqrt(a*x + b*x**4)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(b**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) - 3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(6*b**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) + x*(1 + sqrt(3))*(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_102():
    f = 1/(x*sqrt(a*x + b*x**4))
    F = b**(sympy.S(1)/3)*x*(2 + 2*sqrt(3))*(a + b*x**3)/(a*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x + b*x**4)) - 2*sqrt(a*x + b*x**4)/(a*x) - 2*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4)) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_103():
    f = x**2/sqrt(a*x + b*sqrt(x))
    F = 2*x**2*sqrt(a*x + b*sqrt(x))/(5*a) - 9*b*x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))/(20*a**2) + 21*b**2*x*sqrt(a*x + b*sqrt(x))/(40*a**3) - 21*b**3*sqrt(x)*sqrt(a*x + b*sqrt(x))/(32*a**4) + 63*b**4*sqrt(a*x + b*sqrt(x))/(64*a**5) - 63*b**5*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(64*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_104():
    f = x/sqrt(a*x + b*sqrt(x))
    F = 2*x*sqrt(a*x + b*sqrt(x))/(3*a) - 5*b*sqrt(x)*sqrt(a*x + b*sqrt(x))/(6*a**2) + 5*b**2*sqrt(a*x + b*sqrt(x))/(4*a**3) - 5*b**3*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_105():
    f = 1/sqrt(a*x + b*sqrt(x))
    F = 2*sqrt(a*x + b*sqrt(x))/a - 2*b*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_106():
    f = 1/(x*sqrt(a*x + b*sqrt(x)))
    F = -4*sqrt(a*x + b*sqrt(x))/(b*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_107():
    f = 1/(x**2*sqrt(a*x + b*sqrt(x)))
    F = -32*a**2*sqrt(a*x + b*sqrt(x))/(15*b**3*sqrt(x)) + 16*a*sqrt(a*x + b*sqrt(x))/(15*b**2*x) - 4*sqrt(a*x + b*sqrt(x))/(5*b*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_108():
    f = 1/(x**3*sqrt(a*x + b*sqrt(x)))
    F = -512*a**4*sqrt(a*x + b*sqrt(x))/(315*b**5*sqrt(x)) + 256*a**3*sqrt(a*x + b*sqrt(x))/(315*b**4*x) - 64*a**2*sqrt(a*x + b*sqrt(x))/(105*b**3*x**(sympy.S(3)/2)) + 32*a*sqrt(a*x + b*sqrt(x))/(63*b**2*x**2) - 4*sqrt(a*x + b*sqrt(x))/(9*b*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_109():
    f = 1/(x**4*sqrt(a*x + b*sqrt(x)))
    F = -4096*a**6*sqrt(a*x + b*sqrt(x))/(3003*b**7*sqrt(x)) + 2048*a**5*sqrt(a*x + b*sqrt(x))/(3003*b**6*x) - 512*a**4*sqrt(a*x + b*sqrt(x))/(1001*b**5*x**(sympy.S(3)/2)) + 1280*a**3*sqrt(a*x + b*sqrt(x))/(3003*b**4*x**2) - 160*a**2*sqrt(a*x + b*sqrt(x))/(429*b**3*x**(sympy.S(5)/2)) + 48*a*sqrt(a*x + b*sqrt(x))/(143*b**2*x**3) - 4*sqrt(a*x + b*sqrt(x))/(13*b*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_110():
    f = x**3/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*x**3/(a*sqrt(a*x + b*sqrt(x))) + 22*x**2*sqrt(a*x + b*sqrt(x))/(5*a**2) - 99*b*x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))/(20*a**3) + 231*b**2*x*sqrt(a*x + b*sqrt(x))/(40*a**4) - 231*b**3*sqrt(x)*sqrt(a*x + b*sqrt(x))/(32*a**5) + 693*b**4*sqrt(a*x + b*sqrt(x))/(64*a**6) - 693*b**5*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(64*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_111():
    f = x**2/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*x**2/(a*sqrt(a*x + b*sqrt(x))) + 14*x*sqrt(a*x + b*sqrt(x))/(3*a**2) - 35*b*sqrt(x)*sqrt(a*x + b*sqrt(x))/(6*a**3) + 35*b**2*sqrt(a*x + b*sqrt(x))/(4*a**4) - 35*b**3*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(4*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_112():
    f = x/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*x/(a*sqrt(a*x + b*sqrt(x))) + 6*sqrt(a*x + b*sqrt(x))/a**2 - 6*b*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_113():
    f = (a*x + b*sqrt(x))**(sympy.S(-3)/2)
    F = 4*sqrt(x)/(b*sqrt(a*x + b*sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_114():
    f = 1/(x*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = 32*a*sqrt(a*x + b*sqrt(x))/(3*b**3*sqrt(x)) + 4/(b*sqrt(x)*sqrt(a*x + b*sqrt(x))) - 16*sqrt(a*x + b*sqrt(x))/(3*b**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_115():
    f = 1/(x**2*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = 512*a**3*sqrt(a*x + b*sqrt(x))/(35*b**5*sqrt(x)) - 256*a**2*sqrt(a*x + b*sqrt(x))/(35*b**4*x) + 192*a*sqrt(a*x + b*sqrt(x))/(35*b**3*x**(sympy.S(3)/2)) + 4/(b*x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))) - 32*sqrt(a*x + b*sqrt(x))/(7*b**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_116():
    f = 1/(x**3*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = 4096*a**5*sqrt(a*x + b*sqrt(x))/(231*b**7*sqrt(x)) - 2048*a**4*sqrt(a*x + b*sqrt(x))/(231*b**6*x) + 512*a**3*sqrt(a*x + b*sqrt(x))/(77*b**5*x**(sympy.S(3)/2)) - 1280*a**2*sqrt(a*x + b*sqrt(x))/(231*b**4*x**2) + 160*a*sqrt(a*x + b*sqrt(x))/(33*b**3*x**(sympy.S(5)/2)) + 4/(b*x**(sympy.S(5)/2)*sqrt(a*x + b*sqrt(x))) - 48*sqrt(a*x + b*sqrt(x))/(11*b**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_117():
    f = x**(sympy.S(5)/2)/sqrt(a*x + b*sqrt(x))
    F = x**(sympy.S(5)/2)*sqrt(a*x + b*sqrt(x))/(3*a) - 11*b*x**2*sqrt(a*x + b*sqrt(x))/(30*a**2) + 33*b**2*x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))/(80*a**3) - 77*b**3*x*sqrt(a*x + b*sqrt(x))/(160*a**4) + 77*b**4*sqrt(x)*sqrt(a*x + b*sqrt(x))/(128*a**5) - 231*b**5*sqrt(a*x + b*sqrt(x))/(256*a**6) + 231*b**6*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(256*a**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_118():
    f = x**(sympy.S(3)/2)/sqrt(a*x + b*sqrt(x))
    F = x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))/(2*a) - 7*b*x*sqrt(a*x + b*sqrt(x))/(12*a**2) + 35*b**2*sqrt(x)*sqrt(a*x + b*sqrt(x))/(48*a**3) - 35*b**3*sqrt(a*x + b*sqrt(x))/(32*a**4) + 35*b**4*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(32*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_119():
    f = sqrt(x)/sqrt(a*x + b*sqrt(x))
    F = sqrt(x)*sqrt(a*x + b*sqrt(x))/a - 3*b*sqrt(a*x + b*sqrt(x))/(2*a**2) + 3*b**2*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_120():
    f = 1/(sqrt(x)*sqrt(a*x + b*sqrt(x)))
    F = 4*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_121():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x)))
    F = 8*a*sqrt(a*x + b*sqrt(x))/(3*b**2*sqrt(x)) - 4*sqrt(a*x + b*sqrt(x))/(3*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_122():
    f = 1/(x**(sympy.S(5)/2)*sqrt(a*x + b*sqrt(x)))
    F = 64*a**3*sqrt(a*x + b*sqrt(x))/(35*b**4*sqrt(x)) - 32*a**2*sqrt(a*x + b*sqrt(x))/(35*b**3*x) + 24*a*sqrt(a*x + b*sqrt(x))/(35*b**2*x**(sympy.S(3)/2)) - 4*sqrt(a*x + b*sqrt(x))/(7*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_123():
    f = 1/(x**(sympy.S(7)/2)*sqrt(a*x + b*sqrt(x)))
    F = 1024*a**5*sqrt(a*x + b*sqrt(x))/(693*b**6*sqrt(x)) - 512*a**4*sqrt(a*x + b*sqrt(x))/(693*b**5*x) + 128*a**3*sqrt(a*x + b*sqrt(x))/(231*b**4*x**(sympy.S(3)/2)) - 320*a**2*sqrt(a*x + b*sqrt(x))/(693*b**3*x**2) + 40*a*sqrt(a*x + b*sqrt(x))/(99*b**2*x**(sympy.S(5)/2)) - 4*sqrt(a*x + b*sqrt(x))/(11*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_124():
    f = x**(sympy.S(5)/2)/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*x**(sympy.S(5)/2)/(a*sqrt(a*x + b*sqrt(x))) + 9*x**(sympy.S(3)/2)*sqrt(a*x + b*sqrt(x))/(2*a**2) - 21*b*x*sqrt(a*x + b*sqrt(x))/(4*a**3) + 105*b**2*sqrt(x)*sqrt(a*x + b*sqrt(x))/(16*a**4) - 315*b**3*sqrt(a*x + b*sqrt(x))/(32*a**5) + 315*b**4*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(32*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_125():
    f = x**(sympy.S(3)/2)/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*x**(sympy.S(3)/2)/(a*sqrt(a*x + b*sqrt(x))) + 5*sqrt(x)*sqrt(a*x + b*sqrt(x))/a**2 - 15*b*sqrt(a*x + b*sqrt(x))/(2*a**3) + 15*b**2*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_126():
    f = sqrt(x)/(a*x + b*sqrt(x))**(sympy.S(3)/2)
    F = -4*sqrt(x)/(a*sqrt(a*x + b*sqrt(x))) + 4*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*sqrt(x)))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_127():
    f = 1/(sqrt(x)*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = -(8*a*sqrt(x) + 4*b)/(b**2*sqrt(a*x + b*sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_128():
    f = 1/(x**(sympy.S(3)/2)*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = -64*a**2*sqrt(a*x + b*sqrt(x))/(5*b**4*sqrt(x)) + 32*a*sqrt(a*x + b*sqrt(x))/(5*b**3*x) + 4/(b*x*sqrt(a*x + b*sqrt(x))) - 24*sqrt(a*x + b*sqrt(x))/(5*b**2*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_129():
    f = 1/(x**(sympy.S(5)/2)*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = -1024*a**4*sqrt(a*x + b*sqrt(x))/(63*b**6*sqrt(x)) + 512*a**3*sqrt(a*x + b*sqrt(x))/(63*b**5*x) - 128*a**2*sqrt(a*x + b*sqrt(x))/(21*b**4*x**(sympy.S(3)/2)) + 320*a*sqrt(a*x + b*sqrt(x))/(63*b**3*x**2) + 4/(b*x**2*sqrt(a*x + b*sqrt(x))) - 40*sqrt(a*x + b*sqrt(x))/(9*b**2*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_130():
    f = 1/(x**(sympy.S(7)/2)*(a*x + b*sqrt(x))**(sympy.S(3)/2))
    F = -8192*a**6*sqrt(a*x + b*sqrt(x))/(429*b**8*sqrt(x)) + 4096*a**5*sqrt(a*x + b*sqrt(x))/(429*b**7*x) - 1024*a**4*sqrt(a*x + b*sqrt(x))/(143*b**6*x**(sympy.S(3)/2)) + 2560*a**3*sqrt(a*x + b*sqrt(x))/(429*b**5*x**2) - 2240*a**2*sqrt(a*x + b*sqrt(x))/(429*b**4*x**(sympy.S(5)/2)) + 672*a*sqrt(a*x + b*sqrt(x))/(143*b**3*x**3) + 4/(b*x**3*sqrt(a*x + b*sqrt(x))) - 56*sqrt(a*x + b*sqrt(x))/(13*b**2*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_131():
    f = x**3*sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**4*sqrt(a*x + b*x**(sympy.S(1)/3))/9 + 4*b*x**(sympy.S(10)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(207*a) - 28*b**2*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1311*a**2) + 476*b**3*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(19665*a**3) - 6188*b**4*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(216315*a**4) + 884*b**5*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(24035*a**5) - 884*b**6*sqrt(a*x + b*x**(sympy.S(1)/3))/(14421*a**6) + 442*b**(sympy.S(27)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(14421*a**(sympy.S(25)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_132():
    f = x**2*sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**3*sqrt(a*x + b*x**(sympy.S(1)/3))/7 + 4*b*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(119*a) - 60*b**2*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*a**2) + 220*b**3*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*a**3) - 44*b**4*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(663*a**4) + 44*b**5*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(221*a**(sympy.S(9)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 44*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*a**(sympy.S(19)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 22*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*a**(sympy.S(19)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_133():
    f = x*sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/5 + 4*b*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(55*a) - 36*b**2*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*a**2) + 12*b**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*a**3) - 6*b**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*a**(sympy.S(13)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_134():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x*sqrt(a*x + b*x**(sympy.S(1)/3))/3 + 4*b*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(15*a) - 4*b**2*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*a**(sympy.S(3)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 4*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(7)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 2*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(7)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_135():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))/x
    F = 2*sqrt(a*x + b*x**(sympy.S(1)/3)) + 2*b**(sympy.S(3)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(1)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_136():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))/x**2
    F = -12*a**(sympy.S(5)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 6*a**(sympy.S(5)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 12*a**(sympy.S(3)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*b*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*b*x**(sympy.S(1)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_137():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))/x**3
    F = 10*a**(sympy.S(11)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*b**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 20*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*b**2*x**(sympy.S(2)/3)) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*b*x**(sympy.S(4)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(11*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_138():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))/x**4
    F = 308*a**(sympy.S(17)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 154*a**(sympy.S(17)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 308*a**(sympy.S(9)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(1105*b**4*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 308*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(1105*b**4*x**(sympy.S(1)/3)) - 308*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(3315*b**3*x) + 44*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(663*b**2*x**(sympy.S(5)/3)) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(221*b*x**(sympy.S(7)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(17*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_139():
    f = sqrt(a*x + b*x**(sympy.S(1)/3))/x**5
    F = -1326*a**(sympy.S(23)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(33649*b**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 2652*a**5*sqrt(a*x + b*x**(sympy.S(1)/3))/(33649*b**5*x**(sympy.S(2)/3)) + 7956*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(168245*b**4*x**(sympy.S(4)/3)) - 884*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(24035*b**3*x**2) + 68*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(2185*b**2*x**(sympy.S(8)/3)) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(437*b*x**(sympy.S(10)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(23*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_140():
    f = x**2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = 4*b*x**(sympy.S(10)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/69 + 2*x**3*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/9 + 8*b**2*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1311*a) - 136*b**3*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(19665*a**2) + 1768*b**4*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(216315*a**3) - 1768*b**5*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(168245*a**4) + 1768*b**6*sqrt(a*x + b*x**(sympy.S(1)/3))/(100947*a**5) - 884*b**(sympy.S(27)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(100947*a**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_141():
    f = x*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = 12*b*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/119 + 2*x**2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/7 + 24*b**2*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*a) - 88*b**3*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*a**2) + 88*b**4*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(3315*a**3) - 88*b**5*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(1105*a**(sympy.S(7)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 88*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*a**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 44*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*a**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_142():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = 12*b*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/55 + 2*x*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/5 + 24*b**2*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*a) - 8*b**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*a**2) + 4*b**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*a**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_143():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x
    F = 4*b*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/5 + 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/3 + 8*b**2*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*sqrt(a)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 8*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 4*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_144():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x**2
    F = 4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a*x + b*x**(sympy.S(1)/3)) + 4*a*sqrt(a*x + b*x**(sympy.S(1)/3)) - 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_145():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x**3
    F = -8*a**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 4*a**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 8*a**(sympy.S(5)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*b*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 8*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*b*x**(sympy.S(1)/3)) - 4*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*x) - 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_146():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x**4
    F = 4*a**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*b**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 8*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*b**2*x**(sympy.S(2)/3)) - 24*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*b*x**(sympy.S(4)/3)) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(55*x**2) - 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/(5*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_147():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x**5
    F = 88*a**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 44*a**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1105*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 88*a**(sympy.S(11)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(1105*b**4*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 88*a**5*sqrt(a*x + b*x**(sympy.S(1)/3))/(1105*b**4*x**(sympy.S(1)/3)) - 88*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(3315*b**3*x) + 88*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*b**2*x**(sympy.S(5)/3)) - 24*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*b*x**(sympy.S(7)/3)) - 12*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(119*x**3) - 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/(7*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_148():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/x**6
    F = -884*a**(sympy.S(27)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(100947*b**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 1768*a**6*sqrt(a*x + b*x**(sympy.S(1)/3))/(100947*b**5*x**(sympy.S(2)/3)) + 1768*a**5*sqrt(a*x + b*x**(sympy.S(1)/3))/(168245*b**4*x**(sympy.S(4)/3)) - 1768*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(216315*b**3*x**2) + 136*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(19665*b**2*x**(sympy.S(8)/3)) - 8*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(1311*b*x**(sympy.S(10)/3)) - 4*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(69*x**4) - 2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)/(9*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_149():
    f = x**4/sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(9*a) - 50*b*x**(sympy.S(10)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(207*a**2) + 350*b**2*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1311*a**3) - 1190*b**3*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(3933*a**4) + 15470*b**4*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(43263*a**5) - 2210*b**5*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(4807*a**6) + 11050*b**6*sqrt(a*x + b*x**(sympy.S(1)/3))/(14421*a**7) - 5525*b**(sympy.S(27)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(14421*a**(sympy.S(29)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_150():
    f = x**3/sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(7*a) - 38*b*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(119*a**2) + 570*b**2*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*a**3) - 2090*b**3*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*a**4) + 418*b**4*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(663*a**5) - 418*b**5*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(221*a**(sympy.S(11)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 418*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*a**(sympy.S(23)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 209*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*a**(sympy.S(23)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_151():
    f = x**2/sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*a) - 26*b*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(55*a**2) + 234*b**2*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*a**3) - 78*b**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*a**4) + 39*b**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(77*a**(sympy.S(17)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_152():
    f = x/sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(3*a) - 14*b*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(15*a**2) + 14*b**2*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*a**(sympy.S(5)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 14*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(11)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 7*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(11)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_153():
    f = 1/sqrt(a*x + b*x**(sympy.S(1)/3))
    F = 2*sqrt(a*x + b*x**(sympy.S(1)/3))/a - b**(sympy.S(3)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(5)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_154():
    f = 1/(x*sqrt(a*x + b*x**(sympy.S(1)/3)))
    F = -6*a**(sympy.S(1)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 3*a**(sympy.S(1)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 6*sqrt(a)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(b*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(b*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_155():
    f = 1/(x**2*sqrt(a*x + b*x**(sympy.S(1)/3)))
    F = 5*a**(sympy.S(7)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(7*b**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 10*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(7*b**2*x**(sympy.S(2)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(7*b*x**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_156():
    f = 1/(x**3*sqrt(a*x + b*x**(sympy.S(1)/3)))
    F = 154*a**(sympy.S(13)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 77*a**(sympy.S(13)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 154*a**(sympy.S(7)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(65*b**4*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 154*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(65*b**4*x**(sympy.S(1)/3)) - 154*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(195*b**3*x) + 22*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(39*b**2*x**(sympy.S(5)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(13*b*x**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_157():
    f = 1/(x**4*sqrt(a*x + b*x**(sympy.S(1)/3)))
    F = -663*a**(sympy.S(19)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(1463*b**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 1326*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(1463*b**5*x**(sympy.S(2)/3)) + 3978*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(7315*b**4*x**(sympy.S(4)/3)) - 442*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(1045*b**3*x**2) + 34*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(95*b**2*x**(sympy.S(8)/3)) - 6*sqrt(a*x + b*x**(sympy.S(1)/3))/(19*b*x**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_158():
    f = x**4/(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = -3*x**4/(a*sqrt(a*x + b*x**(sympy.S(1)/3))) + 23*x**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(7*a**2) - 437*b*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(119*a**3) + 6555*b**2*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*a**4) - 24035*b**3*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*a**5) + 4807*b**4*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(663*a**6) - 4807*b**5*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(221*a**(sympy.S(13)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 4807*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*a**(sympy.S(27)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 4807*b**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(442*a**(sympy.S(27)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_159():
    f = x**3/(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = -3*x**3/(a*sqrt(a*x + b*x**(sympy.S(1)/3))) + 17*x**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*a**2) - 221*b*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(55*a**3) + 1989*b**2*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*a**4) - 663*b**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*a**5) + 663*b**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(154*a**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_160():
    f = x**2/(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = -3*x**2/(a*sqrt(a*x + b*x**(sympy.S(1)/3))) + 11*x*sqrt(a*x + b*x**(sympy.S(1)/3))/(3*a**2) - 77*b*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))/(15*a**3) + 77*b**2*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*a**(sympy.S(7)/2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 77*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 77*b**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*a**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_161():
    f = x/(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2)
    F = -3*x/(a*sqrt(a*x + b*x**(sympy.S(1)/3))) + 5*sqrt(a*x + b*x**(sympy.S(1)/3))/a**2 - 5*b**(sympy.S(3)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_162():
    f = (a*x + b*x**(sympy.S(1)/3))**(sympy.S(-3)/2)
    F = 3*x**(sympy.S(2)/3)/(b*sqrt(a*x + b*x**(sympy.S(1)/3))) - x**(sympy.S(1)/3)*(3*a*x**(sympy.S(2)/3) + 3*b)/(sqrt(a)*b*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(3*sqrt(a)*x**(sympy.S(1)/3) + 3*sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(3*sqrt(a)*x**(sympy.S(1)/3) + 3*sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a*x + b*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_163():
    f = 1/(x*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2))
    F = -5*a**(sympy.S(3)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(9)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 3/(b*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 5*sqrt(a*x + b*x**(sympy.S(1)/3))/(b**2*x**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_164():
    f = 1/(x**2*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2))
    F = -77*a**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 77*a**(sympy.S(9)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(10*b**(sympy.S(15)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 77*a**(sympy.S(5)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(5*b**4*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) - 77*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*b**4*x**(sympy.S(1)/3)) + 77*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(15*b**3*x) + 3/(b*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 11*sqrt(a*x + b*x**(sympy.S(1)/3))/(3*b**2*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_165():
    f = 1/(x**3*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2))
    F = 663*a**(sympy.S(15)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(154*b**(sympy.S(21)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) + 663*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(77*b**5*x**(sympy.S(2)/3)) - 1989*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(385*b**4*x**(sympy.S(4)/3)) + 221*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(55*b**3*x**2) + 3/(b*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 17*sqrt(a*x + b*x**(sympy.S(1)/3))/(5*b**2*x**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_166():
    f = 1/(x**4*(a*x + b*x**(sympy.S(1)/3))**(sympy.S(3)/2))
    F = 4807*a**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_e(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(221*b**(sympy.S(27)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 4807*a**(sympy.S(21)/4)*x**(sympy.S(1)/6)*sqrt((a*x**(sympy.S(2)/3) + b)/(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))**2)*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*elliptic_f(2*atan(a**(sympy.S(1)/4)*x**(sympy.S(1)/6)/b**(sympy.S(1)/4)), sympy.S.Half)/(442*b**(sympy.S(27)/4)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 4807*a**(sympy.S(11)/2)*x**(sympy.S(1)/3)*(a*x**(sympy.S(2)/3) + b)/(221*b**7*(sqrt(a)*x**(sympy.S(1)/3) + sqrt(b))*sqrt(a*x + b*x**(sympy.S(1)/3))) + 4807*a**5*sqrt(a*x + b*x**(sympy.S(1)/3))/(221*b**7*x**(sympy.S(1)/3)) - 4807*a**4*sqrt(a*x + b*x**(sympy.S(1)/3))/(663*b**6*x) + 24035*a**3*sqrt(a*x + b*x**(sympy.S(1)/3))/(4641*b**5*x**(sympy.S(5)/3)) - 6555*a**2*sqrt(a*x + b*x**(sympy.S(1)/3))/(1547*b**4*x**(sympy.S(7)/3)) + 437*a*sqrt(a*x + b*x**(sympy.S(1)/3))/(119*b**3*x**3) + 3/(b*x**(sympy.S(10)/3)*sqrt(a*x + b*x**(sympy.S(1)/3))) - 23*sqrt(a*x + b*x**(sympy.S(1)/3))/(7*b**2*x**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_167():
    f = x**3*sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(9*a) - 16*b*x**(sympy.S(8)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(75*a**2) + 352*b**2*x**(sympy.S(7)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(1725*a**3) - 1408*b**3*x**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(7245*a**4) + 2816*b**4*x**(sympy.S(5)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(15295*a**5) - 45056*b**5*x**(sympy.S(4)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(260015*a**6) + 90112*b**6*x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(557175*a**7) - 360448*b**7*x**(sympy.S(2)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(2414425*a**8) + 65536*b**8*x**(sympy.S(1)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(482885*a**9) - 524288*b**9*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(4345965*a**10) + 1048576*b**10*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(10140585*a**11*x**(sympy.S(1)/3)) - 4194304*b**11*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(50702925*a**12*x**(sympy.S(2)/3)) + 8388608*b**12*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(152108775*a**13*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_168():
    f = x**2*sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(7*a) - 36*b*x**(sympy.S(5)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(133*a**2) + 576*b**2*x**(sympy.S(4)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(2261*a**3) - 384*b**3*x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(1615*a**4) + 4608*b**4*x**(sympy.S(2)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(20995*a**5) - 9216*b**5*x**(sympy.S(1)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(46189*a**6) + 8192*b**6*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(46189*a**7) - 49152*b**7*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(323323*a**8*x**(sympy.S(1)/3)) + 196608*b**8*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(1616615*a**9*x**(sympy.S(2)/3)) - 131072*b**9*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(1616615*a**10*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_169():
    f = x*sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(5*a) - 24*b*x**(sympy.S(2)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(65*a**2) + 48*b**2*x**(sympy.S(1)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(143*a**3) - 128*b**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(429*a**4) + 256*b**4*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(1001*a**5*x**(sympy.S(1)/3)) - 1024*b**5*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(5005*a**6*x**(sympy.S(2)/3)) + 2048*b**6*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(15015*a**7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_170():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(3*a) - 4*b*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(7*a**2*x**(sympy.S(1)/3)) + 16*b**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(35*a**3*x**(sympy.S(2)/3)) - 32*b**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(105*a**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_171():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))/x
    F = 2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_172():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))/x**2
    F = 3*a**2*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(4*b**(sympy.S(3)/2)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(4*b*x**(sympy.S(2)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_173():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))/x**3
    F = -21*a**5*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(128*b**(sympy.S(9)/2)) + 21*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(128*b**4*x**(sympy.S(2)/3)) - 7*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(64*b**3*x) + 7*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(80*b**2*x**(sympy.S(4)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(40*b*x**(sympy.S(5)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(5*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_174():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))/x**4
    F = 1287*a**8*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(16384*b**(sympy.S(15)/2)) - 1287*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(16384*b**7*x**(sympy.S(2)/3)) + 429*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(8192*b**6*x) - 429*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(10240*b**5*x**(sympy.S(4)/3)) + 1287*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(35840*b**4*x**(sympy.S(5)/3)) - 143*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(4480*b**3*x**2) + 13*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(448*b**2*x**(sympy.S(7)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(112*b*x**(sympy.S(8)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(8*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_175():
    f = sqrt(a*x + b*x**(sympy.S(2)/3))/x**5
    F = -12597*a**11*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(262144*b**(sympy.S(21)/2)) + 12597*a**10*sqrt(a*x + b*x**(sympy.S(2)/3))/(262144*b**10*x**(sympy.S(2)/3)) - 4199*a**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(131072*b**9*x) + 4199*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(163840*b**8*x**(sympy.S(4)/3)) - 12597*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(573440*b**7*x**(sympy.S(5)/3)) + 4199*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(215040*b**6*x**2) - 4199*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(236544*b**5*x**(sympy.S(7)/3)) + 323*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(19712*b**4*x**(sympy.S(8)/3)) - 323*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(21120*b**3*x**3) + 19*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(1320*b**2*x**(sympy.S(10)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(220*b*x**(sympy.S(11)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(11*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_176():
    f = x**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = 2*x**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(9*a) - 44*b*x**(sympy.S(5)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(225*a**2) + 176*b**2*x**(sympy.S(4)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(1035*a**3) - 352*b**3*x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(2415*a**4) + 5632*b**4*x**(sympy.S(2)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(45885*a**5) - 11264*b**5*x**(sympy.S(1)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(111435*a**6) + 45056*b**6*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(557175*a**7) - 90112*b**7*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(1448655*a**8*x**(sympy.S(1)/3)) + 65536*b**8*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(1448655*a**9*x**(sympy.S(2)/3)) - 131072*b**9*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(4345965*a**10*x) + 524288*b**10*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(30421755*a**11*x**(sympy.S(4)/3)) - 1048576*b**11*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(152108775*a**12*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_177():
    f = x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = 2*x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(7*a) - 32*b*x**(sympy.S(2)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(133*a**2) + 64*b**2*x**(sympy.S(1)/3)*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(323*a**3) - 256*b**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(1615*a**4) + 512*b**4*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(4199*a**5*x**(sympy.S(1)/3)) - 4096*b**5*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(46189*a**6*x**(sympy.S(2)/3)) + 8192*b**6*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(138567*a**7*x) - 32768*b**7*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(969969*a**8*x**(sympy.S(4)/3)) + 65536*b**8*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(4849845*a**9*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_178():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = 2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(5*a) - 4*b*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(13*a**2*x**(sympy.S(1)/3)) + 32*b**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(143*a**3*x**(sympy.S(2)/3)) - 64*b**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(429*a**4*x) + 256*b**4*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(3003*a**5*x**(sympy.S(4)/3)) - 512*b**5*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(15015*a**6*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_179():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x
    F = 2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(3*a*x) - 8*b*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(21*a**2*x**(sympy.S(4)/3)) + 16*b**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(5)/2)/(105*a**3*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_180():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**2
    F = -6*b**(sympy.S(3)/2)*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3))) + 6*b*sqrt(a*x + b*x**(sympy.S(2)/3))/x**(sympy.S(1)/3) + 2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_181():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**3
    F = 3*a**3*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(8*b**(sympy.S(3)/2)) - 3*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(8*b*x**(sympy.S(2)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(4*x) - (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_182():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**4
    F = -21*a**6*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(512*b**(sympy.S(9)/2)) + 21*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(512*b**4*x**(sympy.S(2)/3)) - 7*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(256*b**3*x) + 7*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(320*b**2*x**(sympy.S(4)/3)) - 3*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(160*b*x**(sympy.S(5)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(20*x**2) - (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_183():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**5
    F = 429*a**9*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(32768*b**(sympy.S(15)/2)) - 429*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(32768*b**7*x**(sympy.S(2)/3)) + 143*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(16384*b**6*x) - 143*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(20480*b**5*x**(sympy.S(4)/3)) + 429*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(71680*b**4*x**(sympy.S(5)/3)) - 143*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(26880*b**3*x**2) + 13*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(2688*b**2*x**(sympy.S(7)/3)) - a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(224*b*x**(sympy.S(8)/3)) - a*sqrt(a*x + b*x**(sympy.S(2)/3))/(16*x**3) - (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(3*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_184():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/x**6
    F = -12597*a**12*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(2097152*b**(sympy.S(21)/2)) + 12597*a**11*sqrt(a*x + b*x**(sympy.S(2)/3))/(2097152*b**10*x**(sympy.S(2)/3)) - 4199*a**10*sqrt(a*x + b*x**(sympy.S(2)/3))/(1048576*b**9*x) + 4199*a**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(1310720*b**8*x**(sympy.S(4)/3)) - 12597*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(4587520*b**7*x**(sympy.S(5)/3)) + 4199*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(1720320*b**6*x**2) - 4199*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(1892352*b**5*x**(sympy.S(7)/3)) + 323*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(157696*b**4*x**(sympy.S(8)/3)) - 323*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(168960*b**3*x**3) + 19*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(10560*b**2*x**(sympy.S(10)/3)) - 3*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(1760*b*x**(sympy.S(11)/3)) - 3*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(88*x**4) - (a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)/(4*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_185():
    f = x**4/sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(9*a) - 52*b*x**(sympy.S(11)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(225*a**2) + 416*b**2*x**(sympy.S(10)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(1725*a**3) - 9152*b**3*x**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(36225*a**4) + 36608*b**4*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(137655*a**5) - 73216*b**5*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(260015*a**6) + 1171456*b**6*x**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(3900225*a**7) - 180224*b**7*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(557175*a**8) + 65536*b**8*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(185725*a**9) - 131072*b**9*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(334305*a**10) + 1048576*b**10*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(2340135*a**11) - 2097152*b**11*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(3900225*a**12) + 8388608*b**12*sqrt(a*x + b*x**(sympy.S(2)/3))/(11700675*a**13) - 16777216*b**13*sqrt(a*x + b*x**(sympy.S(2)/3))/(11700675*a**14*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_186():
    f = x**3/sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(7*a) - 40*b*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(133*a**2) + 720*b**2*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(2261*a**3) - 768*b**3*x**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(2261*a**4) + 1536*b**4*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(4199*a**5) - 18432*b**5*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(46189*a**6) + 20480*b**6*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(46189*a**7) - 163840*b**7*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(323323*a**8) + 196608*b**8*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(323323*a**9) - 262144*b**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(323323*a**10) + 524288*b**10*sqrt(a*x + b*x**(sympy.S(2)/3))/(323323*a**11*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_187():
    f = x**2/sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(5*a) - 28*b*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(65*a**2) + 336*b**2*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(715*a**3) - 224*b**3*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(429*a**4) + 256*b**4*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(429*a**5) - 512*b**5*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(715*a**6) + 2048*b**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(2145*a**7) - 4096*b**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(2145*a**8*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_188():
    f = x/sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(3*a) - 16*b*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(21*a**2) + 32*b**2*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(35*a**3) - 128*b**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(105*a**4) + 256*b**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(105*a**5*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_189():
    f = 1/sqrt(a*x + b*x**(sympy.S(2)/3))
    F = 2*sqrt(a*x + b*x**(sympy.S(2)/3))/a - 4*b*sqrt(a*x + b*x**(sympy.S(2)/3))/(a**2*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_190():
    f = 1/(x*sqrt(a*x + b*x**(sympy.S(2)/3)))
    F = 3*a*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/b**(sympy.S(3)/2) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(b*x**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_191():
    f = 1/(x**2*sqrt(a*x + b*x**(sympy.S(2)/3)))
    F = -105*a**4*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(64*b**(sympy.S(9)/2)) + 105*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(64*b**4*x**(sympy.S(2)/3)) - 35*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(32*b**3*x) + 7*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(8*b**2*x**(sympy.S(4)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(4*b*x**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_192():
    f = 1/(x**3*sqrt(a*x + b*x**(sympy.S(2)/3)))
    F = 1287*a**7*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(1024*b**(sympy.S(15)/2)) - 1287*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(1024*b**7*x**(sympy.S(2)/3)) + 429*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(512*b**6*x) - 429*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(640*b**5*x**(sympy.S(4)/3)) + 1287*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(2240*b**4*x**(sympy.S(5)/3)) - 143*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(280*b**3*x**2) + 13*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(28*b**2*x**(sympy.S(7)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(7*b*x**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_193():
    f = 1/(x**4*sqrt(a*x + b*x**(sympy.S(2)/3)))
    F = -138567*a**10*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(131072*b**(sympy.S(21)/2)) + 138567*a**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(131072*b**10*x**(sympy.S(2)/3)) - 46189*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(65536*b**9*x) + 46189*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(81920*b**8*x**(sympy.S(4)/3)) - 138567*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(286720*b**7*x**(sympy.S(5)/3)) + 46189*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(107520*b**6*x**2) - 4199*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(10752*b**5*x**(sympy.S(7)/3)) + 323*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(896*b**4*x**(sympy.S(8)/3)) - 323*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(960*b**3*x**3) + 19*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(60*b**2*x**(sympy.S(10)/3)) - 3*sqrt(a*x + b*x**(sympy.S(2)/3))/(10*b*x**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_194():
    f = x**4/(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = -6*x**4/(a*sqrt(a*x + b*x**(sympy.S(2)/3))) + 44*x**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(7*a**2) - 880*b*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(133*a**3) + 15840*b**2*x**(sympy.S(7)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(2261*a**4) - 16896*b**3*x**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(2261*a**5) + 33792*b**4*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(4199*a**6) - 36864*b**5*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(4199*a**7) + 40960*b**6*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(4199*a**8) - 327680*b**7*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(29393*a**9) + 393216*b**8*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(29393*a**10) - 524288*b**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(29393*a**11) + 1048576*b**10*sqrt(a*x + b*x**(sympy.S(2)/3))/(29393*a**12*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_195():
    f = x**3/(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = -6*x**3/(a*sqrt(a*x + b*x**(sympy.S(2)/3))) + 32*x**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(5*a**2) - 448*b*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(65*a**3) + 5376*b**2*x**(sympy.S(4)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(715*a**4) - 3584*b**3*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(429*a**5) + 4096*b**4*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(429*a**6) - 8192*b**5*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(715*a**7) + 32768*b**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(2145*a**8) - 65536*b**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(2145*a**9*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_196():
    f = x**2/(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = -6*x**2/(a*sqrt(a*x + b*x**(sympy.S(2)/3))) + 20*x*sqrt(a*x + b*x**(sympy.S(2)/3))/(3*a**2) - 160*b*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(21*a**3) + 64*b**2*x**(sympy.S(1)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))/(7*a**4) - 256*b**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(21*a**5) + 512*b**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(21*a**6*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_197():
    f = x/(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = -6*x/(a*sqrt(a*x + b*x**(sympy.S(2)/3))) + 8*sqrt(a*x + b*x**(sympy.S(2)/3))/a**2 - 16*b*sqrt(a*x + b*x**(sympy.S(2)/3))/(a**3*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_198():
    f = (a*x + b*x**(sympy.S(2)/3))**(sympy.S(-3)/2)
    F = 6*x**(sympy.S(1)/3)/(b*sqrt(a*x + b*x**(sympy.S(2)/3))) - 6*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_199():
    f = 1/(x*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2))
    F = 105*a**3*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(8*b**(sympy.S(9)/2)) - 105*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(8*b**4*x**(sympy.S(2)/3)) + 35*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(4*b**3*x) + 6/(b*x**(sympy.S(2)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))) - 7*sqrt(a*x + b*x**(sympy.S(2)/3))/(b**2*x**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_200():
    f = 1/(x**2*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2))
    F = -9009*a**6*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(512*b**(sympy.S(15)/2)) + 9009*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(512*b**7*x**(sympy.S(2)/3)) - 3003*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(256*b**6*x) + 3003*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(320*b**5*x**(sympy.S(4)/3)) - 1287*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(160*b**4*x**(sympy.S(5)/3)) + 143*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(20*b**3*x**2) + 6/(b*x**(sympy.S(5)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))) - 13*sqrt(a*x + b*x**(sympy.S(2)/3))/(2*b**2*x**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_201():
    f = 1/(x**3*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2))
    F = 692835*a**9*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(32768*b**(sympy.S(21)/2)) - 692835*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(32768*b**10*x**(sympy.S(2)/3)) + 230945*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(16384*b**9*x) - 46189*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(4096*b**8*x**(sympy.S(4)/3)) + 138567*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(14336*b**7*x**(sympy.S(5)/3)) - 46189*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(5376*b**6*x**2) + 20995*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(2688*b**5*x**(sympy.S(7)/3)) - 1615*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(224*b**4*x**(sympy.S(8)/3)) + 323*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(48*b**3*x**3) + 6/(b*x**(sympy.S(8)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))) - 19*sqrt(a*x + b*x**(sympy.S(2)/3))/(3*b**2*x**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_202():
    f = 1/(x**4*(a*x + b*x**(sympy.S(2)/3))**(sympy.S(3)/2))
    F = -50702925*a**12*atanh(sqrt(b)*x**(sympy.S(1)/3)/sqrt(a*x + b*x**(sympy.S(2)/3)))/(2097152*b**(sympy.S(27)/2)) + 50702925*a**11*sqrt(a*x + b*x**(sympy.S(2)/3))/(2097152*b**13*x**(sympy.S(2)/3)) - 16900975*a**10*sqrt(a*x + b*x**(sympy.S(2)/3))/(1048576*b**12*x) + 3380195*a**9*sqrt(a*x + b*x**(sympy.S(2)/3))/(262144*b**11*x**(sympy.S(4)/3)) - 1448655*a**8*sqrt(a*x + b*x**(sympy.S(2)/3))/(131072*b**10*x**(sympy.S(5)/3)) + 482885*a**7*sqrt(a*x + b*x**(sympy.S(2)/3))/(49152*b**9*x**2) - 2414425*a**6*sqrt(a*x + b*x**(sympy.S(2)/3))/(270336*b**8*x**(sympy.S(7)/3)) + 185725*a**5*sqrt(a*x + b*x**(sympy.S(2)/3))/(22528*b**7*x**(sympy.S(8)/3)) - 260015*a**4*sqrt(a*x + b*x**(sympy.S(2)/3))/(33792*b**6*x**3) + 15295*a**3*sqrt(a*x + b*x**(sympy.S(2)/3))/(2112*b**5*x**(sympy.S(10)/3)) - 2415*a**2*sqrt(a*x + b*x**(sympy.S(2)/3))/(352*b**4*x**(sympy.S(11)/3)) + 575*a*sqrt(a*x + b*x**(sympy.S(2)/3))/(88*b**3*x**4) + 6/(b*x**(sympy.S(11)/3)*sqrt(a*x + b*x**(sympy.S(2)/3))) - 25*sqrt(a*x + b*x**(sympy.S(2)/3))/(4*b**2*x**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_203():
    f = x**2*(a*x**2 + b*x**3)
    F = a*x**5/5 + b*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_204():
    f = x*(a*x**2 + b*x**3)
    F = a*x**4/4 + b*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_205():
    f = a*x**2 + b*x**3
    F = a*x**3/3 + b*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_206():
    f = (a*x**2 + b*x**3)/x
    F = a*x**2/2 + b*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_207():
    f = (a*x**2 + b*x**3)/x**2
    F = a*x + b*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_208():
    f = x**2*(a*x**2 + b*x**3)**2
    F = a**2*x**7/7 + a*b*x**8/4 + b**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_209():
    f = x*(a*x**2 + b*x**3)**2
    F = a**2*x**6/6 + 2*a*b*x**7/7 + b**2*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_210():
    f = (a*x**2 + b*x**3)**2
    F = a**2*x**5/5 + a*b*x**6/3 + b**2*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_211():
    f = (a*x**2 + b*x**3)**2/x
    F = a**2*x**4/4 + 2*a*b*x**5/5 + b**2*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_212():
    f = (a*x**2 + b*x**3)**2/x**2
    F = a**2*x**3/3 + a*b*x**4/2 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_213():
    f = x**6/(a*x**2 + b*x**3)
    F = a**4*log(a + b*x)/b**5 - a**3*x/b**4 + a**2*x**2/(2*b**3) - a*x**3/(3*b**2) + x**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_214():
    f = x**5/(a*x**2 + b*x**3)
    F = -a**3*log(a + b*x)/b**4 + a**2*x/b**3 - a*x**2/(2*b**2) + x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_215():
    f = x**4/(a*x**2 + b*x**3)
    F = a**2*log(a + b*x)/b**3 - a*x/b**2 + x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_216():
    f = x**3/(a*x**2 + b*x**3)
    F = -a*log(a + b*x)/b**2 + x/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_217():
    f = x**2/(a*x**2 + b*x**3)
    F = log(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_218():
    f = x/(a*x**2 + b*x**3)
    F = log(x)/a - log(a + b*x)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_219():
    f = 1/(a*x**2 + b*x**3)
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_220():
    f = 1/(x*(a*x**2 + b*x**3))
    F = -1/(2*a*x**2) + b/(a**2*x) + b**2*log(x)/a**3 - b**2*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_221():
    f = 1/(x**2*(a*x**2 + b*x**3))
    F = -1/(3*a*x**3) + b/(2*a**2*x**2) - b**2/(a**3*x) - b**3*log(x)/a**4 + b**3*log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_222():
    f = x**8/(a*x**2 + b*x**3)**2
    F = -a**4/(b**5*(a + b*x)) - 4*a**3*log(a + b*x)/b**5 + 3*a**2*x/b**4 - a*x**2/b**3 + x**3/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_223():
    f = x**7/(a*x**2 + b*x**3)**2
    F = a**3/(b**4*(a + b*x)) + 3*a**2*log(a + b*x)/b**4 - 2*a*x/b**3 + x**2/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_224():
    f = x**6/(a*x**2 + b*x**3)**2
    F = -a**2/(b**3*(a + b*x)) - 2*a*log(a + b*x)/b**3 + x/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_225():
    f = x**5/(a*x**2 + b*x**3)**2
    F = a/(b**2*(a + b*x)) + log(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_226():
    f = x**4/(a*x**2 + b*x**3)**2
    F = -1/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_227():
    f = x**3/(a*x**2 + b*x**3)**2
    F = 1/(a*(a + b*x)) + log(x)/a**2 - log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_228():
    f = x**2/(a*x**2 + b*x**3)**2
    F = -b/(a**2*(a + b*x)) - 1/(a**2*x) - 2*b*log(x)/a**3 + 2*b*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_229():
    f = x/(a*x**2 + b*x**3)**2
    F = -1/(2*a**2*x**2) + b**2/(a**3*(a + b*x)) + 2*b/(a**3*x) + 3*b**2*log(x)/a**4 - 3*b**2*log(a + b*x)/a**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_230():
    f = (a*x**2 + b*x**3)**(-2)
    F = -1/(3*a**2*x**3) + b/(a**3*x**2) - b**3/(a**4*(a + b*x)) - 3*b**2/(a**4*x) - 4*b**3*log(x)/a**5 + 4*b**3*log(a + b*x)/a**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_231():
    f = 1/(x*(a*x**2 + b*x**3)**2)
    F = -1/(4*a**2*x**4) + 2*b/(3*a**3*x**3) - 3*b**2/(2*a**4*x**2) + b**4/(a**5*(a + b*x)) + 4*b**3/(a**5*x) + 5*b**4*log(x)/a**6 - 5*b**4*log(a + b*x)/a**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_232():
    f = x**2*sqrt(a*x**2 + b*x**3)
    F = -32*a**3*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(315*b**4*x**3) + 16*a**2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(105*b**3*x**2) - 4*a*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(21*b**2*x) + 2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_233():
    f = x*sqrt(a*x**2 + b*x**3)
    F = 16*a**2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(105*b**3*x**3) - 8*a*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(35*b**2*x**2) + 2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(7*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_234():
    f = sqrt(a*x**2 + b*x**3)
    F = -4*a*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(15*b**2*x**3) + 2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(5*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_235():
    f = sqrt(a*x**2 + b*x**3)/x
    F = 2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(3*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_236():
    f = sqrt(a*x**2 + b*x**3)/x**2
    F = -2*sqrt(a)*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3)) + 2*sqrt(a*x**2 + b*x**3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_237():
    f = sqrt(a*x**2 + b*x**3)/x**3
    F = -sqrt(a*x**2 + b*x**3)/x**2 - b*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_238():
    f = sqrt(a*x**2 + b*x**3)/x**4
    F = -sqrt(a*x**2 + b*x**3)/(2*x**3) - b*sqrt(a*x**2 + b*x**3)/(4*a*x**2) + b**2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_239():
    f = sqrt(a*x**2 + b*x**3)/x**5
    F = -sqrt(a*x**2 + b*x**3)/(3*x**4) - b*sqrt(a*x**2 + b*x**3)/(12*a*x**3) + b**2*sqrt(a*x**2 + b*x**3)/(8*a**2*x**2) - b**3*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_240():
    f = x**2*(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = -512*a**5*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(45045*b**6*x**5) + 256*a**4*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(9009*b**5*x**4) - 64*a**3*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(1287*b**4*x**3) + 32*a**2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(429*b**3*x**2) - 4*a*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(39*b**2*x) + 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_241():
    f = x*(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = 256*a**4*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(15015*b**5*x**5) - 128*a**3*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(3003*b**4*x**4) + 32*a**2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(429*b**3*x**3) - 16*a*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(143*b**2*x**2) + 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(13*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_242():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = -32*a**3*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(1155*b**4*x**5) + 16*a**2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(231*b**3*x**4) - 4*a*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(33*b**2*x**3) + 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(11*b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_243():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x
    F = 16*a**2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(315*b**3*x**5) - 8*a*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(63*b**2*x**4) + 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_244():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**2
    F = -4*a*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(35*b**2*x**5) + 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(7*b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_245():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**3
    F = 2*(a*x**2 + b*x**3)**(sympy.S(5)/2)/(5*b*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_246():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**4
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3)) + 2*a*sqrt(a*x**2 + b*x**3)/x + 2*(a*x**2 + b*x**3)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_247():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**5
    F = -3*sqrt(a)*b*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3)) + 3*b*sqrt(a*x**2 + b*x**3)/x - (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_248():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**6
    F = -3*b*sqrt(a*x**2 + b*x**3)/(4*x**2) - (a*x**2 + b*x**3)**(sympy.S(3)/2)/(2*x**5) - 3*b**2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_249():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**7
    F = -b*sqrt(a*x**2 + b*x**3)/(4*x**3) - (a*x**2 + b*x**3)**(sympy.S(3)/2)/(3*x**6) - b**2*sqrt(a*x**2 + b*x**3)/(8*a*x**2) + b**3*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_250():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**8
    F = -b*sqrt(a*x**2 + b*x**3)/(8*x**4) - (a*x**2 + b*x**3)**(sympy.S(3)/2)/(4*x**7) - b**2*sqrt(a*x**2 + b*x**3)/(32*a*x**3) + 3*b**3*sqrt(a*x**2 + b*x**3)/(64*a**2*x**2) - 3*b**4*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(64*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_251():
    f = (a*x**2 + b*x**3)**(sympy.S(3)/2)/x**9
    F = -3*b*sqrt(a*x**2 + b*x**3)/(40*x**5) - (a*x**2 + b*x**3)**(sympy.S(3)/2)/(5*x**8) - b**2*sqrt(a*x**2 + b*x**3)/(80*a*x**4) + b**3*sqrt(a*x**2 + b*x**3)/(64*a**2*x**3) - 3*b**4*sqrt(a*x**2 + b*x**3)/(128*a**3*x**2) + 3*b**5*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(128*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_252():
    f = x**4/sqrt(a*x**2 + b*x**3)
    F = -32*a**3*sqrt(a*x**2 + b*x**3)/(35*b**4*x) + 16*a**2*sqrt(a*x**2 + b*x**3)/(35*b**3) - 12*a*x*sqrt(a*x**2 + b*x**3)/(35*b**2) + 2*x**2*sqrt(a*x**2 + b*x**3)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_253():
    f = x**3/sqrt(a*x**2 + b*x**3)
    F = 16*a**2*sqrt(a*x**2 + b*x**3)/(15*b**3*x) - 8*a*sqrt(a*x**2 + b*x**3)/(15*b**2) + 2*x*sqrt(a*x**2 + b*x**3)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_254():
    f = x**2/sqrt(a*x**2 + b*x**3)
    F = -4*a*sqrt(a*x**2 + b*x**3)/(3*b**2*x) + 2*sqrt(a*x**2 + b*x**3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_255():
    f = x/sqrt(a*x**2 + b*x**3)
    F = 2*sqrt(a*x**2 + b*x**3)/(b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_256():
    f = 1/sqrt(a*x**2 + b*x**3)
    F = -2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_257():
    f = 1/(x*sqrt(a*x**2 + b*x**3))
    F = -sqrt(a*x**2 + b*x**3)/(a*x**2) + b*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_258():
    f = 1/(x**2*sqrt(a*x**2 + b*x**3))
    F = -sqrt(a*x**2 + b*x**3)/(2*a*x**3) + 3*b*sqrt(a*x**2 + b*x**3)/(4*a**2*x**2) - 3*b**2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(4*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_259():
    f = 1/(x**3*sqrt(a*x**2 + b*x**3))
    F = -sqrt(a*x**2 + b*x**3)/(3*a*x**4) + 5*b*sqrt(a*x**2 + b*x**3)/(12*a**2*x**3) - 5*b**2*sqrt(a*x**2 + b*x**3)/(8*a**3*x**2) + 5*b**3*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_260():
    f = x**6/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = 32*a**2*sqrt(a*x**2 + b*x**3)/(5*b**4*x) - 16*a*sqrt(a*x**2 + b*x**3)/(5*b**3) - 2*x**4/(b*sqrt(a*x**2 + b*x**3)) + 12*x*sqrt(a*x**2 + b*x**3)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_261():
    f = x**5/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = -16*a*sqrt(a*x**2 + b*x**3)/(3*b**3*x) - 2*x**3/(b*sqrt(a*x**2 + b*x**3)) + 8*sqrt(a*x**2 + b*x**3)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_262():
    f = x**4/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = -2*x**2/(b*sqrt(a*x**2 + b*x**3)) + 4*sqrt(a*x**2 + b*x**3)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_263():
    f = x**3/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = -2*x/(b*sqrt(a*x**2 + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_264():
    f = x**2/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = 2*x/(a*sqrt(a*x**2 + b*x**3)) - 2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_265():
    f = x/(a*x**2 + b*x**3)**(sympy.S(3)/2)
    F = 2/(a*sqrt(a*x**2 + b*x**3)) - 3*sqrt(a*x**2 + b*x**3)/(a**2*x**2) + 3*b*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_266():
    f = (a*x**2 + b*x**3)**(sympy.S(-3)/2)
    F = 2/(a*x*sqrt(a*x**2 + b*x**3)) - 5*sqrt(a*x**2 + b*x**3)/(2*a**2*x**3) + 15*b*sqrt(a*x**2 + b*x**3)/(4*a**3*x**2) - 15*b**2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(4*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_267():
    f = 1/(x*(a*x**2 + b*x**3)**(sympy.S(3)/2))
    F = 2/(a*x**2*sqrt(a*x**2 + b*x**3)) - 7*sqrt(a*x**2 + b*x**3)/(3*a**2*x**4) + 35*b*sqrt(a*x**2 + b*x**3)/(12*a**3*x**3) - 35*b**2*sqrt(a*x**2 + b*x**3)/(8*a**4*x**2) + 35*b**3*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(8*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_268():
    f = 1/(x**2*(a*x**2 + b*x**3)**(sympy.S(3)/2))
    F = 2/(a*x**3*sqrt(a*x**2 + b*x**3)) - 9*sqrt(a*x**2 + b*x**3)/(4*a**2*x**5) + 21*b*sqrt(a*x**2 + b*x**3)/(8*a**3*x**4) - 105*b**2*sqrt(a*x**2 + b*x**3)/(32*a**4*x**3) + 315*b**3*sqrt(a*x**2 + b*x**3)/(64*a**5*x**2) - 315*b**4*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**3))/(64*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_269():
    f = x**(sympy.S(7)/2)/sqrt(a*x**2 + b*x**3)
    F = -5*a**3*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**3))/(8*b**(sympy.S(7)/2)) + 5*a**2*sqrt(a*x**2 + b*x**3)/(8*b**3*sqrt(x)) - 5*a*sqrt(x)*sqrt(a*x**2 + b*x**3)/(12*b**2) + x**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_270():
    f = x**(sympy.S(5)/2)/sqrt(a*x**2 + b*x**3)
    F = 3*a**2*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**3))/(4*b**(sympy.S(5)/2)) - 3*a*sqrt(a*x**2 + b*x**3)/(4*b**2*sqrt(x)) + sqrt(x)*sqrt(a*x**2 + b*x**3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_271():
    f = x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**3)
    F = -a*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**3))/b**(sympy.S(3)/2) + sqrt(a*x**2 + b*x**3)/(b*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_272():
    f = sqrt(x)/sqrt(a*x**2 + b*x**3)
    F = 2*atanh(sqrt(b)*x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**3))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_273():
    f = 1/(sqrt(x)*sqrt(a*x**2 + b*x**3))
    F = -2*sqrt(a*x**2 + b*x**3)/(a*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_274():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3))
    F = -2*sqrt(a*x**2 + b*x**3)/(3*a*x**(sympy.S(5)/2)) + 4*b*sqrt(a*x**2 + b*x**3)/(3*a**2*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_275():
    f = 1/(x**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**3))
    F = -2*sqrt(a*x**2 + b*x**3)/(5*a*x**(sympy.S(7)/2)) + 8*b*sqrt(a*x**2 + b*x**3)/(15*a**2*x**(sympy.S(5)/2)) - 16*b**2*sqrt(a*x**2 + b*x**3)/(15*a**3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_276():
    f = 1/(x**(sympy.S(7)/2)*sqrt(a*x**2 + b*x**3))
    F = -2*sqrt(a*x**2 + b*x**3)/(7*a*x**(sympy.S(9)/2)) + 12*b*sqrt(a*x**2 + b*x**3)/(35*a**2*x**(sympy.S(7)/2)) - 16*b**2*sqrt(a*x**2 + b*x**3)/(35*a**3*x**(sympy.S(5)/2)) + 32*b**3*sqrt(a*x**2 + b*x**3)/(35*a**4*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_277():
    f = x**(1 - 3*n)*(a*x**2 + b*x**3)**n
    F = x**(2 - 3*n)*(a*x**2 + b*x**3)**n*hyper((-n, 2 - n), (3 - n,), -b*x/a)/((1 + b*x/a)**n*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_278():
    f = x**(-3*n - 1)*(a*x**2 + b*x**3)**n
    F = -(a*x**2 + b*x**3)**n*hyper((-n, -n), (1 - n,), -b*x/a)/(n*x**(3*n)*(1 + b*x/a)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_279():
    f = x**(-3*n - 2)*(a*x**2 + b*x**3)**n
    F = -x**(-3*n - 3)*(a*x**2 + b*x**3)**(n + 1)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_280():
    f = x**(-3*n - 3)*(a*x**2 + b*x**3)**n
    F = -x**(-3*n - 4)*(a*x**2 + b*x**3)**(n + 1)/(a*(n + 2)) + b*x**(-3*n - 3)*(a*x**2 + b*x**3)**(n + 1)/(a**2*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_281():
    f = x**(-3*n - 4)*(a*x**2 + b*x**3)**n
    F = -x**(-3*n - 5)*(a*x**2 + b*x**3)**(n + 1)/(a*(n + 3)) + 2*b*x**(-3*n - 4)*(a*x**2 + b*x**3)**(n + 1)/(a**2*(n + 2)*(n + 3)) - 2*b**2*x**(-3*n - 3)*(a*x**2 + b*x**3)**(n + 1)/(a**3*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_282():
    f = x**11/(a*x**2 + b*x**5)**3
    F = x**6/(6*a*(a + b*x**3)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_283():
    f = x**9/sqrt(a*x**2 + b*x**5)
    F = 16*a**2*sqrt(a*x**2 + b*x**5)/(45*b**3*x) - 8*a*x**2*sqrt(a*x**2 + b*x**5)/(45*b**2) + 2*x**5*sqrt(a*x**2 + b*x**5)/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_284():
    f = x**6/sqrt(a*x**2 + b*x**5)
    F = -4*a*sqrt(a*x**2 + b*x**5)/(9*b**2*x) + 2*x**2*sqrt(a*x**2 + b*x**5)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_285():
    f = x**3/sqrt(a*x**2 + b*x**5)
    F = 2*sqrt(a*x**2 + b*x**5)/(3*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_286():
    f = 1/sqrt(a*x**2 + b*x**5)
    F = -2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**5))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_287():
    f = 1/(x**3*sqrt(a*x**2 + b*x**5))
    F = -sqrt(a*x**2 + b*x**5)/(3*a*x**4) + b*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**5))/(3*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_288():
    f = x**4/sqrt(a*x**2 + b*x**5)
    F = -4*3**(sympy.S(3)/4)*a*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(15*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) + 2*sqrt(a*x**2 + b*x**5)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_289():
    f = x/sqrt(a*x**2 + b*x**5)
    F = 2*3**(sympy.S(3)/4)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_290():
    f = 1/(x**2*sqrt(a*x**2 + b*x**5))
    F = -3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) - sqrt(a*x**2 + b*x**5)/(2*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_291():
    f = x**5/sqrt(a*x**2 + b*x**5)
    F = 4*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) - 8*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) - 8*a*x*(a + b*x**3)/(7*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)*sqrt(a*x**2 + b*x**5)) + 2*x*sqrt(a*x**2 + b*x**5)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_292():
    f = x**2/sqrt(a*x**2 + b*x**5)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) + 2*x*(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_293():
    f = 1/(x*sqrt(a*x**2 + b*x**5))
    F = b**(sympy.S(1)/3)*x*(a + b*x**3)/(a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)*sqrt(a*x**2 + b*x**5)) - sqrt(a*x**2 + b*x**5)/(a*x**2) - 3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5)) + sqrt(2)*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*x*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_294():
    f = x**(sympy.S(13)/2)/sqrt(a*x**2 + b*x**5)
    F = 7*3**(sympy.S(3)/4)*a**(sympy.S(5)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(120*b**2*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) - 7*a*sqrt(a*x**2 + b*x**5)/(20*b**2*sqrt(x)) + x**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**5)/(5*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_295():
    f = x**(sympy.S(11)/2)/sqrt(a*x**2 + b*x**5)
    F = 5*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(8*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) + 3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(5 - 5*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(48*b**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) - a*x**(sympy.S(3)/2)*(5 + 5*sqrt(3))*(a + b*x**3)/(8*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x**2 + b*x**5)) + x**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**5)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_296():
    f = x**(sympy.S(9)/2)/sqrt(a*x**2 + b*x**5)
    F = -a*atanh(sqrt(b)*x**(sympy.S(5)/2)/sqrt(a*x**2 + b*x**5))/(3*b**(sympy.S(3)/2)) + sqrt(x)*sqrt(a*x**2 + b*x**5)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_297():
    f = x**(sympy.S(7)/2)/sqrt(a*x**2 + b*x**5)
    F = -3**(sympy.S(3)/4)*a**(sympy.S(2)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(12*b*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) + sqrt(a*x**2 + b*x**5)/(2*b*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_298():
    f = x**(sympy.S(5)/2)/sqrt(a*x**2 + b*x**5)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(b**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) - 3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(6*b**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) + x**(sympy.S(3)/2)*(1 + sqrt(3))*(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_299():
    f = x**(sympy.S(3)/2)/sqrt(a*x**2 + b*x**5)
    F = 2*atanh(sqrt(b)*x**(sympy.S(5)/2)/sqrt(a*x**2 + b*x**5))/(3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_300():
    f = sqrt(x)/sqrt(a*x**2 + b*x**5)
    F = 3**(sympy.S(3)/4)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(1)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_301():
    f = 1/(sqrt(x)*sqrt(a*x**2 + b*x**5))
    F = b**(sympy.S(1)/3)*x**(sympy.S(3)/2)*(2 + 2*sqrt(3))*(a + b*x**3)/(a*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x**2 + b*x**5)) - 2*sqrt(a*x**2 + b*x**5)/(a*x**(sympy.S(3)/2)) - 2*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(1 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(3*a**(sympy.S(2)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_302():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**5))
    F = -2*sqrt(a*x**2 + b*x**5)/(3*a*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_303():
    f = 1/(x**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**5))
    F = -2*sqrt(a*x**2 + b*x**5)/(5*a*x**(sympy.S(7)/2)) - 2*3**(sympy.S(3)/4)*b*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(15*a**(sympy.S(4)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_304():
    f = 1/(x**(sympy.S(7)/2)*sqrt(a*x**2 + b*x**5))
    F = -2*sqrt(a*x**2 + b*x**5)/(7*a*x**(sympy.S(9)/2)) - b**(sympy.S(4)/3)*x**(sympy.S(3)/2)*(8 + 8*sqrt(3))*(a + b*x**3)/(7*a**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))*sqrt(a*x**2 + b*x**5)) + 8*b*sqrt(a*x**2 + b*x**5)/(7*a**2*x**(sympy.S(3)/2)) + 8*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(7*a**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5)) + 3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(4 - 4*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(21*a**(sympy.S(5)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_305():
    f = 1/(x**(sympy.S(9)/2)*sqrt(a*x**2 + b*x**5))
    F = -2*sqrt(a*x**2 + b*x**5)/(9*a*x**(sympy.S(11)/2)) + 4*b*sqrt(a*x**2 + b*x**5)/(9*a**2*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_306():
    f = 1/(x**(sympy.S(11)/2)*sqrt(a*x**2 + b*x**5))
    F = -2*sqrt(a*x**2 + b*x**5)/(11*a*x**(sympy.S(13)/2)) + 16*b*sqrt(a*x**2 + b*x**5)/(55*a**2*x**(sympy.S(7)/2)) + 16*3**(sympy.S(3)/4)*b**2*x**(sympy.S(3)/2)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(acos((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 - sqrt(3)))/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))), sqrt(3)/4 + sympy.S.Half)/(165*a**(sympy.S(7)/3)*sqrt(b**(sympy.S(1)/3)*x*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x*(1 + sqrt(3)))**2)*sqrt(a*x**2 + b*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_307():
    f = x/(a*x**3 + b*x**4)
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_308():
    f = 1/(a*x**3 + b*x**4)
    F = -1/(2*a*x**2) + b/(a**2*x) + b**2*log(x)/a**3 - b**2*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_309():
    f = x**4/sqrt(a*x**3 + b*x**4)
    F = -5*a**3*atanh(sqrt(b)*x**2/sqrt(a*x**3 + b*x**4))/(8*b**(sympy.S(7)/2)) + 5*a**2*sqrt(a*x**3 + b*x**4)/(8*b**3*x) - 5*a*sqrt(a*x**3 + b*x**4)/(12*b**2) + x*sqrt(a*x**3 + b*x**4)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_310():
    f = x**3/sqrt(a*x**3 + b*x**4)
    F = 3*a**2*atanh(sqrt(b)*x**2/sqrt(a*x**3 + b*x**4))/(4*b**(sympy.S(5)/2)) - 3*a*sqrt(a*x**3 + b*x**4)/(4*b**2*x) + sqrt(a*x**3 + b*x**4)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_311():
    f = x**2/sqrt(a*x**3 + b*x**4)
    F = -a*atanh(sqrt(b)*x**2/sqrt(a*x**3 + b*x**4))/b**(sympy.S(3)/2) + sqrt(a*x**3 + b*x**4)/(b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_312():
    f = x/sqrt(a*x**3 + b*x**4)
    F = 2*atanh(sqrt(b)*x**2/sqrt(a*x**3 + b*x**4))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_313():
    f = 1/sqrt(a*x**3 + b*x**4)
    F = -2*sqrt(a*x**3 + b*x**4)/(a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_314():
    f = 1/(x*sqrt(a*x**3 + b*x**4))
    F = -2*sqrt(a*x**3 + b*x**4)/(3*a*x**3) + 4*b*sqrt(a*x**3 + b*x**4)/(3*a**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_315():
    f = 1/(x**2*sqrt(a*x**3 + b*x**4))
    F = -2*sqrt(a*x**3 + b*x**4)/(5*a*x**4) + 8*b*sqrt(a*x**3 + b*x**4)/(15*a**2*x**3) - 16*b**2*sqrt(a*x**3 + b*x**4)/(15*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_316():
    f = 1/(x**3*sqrt(a*x**3 + b*x**4))
    F = -2*sqrt(a*x**3 + b*x**4)/(7*a*x**5) + 12*b*sqrt(a*x**3 + b*x**4)/(35*a**2*x**4) - 16*b**2*sqrt(a*x**3 + b*x**4)/(35*a**3*x**3) + 32*b**3*sqrt(a*x**3 + b*x**4)/(35*a**4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_317():
    f = 1/(x**4*sqrt(a*x**3 + b*x**4))
    F = -2*sqrt(a*x**3 + b*x**4)/(9*a*x**6) + 16*b*sqrt(a*x**3 + b*x**4)/(63*a**2*x**5) - 32*b**2*sqrt(a*x**3 + b*x**4)/(105*a**3*x**4) + 128*b**3*sqrt(a*x**3 + b*x**4)/(315*a**4*x**3) - 256*b**4*sqrt(a*x**3 + b*x**4)/(315*a**5*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_318():
    f = 1/(b*x**5 + x**3)
    F = -b*log(x) + b*log(b*x**2 + 1)/2 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_319():
    f = 1/(b*x**5 - x**3)
    F = -b*log(x) + b*log(-b*x**2 + 1)/2 + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_320():
    f = 1/(a*x + b*x)
    F = log(x)/(a + b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_321():
    f = (a*x + b*x)**(-2)
    F = -1/(x*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_322():
    f = (a*x + b*x)**(-3)
    F = -1/(2*x**2*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_323():
    f = 1/(a*x**2 + b*x**2)
    F = -1/(x*(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_324():
    f = 1/(a*x**n + b*x**n)
    F = x**(1 - n)/((1 - n)*(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_325():
    f = (a*x**n + b*x**n)**(-2)
    F = x**(1 - 2*n)/((1 - 2*n)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_326():
    f = (a*x**n + b*x**n)**(-3)
    F = x**(1 - 3*n)/((1 - 3*n)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_327():
    f = (a*x + b*x**14)**12
    F = (a + b*x**13)**13/(169*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_328():
    f = x**12*(a*x + b*x**26)**12
    F = (a + b*x**25)**13/(325*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_329():
    f = x**24*(a*x + b*x**38)**12
    F = (a + b*x**37)**13/(481*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_330():
    f = x**(12*m - 12)*(a*x + b*x**(12*m + 2))**12
    F = (a + b*x**(12*m + 1))**13/(13*b*(12*m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_331():
    f = (a*x + b*x**14)**12
    F = (a + b*x**13)**13/(169*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_332():
    f = (a*x**2 + b*x**27)**12
    F = (a + b*x**25)**13/(325*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_333():
    f = (a*x**3 + b*x**40)**12
    F = (a + b*x**37)**13/(481*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_334():
    f = (a*x**m + b*x**(13*m + 1))**12
    F = (a + b*x**(12*m + 1))**13/(13*b*(12*m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_335():
    f = (a*x**m + b*x**(6*m + 1))**5
    F = (a + b*x**(5*m + 1))**6/(6*b*(5*m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_336():
    f = (a*x**m + b*x**(1 - 2*m))**(-3)
    F = -1/(2*b*(1 - 3*m)*(a + b*x**(1 - 3*m))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_337():
    f = 1/(a*x + b/x)
    F = log(a*x**2 + b)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_338():
    f = 1/(a*x + b/x**2)
    F = log(a*x**3 + b)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_339():
    f = 1/(a*x + b/x**3)
    F = log(a*x**4 + b)/(4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_340():
    f = (a*x + b/x)**(-3)
    F = x**4/(4*b*(a*x**2 + b)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_341():
    f = (a*x**2 + b/x**3)**(-3)
    F = x**10/(10*b*(a*x**5 + b)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_342():
    f = (a*x**3 + b/x**5)**(-3)
    F = x**16/(16*b*(a*x**8 + b)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_343():
    f = (a/x + b*x)**2
    F = -a**2/x + 2*a*b*x + b**2*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_344():
    f = (a/x + b*x)**3
    F = -a**3/(2*x**2) + 3*a**2*b*log(x) + 3*a*b**2*x**2/2 + b**3*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_345():
    f = (a/x + b*x)**4
    F = -a**4/(3*x**3) - 4*a**3*b/x + 6*a**2*b**2*x + 4*a*b**3*x**3/3 + b**4*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_346():
    f = 1/(x**3 + x**(-2))
    F = log(x + 1)/5 - (sympy.S(1)/20 + sqrt(5)/20)*log(x**2 - x*(sympy.S.Half - sqrt(5)/2) + 1) - (sympy.S(1)/20 - sqrt(5)/20)*log(x**2 - x*(sympy.S.Half + sqrt(5)/2) + 1) + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(x*sqrt(2*sqrt(5)/5 + 2) - sqrt(2*sqrt(5)/5 + 1))/5 - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(2*sqrt(2)*x/sqrt(sqrt(5) + 5) + sqrt(1 - 2*sqrt(5)/5))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_347():
    f = x**p*(a*x**n + b*x**(13*n + p + 1))**12
    F = (a + b*x**(12*n + p + 1))**13/(13*b*(12*n + p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_348():
    f = x**12*(a + b*x**13)**12
    F = (a + b*x**13)**13/(169*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_349():
    f = x**12*(a*x + b*x**26)**12
    F = (a + b*x**25)**13/(325*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_350():
    f = x**12*(a*x**2 + b*x**39)**12
    F = (a + b*x**37)**13/(481*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_351():
    f = x**24*(a + b*x**25)**12
    F = (a + b*x**25)**13/(325*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_352():
    f = x**24*(a*x + b*x**38)**12
    F = (a + b*x**37)**13/(481*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_353():
    f = x**36*(a + b*x**37)**12
    F = (a + b*x**37)**13/(481*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_354():
    f = 1/(a*x + b*x**n)
    F = log(a*x**(1 - n) + b)/(a*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_355():
    f = 1/(a*x + b*x**(n + 1))
    F = log(x)/a - log(a + b*x**n)/(a*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_356():
    f = 1/(a*x + b*x**(1 - n))
    F = log(a*x**n + b)/(a*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_357():
    f = 1/(2*x + 3*x**(n + 1))
    F = log(x)/2 - log(3*x**n + 2)/(2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_358():
    f = 1/(2*x + 3*x**(1 - n))
    F = log(2*x**n + 3)/(2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_359():
    f = 1/(-sqrt(x) + x)
    F = 2*log(1 - sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_360():
    f = 1/(-x**(sympy.S(3)/5) + x)
    F = 5*log(1 - x**(sympy.S(2)/5))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_361():
    f = 1/(x + x**(sympy.S(-1)/3))
    F = 3*log(x**(sympy.S(4)/3) + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_362():
    f = 1/(x + x**(sqrt(2)))
    F = log(x) - (1 + sqrt(2))*log(x**(-1 + sqrt(2)) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_363():
    f = x**(-j/2 - 1)*sqrt(a*x**j + b*x**n)
    F = 2*sqrt(a)*atanh(sqrt(a)*x**(j/2)/sqrt(a*x**j + b*x**n))/(j - n) - 2*sqrt(a*x**j + b*x**n)/(x**(j/2)*(j - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_364():
    f = (c*x)**(-j/2 - 1)*sqrt(a*x**j + b*x**n)
    F = 2*sqrt(a)*x**(j/2)*atanh(sqrt(a)*x**(j/2)/sqrt(a*x**j + b*x**n))/(c*(c*x)**(j/2)*(j - n)) - 2*sqrt(a*x**j + b*x**n)/(c*(c*x)**(j/2)*(j - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_365():
    f = sqrt(a*x**3 + b*x**n)/(c*x)**(sympy.S(5)/2)
    F = 2*sqrt(a)*sqrt(c*x)*atanh(sqrt(a)*x**(sympy.S(3)/2)/sqrt(a*x**3 + b*x**n))/(c**3*sqrt(x)*(3 - n)) - 2*sqrt(a*x**3 + b*x**n)/(c*(c*x)**(sympy.S(3)/2)*(3 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_366():
    f = sqrt(a*x**2 + b*x**n)/(c**2*x**2)
    F = 2*sqrt(a)*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**n))/(c**2*(2 - n)) - 2*sqrt(a*x**2 + b*x**n)/(c**2*x*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_367():
    f = sqrt(a*x + b*x**n)/(c*x)**(sympy.S(3)/2)
    F = 2*sqrt(a)*sqrt(x)*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**n))/(c*sqrt(c*x)*(1 - n)) - 2*sqrt(a*x + b*x**n)/(c*sqrt(c*x)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_368():
    f = sqrt(a + b*x**n)/(c*x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*x**n)/sqrt(a))/(c*n) + 2*sqrt(a + b*x**n)/(c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_369():
    f = sqrt(a/x + b*x**n)/sqrt(c*x)
    F = -2*sqrt(a)*sqrt(x)*atanh(sqrt(a)/(sqrt(x)*sqrt(a/x + b*x**n)))/(sqrt(c*x)*(n + 1)) + 2*sqrt(c*x)*sqrt(a/x + b*x**n)/(c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_370():
    f = sqrt(a/x**2 + b*x**n)
    F = -2*sqrt(a)*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x**n)))/(n + 2) + 2*x*sqrt(a/x**2 + b*x**n)/(n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_371():
    f = sqrt(c*x)*sqrt(a/x**3 + b*x**n)
    F = -2*sqrt(a)*c*sqrt(x)*atanh(sqrt(a)/(x**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)))/(sqrt(c*x)*(n + 3)) + 2*(c*x)**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)/(c*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_372():
    f = (c*x)**(-3*j/2 - 1)*(a*x**j + b*x**n)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*x**(3*j/2)*atanh(sqrt(a)*x**(j/2)/sqrt(a*x**j + b*x**n))/(c*(c*x)**(3*j/2)*(j - n)) - 2*a*x**j*sqrt(a*x**j + b*x**n)/(c*(c*x)**(3*j/2)*(j - n)) - 2*(a*x**j + b*x**n)**(sympy.S(3)/2)/(3*c*(c*x)**(3*j/2)*(j - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_373():
    f = (a*x**3 + b*x**n)**(sympy.S(3)/2)/(c*x)**(sympy.S(11)/2)
    F = 2*a**(sympy.S(3)/2)*sqrt(c*x)*atanh(sqrt(a)*x**(sympy.S(3)/2)/sqrt(a*x**3 + b*x**n))/(c**6*sqrt(x)*(3 - n)) - 2*a*sqrt(a*x**3 + b*x**n)/(c**4*(c*x)**(sympy.S(3)/2)*(3 - n)) - 2*(a*x**3 + b*x**n)**(sympy.S(3)/2)/(3*c*(c*x)**(sympy.S(9)/2)*(3 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_374():
    f = (a*x**2 + b*x**n)**(sympy.S(3)/2)/(c**4*x**4)
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**n))/(c**4*(2 - n)) - 2*a*sqrt(a*x**2 + b*x**n)/(c**4*x*(2 - n)) - 2*(a*x**2 + b*x**n)**(sympy.S(3)/2)/(3*c**4*x**3*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_375():
    f = (a*x + b*x**n)**(sympy.S(3)/2)/(c*x)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(3)/2)*sqrt(x)*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**n))/(c**2*sqrt(c*x)*(1 - n)) - 2*a*sqrt(a*x + b*x**n)/(c**2*sqrt(c*x)*(1 - n)) - 2*(a*x + b*x**n)**(sympy.S(3)/2)/(3*c*(c*x)**(sympy.S(3)/2)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_376():
    f = (a + b*x**n)**(sympy.S(3)/2)/(c*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x**n)/sqrt(a))/(c*n) + 2*a*sqrt(a + b*x**n)/(c*n) + 2*(a + b*x**n)**(sympy.S(3)/2)/(3*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_377():
    f = sqrt(c*x)*(a/x + b*x**n)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*c*sqrt(x)*atanh(sqrt(a)/(sqrt(x)*sqrt(a/x + b*x**n)))/(sqrt(c*x)*(n + 1)) + 2*a*sqrt(c*x)*sqrt(a/x + b*x**n)/(n + 1) + 2*(c*x)**(sympy.S(3)/2)*(a/x + b*x**n)**(sympy.S(3)/2)/(3*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_378():
    f = c**2*x**2*(a/x**2 + b*x**n)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*c**2*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x**n)))/(n + 2) + 2*a*c**2*x*sqrt(a/x**2 + b*x**n)/(n + 2) + 2*c**2*x**3*(a/x**2 + b*x**n)**(sympy.S(3)/2)/(3*n + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_379():
    f = (c*x)**(sympy.S(7)/2)*(a/x**3 + b*x**n)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*c**4*sqrt(x)*atanh(sqrt(a)/(x**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)))/(sqrt(c*x)*(n + 3)) + 2*a*c**2*(c*x)**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)/(n + 3) + 2*(c*x)**(sympy.S(9)/2)*(a/x**3 + b*x**n)**(sympy.S(3)/2)/(3*c*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_380():
    f = c**5*x**5*(a/x**4 + b*x**n)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*c**5*atanh(sqrt(a)/(x**2*sqrt(a/x**4 + b*x**n)))/(n + 4) + 2*a*c**5*x**2*sqrt(a/x**4 + b*x**n)/(n + 4) + 2*c**5*x**6*(a/x**4 + b*x**n)**(sympy.S(3)/2)/(3*n + 12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_381():
    f = sqrt((a + b*x)/x**2)
    F = -2*sqrt(a)*atanh(sqrt(a)/(x*sqrt(a/x**2 + b/x))) + 2*x*sqrt(a/x**2 + b/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_382():
    f = sqrt((a + b*x**2)/x**2)
    F = -sqrt(a)*atanh(sqrt(a)/(x*sqrt(a/x**2 + b))) + x*sqrt(a/x**2 + b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_383():
    f = sqrt((a + b*x**3)/x**2)
    F = -2*sqrt(a)*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x)))/3 + 2*x*sqrt(a/x**2 + b*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_384():
    f = sqrt((a + b*x**n)/x**2)
    F = -2*sqrt(a)*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x**(n - 2))))/n + 2*x*sqrt(a/x**2 + b*x**(n - 2))/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_385():
    f = sqrt((-a + b*x)/x**2)
    F = 2*sqrt(a)*atan(sqrt(a)/(x*sqrt(-a/x**2 + b/x))) + 2*x*sqrt(-a/x**2 + b/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_386():
    f = sqrt((-a + b*x**2)/x**2)
    F = sqrt(a)*atan(sqrt(a)/(x*sqrt(-a/x**2 + b))) + x*sqrt(-a/x**2 + b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_387():
    f = sqrt((-a + b*x**3)/x**2)
    F = 2*sqrt(a)*atan(sqrt(a)/(x*sqrt(-a/x**2 + b*x)))/3 + 2*x*sqrt(-a/x**2 + b*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_388():
    f = sqrt((-a + b*x**n)/x**2)
    F = 2*sqrt(a)*atan(sqrt(a)/(x*sqrt(-a/x**2 + b*x**(n - 2))))/n + 2*x*sqrt(-a/x**2 + b*x**(n - 2))/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_389():
    f = (c*x)**(j/2 - 1)/sqrt(a*x**j + b*x**n)
    F = 2*(c*x)**(j/2)*atanh(sqrt(a)*x**(j/2)/sqrt(a*x**j + b*x**n))/(sqrt(a)*c*x**(j/2)*(j - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_390():
    f = sqrt(c*x)/sqrt(a*x**3 + b*x**n)
    F = 2*sqrt(c*x)*atanh(sqrt(a)*x**(sympy.S(3)/2)/sqrt(a*x**3 + b*x**n))/(sqrt(a)*sqrt(x)*(3 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_391():
    f = 1/sqrt(a*x**2 + b*x**n)
    F = 2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**n))/(sqrt(a)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_392():
    f = 1/(sqrt(c*x)*sqrt(a*x + b*x**n))
    F = 2*sqrt(x)*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**n))/(sqrt(a)*sqrt(c*x)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_393():
    f = 1/(c*x*sqrt(a + b*x**n))
    F = -2*atanh(sqrt(a + b*x**n)/sqrt(a))/(sqrt(a)*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_394():
    f = 1/((c*x)**(sympy.S(3)/2)*sqrt(a/x + b*x**n))
    F = -2*sqrt(x)*atanh(sqrt(a)/(sqrt(x)*sqrt(a/x + b*x**n)))/(sqrt(a)*c*sqrt(c*x)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_395():
    f = 1/(c**2*x**2*sqrt(a/x**2 + b*x**n))
    F = -2*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x**n)))/(sqrt(a)*c**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_396():
    f = 1/((c*x)**(sympy.S(5)/2)*sqrt(a/x**3 + b*x**n))
    F = -2*sqrt(x)*atanh(sqrt(a)/(x**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)))/(sqrt(a)*c**2*sqrt(c*x)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_397():
    f = (c*x)**(3*j/2 - 1)/(a*x**j + b*x**n)**(sympy.S(3)/2)
    F = -2*(c*x)**(3*j/2)/(a*c*x**j*(j - n)*sqrt(a*x**j + b*x**n)) + 2*(c*x)**(3*j/2)*atanh(sqrt(a)*x**(j/2)/sqrt(a*x**j + b*x**n))/(a**(sympy.S(3)/2)*c*x**(3*j/2)*(j - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_398():
    f = (c*x)**(sympy.S(7)/2)/(a*x**3 + b*x**n)**(sympy.S(3)/2)
    F = -2*c**2*(c*x)**(sympy.S(3)/2)/(a*(3 - n)*sqrt(a*x**3 + b*x**n)) + 2*c**3*sqrt(c*x)*atanh(sqrt(a)*x**(sympy.S(3)/2)/sqrt(a*x**3 + b*x**n))/(a**(sympy.S(3)/2)*sqrt(x)*(3 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_399():
    f = c**2*x**2/(a*x**2 + b*x**n)**(sympy.S(3)/2)
    F = -2*c**2*x/(a*(2 - n)*sqrt(a*x**2 + b*x**n)) + 2*c**2*atanh(sqrt(a)*x/sqrt(a*x**2 + b*x**n))/(a**(sympy.S(3)/2)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_400():
    f = sqrt(c*x)/(a*x + b*x**n)**(sympy.S(3)/2)
    F = -2*sqrt(c*x)/(a*(1 - n)*sqrt(a*x + b*x**n)) + 2*c*sqrt(x)*atanh(sqrt(a)*sqrt(x)/sqrt(a*x + b*x**n))/(a**(sympy.S(3)/2)*sqrt(c*x)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_401():
    f = 1/(c*x*(a + b*x**n)**(sympy.S(3)/2))
    F = 2/(a*c*n*sqrt(a + b*x**n)) - 2*atanh(sqrt(a + b*x**n)/sqrt(a))/(a**(sympy.S(3)/2)*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_402():
    f = 1/((c*x)**(sympy.S(5)/2)*(a/x + b*x**n)**(sympy.S(3)/2))
    F = 2/(a*c**2*sqrt(c*x)*(n + 1)*sqrt(a/x + b*x**n)) - 2*sqrt(x)*atanh(sqrt(a)/(sqrt(x)*sqrt(a/x + b*x**n)))/(a**(sympy.S(3)/2)*c**2*sqrt(c*x)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_403():
    f = 1/(c**4*x**4*(a/x**2 + b*x**n)**(sympy.S(3)/2))
    F = 2/(a*c**4*x*(n + 2)*sqrt(a/x**2 + b*x**n)) - 2*atanh(sqrt(a)/(x*sqrt(a/x**2 + b*x**n)))/(a**(sympy.S(3)/2)*c**4*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_404():
    f = 1/((c*x)**(sympy.S(11)/2)*(a/x**3 + b*x**n)**(sympy.S(3)/2))
    F = 2/(a*c**4*(c*x)**(sympy.S(3)/2)*(n + 3)*sqrt(a/x**3 + b*x**n)) - 2*sqrt(x)*atanh(sqrt(a)/(x**(sympy.S(3)/2)*sqrt(a/x**3 + b*x**n)))/(a**(sympy.S(3)/2)*c**5*sqrt(c*x)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_405():
    f = 1/(c**7*x**7*(a/x**4 + b*x**n)**(sympy.S(3)/2))
    F = 2/(a*c**7*x**2*(n + 4)*sqrt(a/x**4 + b*x**n)) - 2*atanh(sqrt(a)/(x**2*sqrt(a/x**4 + b*x**n)))/(a**(sympy.S(3)/2)*c**7*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_406():
    f = 1/sqrt((a + b*x**3)/x)
    F = 2*atanh(sqrt(b)*x/sqrt(a/x + b*x**2))/(3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_407():
    f = 1/sqrt((a + b*x**4)/x**2)
    F = atanh(sqrt(b)*x/sqrt(a/x**2 + b*x**2))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_408():
    f = 1/sqrt((a + b*x**5)/x**3)
    F = 2*atanh(sqrt(b)*x/sqrt(a/x**3 + b*x**2))/(5*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_409():
    f = 1/sqrt(x**(2 - n)*(a + b*x**n))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x**(2 - n) + b*x**2))/(sqrt(b)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_410():
    f = 1/sqrt((a - b*x**3)/x)
    F = 2*atan(sqrt(b)*x/sqrt(a/x - b*x**2))/(3*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_411():
    f = 1/sqrt((a - b*x**4)/x**2)
    F = atan(sqrt(b)*x/sqrt(a/x**2 - b*x**2))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_412():
    f = 1/sqrt((a - b*x**5)/x**3)
    F = 2*atan(sqrt(b)*x/sqrt(a/x**3 - b*x**2))/(5*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_413():
    f = 1/sqrt(x**(2 - n)*(a - b*x**n))
    F = 2*atan(sqrt(b)*x/sqrt(a*x**(2 - n) - b*x**2))/(sqrt(b)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_414():
    f = 1/sqrt(x**n*(a + b*x**(2 - n)))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x**n + b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_415():
    f = 1/sqrt(x**2*(a*x**(n - 2) + b))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x**n + b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_416():
    f = 1/sqrt(x*(a*x**(n - 1) + b*x))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x**n + b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_417():
    f = 1/sqrt(x**n*(a - b*x**(2 - n)))
    F = 2*atan(sqrt(b)*x/sqrt(a*x**n - b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_418():
    f = 1/sqrt(x**2*(a*x**(n - 2) - b))
    F = 2*atan(sqrt(b)*x/sqrt(a*x**n - b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_419():
    f = 1/sqrt(x*(a*x**(n - 1) - b*x))
    F = 2*atan(sqrt(b)*x/sqrt(a*x**n - b*x**2))/(sqrt(b)*(2 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_420():
    f = (c*x)**m*(a*x**j + b*x**n)**(sympy.S(3)/2)
    F = 2*b*x**(n + 1)*(c*x)**m*sqrt(a*x**j + b*x**n)*hyper((sympy.S(-3)/2, (m + 3*n/2 + 1)/(j - n)), (1 + (m + 3*n/2 + 1)/(j - n),), -a*x**(j - n)/b)/(sqrt(a*x**(j - n)/b + 1)*(2*m + 3*n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_421():
    f = (c*x)**m*sqrt(a*x**j + b*x**n)
    F = 2*x*(c*x)**m*sqrt(a*x**j + b*x**n)*hyper((sympy.S(-1)/2, (m + n/2 + 1)/(j - n)), (1 + (2*m + n + 2)/(2*j - 2*n),), -a*x**(j - n)/b)/(sqrt(a*x**(j - n)/b + 1)*(2*m + n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_422():
    f = (c*x)**m/sqrt(a*x**j + b*x**n)
    F = 2*x*(c*x)**m*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S.Half, (m - n/2 + 1)/(j - n)), (1 + (m - n/2 + 1)/(j - n),), -a*x**(j - n)/b)/(sqrt(a*x**j + b*x**n)*(2*m - n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_423():
    f = (c*x)**m/(a*x**j + b*x**n)**(sympy.S(3)/2)
    F = 2*x**(1 - n)*(c*x)**m*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S(3)/2, (m - 3*n/2 + 1)/(j - n)), (1 + (m - 3*n/2 + 1)/(j - n),), -a*x**(j - n)/b)/(b*sqrt(a*x**j + b*x**n)*(2*m - 3*n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_424():
    f = (c*x)**m/(a*x**j + b*x**n)**(sympy.S(5)/2)
    F = 2*x**(1 - 2*n)*(c*x)**m*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S(5)/2, (m - 5*n/2 + 1)/(j - n)), (1 + (m - 5*n/2 + 1)/(j - n),), -a*x**(j - n)/b)/(b**2*sqrt(a*x**j + b*x**n)*(2*m - 5*n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_425():
    f = (a*x**j + b*x**n)**(sympy.S(3)/2)
    F = 2*b*x**(n + 1)*sqrt(a*x**j + b*x**n)*hyper((sympy.S(-3)/2, (3*n/2 + 1)/(j - n)), ((2*j + n + 2)/(2*j - 2*n),), -a*x**(j - n)/b)/((3*n + 2)*sqrt(a*x**(j - n)/b + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_426():
    f = sqrt(a*x**j + b*x**n)
    F = 2*x*sqrt(a*x**j + b*x**n)*hyper((sympy.S(-1)/2, (n + 2)/(2*j - 2*n)), (1 + (n + 2)/(2*j - 2*n),), -a*x**(j - n)/b)/((n + 2)*sqrt(a*x**(j - n)/b + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_427():
    f = 1/sqrt(a*x**j + b*x**n)
    F = 2*x*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S.Half, (2 - n)/(2*j - 2*n)), ((2 - n)/(2*(j - n)) + 1,), -a*x**(j - n)/b)/((2 - n)*sqrt(a*x**j + b*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_428():
    f = (a*x**j + b*x**n)**(sympy.S(-3)/2)
    F = 2*x**(1 - n)*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S(3)/2, (1 - 3*n/2)/(j - n)), ((2 - 3*n)/(2*j - 2*n) + 1,), -a*x**(j - n)/b)/(b*(2 - 3*n)*sqrt(a*x**j + b*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_429():
    f = (a*x**j + b*x**n)**(sympy.S(-5)/2)
    F = 2*x**(1 - 2*n)*sqrt(a*x**(j - n)/b + 1)*hyper((sympy.S(5)/2, (1 - 5*n/2)/(j - n)), ((2 - 5*n)/(2*j - 2*n) + 1,), -a*x**(j - n)/b)/(b**2*(2 - 5*n)*sqrt(a*x**j + b*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_430():
    f = sqrt((x + 1)/x**5)
    F = -2*x**6*(x**(-4) + x**(-5))**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_431():
    f = sqrt(x**(sympy.S(5)/2) + x)
    F = 4*(x**(sympy.S(5)/2) + x)**(sympy.S(3)/2)/(9*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_432():
    f = 1/(x**(sympy.S(3)/2) + sqrt(x))
    F = 2*atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_433():
    f = x*sqrt(x**2*(a + b*x**3))
    F = 2*(x**2*(a + b*x**3))**(sympy.S(3)/2)/(9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_434():
    f = x*sqrt(a*x**2 + b*x**5)
    F = 2*(a*x**2 + b*x**5)**(sympy.S(3)/2)/(9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_435():
    f = sqrt(x**4*(a + b*x**3))
    F = 2*(a*x**4 + b*x**7)**(sympy.S(3)/2)/(9*b*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_436():
    f = (a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(-1)/3)
    F = -45*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*a**4*sqrt((2*2**(sympy.S(1)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))/a**2)**(sympy.S(1)/3)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(56*b**3*sqrt(-(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + 2*b*x**(sympy.S(1)/3))*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(1)/3)) + 15*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*a**4*sqrt((2*2**(sympy.S(1)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))/a**2)**(sympy.S(1)/3)*(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(14*b**3*sqrt(-(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + 2*b*x**(sympy.S(1)/3))*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(1)/3)) - 45*2**(sympy.S(2)/3)*a**2*(-b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))/a**2)**(sympy.S(1)/3)*(a + 2*b*x**(sympy.S(1)/3))/(28*b**3*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(1)/3)*(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)) - 45*a*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/(28*b**2*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(1)/3)) + x**(sympy.S(2)/3)*(9*a + 9*b*x**(sympy.S(1)/3))/(7*b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_437():
    f = (a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(-2)/3)
    F = 6*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*a**4*sqrt((2*2**(sympy.S(1)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))/a**2)**(sympy.S(2)/3)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(5*b**3*sqrt(-(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-b*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/a**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(a + 2*b*x**(sympy.S(1)/3))*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(2)/3)) - 18*a*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/(5*b**2*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(2)/3)) + x**(sympy.S(2)/3)*(9*a + 9*b*x**(sympy.S(1)/3))/(5*b*(a*x**(sympy.S(1)/3) + b*x**(sympy.S(2)/3))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_438():
    f = x**(-p*q - 1)*(a*x**q + b*x**n)**p
    F = -(a + b*x**(n - q))*(a*x**q + b*x**n)**p*hyper((1, p + 1), (p + 2,), 1 + b*x**(n - q)/a)/(a*x**(p*q)*(n - q)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_439():
    f = x**(-n - q*(p - 1) - 1)*(a*x**q + b*x**n)**p
    F = b*(a + b*x**(n - q))*(a*x**q + b*x**n)**p*hyper((2, p + 1), (p + 2,), 1 + b*x**(n - q)/a)/(a**2*x**(p*q)*(n - q)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_440():
    f = x**(-n*(p - 1) - q - 1)*(a*x**q + b*x**n)**p
    F = x**(-n*p + n - q)*(a*x**q + b*x**n)**p*hyper((-p, 1 - p), (2 - p,), -b*x**(n - q)/a)/((1 - p)*(1 + b*x**(n - q)/a)**p*(n - q))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_441():
    f = (a*x**m + b*x**(m*p + m + 1))**p
    F = (a*x**m + b*x**(m*p + m + 1))**(p + 1)/(b*x**(m*(p + 1))*(p + 1)*(m*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_442():
    f = (x**m*(a + b*x**(m*p + 1)))**p
    F = (a*x**m + b*x**(m*p + m + 1))**(p + 1)/(b*x**(m*(p + 1))*(p + 1)*(m*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_443():
    f = x**n*(x**m*(a + b*x**(m*p + n + 1)))**p
    F = (a*x**m + b*x**(m*p + m + n + 1))**(p + 1)/(b*x**(m*(p + 1))*(p + 1)*(m*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_444():
    f = x**n*(a*x**m + b*x**(m*p + m + n + 1))**p
    F = (a*x**m + b*x**(m*p + m + n + 1))**(p + 1)/(b*x**(m*(p + 1))*(p + 1)*(m*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_445():
    f = sqrt(x**(2*n - 2)*(a + b*x**n))
    F = 2*x**(3 - 3*n)*(a*x**(2*n - 2) + b*x**(3*n - 2))**(sympy.S(3)/2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_446():
    f = (x**(3*n - 3)*(a + b*x**n))**(sympy.S(1)/3)
    F = 3*x**(4 - 4*n)*(a*x**(3*n - 3) + b*x**(4*n - 3))**(sympy.S(4)/3)/(4*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_447():
    f = (x**(4*n - 4)*(a + b*x**n))**(sympy.S(1)/4)
    F = 4*x**(5 - 5*n)*(a*x**(4*n - 4) + b*x**(5*n - 4))**(sympy.S(5)/4)/(5*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_448():
    f = (x**(p*(n - 1))*(a + b*x**n))**(1/p)
    F = p*x**((1 - n)*(p + 1))*(a/x**(p*(1 - n)) + b*x**(n - p*(1 - n)))**(1 + 1/p)/(b*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_449():
    f = (x**((n - 1)/p)*(a + b*x**n))**p
    F = x**((1 - n)*(p + 1)/p)*(a/x**((1 - n)/p) + b*x**(n - (1 - n)/p))**(p + 1)/(b*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_450():
    f = x**(n - p*(q + 1) - 1)*(a*x**n + b*x**p)**q
    F = (a*x**n + b*x**p)**(q + 1)/(a*x**(p*(q + 1))*(n - p)*(q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_4_Improper_1_1_4_2_c_x_pow_m_a_x_pow_j_plus_b_x_pow_n_pow_p_451():
    f = x**(-n*q - p*(q + 1) - 1)*(x**n*(a + b*x**p))**q
    F = -(a*x**n + b*x**(n + p))**(q + 1)/(a*p*x**((n + p)*(q + 1))*(q + 1))
    assert integrate(f, x) == F

