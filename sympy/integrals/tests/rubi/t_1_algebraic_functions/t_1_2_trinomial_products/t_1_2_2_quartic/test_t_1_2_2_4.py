"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.4 (f x)^m (d+e x^2)^q (a+b x^2+c x^4)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, p, q = symbols('A B a b c d e f m p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1():
    f = x**3*(a + c*x**4)**5*(d + e*x**2)
    F = a**5*d*x**4/4 + a**5*e*x**6/6 + 5*a**4*c*d*x**8/8 + a**4*c*e*x**10/2 + 5*a**3*c**2*d*x**12/6 + 5*a**3*c**2*e*x**14/7 + 5*a**2*c**3*d*x**16/8 + 5*a**2*c**3*e*x**18/9 + a*c**4*d*x**20/4 + 5*a*c**4*e*x**22/22 + c**5*d*x**24/24 + c**5*e*x**26/26
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_2():
    f = x**2*(a + c*x**4)**5*(d + e*x**2)
    F = a**5*d*x**3/3 + a**5*e*x**5/5 + 5*a**4*c*d*x**7/7 + 5*a**4*c*e*x**9/9 + 10*a**3*c**2*d*x**11/11 + 10*a**3*c**2*e*x**13/13 + 2*a**2*c**3*d*x**15/3 + 10*a**2*c**3*e*x**17/17 + 5*a*c**4*d*x**19/19 + 5*a*c**4*e*x**21/21 + c**5*d*x**23/23 + c**5*e*x**25/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_3():
    f = x*(a + c*x**4)**5*(d + e*x**2)
    F = a**5*d*x**2/2 + 5*a**4*c*d*x**6/6 + a**3*c**2*d*x**10 + 5*a**2*c**3*d*x**14/7 + 5*a*c**4*d*x**18/18 + c**5*d*x**22/22 + e*(a + c*x**4)**6/(24*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_4():
    f = (a + c*x**4)**5*(d + e*x**2)
    F = a**5*d*x + a**5*e*x**3/3 + a**4*c*d*x**5 + 5*a**4*c*e*x**7/7 + 10*a**3*c**2*d*x**9/9 + 10*a**3*c**2*e*x**11/11 + 10*a**2*c**3*d*x**13/13 + 2*a**2*c**3*e*x**15/3 + 5*a*c**4*d*x**17/17 + 5*a*c**4*e*x**19/19 + c**5*d*x**21/21 + c**5*e*x**23/23
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_5():
    f = (a + c*x**4)**5*(d + e*x**2)/x
    F = a**5*d*log(x) + a**5*e*x**2/2 + 5*a**4*c*d*x**4/4 + 5*a**4*c*e*x**6/6 + 5*a**3*c**2*d*x**8/4 + a**3*c**2*e*x**10 + 5*a**2*c**3*d*x**12/6 + 5*a**2*c**3*e*x**14/7 + 5*a*c**4*d*x**16/16 + 5*a*c**4*e*x**18/18 + c**5*d*x**20/20 + c**5*e*x**22/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_6():
    f = (a + c*x**4)**5*(d + e*x**2)/x**2
    F = -a**5*d/x + a**5*e*x + 5*a**4*c*d*x**3/3 + a**4*c*e*x**5 + 10*a**3*c**2*d*x**7/7 + 10*a**3*c**2*e*x**9/9 + 10*a**2*c**3*d*x**11/11 + 10*a**2*c**3*e*x**13/13 + a*c**4*d*x**15/3 + 5*a*c**4*e*x**17/17 + c**5*d*x**19/19 + c**5*e*x**21/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_7():
    f = (a + c*x**4)**5*(d + e*x**2)/x**3
    F = -a**5*d/(2*x**2) + a**5*e*log(x) + 5*a**4*c*d*x**2/2 + 5*a**4*c*e*x**4/4 + 5*a**3*c**2*d*x**6/3 + 5*a**3*c**2*e*x**8/4 + a**2*c**3*d*x**10 + 5*a**2*c**3*e*x**12/6 + 5*a*c**4*d*x**14/14 + 5*a*c**4*e*x**16/16 + c**5*d*x**18/18 + c**5*e*x**20/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_8():
    f = x**5*(3*x**2 + 2)*sqrt(x**4 + 5)
    F = 3*x**4*(x**4 + 5)**(sympy.S(3)/2)/10 - 5*x**2*sqrt(x**4 + 5)/8 - (4 - x**2)*(x**4 + 5)**(sympy.S(3)/2)/4 - 25*asinh(sqrt(5)*x**2/5)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_9():
    f = x**3*(3*x**2 + 2)*sqrt(x**4 + 5)
    F = -15*x**2*sqrt(x**4 + 5)/16 + (9*x**2 + 8)*(x**4 + 5)**(sympy.S(3)/2)/24 - 75*asinh(sqrt(5)*x**2/5)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_10():
    f = x*(3*x**2 + 2)*sqrt(x**4 + 5)
    F = x**2*sqrt(x**4 + 5)/2 + (x**4 + 5)**(sympy.S(3)/2)/2 + 5*asinh(sqrt(5)*x**2/5)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_11():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x
    F = (3*x**2 + 4)*sqrt(x**4 + 5)/4 + 15*asinh(sqrt(5)*x**2/5)/4 - sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_12():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x**3
    F = asinh(sqrt(5)*x**2/5) - 3*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/2 - (2 - 3*x**2)*sqrt(x**4 + 5)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_13():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x**5
    F = 3*asinh(sqrt(5)*x**2/5)/2 - sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/10 - (3*x**2 + 1)*sqrt(x**4 + 5)/(2*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_14():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x**7
    F = -3*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/20 - 3*sqrt(x**4 + 5)/(4*x**4) - (x**4 + 5)**(sympy.S(3)/2)/(15*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_15():
    f = x**4*(3*x**2 + 2)*sqrt(x**4 + 5)
    F = x**5*(7*x**2 + 6)*sqrt(x**4 + 5)/21 + 2*x**3*sqrt(x**4 + 5)/3 + 20*x*sqrt(x**4 + 5)/21 - 10*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) + 10*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) - 5*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*sqrt(5) + 21)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(21*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_16():
    f = x**2*(3*x**2 + 2)*sqrt(x**4 + 5)
    F = x**3*(15*x**2 + 14)*sqrt(x**4 + 5)/35 + 10*x*sqrt(x**4 + 5)/7 + 4*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 4*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(14 - 5*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(7*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_17():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)
    F = x*(9*x**2 + 10)*sqrt(x**4 + 5)/15 + 6*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 6*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*sqrt(5) + 9)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(3*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_18():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x**2
    F = 4*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 4*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 + sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) - (2 - x**2)*sqrt(x**4 + 5)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_19():
    f = (3*x**2 + 2)*sqrt(x**4 + 5)/x**4
    F = 6*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 6*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(3)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 + 9*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(15*sqrt(x**4 + 5)) - 6*sqrt(x**4 + 5)/x - (2 - 9*x**2)*sqrt(x**4 + 5)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_20():
    f = x**5*(3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = 3*x**4*(x**4 + 5)**(sympy.S(5)/2)/14 - 5*x**2*(x**4 + 5)**(sympy.S(3)/2)/24 - 25*x**2*sqrt(x**4 + 5)/16 - (sympy.S(3)/7 - x**2/6)*(x**4 + 5)**(sympy.S(5)/2) - 125*asinh(sqrt(5)*x**2/5)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_21():
    f = x**3*(3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = -5*x**2*(x**4 + 5)**(sympy.S(3)/2)/16 - 75*x**2*sqrt(x**4 + 5)/32 + (x**2/4 + sympy.S(1)/5)*(x**4 + 5)**(sympy.S(5)/2) - 375*asinh(sqrt(5)*x**2/5)/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_22():
    f = x*(3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = x**2*(x**4 + 5)**(sympy.S(3)/2)/4 + 15*x**2*sqrt(x**4 + 5)/8 + 3*(x**4 + 5)**(sympy.S(5)/2)/10 + 75*asinh(sqrt(5)*x**2/5)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_23():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x
    F = (3*x**2/8 + sympy.S(1)/3)*(x**4 + 5)**(sympy.S(3)/2) + (45*x**2/16 + 5)*sqrt(x**4 + 5) + 225*asinh(sqrt(5)*x**2/5)/16 - 5*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_24():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x**3
    F = (3*x**2/2 + sympy.S(15)/2)*sqrt(x**4 + 5) + 15*asinh(sqrt(5)*x**2/5)/2 - 15*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/2 - (2 - x**2)*(x**4 + 5)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_25():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x**5
    F = 45*asinh(sqrt(5)*x**2/5)/4 - 3*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/2 - (45 - 6*x**2)*sqrt(x**4 + 5)/(4*x**2) - (2 - 3*x**2)*(x**4 + 5)**(sympy.S(3)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_26():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x**7
    F = asinh(sqrt(5)*x**2/5) - 9*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/4 - (4 - 9*x**2)*sqrt(x**4 + 5)/(4*x**2) - (9*x**2 + 4)*(x**4 + 5)**(sympy.S(3)/2)/(12*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_27():
    f = x**4*(3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = x**5*(33*x**2 + 26)*(x**4 + 5)**(sympy.S(3)/2)/143 + 10*x**5*(77*x**2 + 78)*sqrt(x**4 + 5)/1001 + 20*x**3*sqrt(x**4 + 5)/13 + 200*x*sqrt(x**4 + 5)/77 - 300*x*sqrt(x**4 + 5)/(13*x**2 + 13*sqrt(5)) + 300*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(13*sqrt(x**4 + 5)) - 50*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(26*sqrt(5) + 231)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(1001*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_28():
    f = x**2*(3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = x**3*(27*x**2 + 22)*(x**4 + 5)**(sympy.S(3)/2)/99 + 2*x**3*(135*x**2 + 154)*sqrt(x**4 + 5)/231 + 300*x*sqrt(x**4 + 5)/77 + 40*x*sqrt(x**4 + 5)/(3*x**2 + 3*sqrt(5)) - 40*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(3*sqrt(x**4 + 5)) + 10*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(154 - 45*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(231*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_29():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)
    F = x*(7*x**2 + 6)*(x**4 + 5)**(sympy.S(3)/2)/21 + 2*x*(7*x**2 + 10)*sqrt(x**4 + 5)/7 + 20*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 20*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 10*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*sqrt(5) + 7)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(7*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_30():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x**2
    F = 6*x*(14*x**2 + 25)*sqrt(x**4 + 5)/35 + 24*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 24*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 6*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(5*sqrt(5) + 14)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(7*sqrt(x**4 + 5)) - (14 - 3*x**2)*(x**4 + 5)**(sympy.S(3)/2)/(7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_31():
    f = (3*x**2 + 2)*(x**4 + 5)**(sympy.S(3)/2)/x**4
    F = 36*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 36*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 2*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*sqrt(5) + 27)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(3*sqrt(x**4 + 5)) - (54 - 4*x**2)*sqrt(x**4 + 5)/(3*x) - (10 - 9*x**2)*(x**4 + 5)**(sympy.S(3)/2)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_32():
    f = x**7*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = 3*x**6*sqrt(x**4 + 5)/8 + x**4*sqrt(x**4 + 5)/3 - (45*x**2/16 + sympy.S(10)/3)*sqrt(x**4 + 5) + 225*asinh(sqrt(5)*x**2/5)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_33():
    f = x**5*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = x**4*sqrt(x**4 + 5)/2 - (10 - x**2)*sqrt(x**4 + 5)/2 - 5*asinh(sqrt(5)*x**2/5)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_34():
    f = x**3*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = (3*x**2 + 4)*sqrt(x**4 + 5)/4 - 15*asinh(sqrt(5)*x**2/5)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_35():
    f = x*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = 3*sqrt(x**4 + 5)/2 + asinh(sqrt(5)*x**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_36():
    f = (3*x**2 + 2)/(x*sqrt(x**4 + 5))
    F = 3*asinh(sqrt(5)*x**2/5)/2 - sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_37():
    f = (3*x**2 + 2)/(x**3*sqrt(x**4 + 5))
    F = -3*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/10 - sqrt(x**4 + 5)/(5*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_38():
    f = (3*x**2 + 2)/(x**5*sqrt(x**4 + 5))
    F = sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/50 - 3*sqrt(x**4 + 5)/(10*x**2) - sqrt(x**4 + 5)/(10*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_39():
    f = x**4*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = 3*x**3*sqrt(x**4 + 5)/5 + 2*x*sqrt(x**4 + 5)/3 - 9*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) + 9*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*sqrt(5) + 27)*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(6*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_40():
    f = x**2*(3*x**2 + 2)/sqrt(x**4 + 5)
    F = x*sqrt(x**4 + 5) + 2*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 2*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 - sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(2*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_41():
    f = (3*x**2 + 2)/sqrt(x**4 + 5)
    F = 3*x*sqrt(x**4 + 5)/(x**2 + sqrt(5)) - 3*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/sqrt(x**4 + 5) + 5**(sympy.S(3)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 + 3*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(10*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_42():
    f = (3*x**2 + 2)/(x**2*sqrt(x**4 + 5))
    F = 2*x*sqrt(x**4 + 5)/(5*x**2 + 5*sqrt(5)) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 + 3*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(10*sqrt(x**4 + 5)) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2*x**2 + 2*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(5*sqrt(x**4 + 5)) - 2*sqrt(x**4 + 5)/(5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_43():
    f = (3*x**2 + 2)/(x**4*sqrt(x**4 + 5))
    F = 3*x*sqrt(x**4 + 5)/(5*x**2 + 5*sqrt(5)) - 5**(sympy.S(3)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 - 9*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(150*sqrt(x**4 + 5)) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(3*x**2 + 3*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(5*sqrt(x**4 + 5)) - 3*sqrt(x**4 + 5)/(5*x) - 2*sqrt(x**4 + 5)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_44():
    f = x**7*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = -x**4*(3*x**2 + 2)/(2*sqrt(x**4 + 5)) + (9*x**2/4 + 2)*sqrt(x**4 + 5) - 45*asinh(sqrt(5)*x**2/5)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_45():
    f = x**5*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = -x**2*(3*x**2 + 2)/(2*sqrt(x**4 + 5)) + 3*sqrt(x**4 + 5) + asinh(sqrt(5)*x**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_46():
    f = x**3*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = (-3*x**2 - 2)/(2*sqrt(x**4 + 5)) + 3*asinh(sqrt(5)*x**2/5)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_47():
    f = x*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = (2*x**2 - 15)/(10*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_48():
    f = (3*x**2 + 2)/(x*(x**4 + 5)**(sympy.S(3)/2))
    F = (3*x**2 + 2)/(10*sqrt(x**4 + 5)) - sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_49():
    f = (3*x**2 + 2)/(x**3*(x**4 + 5)**(sympy.S(3)/2))
    F = -3*sqrt(5)*atanh(sqrt(5)*sqrt(x**4 + 5)/5)/50 + (3*x**2 + 2)/(10*x**2*sqrt(x**4 + 5)) - 2*sqrt(x**4 + 5)/(25*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_50():
    f = x**4*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = -x**3*(15 - 2*x**2)/(10*sqrt(x**4 + 5)) - x*sqrt(x**4 + 5)/5 + 9*x*sqrt(x**4 + 5)/(2*x**2 + 2*sqrt(5)) - 9*5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(2*sqrt(x**4 + 5)) + 5**(sympy.S(3)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 + 9*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(20*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_51():
    f = x**2*(3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = -x*(15 - 2*x**2)/(10*sqrt(x**4 + 5)) - x*sqrt(x**4 + 5)/(5*x**2 + 5*sqrt(5)) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(x**2 + sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(5*sqrt(x**4 + 5)) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 - 3*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(20*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_52():
    f = (3*x**2 + 2)/(x**4 + 5)**(sympy.S(3)/2)
    F = x*(3*x**2 + 2)/(10*sqrt(x**4 + 5)) - 3*x*sqrt(x**4 + 5)/(10*x**2 + 10*sqrt(5)) + 5**(sympy.S(3)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(2 - 3*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(100*sqrt(x**4 + 5)) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(3*x**2 + 3*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(10*sqrt(x**4 + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_53():
    f = (3*x**2 + 2)/(x**2*(x**4 + 5)**(sympy.S(3)/2))
    F = 3*x*sqrt(x**4 + 5)/(25*x**2 + 25*sqrt(5)) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(6 + 3*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(100*sqrt(x**4 + 5)) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(3*x**2 + 3*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(25*sqrt(x**4 + 5)) + (3*x**2 + 2)/(10*x*sqrt(x**4 + 5)) - 3*sqrt(x**4 + 5)/(25*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_54():
    f = (3*x**2 + 2)/(x**4*(x**4 + 5)**(sympy.S(3)/2))
    F = 9*x*sqrt(x**4 + 5)/(50*x**2 + 50*sqrt(5)) + 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(27 - 2*sqrt(5))*(x**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(300*sqrt(x**4 + 5)) - 5**(sympy.S(1)/4)*sqrt((x**4 + 5)/(x**2 + sqrt(5))**2)*(9*x**2 + 9*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*x/5), sympy.S.Half)/(50*sqrt(x**4 + 5)) - 9*sqrt(x**4 + 5)/(50*x) + (3*x**2 + 2)/(10*x**3*sqrt(x**4 + 5)) - sqrt(x**4 + 5)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_55():
    f = (f*x)**m*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = d*(f*x)**(m + 1)/(f*(m + 1)) + e*(f*x)**(m + 23)/(f**23*(m + 23)) + (f*x)**(m + 3)*(10*d + e)/(f**3*(m + 3)) + (f*x)**(m + 5)*(45*d + 10*e)/(f**5*(m + 5)) + (f*x)**(m + 7)*(120*d + 45*e)/(f**7*(m + 7)) + (f*x)**(m + 9)*(210*d + 120*e)/(f**9*(m + 9)) + (f*x)**(m + 11)*(252*d + 210*e)/(f**11*(m + 11)) + (f*x)**(m + 13)*(210*d + 252*e)/(f**13*(m + 13)) + (f*x)**(m + 15)*(120*d + 210*e)/(f**15*(m + 15)) + (f*x)**(m + 17)*(45*d + 120*e)/(f**17*(m + 17)) + (f*x)**(m + 19)*(10*d + 45*e)/(f**19*(m + 19)) + (f*x)**(m + 21)*(d + 10*e)/(f**21*(m + 21))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_56():
    f = x**5*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = e*(x**2 + 1)**14/28 + (d/26 - 3*e/26)*(x**2 + 1)**13 + (d/22 - e/22)*(x**2 + 1)**11 - (d/12 - e/8)*(x**2 + 1)**12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_57():
    f = x**4*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = d*x**5/5 + e*x**27/27 + x**25*(d/25 + 2*e/5) + x**23*(10*d/23 + 45*e/23) + x**21*(15*d/7 + 40*e/7) + x**19*(120*d/19 + 210*e/19) + x**17*(210*d/17 + 252*e/17) + x**15*(84*d/5 + 14*e) + x**13*(210*d/13 + 120*e/13) + x**11*(120*d/11 + 45*e/11) + x**9*(5*d + 10*e/9) + x**7*(10*d/7 + e/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_58():
    f = x**3*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = e*(x**2 + 1)**13/26 + (-d/22 + e/22)*(x**2 + 1)**11 + (d/24 - e/12)*(x**2 + 1)**12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_59():
    f = x**2*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = d*x**3/3 + e*x**25/25 + x**23*(d/23 + 10*e/23) + x**21*(10*d/21 + 15*e/7) + x**19*(45*d/19 + 120*e/19) + x**17*(120*d/17 + 210*e/17) + x**15*(14*d + 84*e/5) + x**13*(252*d/13 + 210*e/13) + x**11*(210*d/11 + 120*e/11) + x**9*(40*d/3 + 5*e) + x**7*(45*d/7 + 10*e/7) + x**5*(2*d + e/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_60():
    f = x*(d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = e*(x**2 + 1)**12/24 + (d/22 - e/22)*(x**2 + 1)**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_61():
    f = (d + e*x**2)*(x**4 + 2*x**2 + 1)**5
    F = d*x + e*x**23/23 + x**21*(d/21 + 10*e/21) + x**19*(10*d/19 + 45*e/19) + x**17*(45*d/17 + 120*e/17) + x**15*(8*d + 14*e) + x**13*(210*d/13 + 252*e/13) + x**11*(252*d/11 + 210*e/11) + x**9*(70*d/3 + 40*e/3) + x**7*(120*d/7 + 45*e/7) + x**5*(9*d + 2*e) + x**3*(10*d/3 + e/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_62():
    f = (d + e*x**2)*(x**4 + 2*x**2 + 1)**5/x
    F = d*x**20/20 + 5*d*x**18/9 + 45*d*x**16/16 + 60*d*x**14/7 + 35*d*x**12/2 + 126*d*x**10/5 + 105*d*x**8/4 + 20*d*x**6 + 45*d*x**4/4 + 5*d*x**2 + d*log(x) + e*(x**2 + 1)**11/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_63():
    f = (d + e*x**2)*(x**4 + 2*x**2 + 1)**5/x**2
    F = -d/x + e*x**21/21 + x**19*(d/19 + 10*e/19) + x**17*(10*d/17 + 45*e/17) + x**15*(3*d + 8*e) + x**13*(120*d/13 + 210*e/13) + x**11*(210*d/11 + 252*e/11) + x**9*(28*d + 70*e/3) + x**7*(30*d + 120*e/7) + x**5*(24*d + 9*e) + x**3*(15*d + 10*e/3) + x*(10*d + e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_64():
    f = (d + e*x**2)*(x**4 + 2*x**2 + 1)**5/x**3
    F = -d/(2*x**2) + e*x**20/20 + x**18*(d/18 + 5*e/9) + x**16*(5*d/8 + 45*e/16) + x**14*(45*d/14 + 60*e/7) + x**12*(10*d + 35*e/2) + x**10*(21*d + 126*e/5) + x**8*(63*d/2 + 105*e/4) + x**6*(35*d + 20*e) + x**4*(30*d + 45*e/4) + x**2*(45*d/2 + 5*e) + (10*d + e)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_65():
    f = (f*x)**m*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = (f*x)**(m + 1)/(f*(m + 1)) + 11*(f*x)**(m + 3)/(f**3*(m + 3)) + 55*(f*x)**(m + 5)/(f**5*(m + 5)) + 165*(f*x)**(m + 7)/(f**7*(m + 7)) + 330*(f*x)**(m + 9)/(f**9*(m + 9)) + 462*(f*x)**(m + 11)/(f**11*(m + 11)) + 462*(f*x)**(m + 13)/(f**13*(m + 13)) + 330*(f*x)**(m + 15)/(f**15*(m + 15)) + 165*(f*x)**(m + 17)/(f**17*(m + 17)) + 55*(f*x)**(m + 19)/(f**19*(m + 19)) + 11*(f*x)**(m + 21)/(f**21*(m + 21)) + (f*x)**(m + 23)/(f**23*(m + 23))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_66():
    f = x**5*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = (x**2 + 1)**14/28 - (x**2 + 1)**13/13 + (x**2 + 1)**12/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_67():
    f = x**4*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = x**27/27 + 11*x**25/25 + 55*x**23/23 + 55*x**21/7 + 330*x**19/19 + 462*x**17/17 + 154*x**15/5 + 330*x**13/13 + 15*x**11 + 55*x**9/9 + 11*x**7/7 + x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_68():
    f = x**3*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = (x**2 + 1)**13/26 - (x**2 + 1)**12/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_69():
    f = x**2*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = x**25/25 + 11*x**23/23 + 55*x**21/21 + 165*x**19/19 + 330*x**17/17 + 154*x**15/5 + 462*x**13/13 + 30*x**11 + 55*x**9/3 + 55*x**7/7 + 11*x**5/5 + x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_70():
    f = x*(x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = (x**2 + 1)**12/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_71():
    f = (x**2 + 1)*(x**4 + 2*x**2 + 1)**5
    F = x**23/23 + 11*x**21/21 + 55*x**19/19 + 165*x**17/17 + 22*x**15 + 462*x**13/13 + 42*x**11 + 110*x**9/3 + 165*x**7/7 + 11*x**5 + 11*x**3/3 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_72():
    f = (x**2 + 1)*(x**4 + 2*x**2 + 1)**5/x
    F = x**22/22 + 11*x**20/20 + 55*x**18/18 + 165*x**16/16 + 165*x**14/7 + 77*x**12/2 + 231*x**10/5 + 165*x**8/4 + 55*x**6/2 + 55*x**4/4 + 11*x**2/2 + log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_73():
    f = (x**2 + 1)*(x**4 + 2*x**2 + 1)**5/x**2
    F = x**21/21 + 11*x**19/19 + 55*x**17/17 + 11*x**15 + 330*x**13/13 + 42*x**11 + 154*x**9/3 + 330*x**7/7 + 33*x**5 + 55*x**3/3 + 11*x - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_74():
    f = (x**2 + 1)*(x**4 + 2*x**2 + 1)**5/x**3
    F = x**20/20 + 11*x**18/18 + 55*x**16/16 + 165*x**14/14 + 55*x**12/2 + 231*x**10/5 + 231*x**8/4 + 55*x**6 + 165*x**4/4 + 55*x**2/2 + 11*log(x) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_75():
    f = x**2*(d + e*x**2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = -sqrt(a)*(a + b*x**2)*(-a*e + b*d)*atan(sqrt(b)*x/sqrt(a))/(b**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + e*x**3*(a + b*x**2)/(3*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x*(a + b*x**2)*(-a*e + b*d)/(b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_76():
    f = x*(d + e*x**2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = e*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*b**2) + (a + b*x**2)*(-a*e + b*d)*log(a + b*x**2)/(2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_77():
    f = (d + e*x**2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = e*x*(a + b*x**2)/(b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*(-a*e + b*d)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_78():
    f = (d + e*x**2)/(x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = d*(a + b*x**2)*log(x)/(a*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*(-a*e + b*d)*log(a + b*x**2)/(2*a*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_79():
    f = (d + e*x**2)/(x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -d*(a + b*x**2)/(a*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*(-a*e + b*d)*atan(sqrt(b)*x/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(b)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_80():
    f = (d + e*x**2)/(x**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    F = -d*(a + b*x**2)/(2*a*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*(-a*e + b*d)*log(x)/(a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*(-a*e + b*d)*log(a + b*x**2)/(2*a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_81():
    f = x**2*(d + e*x**2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -x*(-a*e + b*d)/(4*b**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x*(-5*a*e + b*d)/(8*a*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*(3*a*e + b*d)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(3)/2)*b**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_82():
    f = x*(d + e*x**2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = -e/(2*b**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (-a*e + b*d)/(4*b**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_83():
    f = (d + e*x**2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = x*(-a*e + b*d)/(4*a*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + x*(a*e + 3*b*d)/(8*a**2*b*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*(a*e + 3*b*d)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_84():
    f = (d + e*x**2)/(x*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = (-a*e + b*d)/(4*a*b*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + d/(2*a**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + d*(a + b*x**2)*log(x)/(a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(a + b*x**2)*log(a + b*x**2)/(2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_85():
    f = (d + e*x**2)/(x**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = -x*(-a*e + b*d)/(4*a**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(a + b*x**2)/(a**3*x*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - x*(-3*a*e + 7*b*d)/(8*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*(-3*a*e + 15*b*d)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(7)/2)*sqrt(b)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_86():
    f = (d + e*x**2)/(x**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2))
    F = -(-a*e + b*d)/(4*a**2*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - d*(a + b*x**2)/(2*a**3*x**2*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (-a*e + 2*b*d)/(2*a**3*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) - (a + b*x**2)*(-a*e + 3*b*d)*log(x)/(a**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (a + b*x**2)*(-a*e + 3*b*d)*log(a + b*x**2)/(2*a**4*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_87():
    f = (f*x)**m*(d + e*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(5)/2)
    F = a**5*d*(f*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f*(a + b*x**2)*(m + 1)) + a**4*(f*x)**(m + 3)*(a*e + 5*b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**3*(a + b*x**2)*(m + 3)) + 5*a**3*b*(f*x)**(m + 5)*(a*e + 2*b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**5*(a + b*x**2)*(m + 5)) + 10*a**2*b**2*(f*x)**(m + 7)*(a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**7*(a + b*x**2)*(m + 7)) + 5*a*b**3*(f*x)**(m + 9)*(2*a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**9*(a + b*x**2)*(m + 9)) + b**5*e*(f*x)**(m + 13)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**13*(a + b*x**2)*(m + 13)) + b**4*(f*x)**(m + 11)*(5*a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**11*(a + b*x**2)*(m + 11))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_88():
    f = (f*x)**m*(d + e*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = a**3*d*(f*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f*(a + b*x**2)*(m + 1)) + a**2*(f*x)**(m + 3)*(a*e + 3*b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**3*(a + b*x**2)*(m + 3)) + 3*a*b*(f*x)**(m + 5)*(a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**5*(a + b*x**2)*(m + 5)) + b**3*e*(f*x)**(m + 9)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**9*(a + b*x**2)*(m + 9)) + b**2*(f*x)**(m + 7)*(3*a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**7*(a + b*x**2)*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_89():
    f = (f*x)**m*(d + e*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = a*d*(f*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f*(a + b*x**2)*(m + 1)) + b*e*(f*x)**(m + 5)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**5*(a + b*x**2)*(m + 5)) + (f*x)**(m + 3)*(a*e + b*d)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(f**3*(a + b*x**2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_90():
    f = (f*x)**m*(d + e*x**2)/sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = e*(f*x)**(m + 1)*(a + b*x**2)/(b*f*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (f*x)**(m + 1)*(a + b*x**2)*(-a*e + b*d)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b*f*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_91():
    f = (f*x)**m*(d + e*x**2)/(a**2 + 2*a*b*x**2 + b**2*x**4)**(sympy.S(3)/2)
    F = (f*x)**(m + 1)*(-a*e + b*d)/(4*a*b*f*(a + b*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)) + (f*x)**(m + 1)*(a + b*x**2)*(a*e*(m + 1) + b*d*(3 - m))*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(4*a**3*b*f*(m + 1)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_92():
    f = x*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = (a**2 + 2*a*b*x**2 + b**2*x**4)**(p + 1)/(4*b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_93():
    f = x**3*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = -a*(a + b*x**2)**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**2*(p + 1)) + (a + b*x**2)**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(2*b**2*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_94():
    f = x**5*(a + b*x**2)*(a**2 + 2*a*b*x**2 + b**2*x**4)**p
    F = a**2*(a + b*x**2)**2*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**3*(p + 1)) - a*(a + b*x**2)**3*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(b**3*(2*p + 3)) + (a + b*x**2)**4*(a**2 + 2*a*b*x**2 + b**2*x**4)**p/(4*b**3*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_95():
    f = x**3*(A + B*x**2)*(a + b*x**2 + c*x**4)**3
    F = A*a**3*x**4/4 + B*c**3*x**18/18 + a**2*x**6*(3*A*b + B*a)/6 + 3*a*x**8*(A*(a*c + b**2) + B*a*b)/8 + c**2*x**16*(A*c + 3*B*b)/16 + 3*c*x**14*(A*b*c + B*a*c + B*b**2)/14 + x**12*(A*a*c**2/4 + A*b**2*c/4 + B*a*b*c/2 + B*b**3/12) + x**10*(A*(6*a*b*c + b**3)/10 + 3*B*a*(a*c + b**2)/10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_96():
    f = x**2*(A + B*x**2)*(a + b*x**2 + c*x**4)**3
    F = A*a**3*x**3/3 + B*c**3*x**17/17 + a**2*x**5*(3*A*b + B*a)/5 + 3*a*x**7*(A*(a*c + b**2) + B*a*b)/7 + c**2*x**15*(A*c + 3*B*b)/15 + 3*c*x**13*(A*b*c + B*a*c + B*b**2)/13 + x**11*(3*A*a*c**2/11 + 3*A*b**2*c/11 + 6*B*a*b*c/11 + B*b**3/11) + x**9*(A*(6*a*b*c + b**3)/9 + B*a*(a*c + b**2)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_97():
    f = x*(A + B*x**2)*(a + b*x**2 + c*x**4)**3
    F = A*a**3*x**2/2 + B*c**3*x**16/16 + a**2*x**4*(3*A*b + B*a)/4 + a*x**6*(A*(a*c + b**2) + B*a*b)/2 + c**2*x**14*(A*c + 3*B*b)/14 + c*x**12*(A*b*c + B*a*c + B*b**2)/4 + x**10*(3*A*a*c**2/10 + 3*A*b**2*c/10 + 3*B*a*b*c/5 + B*b**3/10) + x**8*(A*(6*a*b*c + b**3)/8 + 3*B*a*(a*c + b**2)/8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_98():
    f = (A + B*x**2)*(a + b*x**2 + c*x**4)**3
    F = A*a**3*x + B*c**3*x**15/15 + a**2*x**3*(3*A*b + B*a)/3 + 3*a*x**5*(A*(a*c + b**2) + B*a*b)/5 + c**2*x**13*(A*c + 3*B*b)/13 + 3*c*x**11*(A*b*c + B*a*c + B*b**2)/11 + x**9*(A*a*c**2/3 + A*b**2*c/3 + 2*B*a*b*c/3 + B*b**3/9) + x**7*(A*(6*a*b*c + b**3)/7 + 3*B*a*(a*c + b**2)/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_99():
    f = (A + B*x**2)*(a + b*x**2 + c*x**4)**3/x
    F = A*a**3*log(x) + B*c**3*x**14/14 + a**2*x**2*(3*A*b + B*a)/2 + 3*a*x**4*(A*(a*c + b**2) + B*a*b)/4 + c**2*x**12*(A*c + 3*B*b)/12 + 3*c*x**10*(A*b*c + B*a*c + B*b**2)/10 + x**8*(3*A*a*c**2/8 + 3*A*b**2*c/8 + 3*B*a*b*c/4 + B*b**3/8) + x**6*(A*(6*a*b*c + b**3)/6 + B*a*(a*c + b**2)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_100():
    f = (A + B*x**2)*(a + b*x**2 + c*x**4)**3/x**2
    F = -A*a**3/x + B*c**3*x**13/13 + a**2*x*(3*A*b + B*a) + a*x**3*(A*(a*c + b**2) + B*a*b) + c**2*x**11*(A*c + 3*B*b)/11 + c*x**9*(A*b*c + B*a*c + B*b**2)/3 + x**7*(3*A*a*c**2/7 + 3*A*b**2*c/7 + 6*B*a*b*c/7 + B*b**3/7) + x**5*(A*(6*a*b*c + b**3)/5 + 3*B*a*(a*c + b**2)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_101():
    f = (A + B*x**2)*(a + b*x**2 + c*x**4)**3/x**3
    F = -A*a**3/(2*x**2) + B*c**3*x**12/12 + a**2*(3*A*b + B*a)*log(x) + 3*a*x**2*(A*(a*c + b**2) + B*a*b)/2 + c**2*x**10*(A*c + 3*B*b)/10 + 3*c*x**8*(A*b*c + B*a*c + B*b**2)/8 + x**6*(A*a*c**2/2 + A*b**2*c/2 + B*a*b*c + B*b**3/6) + x**4*(A*(6*a*b*c + b**3)/4 + 3*B*a*(a*c + b**2)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_102():
    f = x**5*(A + B*x**2)/(a + b*x**2 + c*x**4)
    F = B*x**4/(4*c) - x**2*(-A*c + B*b)/(2*c**2) + (-A*b*c - B*a*c + B*b**2)*log(a + b*x**2 + c*x**4)/(4*c**3) + (2*A*a*c**2 - A*b**2*c - 3*B*a*b*c + B*b**3)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_103():
    f = x**3*(A + B*x**2)/(a + b*x**2 + c*x**4)
    F = B*x**2/(2*c) - (-A*c + B*b)*log(a + b*x**2 + c*x**4)/(4*c**2) - (-A*b*c - 2*B*a*c + B*b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_104():
    f = x*(A + B*x**2)/(a + b*x**2 + c*x**4)
    F = B*log(a + b*x**2 + c*x**4)/(4*c) + (-2*A*c + B*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_105():
    f = (A + B*x**2)/(x*(a + b*x**2 + c*x**4))
    F = A*log(x)/a - A*log(a + b*x**2 + c*x**4)/(4*a) + (A*b - 2*B*a)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_106():
    f = (A + B*x**2)/(x**3*(a + b*x**2 + c*x**4))
    F = -A/(2*a*x**2) - (A*b - B*a)*log(x)/a**2 + (A*b - B*a)*log(a + b*x**2 + c*x**4)/(4*a**2) - (-2*A*a*c + A*b**2 - B*a*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_107():
    f = x**4*(A + B*x**2)/(a + b*x**2 + c*x**4)
    F = B*x**3/(3*c) - x*(-A*c + B*b)/c**2 + sqrt(2)*(-A*b*c - B*a*c + B*b**2 + (2*A*a*c**2 - A*b**2*c - 3*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-A*b*c - B*a*c + B*b**2 - (2*A*a*c**2 - A*b**2*c - 3*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_108():
    f = x**2*(A + B*x**2)/(a + b*x**2 + c*x**4)
    F = B*x/c - sqrt(2)*(-A*c + B*b + (-A*b*c - 2*B*a*c + B*b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*(-A*c + B*b - (-A*b*c - 2*B*a*c + B*b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_109():
    f = (A + B*x**2)/(a + b*x**2 + c*x**4)
    F = sqrt(2)*(B - (-2*A*c + B*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))) + sqrt(2)*(B + (-2*A*c + B*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_110():
    f = (A + B*x**2)/(x**2*(a + b*x**2 + c*x**4))
    F = -A/(a*x) - sqrt(2)*sqrt(c)*(A - (A*b - 2*B*a)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(A + (A*b - 2*B*a)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_111():
    f = (A + B*x**2)/(x**4*(a + b*x**2 + c*x**4))
    F = -A/(3*a*x**3) + sqrt(2)*sqrt(c)*(-A*(-2*a*c + b**2 - b*sqrt(-4*a*c + b**2)) + B*a*(b - sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(-A*(-2*a*c + b**2 + b*sqrt(-4*a*c + b**2)) + B*a*(b + sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + (A*b - B*a)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_112():
    f = x**7*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -x**4*(a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + x**2*(-A*b*c - 6*B*a*c + 2*B*b**2)/(2*c**2*(-4*a*c + b**2)) - (-A*c + 2*B*b)*log(a + b*x**2 + c*x**4)/(4*c**3) - (6*A*a*b*c**2 - A*b**3*c + 12*B*a**2*c**2 - 12*B*a*b**2*c + 2*B*b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_113():
    f = x**5*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = B*log(a + b*x**2 + c*x**4)/(4*c**2) - x**2*(a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (4*A*a*c**2 - 6*B*a*b*c + B*b**3)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_114():
    f = x**3*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -(A*b - 2*B*a)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_115():
    f = x*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -(-2*A*c + B*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_116():
    f = (A + B*x**2)/(x*(a + b*x**2 + c*x**4)**2)
    F = A*log(x)/a**2 - A*log(a + b*x**2 + c*x**4)/(4*a**2) - (-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (A*(-6*a*b*c + b**3) + 4*B*a**2*c)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_117():
    f = (A + B*x**2)/(x**3*(a + b*x**2 + c*x**4)**2)
    F = -(-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(2*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-6*A*a*c + 2*A*b**2 - B*a*b)/(2*a**2*x**2*(-4*a*c + b**2)) - (2*A*b - B*a)*log(x)/a**3 + (2*A*b - B*a)*log(a + b*x**2 + c*x**4)/(4*a**3) + (-2*A*(6*a**2*c**2 - 6*a*b**2*c + b**4) + B*a*b*(-6*a*c + b**2))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_118():
    f = x**6*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -x**5*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x**3*(-2*A*c + B*b)/(2*c*(-4*a*c + b**2)) + x*(-A*b*c - 10*B*a*c + 3*B*b**2)/(2*c**2*(-4*a*c + b**2)) - sqrt(2)*(6*A*a*c**2 - A*b**2*c - 13*B*a*b*c + 3*B*b**3 + (8*A*a*b*c**2 - A*b**3*c + 20*B*a**2*c**2 - 19*B*a*b**2*c + 3*B*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(6*A*a*c**2 - A*b**2*c - 13*B*a*b*c + 3*B*b**3 - (8*A*a*b*c**2 - A*b**3*c + 20*B*a**2*c**2 - 19*B*a*b**2*c + 3*B*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_119():
    f = x**4*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -x**3*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - x*(-2*A*c + B*b)/(2*c*(-4*a*c + b**2)) + sqrt(2)*(A*b*c - 6*B*a*c + B*b**2 + (4*A*a*c**2 + A*b**2*c - 8*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(A*b*c - 6*B*a*c + B*b**2 - (4*A*a*c**2 + A*b**2*c - 8*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_120():
    f = x**2*(A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = -x*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(-2*A*c + B*b + (-4*A*b*c + 4*B*a*c + B*b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(-2*A*c + B*b - (-4*A*b*c + 4*B*a*c + B*b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_121():
    f = (A + B*x**2)/(a + b*x**2 + c*x**4)**2
    F = sqrt(2)*sqrt(c)*(A*b - 2*B*a - (-12*A*a*c + A*b**2 + 4*B*a*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(A*b - 2*B*a + (A*(-12*a*c + b**2) + 4*B*a*b)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-2*A*a*c + A*b**2 - B*a*b + c*x**2*(A*b - 2*B*a))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_122():
    f = (A + B*x**2)/(x**2*(a + b*x**2 + c*x**4)**2)
    F = -(-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(2*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*sqrt(c)*(-10*A*a*c + 3*A*b**2 - B*a*b + (-A*(-16*a*b*c + 3*b**3) + B*a*(-12*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-A*(-16*a*b*c - 10*a*c*sqrt(-4*a*c + b**2) + 3*b**3 + 3*b**2*sqrt(-4*a*c + b**2)) + B*a*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*A*a*c + 3*A*b**2 - B*a*b)/(2*a**2*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_123():
    f = (A + B*x**2)/(x**4*(a + b*x**2 + c*x**4)**2)
    F = -(-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(2*a*x**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-14*A*a*c + 5*A*b**2 - 3*B*a*b)/(6*a**2*x**3*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-A*(28*a**2*c**2 - 29*a*b**2*c + 19*a*b*c*sqrt(-4*a*c + b**2) + 5*b**4 - 5*b**3*sqrt(-4*a*c + b**2)) + B*a*(-16*a*b*c + 10*a*c*sqrt(-4*a*c + b**2) + 3*b**3 - 3*b**2*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-A*(28*a**2*c**2 - 29*a*b**2*c - 19*a*b*c*sqrt(-4*a*c + b**2) + 5*b**4 + 5*b**3*sqrt(-4*a*c + b**2)) + B*a*(-16*a*b*c - 10*a*c*sqrt(-4*a*c + b**2) + 3*b**3 + 3*b**2*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-A*(-19*a*b*c + 5*b**3) + B*a*(-10*a*c + 3*b**2))/(2*a**3*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_124():
    f = x**11*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x**8*(a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(4*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - x**4*(a*(16*A*a*c**2 - A*b**2*c - 18*B*a*b*c + 3*B*b**3) + x**2*(10*A*a*b*c**2 - A*b**3*c + 20*B*a**2*c**2 - 20*B*a*b**2*c + 3*B*b**4))/(4*c**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + x**2*(7*A*a*b*c**2 - A*b**3*c + 30*B*a**2*c**2 - 21*B*a*b**2*c + 3*B*b**4)/(2*c**3*(-4*a*c + b**2)**2) - (-A*c + 3*B*b)*log(a + b*x**2 + c*x**4)/(4*c**4) - (-30*A*a**2*b*c**3 + 10*A*a*b**3*c**2 - A*b**5*c - 60*B*a**3*c**3 + 90*B*a**2*b**2*c**2 - 30*B*a*b**4*c + 3*B*b**6)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**4*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_125():
    f = x**9*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = B*log(a + b*x**2 + c*x**4)/(4*c**3) - x**6*(a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(4*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - x**2*(2*a*(6*A*a*c**2 - 7*B*a*b*c + B*b**3) + x**2*(6*A*a*b*c**2 + 16*B*a**2*c**2 - 15*B*a*b**2*c + 2*B*b**4))/(4*c**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + (-12*A*a**2*c**3 + 30*B*a**2*b*c**2 - 10*B*a*b**3*c + B*b**5)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_126():
    f = x**7*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = 3*a*(A*b - 2*B*a)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - x**6*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x**2*(2*a + b*x**2)*(3*A*b - 6*B*a)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_127():
    f = x**5*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x**4*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + (-A*(2*a*c + b**2) + 3*B*a*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - (a*(-6*A*b*c + 8*B*a*c + B*b**2) + x**2*(4*A*a*c**2 - 4*A*b**2*c + 2*B*a*b*c + B*b**3))/(4*c*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_128():
    f = x**3*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -(-3*A*b*c + 2*B*a*c + B*b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + (b + 2*c*x**2)*(-3*A*b*c + 2*B*a*c + B*b**2)/(4*c*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (a*(-2*A*c + B*b) + x**2*(-A*b*c - 2*B*a*c + B*b**2))/(4*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_129():
    f = x*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = 3*c*(-2*A*c + B*b)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - (b + 2*c*x**2)*(-6*A*c + 3*B*b)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_130():
    f = (A + B*x**2)/(x*(a + b*x**2 + c*x**4)**3)
    F = A*log(x)/a**3 - A*log(a + b*x**2 + c*x**4)/(4*a**3) - (-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + (A*(16*a**2*c**2 - 15*a*b**2*c + 2*b**4) + 6*B*a**2*b*c + 2*c*x**2*(A*(-7*a*b*c + b**3) + 6*B*a**2*c))/(4*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-A*(30*a**2*b*c**2 - 10*a*b**3*c + b**5) + 12*B*a**3*c**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_131():
    f = (A + B*x**2)/(x**3*(a + b*x**2 + c*x**4)**3)
    F = -(-A*(-2*a*c + b**2) + B*a*b - c*x**2*(A*b - 2*B*a))/(4*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - (-A*(20*a**2*c**2 - 20*a*b**2*c + 3*b**4) + B*a*b*(-10*a*c + b**2) + c*x**2*(-3*A*(-6*a*b*c + b**3) + B*a*(-16*a*c + b**2)))/(4*a**2*x**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + (-3*A*(10*a**2*c**2 - 7*a*b**2*c + b**4) + B*a*b*(-7*a*c + b**2))/(2*a**3*x**2*(-4*a*c + b**2)**2) - (3*A*b - B*a)*log(x)/a**4 + (3*A*b - B*a)*log(a + b*x**2 + c*x**4)/(4*a**4) + (-3*A*(-20*a**3*c**3 + 30*a**2*b**2*c**2 - 10*a*b**4*c + b**6) + B*a*b*(30*a**2*c**2 - 10*a*b**2*c + b**4))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**4*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_132():
    f = x**8*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x**7*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - x**5*(-4*A*a*c + 7*A*b**2 - 12*B*a*b + x**2*(12*A*b*c - 28*B*a*c + B*b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + x**3*(12*A*b*c - 28*B*a*c + B*b**2)/(8*c*(-4*a*c + b**2)**2) - x*(20*A*a*c**2 + A*b**2*c - 24*B*a*b*c + 3*B*b**3)/(8*c**2*(-4*a*c + b**2)**2) + sqrt(2)*(-16*A*a*b*c**2 + A*b**3*c + 84*B*a**2*c**2 - 27*B*a*b**2*c + 3*B*b**4 + (-40*A*a**2*c**3 - 18*A*a*b**2*c**2 + A*b**4*c + 132*B*a**2*b*c**2 - 33*B*a*b**3*c + 3*B*b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(-16*A*a*b*c**2 + A*b**3*c + 84*B*a**2*c**2 - 27*B*a*b**2*c + 3*B*b**4 - (-40*A*a**2*c**3 - 18*A*a*b**2*c**2 + A*b**4*c + 132*B*a**2*b*c**2 - 33*B*a*b**3*c + 3*B*b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_133():
    f = x**6*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x**5*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - x**3*(4*A*a*c + 5*A*b**2 - 12*B*a*b - x**2*(-12*A*b*c + 20*B*a*c + B*b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - x*(-12*A*b*c + 20*B*a*c + B*b**2)/(8*c*(-4*a*c + b**2)**2) + sqrt(2)*(12*A*a*c**2 + 3*A*b**2*c - 16*B*a*b*c + B*b**3 + (36*A*a*b*c**2 + 3*A*b**3*c - 40*B*a**2*c**2 - 18*B*a*b**2*c + B*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(12*A*a*c**2 + 3*A*b**2*c - 16*B*a*b*c + B*b**3 - (36*A*a*b*c**2 + 3*A*b**3*c - 40*B*a**2*c**2 - 18*B*a*b**2*c + B*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_134():
    f = x**4*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x**3*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + 3*x*(-A*(4*a*c + b**2) + 4*B*a*b + x**2*(-4*A*b*c + 4*B*a*c + B*b**2))/(8*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + sqrt(2)*(-12*A*b*c + 12*B*a*c + 3*B*b**2 + 3*(-8*A*a*c**2 - 6*A*b**2*c + 12*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(-12*A*b*c + 12*B*a*c + 3*B*b**2 - 3*(-8*A*a*c**2 - 6*A*b**2*c + 12*B*a*b*c + B*b**3)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_135():
    f = x**2*(A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = -x*(A*b - 2*B*a - x**2*(-2*A*c + B*b))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) - sqrt(2)*sqrt(c)*(A*(-52*a*b*c - 20*a*c*sqrt(-4*a*c + b**2) + b**3 - b**2*sqrt(-4*a*c + b**2)) + 6*B*a*(4*a*c + 3*b**2 + 2*b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + sqrt(2)*sqrt(c)*(A*(-52*a*b*c + 20*a*c*sqrt(-4*a*c + b**2) + b**3 + b**2*sqrt(-4*a*c + b**2)) + 6*B*a*(4*a*c + 3*b**2 - 2*b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) - x*(-A*(8*a*b*c + b**3) + B*a*(-4*a*c + 7*b**2) + c*x**2*(-A*(20*a*c + b**2) + 12*B*a*b))/(8*a*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_136():
    f = (A + B*x**2)/(a + b*x**2 + c*x**4)**3
    F = x*(-2*A*a*c + A*b**2 - B*a*b + c*x**2*(A*b - 2*B*a))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(3*A*(-8*a*b*c + b**3) + B*a*(20*a*c + b**2) - (3*A*(56*a**2*c**2 - 10*a*b**2*c + b**4) + B*a*b*(-52*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(3*A*(-8*a*b*c + b**3) + B*a*(20*a*c + b**2) + (3*A*(56*a**2*c**2 - 10*a*b**2*c + b**4) + B*a*b*(-52*a*c + b**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*(A*(28*a**2*c**2 - 25*a*b**2*c + 3*b**4) + B*a*b*(8*a*c + b**2) + c*x**2*(3*A*(-8*a*b*c + b**3) + B*a*(20*a*c + b**2)))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_137():
    f = x*(4*x**2 - 7)/(x**4 - 5*x**2 + 4)
    F = log(1 - x**2)/2 + 3*log(4 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_138():
    f = (4*x**3 - 7*x)/(x**4 - 5*x**2 + 4)
    F = log(1 - x**2)/2 + 3*log(4 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_139():
    f = x*(x**2 + 2)/(x**4 + x**2 + 1)
    F = log(x**4 + x**2 + 1)/4 + sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_140():
    f = (x**3 + 2*x)/(x**4 + x**2 + 1)
    F = log(x**4 + x**2 + 1)/4 + sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_141():
    f = (2*x**3 + 11*x)/(x**4 + 2*x**2 + 3)**2
    F = (9*x**2 + 5)/(8*x**4 + 16*x**2 + 24) + 9*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_142():
    f = x**5*(3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = 3*x**4*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/10 + (sympy.S(1837)/480 - 17*x**2/16)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2) + (-1633*x**2/128 + sympy.S(-8165)/256)*sqrt(x**4 + 5*x**2 + 3) + 21229*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_143():
    f = x**3*(3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = -(59 - 18*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/48 + (518*x**2 + 1295)*sqrt(x**4 + 5*x**2 + 3)/128 - 3367*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_144():
    f = x*(3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = (-22*x**2 - 55)*sqrt(x**4 + 5*x**2 + 3)/16 + (x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/2 + 143*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_145():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x
    F = (6*x**2 + 23)*sqrt(x**4 + 5*x**2 + 3)/8 + atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/16 - sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_146():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**3
    F = 19*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/4 - 7*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/3 - (2 - 3*x**2)*sqrt(x**4 + 5*x**2 + 3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_147():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**5
    F = 3*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/2 - 77*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/72 - (23*x**2 + 6)*sqrt(x**4 + 5*x**2 + 3)/(12*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_148():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**7
    F = 13*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/108 - (5*x**2 + 6)*sqrt(x**4 + 5*x**2 + 3)/(18*x**4) - (x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(9*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_149():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**9
    F = -871*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/10368 + (335*x**2 + 402)*sqrt(x**4 + 5*x**2 + 3)/(1728*x**4) - 11*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(216*x**6) - (x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(12*x**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_150():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**11
    F = 2093*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/31104 - (805*x**2 + 966)*sqrt(x**4 + 5*x**2 + 3)/(5184*x**4) + 173*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(3240*x**6) - (x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(36*x**8) - (x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(15*x**10)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_151():
    f = x**4*(3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = x**5*(7*x**2 + 11)*sqrt(x**4 + 5*x**2 + 3)/21 - 26*x**3*sqrt(x**4 + 5*x**2 + 3)/35 - 1924*x*(2*x**2 + sqrt(13) + 5)/(105*sqrt(x**4 + 5*x**2 + 3)) + 13*x*sqrt(x**4 + 5*x**2 + 3)/3 + 962*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(2*sqrt(13)/3 + sympy.S(10)/3)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(105*sqrt(x**4 + 5*x**2 + 3)) - 13*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_152():
    f = x**2*(3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = x**3*(15*x**2 + 29)*sqrt(x**4 + 5*x**2 + 3)/35 + 1247*x*(2*x**2 + sqrt(13) + 5)/(210*sqrt(x**4 + 5*x**2 + 3)) - 4*x*sqrt(x**4 + 5*x**2 + 3)/3 - 1247*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(210*sqrt(x**4 + 5*x**2 + 3)) + 2*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_153():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)
    F = x*(9*x**2 + 25)*sqrt(x**4 + 5*x**2 + 3)/15 - 23*x*(2*x**2 + sqrt(13) + 5)/(15*sqrt(x**4 + 5*x**2 + 3)) + 23*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(15*sqrt(x**4 + 5*x**2 + 3)) + sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_154():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**2
    F = 9*x*(2*x**2 + sqrt(13) + 5)/(2*sqrt(x**4 + 5*x**2 + 3)) - 3*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(3*sqrt(13)/2 + sympy.S(15)/2)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(2*sqrt(x**4 + 5*x**2 + 3)) + 8*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3)) - (2 - x**2)*sqrt(x**4 + 5*x**2 + 3)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_155():
    f = (3*x**2 + 2)*sqrt(x**4 + 5*x**2 + 3)/x**4
    F = 32*x*(2*x**2 + sqrt(13) + 5)/(9*sqrt(x**4 + 5*x**2 + 3)) - 16*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(2*sqrt(13)/3 + sympy.S(10)/3)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(9*sqrt(x**4 + 5*x**2 + 3)) + 49*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(3*sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3)) - 64*sqrt(x**4 + 5*x**2 + 3)/(9*x) - (2 - 9*x**2)*sqrt(x**4 + 5*x**2 + 3)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_156():
    f = x**5*(3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = 3*x**4*(x**4 + 5*x**2 + 3)**(sympy.S(5)/2)/14 + (3313 - 1070*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(5)/2)/1680 - (2183*x**2/384 + sympy.S(10915)/768)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2) + (56758*x**2 + 141895)*sqrt(x**4 + 5*x**2 + 3)/2048 - 368927*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/4096
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_157():
    f = x**3*(3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = -(sympy.S(27)/40 - x**2/4)*(x**4 + 5*x**2 + 3)**(sympy.S(5)/2) + (123*x**2/64 + sympy.S(615)/128)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2) - (9594*x**2 + 23985)*sqrt(x**4 + 5*x**2 + 3)/1024 + 62361*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/2048
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_158():
    f = x*(3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = -(11*x**2/16 + sympy.S(55)/32)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2) + (429*x**2/128 + sympy.S(2145)/256)*sqrt(x**4 + 5*x**2 + 3) + 3*(x**4 + 5*x**2 + 3)**(sympy.S(5)/2)/10 - 5577*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/512
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_159():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x
    F = (sympy.S(199)/128 - 37*x**2/64)*sqrt(x**4 + 5*x**2 + 3) + (3*x**2/8 + sympy.S(61)/48)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2) + 2401*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/256 - 3*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_160():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**3
    F = (27*x**2/8 + sympy.S(327)/16)*sqrt(x**4 + 5*x**2 + 3) + 609*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/32 - 12*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3))) - (2 - x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_161():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**5
    F = 453*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/16 - 127*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/8 - (84 - 57*x**2)*sqrt(x**4 + 5*x**2 + 3)/(8*x**2) - (2 - 3*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_162():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**7
    F = 49*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/4 - 527*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/72 - (67 - 32*x**2)*sqrt(x**4 + 5*x**2 + 3)/(12*x**2) - (7*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_163():
    f = x**4*(3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = x**5*(33*x**2 + 71)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/143 - x**5*(272*x**2 + 283)*sqrt(x**4 + 5*x**2 + 3)/429 + 1251*x**3*sqrt(x**4 + 5*x**2 + 3)/715 + 176723*x*(2*x**2 + sqrt(13) + 5)/(4290*sqrt(x**4 + 5*x**2 + 3)) - 4210*x*sqrt(x**4 + 5*x**2 + 3)/429 - 176723*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(4290*sqrt(x**4 + 5*x**2 + 3)) + 2105*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(143*sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_164():
    f = x**2*(3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = x**3*(27*x**2 + 67)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/99 - x**3*(890*x**2 + 911)*sqrt(x**4 + 5*x**2 + 3)/1155 - 49949*x*(2*x**2 + sqrt(13) + 5)/(3465*sqrt(x**4 + 5*x**2 + 3)) + 353*x*sqrt(x**4 + 5*x**2 + 3)/99 + 49949*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(3465*sqrt(x**4 + 5*x**2 + 3)) - 353*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(33*sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_165():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = x*(x**2 + 3)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/3 - x*(12*x**2 + 5)*sqrt(x**4 + 5*x**2 + 3)/15 + 203*x*(2*x**2 + sqrt(13) + 5)/(30*sqrt(x**4 + 5*x**2 + 3)) - 203*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(30*sqrt(x**4 + 5*x**2 + 3)) + 5*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_166():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**2
    F = x*(129*x**2 + 655)*sqrt(x**4 + 5*x**2 + 3)/35 + 412*x*(2*x**2 + sqrt(13) + 5)/(35*sqrt(x**4 + 5*x**2 + 3)) - 206*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(2*sqrt(13)/3 + sympy.S(10)/3)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(35*sqrt(x**4 + 5*x**2 + 3)) + 19*sqrt(3)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(2*sqrt(13) + 10)*sqrt(x**4 + 5*x**2 + 3)) - (14 - 3*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_167():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**4
    F = 949*x*(2*x**2 + sqrt(13) + 5)/(30*sqrt(x**4 + 5*x**2 + 3)) - 949*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(30*sqrt(x**4 + 5*x**2 + 3)) + 65*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3)) - (312 - 65*x**2)*sqrt(x**4 + 5*x**2 + 3)/(15*x) - (10 - 9*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_168():
    f = (3*x**2 + 2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/x**6
    F = 361*x*(2*x**2 + sqrt(13) + 5)/(15*sqrt(x**4 + 5*x**2 + 3)) - 361*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(15*sqrt(x**4 + 5*x**2 + 3)) + 103*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3)) - 722*sqrt(x**4 + 5*x**2 + 3)/(15*x) - (40 - 87*x**2)*sqrt(x**4 + 5*x**2 + 3)/(5*x**3) - (2 - 5*x**2)*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_169():
    f = x**5*(A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = B*x**4*sqrt(a + b*x**2 + c*x**4)/(6*c) + sqrt(a + b*x**2 + c*x**4)*(-18*A*b*c - 16*B*a*c + 15*B*b**2 - 2*c*x**2*(-6*A*c + 5*B*b))/(48*c**3) - (8*A*a*c**2 - 6*A*b**2*c - 12*B*a*b*c + 5*B*b**3)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_170():
    f = x**3*(A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = -sqrt(a + b*x**2 + c*x**4)*(-4*A*c + 3*B*b - 2*B*c*x**2)/(8*c**2) + (-4*A*b*c - 4*B*a*c + 3*B*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_171():
    f = x*(A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = B*sqrt(a + b*x**2 + c*x**4)/(2*c) - (-2*A*c + B*b)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_172():
    f = (A + B*x**2)/(x*sqrt(a + b*x**2 + c*x**4))
    F = -A*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a)) + B*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_173():
    f = (A + B*x**2)/(x**3*sqrt(a + b*x**2 + c*x**4))
    F = -A*sqrt(a + b*x**2 + c*x**4)/(2*a*x**2) + (A*b - 2*B*a)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_174():
    f = (A + B*x**2)/(x**5*sqrt(a + b*x**2 + c*x**4))
    F = -A*sqrt(a + b*x**2 + c*x**4)/(4*a*x**4) + (3*A*b - 4*B*a)*sqrt(a + b*x**2 + c*x**4)/(8*a**2*x**2) - (-4*A*a*c + 3*A*b**2 - 4*B*a*b)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_175():
    f = (A + B*x**2)/(x**7*sqrt(a + b*x**2 + c*x**4))
    F = -A*sqrt(a + b*x**2 + c*x**4)/(6*a*x**6) + (5*A*b - 6*B*a)*sqrt(a + b*x**2 + c*x**4)/(24*a**2*x**4) - sqrt(a + b*x**2 + c*x**4)*(-16*A*a*c + 15*A*b**2 - 18*B*a*b)/(48*a**3*x**2) + (-12*A*a*b*c + 5*A*b**3 + 8*B*a**2*c - 6*B*a*b**2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(32*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_176():
    f = x**4*(A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = B*x**3*sqrt(a + b*x**2 + c*x**4)/(5*c) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-10*A*b*c - 9*B*a*c + 8*B*b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-10*A*b*c - 9*B*a*c + 8*B*b**2 + sqrt(a)*sqrt(c)*(-5*A*c + 4*B*b))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) - x*(-5*A*c + 4*B*b)*sqrt(a + b*x**2 + c*x**4)/(15*c**2) + x*sqrt(a + b*x**2 + c*x**4)*(-10*A*b*c - 9*B*a*c + 8*B*b**2)/(15*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_177():
    f = x**2*(A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = B*x*sqrt(a + b*x**2 + c*x**4)/(3*c) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*A*c + 2*B*b)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*A*c + B*sqrt(a)*sqrt(c) + 2*B*b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - x*(-3*A*c + 2*B*b)*sqrt(a + b*x**2 + c*x**4)/(3*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_178():
    f = (A + B*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = -B*a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + B*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(A*sqrt(c)/sqrt(a) + B)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_179():
    f = (A + B*x**2)/(x**2*sqrt(a + b*x**2 + c*x**4))
    F = A*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)) - A*sqrt(a + b*x**2 + c*x**4)/(a*x) - A*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(A*sqrt(c) + B*sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_180():
    f = (A + B*x**2)/(x**4*sqrt(a + b*x**2 + c*x**4))
    F = -A*sqrt(a + b*x**2 + c*x**4)/(3*a*x**3) - sqrt(c)*x*(2*A*b - 3*B*a)*sqrt(a + b*x**2 + c*x**4)/(3*a**2*(sqrt(a) + sqrt(c)*x**2)) + (2*A*b - 3*B*a)*sqrt(a + b*x**2 + c*x**4)/(3*a**2*x) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*A*b - 3*B*a)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(A*sqrt(a)*sqrt(c) + 2*A*b - 3*B*a)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*a**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_181():
    f = x**7*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = 3*x**6*sqrt(x**4 + 5*x**2 + 3)/8 - 89*x**4*sqrt(x**4 + 5*x**2 + 3)/48 - (sympy.S(8081)/128 - 1901*x**2/192)*sqrt(x**4 + 5*x**2 + 3) + 32801*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_182():
    f = x**5*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = x**4*sqrt(x**4 + 5*x**2 + 3)/2 + (sympy.S(267)/16 - 21*x**2/8)*sqrt(x**4 + 5*x**2 + 3) - 1083*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_183():
    f = x**3*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = -(37 - 6*x**2)*sqrt(x**4 + 5*x**2 + 3)/8 + 149*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_184():
    f = x*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = 3*sqrt(x**4 + 5*x**2 + 3)/2 - 11*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_185():
    f = (3*x**2 + 2)/(x*sqrt(x**4 + 5*x**2 + 3))
    F = 3*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/2 - sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_186():
    f = (3*x**2 + 2)/(x**3*sqrt(x**4 + 5*x**2 + 3))
    F = -2*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/9 - sqrt(x**4 + 5*x**2 + 3)/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_187():
    f = (3*x**2 + 2)/(x**5*sqrt(x**4 + 5*x**2 + 3))
    F = sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/8 - sqrt(x**4 + 5*x**2 + 3)/(12*x**2) - sqrt(x**4 + 5*x**2 + 3)/(6*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_188():
    f = (3*x**2 + 2)/(x**7*sqrt(x**4 + 5*x**2 + 3))
    F = -61*sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/648 + 13*sqrt(x**4 + 5*x**2 + 3)/(108*x**2) - sqrt(x**4 + 5*x**2 + 3)/(54*x**4) - sqrt(x**4 + 5*x**2 + 3)/(9*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_189():
    f = x**4*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = 3*x**3*sqrt(x**4 + 5*x**2 + 3)/5 + 419*x*(2*x**2 + sqrt(13) + 5)/(30*sqrt(x**4 + 5*x**2 + 3)) - 10*x*sqrt(x**4 + 5*x**2 + 3)/3 - 419*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(30*sqrt(x**4 + 5*x**2 + 3)) + 5*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_190():
    f = x**2*(3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = -4*x*(2*x**2 + sqrt(13) + 5)/sqrt(x**4 + 5*x**2 + 3) + x*sqrt(x**4 + 5*x**2 + 3) + 2*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(2*sqrt(13)/3 + sympy.S(10)/3)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/sqrt(x**4 + 5*x**2 + 3) - sqrt(3)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(2*sqrt(13) + 10)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_191():
    f = (3*x**2 + 2)/sqrt(x**4 + 5*x**2 + 3)
    F = 3*x*(2*x**2 + sqrt(13) + 5)/(2*sqrt(x**4 + 5*x**2 + 3)) - sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(3*sqrt(13)/2 + sympy.S(15)/2)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(2*sqrt(x**4 + 5*x**2 + 3)) + sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_192():
    f = (3*x**2 + 2)/(x**2*sqrt(x**4 + 5*x**2 + 3))
    F = x*(2*x**2 + sqrt(13) + 5)/(3*sqrt(x**4 + 5*x**2 + 3)) - sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(3*sqrt(x**4 + 5*x**2 + 3)) + sqrt(3)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(sqrt(2*sqrt(13) + 10)*sqrt(x**4 + 5*x**2 + 3)) - 2*sqrt(x**4 + 5*x**2 + 3)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_193():
    f = (3*x**2 + 2)/(x**4*sqrt(x**4 + 5*x**2 + 3))
    F = 7*x*(2*x**2 + sqrt(13) + 5)/(54*sqrt(x**4 + 5*x**2 + 3)) - 7*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(54*sqrt(x**4 + 5*x**2 + 3)) - sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(9*sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3)) - 7*sqrt(x**4 + 5*x**2 + 3)/(27*x) - 2*sqrt(x**4 + 5*x**2 + 3)/(9*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_194():
    f = x**5*(3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = -x**2*(47*x**2 + 33)/(13*sqrt(x**4 + 5*x**2 + 3)) + 133*sqrt(x**4 + 5*x**2 + 3)/26 - 41*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_195():
    f = x**3*(3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = (-47*x**2 - 33)/(13*sqrt(x**4 + 5*x**2 + 3)) + 3*atanh((2*x**2 + 5)/(2*sqrt(x**4 + 5*x**2 + 3)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_196():
    f = x*(3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = (11*x**2 + 8)/(13*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_197():
    f = (3*x**2 + 2)/(x*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2))
    F = (-8*x**2 - 7)/(39*sqrt(x**4 + 5*x**2 + 3)) - sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_198():
    f = (3*x**2 + 2)/(x**3*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2))
    F = sqrt(3)*atanh(sqrt(3)*(5*x**2 + 6)/(6*sqrt(x**4 + 5*x**2 + 3)))/9 + (-8*x**2 - 7)/(39*x**2*sqrt(x**4 + 5*x**2 + 3)) - 2*sqrt(x**4 + 5*x**2 + 3)/(39*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_199():
    f = x**4*(3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = x**3*(11*x**2 + 8)/(13*sqrt(x**4 + 5*x**2 + 3)) + 43*x*(2*x**2 + sqrt(13) + 5)/(13*sqrt(x**4 + 5*x**2 + 3)) - 11*x*sqrt(x**4 + 5*x**2 + 3)/13 - 43*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(13*sqrt(x**4 + 5*x**2 + 3)) + 11*sqrt(3)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(13*sqrt(2*sqrt(13) + 10)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_200():
    f = x**2*(3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = x*(11*x**2 + 8)/(13*sqrt(x**4 + 5*x**2 + 3)) - 11*x*(2*x**2 + sqrt(13) + 5)/(26*sqrt(x**4 + 5*x**2 + 3)) + 11*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(26*sqrt(x**4 + 5*x**2 + 3)) - 4*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(13*sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_201():
    f = (3*x**2 + 2)/(x**4 + 5*x**2 + 3)**(sympy.S(3)/2)
    F = -x*(8*x**2 + 7)/(39*sqrt(x**4 + 5*x**2 + 3)) + 4*x*(2*x**2 + sqrt(13) + 5)/(39*sqrt(x**4 + 5*x**2 + 3)) - 2*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(2*sqrt(13)/3 + sympy.S(10)/3)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(39*sqrt(x**4 + 5*x**2 + 3)) + 11*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(13*sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_202():
    f = (3*x**2 + 2)/(x**2*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2))
    F = 19*x*(2*x**2 + sqrt(13) + 5)/(234*sqrt(x**4 + 5*x**2 + 3)) - 19*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(234*sqrt(x**4 + 5*x**2 + 3)) - 4*sqrt(2)*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(39*sqrt(3*sqrt(13) + 15)*sqrt(x**4 + 5*x**2 + 3)) - (8*x**2 + 7)/(39*x*sqrt(x**4 + 5*x**2 + 3)) - 19*sqrt(x**4 + 5*x**2 + 3)/(117*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_203():
    f = (3*x**2 + 2)/(x**4*(x**4 + 5*x**2 + 3)**(sympy.S(3)/2))
    F = -133*x*(2*x**2 + sqrt(13) + 5)/(1053*sqrt(x**4 + 5*x**2 + 3)) + 133*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*sqrt(sqrt(13)/6 + sympy.S(5)/6)*(x**2*(sqrt(13) + 5) + 6)*elliptic_e(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(1053*sqrt(x**4 + 5*x**2 + 3)) - 5*sqrt((x**2*(5 - sqrt(13)) + 6)/(x**2*(sqrt(13) + 5) + 6))*(x**2*(sqrt(13) + 5) + 6)*elliptic_f(atan(x*sqrt(sqrt(13)/6 + sympy.S(5)/6)), sympy.S(-13)/6 + 5*sqrt(13)/6)/(351*sqrt(6*sqrt(13) + 30)*sqrt(x**4 + 5*x**2 + 3)) + 266*sqrt(x**4 + 5*x**2 + 3)/(1053*x) - (8*x**2 + 7)/(39*x**3*sqrt(x**4 + 5*x**2 + 3)) - 5*sqrt(x**4 + 5*x**2 + 3)/(351*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_204():
    f = (f*x)**(sympy.S(3)/2)*(d + e*x**2)*sqrt(a + b*x**2 + c*x**4)
    F = 2*d*(f*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*e*(f*x)**(sympy.S(9)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(9)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(13)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(9*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_205():
    f = sqrt(f*x)*(d + e*x**2)*sqrt(a + b*x**2 + c*x**4)
    F = 2*d*(f*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*e*(f*x)**(sympy.S(7)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(7)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(11)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(7*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_206():
    f = (d + e*x**2)*sqrt(a + b*x**2 + c*x**4)/sqrt(f*x)
    F = 2*d*sqrt(f*x)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(1)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*e*(f*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_207():
    f = (d + e*x**2)*sqrt(a + b*x**2 + c*x**4)/(f*x)**(sympy.S(3)/2)
    F = -2*d*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(-1)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(f*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*e*(f*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_208():
    f = (f*x)**(sympy.S(3)/2)*(d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*a*d*(f*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*a*e*(f*x)**(sympy.S(9)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(9)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(13)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(9*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_209():
    f = sqrt(f*x)*(d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*a*d*(f*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*a*e*(f*x)**(sympy.S(7)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(7)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(11)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(7*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_210():
    f = (d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/sqrt(f*x)
    F = 2*a*d*sqrt(f*x)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(1)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*a*e*(f*x)**(sympy.S(5)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(5)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_211():
    f = (d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(f*x)**(sympy.S(3)/2)
    F = -2*a*d*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(-1)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(f*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + 2*a*e*(f*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4)*appellf1(sympy.S(3)/4, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f**3*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_212():
    f = (f*x)**(sympy.S(3)/2)*(d + e*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = 2*d*(f*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S.Half, sympy.S.Half, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(9)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(9)/4, sympy.S.Half, sympy.S.Half, sympy.S(13)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(9*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_213():
    f = sqrt(f*x)*(d + e*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = 2*d*(f*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S.Half, sympy.S.Half, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(7)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(7)/4, sympy.S.Half, sympy.S.Half, sympy.S(11)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(7*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_214():
    f = (d + e*x**2)/(sqrt(f*x)*sqrt(a + b*x**2 + c*x**4))
    F = 2*d*sqrt(f*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/4, sympy.S.Half, sympy.S.Half, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S.Half, sympy.S.Half, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_215():
    f = (d + e*x**2)/((f*x)**(sympy.S(3)/2)*sqrt(a + b*x**2 + c*x**4))
    F = -2*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/4, sympy.S.Half, sympy.S.Half, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*sqrt(f*x)*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S.Half, sympy.S.Half, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_216():
    f = (f*x)**(sympy.S(3)/2)*(d + e*x**2)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*d*(f*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*a*f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(9)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(9)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(13)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(9*a*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_217():
    f = sqrt(f*x)*(d + e*x**2)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = 2*d*(f*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*a*f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(7)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(7)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(11)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(7*a*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_218():
    f = (d + e*x**2)/(sqrt(f*x)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = 2*d*sqrt(f*x)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(5)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*f*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(5)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(5)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(9)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(5*a*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_219():
    f = (d + e*x**2)/((f*x)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -2*d*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(3)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*f*sqrt(f*x)*sqrt(a + b*x**2 + c*x**4)) + 2*e*(f*x)**(sympy.S(3)/2)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S(3)/2, sympy.S(3)/2, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*a*f**3*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_220():
    f = (f*x)**m*(d + e*x**2)*(a + b*x**2 + c*x**4)**3
    F = a**3*d*(f*x)**(m + 1)/(f*(m + 1)) + a**2*(f*x)**(m + 3)*(a*e + 3*b*d)/(f**3*(m + 3)) + 3*a*(f*x)**(m + 5)*(a*b*e + a*c*d + b**2*d)/(f**5*(m + 5)) + c**3*e*(f*x)**(m + 15)/(f**15*(m + 15)) + c**2*(f*x)**(m + 13)*(3*b*e + c*d)/(f**13*(m + 13)) + 3*c*(f*x)**(m + 11)*(a*c*e + b**2*e + b*c*d)/(f**11*(m + 11)) + (f*x)**(m + 7)*(3*a**2*c*e + 3*a*b**2*e + 6*a*b*c*d + b**3*d)/(f**7*(m + 7)) + (f*x)**(m + 9)*(6*a*b*c*e + 3*a*c**2*d + b**3*e + 3*b**2*c*d)/(f**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_221():
    f = (f*x)**m*(d + e*x**2)*(a + b*x**2 + c*x**4)**2
    F = a**2*d*(f*x)**(m + 1)/(f*(m + 1)) + a*(f*x)**(m + 3)*(a*e + 2*b*d)/(f**3*(m + 3)) + c**2*e*(f*x)**(m + 11)/(f**11*(m + 11)) + c*(f*x)**(m + 9)*(2*b*e + c*d)/(f**9*(m + 9)) + (f*x)**(m + 5)*(2*a*b*e + 2*a*c*d + b**2*d)/(f**5*(m + 5)) + (f*x)**(m + 7)*(2*a*c*e + b**2*e + 2*b*c*d)/(f**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_222():
    f = (f*x)**m*(d + e*x**2)*(a + b*x**2 + c*x**4)
    F = a*d*(f*x)**(m + 1)/(f*(m + 1)) + c*e*(f*x)**(m + 7)/(f**7*(m + 7)) + (f*x)**(m + 3)*(a*e + b*d)/(f**3*(m + 3)) + (f*x)**(m + 5)*(b*e + c*d)/(f**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_223():
    f = (f*x)**m*(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = (f*x)**(m + 1)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*(b + sqrt(-4*a*c + b**2))*(m + 1)) + (f*x)**(m + 1)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(f*(b - sqrt(-4*a*c + b**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_224():
    f = (f*x)**m*(d + e*x**2)/(a + b*x**2 + c*x**4)**2
    F = -c*(f*x)**(m + 1)*(2*a*(-2*c*d*(3 - m) + e*(1 - m)*sqrt(-4*a*c + b**2)) + b**2*d*(1 - m) + b*(4*a*e - d*(1 - m)*sqrt(-4*a*c + b**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(2*a*f*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + c*(f*x)**(m + 1)*(-2*a*(2*c*d*(3 - m) + e*(1 - m)*sqrt(-4*a*c + b**2)) + b**2*(-d*m + d) + b*(4*a*e + d*(1 - m)*sqrt(-4*a*c + b**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -2*c*x**2/(b - sqrt(-4*a*c + b**2)))/(2*a*f*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + (f*x)**(m + 1)*(-a*b*e - 2*a*c*d + b**2*d + c*x**2*(-2*a*e + b*d))/(2*a*f*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_225():
    f = (f*x)**m*(d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = a*d*(f*x)**(m + 1)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S.Half, sympy.S(-3)/2, sympy.S(-3)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + a*e*(f*x)**(m + 3)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S(3)/2, sympy.S(-3)/2, sympy.S(-3)/2, m/2 + sympy.S(5)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f**3*(m + 3)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_226():
    f = (f*x)**m*(d + e*x**2)*sqrt(a + b*x**2 + c*x**4)
    F = d*(f*x)**(m + 1)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S.Half, sympy.S(-1)/2, sympy.S(-1)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + e*(f*x)**(m + 3)*sqrt(a + b*x**2 + c*x**4)*appellf1(m/2 + sympy.S(3)/2, sympy.S(-1)/2, sympy.S(-1)/2, m/2 + sympy.S(5)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f**3*(m + 3)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_227():
    f = (f*x)**m*(d + e*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = d*(f*x)**(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S.Half, sympy.S.Half, sympy.S.Half, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*sqrt(a + b*x**2 + c*x**4)) + e*(f*x)**(m + 3)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S(3)/2, sympy.S.Half, sympy.S.Half, m/2 + sympy.S(5)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(f**3*(m + 3)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_228():
    f = (f*x)**m*(d + e*x**2)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = d*(f*x)**(m + 1)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S.Half, sympy.S(3)/2, sympy.S(3)/2, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*f*(m + 1)*sqrt(a + b*x**2 + c*x**4)) + e*(f*x)**(m + 3)*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/2 + sympy.S(3)/2, sympy.S(3)/2, sympy.S(3)/2, m/2 + sympy.S(5)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(a*f**3*(m + 3)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_229():
    f = x**9/((a + c*x**4)*(d + e*x**2))
    F = a**(sympy.S(3)/2)*d*atan(sqrt(c)*x**2/sqrt(a))/(2*c**(sympy.S(3)/2)*(a*e**2 + c*d**2)) - a**2*e*log(a + c*x**4)/(4*c**2*(a*e**2 + c*d**2)) + d**4*log(d + e*x**2)/(2*e**3*(a*e**2 + c*d**2)) - d*x**2/(2*c*e**2) + x**4/(4*c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_230():
    f = x**7/((a + c*x**4)*(d + e*x**2))
    F = -a**(sympy.S(3)/2)*e*atan(sqrt(c)*x**2/sqrt(a))/(2*c**(sympy.S(3)/2)*(a*e**2 + c*d**2)) - a*d*log(a + c*x**4)/(4*c*(a*e**2 + c*d**2)) - d**3*log(d + e*x**2)/(2*e**2*(a*e**2 + c*d**2)) + x**2/(2*c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_231():
    f = x**5/((a + c*x**4)*(d + e*x**2))
    F = -sqrt(a)*d*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(c)*(a*e**2 + c*d**2)) + a*e*log(a + c*x**4)/(4*c*(a*e**2 + c*d**2)) + d**2*log(d + e*x**2)/(2*e*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_232():
    f = x**3/((a + c*x**4)*(d + e*x**2))
    F = sqrt(a)*e*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(c)*(a*e**2 + c*d**2)) + d*log(a + c*x**4)/(4*a*e**2 + 4*c*d**2) - d*log(d + e*x**2)/(2*a*e**2 + 2*c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_233():
    f = x/((a + c*x**4)*(d + e*x**2))
    F = -e*log(a + c*x**4)/(4*a*e**2 + 4*c*d**2) + e*log(d + e*x**2)/(2*a*e**2 + 2*c*d**2) + sqrt(c)*d*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_234():
    f = 1/(x*(a + c*x**4)*(d + e*x**2))
    F = -e**2*log(d + e*x**2)/(2*d*(a*e**2 + c*d**2)) - c*d*log(a + c*x**4)/(4*a*(a*e**2 + c*d**2)) + log(x)/(a*d) - sqrt(c)*e*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_235():
    f = 1/(x**3*(a + c*x**4)*(d + e*x**2))
    F = e**3*log(d + e*x**2)/(2*d**2*(a*e**2 + c*d**2)) + c*e*log(a + c*x**4)/(4*a*(a*e**2 + c*d**2)) - 1/(2*a*d*x**2) - e*log(x)/(a*d**2) - c**(sympy.S(3)/2)*d*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_236():
    f = 1/(x**5*(a + c*x**4)*(d + e*x**2))
    F = -e**4*log(d + e*x**2)/(2*d**3*(a*e**2 + c*d**2)) - 1/(4*a*d*x**4) + e/(2*a*d**2*x**2) + c**2*d*log(a + c*x**4)/(4*a**2*(a*e**2 + c*d**2)) - (-a*e**2 + c*d**2)*log(x)/(a**2*d**3) + c**(sympy.S(3)/2)*e*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_237():
    f = x**8/((a + c*x**4)*(d + e*x**2))
    F = -sqrt(2)*a**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + d**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(e**(sympy.S(5)/2)*(a*e**2 + c*d**2)) - d*x/(c*e**2) + x**3/(3*c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_238():
    f = x**6/((a + c*x**4)*(d + e*x**2))
    F = -sqrt(2)*a**(sympy.S(3)/4)*(-sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(3)/4)*(-sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(3)/4)*(sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(3)/4)*(sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) - d**(sympy.S(5)/2)*atan(sqrt(e)*x/sqrt(d))/(e**(sympy.S(3)/2)*(a*e**2 + c*d**2)) + x/(c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_239():
    f = x**4/((a + c*x**4)*(d + e*x**2))
    F = sqrt(2)*a**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + d**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(e)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_240():
    f = x**2/((a + c*x**4)*(d + e*x**2))
    F = -sqrt(d)*sqrt(e)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2) + sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) + sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_241():
    f = 1/((a + c*x**4)*(d + e*x**2))
    F = e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_242():
    f = 1/(x**2*(a + c*x**4)*(d + e*x**2))
    F = -e**(sympy.S(5)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(3)/2)*(a*e**2 + c*d**2)) - 1/(a*d*x) - sqrt(2)*c**(sympy.S(3)/4)*(-sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(5)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(3)/4)*(-sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(5)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(3)/4)*(sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(5)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(3)/4)*(sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(5)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_243():
    f = 1/(x**4*(a + c*x**4)*(d + e*x**2))
    F = e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(5)/2)*(a*e**2 + c*d**2)) - 1/(3*a*d*x**3) + e/(a*d**2*x) + sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(7)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_244():
    f = x**9/((a + c*x**4)**2*(d + e*x**2))
    F = -sqrt(a)*d*(a*e**2 + 3*c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(4*c**(sympy.S(3)/2)*(a*e**2 + c*d**2)**2) + a*e*(a*e**2 + 2*c*d**2)*log(a + c*x**4)/(4*c**2*(a*e**2 + c*d**2)**2) + a*(a*e + c*d*x**2)/(4*c**2*(a + c*x**4)*(a*e**2 + c*d**2)) + d**4*log(d + e*x**2)/(2*e*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_245():
    f = x**7/((a + c*x**4)**2*(d + e*x**2))
    F = sqrt(a)*e*(a*e**2 + 3*c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(4*c**(sympy.S(3)/2)*(a*e**2 + c*d**2)**2) + a*(d - e*x**2)/(4*c*(a + c*x**4)*(a*e**2 + c*d**2)) + d**3*log(a + c*x**4)/(4*(a*e**2 + c*d**2)**2) - d**3*log(d + e*x**2)/(2*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_246():
    f = x**5/((a + c*x**4)**2*(d + e*x**2))
    F = -d**2*e*log(a + c*x**4)/(4*(a*e**2 + c*d**2)**2) + d**2*e*log(d + e*x**2)/(2*(a*e**2 + c*d**2)**2) + (-a*e - c*d*x**2)/(4*c*(a + c*x**4)*(a*e**2 + c*d**2)) + d*(-a*e**2 + c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(4*sqrt(a)*sqrt(c)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_247():
    f = x**3/((a + c*x**4)**2*(d + e*x**2))
    F = d*e**2*log(a + c*x**4)/(4*(a*e**2 + c*d**2)**2) - d*e**2*log(d + e*x**2)/(2*(a*e**2 + c*d**2)**2) + (-d + e*x**2)/((a + c*x**4)*(4*a*e**2 + 4*c*d**2)) - e*(-a*e**2 + c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(4*sqrt(a)*sqrt(c)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_248():
    f = x/((a + c*x**4)**2*(d + e*x**2))
    F = -e**3*log(a + c*x**4)/(4*(a*e**2 + c*d**2)**2) + e**3*log(d + e*x**2)/(2*(a*e**2 + c*d**2)**2) + (a*e + c*d*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)) + sqrt(c)*d*(3*a*e**2 + c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_249():
    f = 1/(x*(a + c*x**4)**2*(d + e*x**2))
    F = -e**4*log(d + e*x**2)/(2*d*(a*e**2 + c*d**2)**2) + c*(d - e*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)) - c*d*(2*a*e**2 + c*d**2)*log(a + c*x**4)/(4*a**2*(a*e**2 + c*d**2)**2) + log(x)/(a**2*d) - sqrt(c)*e**3*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**2 + c*d**2)**2) - sqrt(c)*e*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_250():
    f = 1/(x**3*(a + c*x**4)**2*(d + e*x**2))
    F = e**5*log(d + e*x**2)/(2*d**2*(a*e**2 + c*d**2)**2) + c*e*(2*a*e**2 + c*d**2)*log(a + c*x**4)/(4*a**2*(a*e**2 + c*d**2)**2) - c*(a*e + c*d*x**2)/(4*a**2*(a + c*x**4)*(a*e**2 + c*d**2)) - 1/(2*a**2*d*x**2) - e*log(x)/(a**2*d**2) - c**(sympy.S(3)/2)*d*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(5)/2)*(a*e**2 + c*d**2)) - c**(sympy.S(3)/2)*d*(2*a*e**2 + c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(5)/2)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_251():
    f = 1/(x**5*(a + c*x**4)**2*(d + e*x**2))
    F = -e**6*log(d + e*x**2)/(2*d**3*(a*e**2 + c*d**2)**2) - c**2*(d - e*x**2)/(4*a**2*(a + c*x**4)*(a*e**2 + c*d**2)) - 1/(4*a**2*d*x**4) + e/(2*a**2*d**2*x**2) + c**2*d*(3*a*e**2 + 2*c*d**2)*log(a + c*x**4)/(4*a**3*(a*e**2 + c*d**2)**2) - (-a*e**2 + 2*c*d**2)*log(x)/(a**3*d**3) + c**(sympy.S(3)/2)*e*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(5)/2)*(a*e**2 + c*d**2)) + c**(sympy.S(3)/2)*e*(2*a*e**2 + c*d**2)*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(5)/2)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_252():
    f = x**8/((a + c*x**4)**2*(d + e*x**2))
    F = sqrt(2)*a**(sympy.S(1)/4)*d**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*a**(sympy.S(1)/4)*d**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*a**(sympy.S(1)/4)*d**2*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*a**(sympy.S(1)/4)*d**2*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*a**(sympy.S(1)/4)*(-3*sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(1)/4)*(-3*sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*a**(sympy.S(1)/4)*(3*sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*a**(sympy.S(1)/4)*(3*sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*c**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + d**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(e)*(a*e**2 + c*d**2)**2) + d*x/(4*c*(a*e**2 + c*d**2)) - x**3*(a*e + c*d*x**2)/(4*c*(a + c*x**4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_253():
    f = x**6/((a + c*x**4)**2*(d + e*x**2))
    F = -d**(sympy.S(5)/2)*sqrt(e)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2)**2 - x*(a*e + c*d*x**2)/(4*c*(a + c*x**4)*(a*e**2 + c*d**2)) + sqrt(2)*d**2*(-sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*d**2*(-sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*d**2*(sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*d**2*(sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)*(a*e**2 + c*d**2)) + sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_254():
    f = x**4/((a + c*x**4)**2*(d + e*x**2))
    F = d**(sympy.S(3)/2)*e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2)**2 - x*(d - e*x**2)/((a + c*x**4)*(4*a*e**2 + 4*c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*d**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*d**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*d**2*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*d**2*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*(sqrt(a)*e + 3*sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(sqrt(a)*e + 3*sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_255():
    f = x**2/((a + c*x**4)**2*(d + e*x**2))
    F = -sqrt(d)*e**(sympy.S(5)/2)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2)**2 + x*(a*e + c*d*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*d*e*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*e*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*d*e*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*e*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*(-3*sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(-3*sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) - sqrt(2)*(3*sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2)) + sqrt(2)*(3*sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_256():
    f = 1/((a + c*x**4)**2*(d + e*x**2))
    F = e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 + c*d**2)**2) + c*x*(d - e*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*e**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*e**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*e**2*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*e**2*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_257():
    f = 1/(x**2*(a + c*x**4)**2*(d + e*x**2))
    F = -e**(sympy.S(9)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(3)/2)*(a*e**2 + c*d**2)**2) - c*x*(a*e + c*d*x**2)/(4*a**2*(a + c*x**4)*(a*e**2 + c*d**2)) - 1/(a**2*d*x) - sqrt(2)*c**(sympy.S(3)/4)*(-3*sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(3)/4)*(-3*sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(3)/4)*(3*sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(3)/4)*(3*sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(3)/4)*(a**(sympy.S(3)/2)*e**3 - sqrt(c)*d*(2*a*e**2 + c*d**2))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(3)/4)*(a**(sympy.S(3)/2)*e**3 - sqrt(c)*d*(2*a*e**2 + c*d**2))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(3)/4)*(a**(sympy.S(3)/2)*e**3 + sqrt(c)*d*(2*a*e**2 + c*d**2))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(3)/4)*(a**(sympy.S(3)/2)*e**3 + sqrt(c)*d*(2*a*e**2 + c*d**2))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(9)/4)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_258():
    f = 1/(x**4*(a + c*x**4)**2*(d + e*x**2))
    F = e**(sympy.S(11)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(5)/2)*(a*e**2 + c*d**2)**2) - c**2*x*(d - e*x**2)/(4*a**2*(a + c*x**4)*(a*e**2 + c*d**2)) - 1/(3*a**2*d*x**3) + e/(a**2*d**2*x) + sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*(2*a*e**2 + c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + sqrt(c)*d)*(2*a*e**2 + c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(5)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*(2*a*e**2 + c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + sqrt(c)*d)*(2*a*e**2 + c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(11)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(5)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(11)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_259():
    f = x**2/((x**2 + 1)*sqrt(x**4 + 1))
    F = sqrt((x**4 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S.Half)/(4*sqrt(x**4 + 1)) - sqrt(2)*atan(sqrt(2)*x/sqrt(x**4 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_260():
    f = x**2/((1 - x**2)*sqrt(x**4 + 1))
    F = -sqrt((x**4 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S.Half)/(4*sqrt(x**4 + 1)) + sqrt(2)*atanh(sqrt(2)*x/sqrt(x**4 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_261():
    f = x**2/(sqrt(1 - x**4)*(x**2 + 1))
    F = -x*(1 - x**2)/(2*sqrt(1 - x**4)) - sqrt(1 - x**2)*sqrt(x**2 + 1)*elliptic_e(asin(x), -1)/(2*sqrt(1 - x**4)) + sqrt(1 - x**2)*sqrt(x**2 + 1)*elliptic_f(asin(x), -1)/sqrt(1 - x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_262():
    f = x**2/((1 - x**2)*sqrt(1 - x**4))
    F = x*(x**2 + 1)/(2*sqrt(1 - x**4)) - sqrt(1 - x**2)*sqrt(x**2 + 1)*elliptic_e(asin(x), -1)/(2*sqrt(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_263():
    f = x**2/((x**2 + 1)*sqrt(x**4 - 1))
    F = -x*(1 - x**2)/(2*sqrt(x**4 - 1)) - sqrt(1 - x**2)*sqrt(x**2 + 1)*elliptic_e(asin(x), -1)/(2*sqrt(x**4 - 1)) + sqrt(2)*sqrt(x**2 - 1)*sqrt(x**2 + 1)*elliptic_f(asin(sqrt(2)*x/sqrt(x**2 - 1)), sympy.S.Half)/(2*sqrt(x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_264():
    f = x**2/((1 - x**2)*sqrt(x**4 - 1))
    F = x*(x**2 + 1)/(2*sqrt(x**4 - 1)) - sqrt(1 - x**2)*sqrt(x**2 + 1)*elliptic_e(asin(x), -1)/(2*sqrt(x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_265():
    f = x**2/((x**2 + 1)*sqrt(-x**4 - 1))
    F = sqrt((x**4 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S.Half)/(4*sqrt(-x**4 - 1)) - sqrt(2)*atanh(sqrt(2)*x/sqrt(-x**4 - 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_266():
    f = x**2/((1 - x**2)*sqrt(-x**4 - 1))
    F = -sqrt((x**4 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S.Half)/(4*sqrt(-x**4 - 1)) + sqrt(2)*atan(sqrt(2)*x/sqrt(-x**4 - 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_267():
    f = x**2*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = b*x**3*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(6*d*(a + b*x**2)) + c**2*(-2*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(16*d**(sympy.S(5)/2)*(a + b*x**2)) - c*x*sqrt(c + d*x**2)*(-2*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(16*d**2*(a + b*x**2)) - x**3*sqrt(c + d*x**2)*(-2*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*d*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_268():
    f = x*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = b*(c + d*x**2)**(sympy.S(5)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(5*d**2*(a + b*x**2)) - (c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d**2*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_269():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)
    F = b*x*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(4*d*(a + b*x**2)) - c*(-4*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(3)/2)*(a + b*x**2)) - x*sqrt(c + d*x**2)*(-4*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(8*d*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_270():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x
    F = -a*sqrt(c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*atanh(sqrt(c + d*x**2)/sqrt(c))/(a + b*x**2) + a*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(a + b*x**2) + b*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(3*d*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_271():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**2
    F = -a*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(c*x*(a + b*x**2)) + (2*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*sqrt(d)*(a + b*x**2)) + x*sqrt(c + d*x**2)*(2*a*d + b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*c*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_272():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/x**3
    F = -a*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*c*x**2*(a + b*x**2)) + sqrt(c + d*x**2)*(a*d + 2*b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)/(2*c*(a + b*x**2)) - (a*d + 2*b*c)*sqrt(a**2 + 2*a*b*x**2 + b**2*x**4)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*sqrt(c)*(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_273():
    f = x**3*(d + e*x**2)**2*(a + b*x**2 + c*x**4)
    F = a*d**2*x**4/4 + c*e**2*x**12/12 + d*x**6*(2*a*e + b*d)/6 + e*x**10*(b*e + 2*c*d)/10 + x**8*(c*d**2/8 + e*(a*e + 2*b*d)/8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_274():
    f = x**2*(d + e*x**2)**2*(a + b*x**2 + c*x**4)
    F = a*d**2*x**3/3 + c*e**2*x**11/11 + d*x**5*(2*a*e + b*d)/5 + e*x**9*(b*e + 2*c*d)/9 + x**7*(c*d**2/7 + e*(a*e + 2*b*d)/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_275():
    f = x*(d + e*x**2)**2*(a + b*x**2 + c*x**4)
    F = c*(d + e*x**2)**5/(10*e**3) - (d + e*x**2)**4*(-b*e + 2*c*d)/(8*e**3) + (d + e*x**2)**3*(a*e**2 - b*d*e + c*d**2)/(6*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_276():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)
    F = a*d**2*x + c*e**2*x**9/9 + d*x**3*(2*a*e + b*d)/3 + e*x**7*(b*e + 2*c*d)/7 + x**5*(c*d**2/5 + e*(a*e + 2*b*d)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_277():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)/x
    F = a*d**2*log(x) + c*e**2*x**8/8 + d*x**2*(2*a*e + b*d)/2 + e*x**6*(b*e + 2*c*d)/6 + x**4*(c*d**2/4 + e*(a*e + 2*b*d)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_278():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)/x**2
    F = -a*d**2/x + c*e**2*x**7/7 + d*x*(2*a*e + b*d) + e*x**5*(b*e + 2*c*d)/5 + x**3*(c*d**2/3 + e*(a*e + 2*b*d)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_279():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)/x**3
    F = -a*d**2/(2*x**2) + c*e**2*x**6/6 + d*(2*a*e + b*d)*log(x) + e*x**4*(b*e + 2*c*d)/4 + x**2*(c*d**2/2 + e*(a*e + 2*b*d)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_280():
    f = x**6*(a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x**7/(7*e**2) + d**(sympy.S(3)/2)*(9*c*d**2 - e*(-5*a*e + 7*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*e**(sympy.S(11)/2)) - d**2*x*(a*e**2 - b*d*e + c*d**2)/(2*e**5*(d + e*x**2)) - d*x*(4*c*d**2 - e*(-2*a*e + 3*b*d))/e**5 - x**5*(-b*e + 2*c*d)/(5*e**3) + x**3*(3*c*d**2 - e*(-a*e + 2*b*d))/(3*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_281():
    f = x**4*(a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x**5/(5*e**2) - sqrt(d)*(7*c*d**2 - e*(-3*a*e + 5*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*e**(sympy.S(9)/2)) + d*x*(a*e**2 - b*d*e + c*d**2)/(2*e**4*(d + e*x**2)) - x**3*(-b*e + 2*c*d)/(3*e**3) + x*(3*c*d**2 - e*(-a*e + 2*b*d))/e**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_282():
    f = x**2*(a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x**3/(3*e**2) - x*(-b*e + 2*c*d)/e**3 - x*(a*e**2 - b*d*e + c*d**2)/(2*e**3*(d + e*x**2)) + (5*c*d**2 - e*(-a*e + 3*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*sqrt(d)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_283():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x/e**2 + x*(a + d*(-b*e + c*d)/e**2)/(2*d*(d + e*x**2)) - (3*c*d**2 - e*(a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_284():
    f = (a + b*x**2 + c*x**4)/(x**2*(d + e*x**2)**2)
    F = -a/(d**2*x) - x*(a*e**2 - b*d*e + c*d**2)/(2*d**2*e*(d + e*x**2)) + (c*d**2 + e*(-3*a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_285():
    f = (a + b*x**2 + c*x**4)/(x**4*(d + e*x**2)**2)
    F = -a/(3*d**2*x**3) + x*(a*e**2 - b*d*e + c*d**2)/(2*d**3*(d + e*x**2)) - (-2*a*e + b*d)/(d**3*x) + (c*d**2 - e*(-5*a*e + 3*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(7)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_286():
    f = (a + b*x**2 + c*x**4)/(x**6*(d + e*x**2)**2)
    F = -a/(5*d**2*x**5) - (-2*a*e + b*d)/(3*d**3*x**3) - e*x*(a*e**2 - b*d*e + c*d**2)/(2*d**4*(d + e*x**2)) - (c*d**2 - e*(-3*a*e + 2*b*d))/(d**4*x) - sqrt(e)*(3*c*d**2 - e*(-7*a*e + 5*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_287():
    f = (a + b*x**2 + c*x**4)/(x**8*(d + e*x**2)**2)
    F = -a/(7*d**2*x**7) - (-2*a*e + b*d)/(5*d**3*x**5) - (c*d**2 - e*(-3*a*e + 2*b*d))/(3*d**4*x**3) + e**2*x*(a*e**2 - b*d*e + c*d**2)/(2*d**5*(d + e*x**2)) + e*(2*c*d**2 - e*(-4*a*e + 3*b*d))/(d**5*x) + e**(sympy.S(3)/2)*(5*c*d**2 - e*(-9*a*e + 7*b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_288():
    f = x**6*(a + b*x**2 + c*x**4)/(d + e*x**2)**3
    F = c*x**5/(5*e**3) - sqrt(d)*(15*a*e**2 - 35*b*d*e + 63*c*d**2)*atan(sqrt(e)*x/sqrt(d))/(8*e**(sympy.S(11)/2)) - d**2*x*(a*e**2 - b*d*e + c*d**2)/(4*e**5*(d + e*x**2)**2) + d*x*(17*c*d**2 - e*(-9*a*e + 13*b*d))/(8*e**5*(d + e*x**2)) - x**3*(-b*e + 3*c*d)/(3*e**4) + x*(6*c*d**2 - e*(-a*e + 3*b*d))/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_289():
    f = x**4*(a + b*x**2 + c*x**4)/(d + e*x**2)**3
    F = c*x**3/(3*e**3) + d*x*(a*e**2 - b*d*e + c*d**2)/(4*e**4*(d + e*x**2)**2) - x*(-b*e + 3*c*d)/e**4 - x*(13*c*d**2 - e*(-5*a*e + 9*b*d))/(8*e**4*(d + e*x**2)) + (35*c*d**2 - 3*e*(-a*e + 5*b*d))*atan(sqrt(e)*x/sqrt(d))/(8*sqrt(d)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_290():
    f = x**2*(a + b*x**2 + c*x**4)/(d + e*x**2)**3
    F = c*x/e**3 - x*(a*e**2 - b*d*e + c*d**2)/(4*e**3*(d + e*x**2)**2) + x*(9*c*d**2 - e*(-a*e + 5*b*d))/(8*d*e**3*(d + e*x**2)) - (15*c*d**2 - e*(a*e + 3*b*d))*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(3)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_291():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**3
    F = x*(a + d*(-b*e + c*d)/e**2)/(4*d*(d + e*x**2)**2) - x*(5*c*d**2 - e*(3*a*e + b*d))/(8*d**2*e**2*(d + e*x**2)) + (3*c*d**2 + e*(3*a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_292():
    f = (a + b*x**2 + c*x**4)/(x**2*(d + e*x**2)**3)
    F = -a/(d**3*x) - x*(a*e**2 - b*d*e + c*d**2)/(4*d**2*e*(d + e*x**2)**2) + x*(c*d**2 + e*(-7*a*e + 3*b*d))/(8*d**3*e*(d + e*x**2)) + (c*d**2 + 3*e*(-5*a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(7)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_293():
    f = (a + b*x**2 + c*x**4)/(x**4*(d + e*x**2)**3)
    F = -a/(3*d**3*x**3) + x*(a*e**2 - b*d*e + c*d**2)/(4*d**3*(d + e*x**2)**2) + x*(3*c*d**2 - e*(-11*a*e + 7*b*d))/(8*d**4*(d + e*x**2)) - (-3*a*e + b*d)/(d**4*x) + (35*a*e**2 - 15*b*d*e + 3*c*d**2)*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(9)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_294():
    f = (a + b*x**2 + c*x**4)/(x**6*(d + e*x**2)**3)
    F = -a/(5*d**3*x**5) - e*x*(a*e**2 - b*d*e + c*d**2)/(4*d**4*(d + e*x**2)**2) - (-3*a*e + b*d)/(3*d**4*x**3) - e*x*(7*c*d**2 - e*(-15*a*e + 11*b*d))/(8*d**5*(d + e*x**2)) - (6*a*e**2 - 3*b*d*e + c*d**2)/(d**5*x) - sqrt(e)*(63*a*e**2 - 35*b*d*e + 15*c*d**2)*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_295():
    f = x**9/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = d**4*log(d + e*x**2)/(2*e**3*(a*e**2 - b*d*e + c*d**2)) + x**4/(4*c*e) - x**2*(b*e + c*d)/(2*c**2*e**2) - (a**2*c*e - a*b**2*e - 2*a*b*c*d + b**3*d)*log(a + b*x**2 + c*x**4)/(4*c**3*(a*e**2 - b*d*e + c*d**2)) - (3*a**2*b*c*e + 2*a**2*c**2*d - a*b**3*e - 4*a*b**2*c*d + b**4*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_296():
    f = x**7/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -d**3*log(d + e*x**2)/(2*e**2*(a*e**2 - b*d*e + c*d**2)) + x**2/(2*c*e) + (-a*b*e - a*c*d + b**2*d)*log(a + b*x**2 + c*x**4)/(4*c**2*(a*e**2 - b*d*e + c*d**2)) + (2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_297():
    f = x**5/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = d**2*log(d + e*x**2)/(2*e*(a*e**2 - b*d*e + c*d**2)) - (-a*e + b*d)*log(a + b*x**2 + c*x**4)/(4*c*(a*e**2 - b*d*e + c*d**2)) - (-a*b*e - 2*a*c*d + b**2*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_298():
    f = x**3/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = d*log(a + b*x**2 + c*x**4)/(4*a*e**2 - 4*b*d*e + 4*c*d**2) - d*log(d + e*x**2)/(2*a*e**2 - 2*b*d*e + 2*c*d**2) + (-2*a*e + b*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_299():
    f = x/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -e*log(a + b*x**2 + c*x**4)/(4*a*e**2 - 4*b*d*e + 4*c*d**2) + e*log(d + e*x**2)/(2*a*e**2 - 2*b*d*e + 2*c*d**2) - (-b*e + 2*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_300():
    f = 1/(x*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -e**2*log(d + e*x**2)/(2*d*(a*e**2 - b*d*e + c*d**2)) - (-b*e + c*d)*log(a + b*x**2 + c*x**4)/(4*a*(a*e**2 - b*d*e + c*d**2)) + (2*a*c*e - b**2*e + b*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + log(x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_301():
    f = 1/(x**3*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = e**3*log(d + e*x**2)/(2*d**2*(a*e**2 - b*d*e + c*d**2)) - 1/(2*a*d*x**2) + (a*c*e - b**2*e + b*c*d)*log(a + b*x**2 + c*x**4)/(4*a**2*(a*e**2 - b*d*e + c*d**2)) - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) - (a*e + b*d)*log(x)/(a**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_302():
    f = 1/(x**5*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -e**4*log(d + e*x**2)/(2*d**3*(a*e**2 - b*d*e + c*d**2)) - 1/(4*a*d*x**4) + (a*e + b*d)/(2*a**2*d**2*x**2) - (2*a*b*c*e - a*c**2*d - b**3*e + b**2*c*d)*log(a + b*x**2 + c*x**4)/(4*a**3*(a*e**2 - b*d*e + c*d**2)) + (-2*a**2*c**2*e + 4*a*b**2*c*e - 3*a*b*c**2*d - b**4*e + b**3*c*d)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**3*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + (a*b*d*e - a*(-a*e**2 + c*d**2) + b**2*d**2)*log(x)/(a**3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_303():
    f = x**8/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = d**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(e**(sympy.S(5)/2)*(a*e**2 - b*d*e + c*d**2)) + x**3/(3*c*e) - x*(b*e + c*d)/(c**2*e**2) - sqrt(2)*(a**2*c*e - a*b**2*e - 2*a*b*c*d + b**3*d + (3*a**2*b*c*e + 2*a**2*c**2*d - a*b**3*e - 4*a*b**2*c*d + b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*(a**2*c*e - a*b**2*e - 2*a*b*c*d + b**3*d - (3*a**2*b*c*e + 2*a**2*c**2*d - a*b**3*e - 4*a*b**2*c*d + b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_304():
    f = x**6/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -d**(sympy.S(5)/2)*atan(sqrt(e)*x/sqrt(d))/(e**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)) + x/(c*e) + sqrt(2)*(-a*b*e - a*c*d + b**2*d + (2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*(-a*b*e - a*c*d + b**2*d - (2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_305():
    f = x**4/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = d**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(e)*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*(-a*e + b*d + (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*(-a*e + b*d - (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_306():
    f = x**2/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = sqrt(2)*sqrt(c)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*sqrt(c)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(d)*sqrt(e)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 - b*d*e + c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_307():
    f = 1/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*sqrt(c)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_308():
    f = 1/(x**2*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -e**(sympy.S(5)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*sqrt(c)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*sqrt(c)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - 1/(a*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_309():
    f = 1/(x**4*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(d**(sympy.S(5)/2)*(a*e**2 - b*d*e + c*d**2)) - 1/(3*a*d*x**3) + sqrt(2)*sqrt(c)*(a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*sqrt(c)*(a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + (a*e + b*d)/(a**2*d**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_310():
    f = 1/(sqrt(f*x)*(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(f*x)/(sqrt(f)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(2*sqrt(f)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) - 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(f*x)/(sqrt(f)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(2*sqrt(f)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(f*x)/(sqrt(f)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(2*sqrt(f)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(f*x)/(sqrt(f)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(2*sqrt(f)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*e**(sympy.S(7)/4)*log(-sqrt(2)*d**(sympy.S(1)/4)*e**(sympy.S(1)/4)*sqrt(f*x) + sqrt(d)*sqrt(f) + sqrt(e)*sqrt(f)*x)/(4*d**(sympy.S(3)/4)*sqrt(f)*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*e**(sympy.S(7)/4)*log(sqrt(2)*d**(sympy.S(1)/4)*e**(sympy.S(1)/4)*sqrt(f*x) + sqrt(d)*sqrt(f) + sqrt(e)*sqrt(f)*x)/(4*d**(sympy.S(3)/4)*sqrt(f)*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*e**(sympy.S(7)/4)*atan(1 - sqrt(2)*e**(sympy.S(1)/4)*sqrt(f*x)/(d**(sympy.S(1)/4)*sqrt(f)))/(2*d**(sympy.S(3)/4)*sqrt(f)*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*e**(sympy.S(7)/4)*atan(1 + sqrt(2)*e**(sympy.S(1)/4)*sqrt(f*x)/(d**(sympy.S(1)/4)*sqrt(f)))/(2*d**(sympy.S(3)/4)*sqrt(f)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_311():
    f = x**5*sqrt(a + b*x**2 + c*x**4)/(d + e*x**2)
    F = d**2*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**4) + (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*c*e) + (-2*c*e*x**2*(b*e + 2*c*d) + (-b*e + 2*c*d)*(b*e + 4*c*d))*sqrt(a + b*x**2 + c*x**4)/(16*c**2*e**3) - (-b**3*e**3 - 2*b*c*e**2*(-2*a*e + b*d) + 16*c**3*d**3 - 8*c**2*d*e*(-a*e + b*d))*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(5)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_312():
    f = x**3*sqrt(a + b*x**2 + c*x**4)/(d + e*x**2)
    F = -d*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**3) - sqrt(a + b*x**2 + c*x**4)*(-b*e + 4*c*d - 2*c*e*x**2)/(8*c*e**2) + (-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-a*e + b*d))*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(3)/2)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_313():
    f = x*sqrt(a + b*x**2 + c*x**4)/(d + e*x**2)
    F = sqrt(a + b*x**2 + c*x**4)/(2*e) + sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**2) - (-b*e + 2*c*d)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_314():
    f = sqrt(a + b*x**2 + c*x**4)/(x*(d + e*x**2))
    F = -sqrt(a)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*d) + sqrt(c)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*e) - sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_315():
    f = sqrt(a + b*x**2 + c*x**4)/(x**3*(d + e*x**2))
    F = sqrt(a)*e*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*d**2) - b*e*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)*d**2) + sqrt(c)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*d) - sqrt(a + b*x**2 + c*x**4)/(2*d*x**2) + sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d**2) - (-b*e + 2*c*d)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)*d**2) - b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_316():
    f = sqrt(2*x**4 + 2*x**2 + 1)/(x**2*(2*x**2 + 3))
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))) * ((Integer(3) * x))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1)))))) + (Integer(-1) * (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(3) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(21) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(5) * ((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(252) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_317():
    f = sqrt(2*x**4 + 2*x**2 + 1)/(x**4*(2*x**2 + 3))
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(9))**(Integer(-1)) * sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) + (Integer(-1) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(9) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(5) * (Integer(3) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(63) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * ((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(378) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_318():
    f = sqrt(2*x**4 + 2*x**2 + 1)/(x**6*(2*x**2 + 3))
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))) * ((Integer(15) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(135) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(45) * x))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(45) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(27))**(Integer(-1))) * sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4) * (Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(45) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(5) * (Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(5) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(189) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(19) + (Integer(-1) * (Integer(2) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(135) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(5) * ((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(567) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_319():
    f = x**5*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(d + e*x**2)
    F = d**2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**6) + (a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*c*e) + (a + b*x**2 + c*x**4)**(sympy.S(3)/2)*(-3*b**2*e**2 - 6*b*c*d*e + 16*c**2*d**2 - 6*c*e*x**2*(b*e + 2*c*d))/(96*c**2*e**3) + sqrt(a + b*x**2 + c*x**4)*(3*b**4*e**4 + 6*b**2*c*e**3*(-2*a*e + b*d) + 8*b*c**2*d*e**2*(-3*a*e + 2*b*d) + 128*c**4*d**4 - 32*c**3*d**2*e*(-4*a*e + 5*b*d) - 2*c*e*x**2*(-3*b**3*e**3 - 6*b*c*e**2*(-2*a*e + b*d) + 32*c**3*d**3 - 8*c**2*d*e*(-3*a*e + 2*b*d)))/(256*c**3*e**5) - (3*b**5*e**5 + 6*b**3*c*e**4*(-4*a*e + b*d) + 16*b*c**2*e**3*(3*a**2*e**2 - 3*a*b*d*e + b**2*d**2) + 256*c**5*d**5 - 384*c**4*d**3*e*(-a*e + b*d) + 96*c**3*d*e**2*(-a*e + b*d)**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(512*c**(sympy.S(7)/2)*e**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_320():
    f = x**3*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(d + e*x**2)
    F = -d*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**5) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)*(-3*b*e + 8*c*d - 6*c*e*x**2)/(48*c*e**2) - sqrt(a + b*x**2 + c*x**4)*(3*b**3*e**3 + 4*b*c*e**2*(-3*a*e + 2*b*d) + 64*c**3*d**3 - 16*c**2*d*e*(-4*a*e + 5*b*d) - 2*c*e*x**2*(-3*b**2*e**2 + 16*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d)))/(128*c**2*e**4) + (3*b**4*e**4 + 8*b**2*c*e**3*(-3*a*e + b*d) + 128*c**4*d**4 - 192*c**3*d**2*e*(-a*e + b*d) + 48*c**2*e**2*(-a*e + b*d)**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(256*c**(sympy.S(5)/2)*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_321():
    f = x*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(d + e*x**2)
    F = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*e) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**4) + sqrt(a + b*x**2 + c*x**4)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x**2*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(16*c*e**3) - (-b*e + 2*c*d)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(3)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_322():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(x*(d + e*x**2))
    F = -a**(sympy.S(3)/2)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*d) + a*b*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)*d) + a*sqrt(a + b*x**2 + c*x**4)/(2*d) - sqrt(a + b*x**2 + c*x**4)*(4*c*d**2 - 2*c*d*e*x**2 - e*(-4*a*e + 5*b*d))/(8*d*e**2) - (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d*e**3) + (b*e**2*(-4*a*e + 3*b*d) + 8*c**2*d**3 - 12*c*d*e*(-a*e + b*d))*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*sqrt(c)*d*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_323():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(x**3*(d + e*x**2))
    F = a**(sympy.S(3)/2)*e*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*d**2) - 3*sqrt(a)*b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*d) + b*e*(-12*a*c + b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(3)/2)*d**2) + (9*b + 6*c*x**2)*sqrt(a + b*x**2 + c*x**4)/(8*d) - (a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(2*d*x**2) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d**2*e**2) - e*sqrt(a + b*x**2 + c*x**4)*(8*a*c + b**2 + 2*b*c*x**2)/(16*c*d**2) + sqrt(a + b*x**2 + c*x**4)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x**2*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(16*c*d**2*e) + (12*a*c + 3*b**2)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*sqrt(c)*d) - (-b*e + 2*c*d)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(3)/2)*d**2*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_324():
    f = (2*x**4 + 2*x**2 + 1)**(sympy.S(3)/2)/(x**2*(3 - 2*x**2))
    F = (Integer(-1) * (((Integer(1) + (x)**(Integer(2))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(17) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * sympy.sqrt(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(17) * (Integer(12))**(Integer(-1))) * sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) + ((Integer(17) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(289) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(84) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(17) * (Integer(5) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(12) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(289) * (Integer(11) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(11) * sympy.sqrt(Integer(2))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(504) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_325():
    f = (2*x**4 + 2*x**2 + 1)**(sympy.S(3)/2)/(x**4*(3 - 2*x**2))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * (x)**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Integer(-1) * (Integer(8) * (x)**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(9) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(17) * (Integer(18))**(Integer(-1))) * sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) + (Integer(-1) * (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(9) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(289) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(126) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(17) * (Integer(5) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(18) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(9) + (Integer(5) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(9) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(289) * (Integer(11) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(11) * sympy.sqrt(Integer(2))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(756) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_326():
    f = (2*x**4 + 2*x**2 + 1)**(sympy.S(3)/2)/(x**6*(3 - 2*x**2))
    F = ((Integer(74) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(135) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(262) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(135) * x))**(Integer(-1)))) + (Integer(-1) * (((Integer(3) + (Integer(40) * (x)**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(45) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(262) * sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(135) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(17) * (Integer(27))**(Integer(-1))) * sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt((Integer(17) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(262) * (Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(135) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(85) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(189) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(37) + (Integer(23) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(135) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(289) * (Integer(11) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(11) * sympy.sqrt(Integer(2))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(1134) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_327():
    f = x**5/((d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = d**2*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e**2*sqrt(a*e**2 - b*d*e + c*d**2)) + sqrt(a + b*x**2 + c*x**4)/(2*c*e) - (b*e + 2*c*d)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*c**(sympy.S(3)/2)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_328():
    f = x**3/((d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = -d*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e*sqrt(a*e**2 - b*d*e + c*d**2)) + atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(c)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_329():
    f = x/((d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*sqrt(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_330():
    f = 1/(x*(d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = -e*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d*sqrt(a*e**2 - b*d*e + c*d**2)) - atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_331():
    f = 1/(x**3*(d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = e**2*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d**2*sqrt(a*e**2 - b*d*e + c*d**2)) - sqrt(a + b*x**2 + c*x**4)/(2*a*d*x**2) + e*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a)*d**2) + b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_332():
    f = x**4/((2*x**2 + 3)*sqrt(2*x**4 + 2*x**2 + 1))
    F = ((x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(3) * (Integer(10))**(Integer(-1)))) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(2)))) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) * ((Integer(4) * (Integer(2) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(2) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(1) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(2) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(2) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(3) * (Integer(3) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(8) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Integer(2) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_333():
    f = x**2/((2*x**2 + 3)*sqrt(2*x**4 + 2*x**2 + 1))
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(3) * (Integer(5))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) + (Integer(-1) * (((Integer(3) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(14) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + ((((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(56) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_334():
    f = 1/((2*x**2 + 3)*sqrt(2*x**4 + 2*x**2 + 1))
    F = (sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(Integer(15))))**(Integer(-1))) + (((Integer(3) + sympy.sqrt(Integer(2))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(14) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(84) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_335():
    f = 1/(x**2*(2*x**2 + 3)*sqrt(2*x**4 + 2*x**2 + 1))
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))) * ((Integer(3) * x))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1)))) * ((Integer(3) * sympy.sqrt(Integer(15))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(5) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(2))))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(21) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + ((((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(126) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_336():
    f = 1/(x**4*(2*x**2 + 3)*sqrt(2*x**4 + 2*x**2 + 1))
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))) * ((Integer(3) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(2) * sympy.atan(((sympy.sqrt((Integer(5) * (Integer(3))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4))))))**(Integer(-1))))) * ((Integer(9) * sympy.sqrt(Integer(15))))**(Integer(-1))) + ((Integer(2) * (Integer(2))**((Integer(4))**(Integer(-1))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(1) + (Integer(19) * sympy.sqrt(Integer(2)))) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(63) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) + sympy.sqrt(Integer(2))))**(Integer(2)) * (Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))) * (((Integer(1) + (sympy.sqrt(Integer(2)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(24))**(Integer(-1)) * (Integer(12) + (Integer(-1) * (Integer(11) * sympy.sqrt(Integer(2)))))), (Integer(2) * sympy.atan(((Integer(2))**((Integer(4))**(Integer(-1))) * x))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2))))))) * ((Integer(189) * (Integer(2))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(2) * (x)**(Integer(2))) + (Integer(2) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_337():
    f = x**7/((d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -d**3*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*e*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) + (a*(-a*b*e - 2*a*c*d + b**2*d) + x**2*(2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2)) + atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_338():
    f = x**5/((d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = d**2*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - (a*(-2*a*e + b*d) + x**2*(-a*b*e - 2*a*c*d + b**2*d))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_339():
    f = x**3/((d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -d*e*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) + (a*(-b*e + 2*c*d) + c*x**2*(-2*a*e + b*d))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_340():
    f = x/((d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = e**2*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - (2*a*c*e - b**2*e + b*c*d + c*x**2*(-b*e + 2*c*d))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_341():
    f = 1/(x*(d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = -e**3*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) + e*(2*a*c*e - b**2*e + b*c*d + c*x**2*(-b*e + 2*c*d))/(d*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2)) + (-2*a*c + b**2 + b*c*x**2)/(a*d*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_342():
    f = 1/(x**3*(d + e*x**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2))
    F = e**4*atanh((-2*a*e + b*d + x**2*(-b*e + 2*c*d))/(2*sqrt(a + b*x**2 + c*x**4)*sqrt(a*e**2 - b*d*e + c*d**2)))/(2*d**2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - e**2*(2*a*c*e - b**2*e + b*c*d + c*x**2*(-b*e + 2*c*d))/(d**2*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2)) + (-2*a*c + b**2 + b*c*x**2)/(a*d*x**2*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - e*(-2*a*c + b**2 + b*c*x**2)/(a*d**2*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - (-8*a*c + 3*b**2)*sqrt(a + b*x**2 + c*x**4)/(2*a**2*d*x**2*(-4*a*c + b**2)) + e*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*a**(sympy.S(3)/2)*d**2) + 3*b*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(4*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_343():
    f = x**7*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = (d + e*x**2)**(sympy.S(5)/2)/(5*c*e**2) - (d + e*x**2)**(sympy.S(3)/2)*(b*e + c*d)/(3*c**2*e**2) + sqrt(d + e*x**2)*(-a*c + b**2)/c**3 - sqrt(2)*(2*a*b*c*e - a*c**2*d - b**3*e + b**2*c*d + (-2*a**2*c**2*e + 4*a*b**2*c*e - 3*a*b*c**2*d - b**4*e + b**3*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(7)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*(2*a*b*c*e - a*c**2*d - b**3*e + b**2*c*d - (-2*a**2*c**2*e + 4*a*b**2*c*e - 3*a*b*c**2*d - b**4*e + b**3*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(7)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_344():
    f = x**5*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = -b*sqrt(d + e*x**2)/c**2 + (d + e*x**2)**(sympy.S(3)/2)/(3*c*e) + sqrt(2)*(a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(5)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(5)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_345():
    f = x**3*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = sqrt(d + e*x**2)/c - sqrt(2)*(2*a*c*e - b**2*e + b*c*d + sqrt(-4*a*c + b**2)*(-b*e + c*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(2*a*c*e - b**2*e + b*c*d - sqrt(-4*a*c + b**2)*(-b*e + c*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_346():
    f = x*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_347():
    f = sqrt(d + e*x**2)/(x*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(-2*a*e + b*d - d*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*sqrt(c)*(-2*a*e + b*d + d*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - sqrt(d)*atanh(sqrt(d + e*x**2)/sqrt(d))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_348():
    f = sqrt(d + e*x**2)/(x**3*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(2*a*x**2) + e*atanh(sqrt(d + e*x**2)/sqrt(d))/(2*a*sqrt(d)) + sqrt(2)*sqrt(c)*(-a*(2*c*d - e*sqrt(-4*a*c + b**2)) + b**2*d - b*(a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*sqrt(c)*(-a*b*e - 2*a*c*d + b**2*d + sqrt(-4*a*c + b**2)*(-a*e + b*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (-a*e + b*d)*atanh(sqrt(d + e*x**2)/sqrt(d))/(a**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_349():
    f = sqrt(d + e*x**2)/(x**5*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(4*a*x**4) + 3*e*sqrt(d + e*x**2)/(8*a*d*x**2) - 3*e**2*atanh(sqrt(d + e*x**2)/sqrt(d))/(8*a*d**(sympy.S(3)/2)) + sqrt(d + e*x**2)*(-a*e + b*d)/(2*a**2*d*x**2) - e*(-a*e + b*d)*atanh(sqrt(d + e*x**2)/sqrt(d))/(2*a**2*d**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-a*b*(3*c*d - e*sqrt(-4*a*c + b**2)) + a*c*(2*a*e + d*sqrt(-4*a*c + b**2)) + b**3*d - b**2*(a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*sqrt(c)*(-a*b*(3*c*d + e*sqrt(-4*a*c + b**2)) - a*c*(-2*a*e + d*sqrt(-4*a*c + b**2)) + b**3*d + b**2*(-a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (-a*b*e - a*c*d + b**2*d)*atanh(sqrt(d + e*x**2)/sqrt(d))/(a**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_350():
    f = x**4*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = x*sqrt(d + e*x**2)/(2*c) - (a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - (a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (-2*b*e + c*d)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**2*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_351():
    f = x**2*sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = sqrt(e)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/c + (-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + (-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_352():
    f = sqrt(d + e*x**2)/(a + b*x**2 + c*x**4)
    F = -sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_353():
    f = sqrt(d + e*x**2)/(x**2*(a + b*x**2 + c*x**4))
    F = -c*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - sqrt(d + e*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_354():
    f = sqrt(d + e*x**2)/(x**4*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(3*a*x**3) + 2*e*sqrt(d + e*x**2)/(3*a*d*x) + c*(-a*e + b*d - (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + c*(-a*e + b*d + (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + sqrt(d + e*x**2)*(-a*e + b*d)/(a**2*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_355():
    f = sqrt(d + e*x**2)/(x**6*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(5*a*x**5) + 4*e*sqrt(d + e*x**2)/(15*a*d*x**3) - 8*e**2*sqrt(d + e*x**2)/(15*a*d**2*x) + sqrt(d + e*x**2)*(-a*e + b*d)/(3*a**2*d*x**3) - 2*e*sqrt(d + e*x**2)*(-a*e + b*d)/(3*a**2*d**2*x) - c*(-a*b*e - a*c*d + b**2*d - (2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**3*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(-a*b*e - a*c*d + b**2*d + (2*a**2*c*e - a*b**2*e - 3*a*b*c*d + b**3*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**3*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - sqrt(d + e*x**2)*(-a*b*e - a*c*d + b**2*d)/(a**3*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_356():
    f = x**3*(d + e*x**2)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = (d + e*x**2)**(sympy.S(3)/2)/(3*c) + sqrt(d + e*x**2)*(-b*e + c*d)/c**2 - sqrt(2)*(b**3*e**2 - b**2*e*(2*c*d - e*sqrt(-4*a*c + b**2)) + b*c*(c*d**2 - e*(3*a*e + 2*d*sqrt(-4*a*c + b**2))) - c*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(4*a*e + d*sqrt(-4*a*c + b**2))))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(5)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(b**3*e**2 - b**2*e*(2*c*d + e*sqrt(-4*a*c + b**2)) + b*c*(c*d**2 + e*(-3*a*e + 2*d*sqrt(-4*a*c + b**2))) + c*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(-4*a*e + d*sqrt(-4*a*c + b**2))))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(5)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_357():
    f = x*(d + e*x**2)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = e*sqrt(d + e*x**2)/c + sqrt(2)*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d - d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_358():
    f = (d + e*x**2)**(sympy.S(3)/2)/(x*(a + b*x**2 + c*x**4))
    F = -d**(sympy.S(3)/2)*atanh(sqrt(d + e*x**2)/sqrt(d))/a - sqrt(2)*(a*e**2*sqrt(-4*a*c + b**2) + b*(a*e**2 + c*d**2) - c*d*(4*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a*sqrt(c)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*(a*e**2*sqrt(-4*a*c + b**2) - b*(a*e**2 + c*d**2) - c*d*(-4*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a*sqrt(c)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_359():
    f = (d + e*x**2)**(sympy.S(3)/2)/(x**3*(a + b*x**2 + c*x**4))
    F = sqrt(d)*e*atanh(sqrt(d + e*x**2)/sqrt(d))/(2*a) - d*sqrt(d + e*x**2)/(2*a*x**2) + sqrt(2)*sqrt(c)*(-2*a*(c*d**2 - e*(a*e + d*sqrt(-4*a*c + b**2))) + b**2*d**2 - b*d*(2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*sqrt(c)*(-2*a*(c*d**2 + e*(-a*e + d*sqrt(-4*a*c + b**2))) + b**2*d**2 + b*d*(-2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x**2)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + sqrt(d)*(-2*a*e + b*d)*atanh(sqrt(d + e*x**2)/sqrt(d))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_360():
    f = x**4*(d + e*x**2)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = x*(d + e*x**2)**(sympy.S(3)/2)/(4*c) + d*(-4*b*e + 3*c*d)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(8*c**2*sqrt(e)) + x*sqrt(d + e*x**2)*(-4*b*e + 3*c*d)/(8*c**2) - sqrt(e)*(a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**3) - sqrt(e)*(a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**3) - sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*c**3*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*c**3*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_361():
    f = x**2*(d + e*x**2)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = d*sqrt(e)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c) + e*x*sqrt(d + e*x**2)/(2*c) + sqrt(e)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**2) + sqrt(e)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**2) + sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*c**2*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*c**2*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_362():
    f = (d + e*x**2)**(sympy.S(3)/2)/(a + b*x**2 + c*x**4)
    F = sqrt(e)*(3*c*d - e*(b - sqrt(-4*a*c + b**2)))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c*sqrt(-4*a*c + b**2)) - sqrt(e)*(3*c*d - e*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c*sqrt(-4*a*c + b**2)) - (b*e**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d + d*sqrt(-4*a*c + b**2)))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + (b*e**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d - d*sqrt(-4*a*c + b**2)))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_363():
    f = (d + e*x**2)**(sympy.S(3)/2)/(x**2*(a + b*x**2 + c*x**4))
    F = (2*c*d - e*(b + sqrt(-4*a*c + b**2)))**(sympy.S(3)/2)*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/((b + sqrt(-4*a*c + b**2))**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)) - (2*c*d - e*(b - sqrt(-4*a*c + b**2)))**(sympy.S(3)/2)*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/((b - sqrt(-4*a*c + b**2))**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)) - d*sqrt(d + e*x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_364():
    f = (d + e*x**2)**(sympy.S(3)/2)/(x**4*(a + b*x**2 + c*x**4))
    F = -(d + e*x**2)**(sympy.S(3)/2)/(3*a*x**3) - sqrt(e)*(-a*e + b*d)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/a**2 + sqrt(e)*(-a*e + b*d - (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*a**2) + sqrt(e)*(-a*e + b*d + (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*a**2) + sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(-a*e + b*d - (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*a**2*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(-a*e + b*d + (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(2*a**2*sqrt(b - sqrt(-4*a*c + b**2))) + sqrt(d + e*x**2)*(-a*e + b*d)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_365():
    f = x**5*sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = -b*sqrt(1 - x**2)/c**2 - (1 - x**2)**(sympy.S(3)/2)/(3*c) + sqrt(2)*(-a*c + b**2 + b*c + (-3*a*b*c - 2*a*c**2 + b**3 + b**2*c)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + sqrt(2)*(-a*c + b**2 + b*c - (-3*a*b*c - 2*a*c**2 + b**3 + b**2*c)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_366():
    f = x**3*sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = sqrt(1 - x**2)/c - sqrt(2)*(b + c - (-2*a*c + b**2 + b*c)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + c + (-2*a*c + b**2 + b*c)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_367():
    f = x*sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(b + 2*c - sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b + 2*c + sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_368():
    f = sqrt(1 - x**2)/(x*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(2*a + b - sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(2*a*sqrt(-4*a*c + b**2)*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + sqrt(2)*sqrt(c)*(2*a + b + sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(2*a*sqrt(-4*a*c + b**2)*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - atanh(sqrt(1 - x**2))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_369():
    f = sqrt(1 - x**2)/(x**3*(a + b*x**2 + c*x**4))
    F = 1/(4*a*(sqrt(1 - x**2) + 1)) - 1/(4*a*(1 - sqrt(1 - x**2))) - sqrt(2)*sqrt(c)*(a + b - (a*(b - 2*c) + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(a + b + (a*(b - 2*c) + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(1 - x**2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(2*a**2*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + (a + 2*b)*atanh(sqrt(1 - x**2))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_370():
    f = x**4*sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = x*sqrt(1 - x**2)/(2*c) + (2*b + c)*asin(x)/(2*c**2) - (-a*c + b**2 + b*c + (-3*a*b*c - 2*a*c**2 + b**3 + b**2*c)/sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c + sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b + sqrt(-4*a*c + b**2))))/(c**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - (-a*c + b**2 + b*c - (-3*a*b*c - 2*a*c**2 + b**3 + b**2*c)/sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c - sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b - sqrt(-4*a*c + b**2))))/(c**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_371():
    f = x**2*sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = -asin(x)/c + (b + c + (-2*a*c + b**2 + b*c)/sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c + sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b + sqrt(-4*a*c + b**2))))/(c*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + (b + c - (-2*a*c + b**2 + b*c)/sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c - sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b - sqrt(-4*a*c + b**2))))/(c*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_372():
    f = sqrt(1 - x**2)/(a + b*x**2 + c*x**4)
    F = -sqrt(b + 2*c + sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c + sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b + sqrt(-4*a*c + b**2))))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(b + 2*c - sqrt(-4*a*c + b**2))*atan(x*sqrt(b + 2*c - sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b - sqrt(-4*a*c + b**2))))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_373():
    f = sqrt(1 - x**2)/(x**2*(a + b*x**2 + c*x**4))
    F = -c*(-(2*a + b)/sqrt(-4*a*c + b**2) + 1)*atan(x*sqrt(b + 2*c + sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b + sqrt(-4*a*c + b**2))))/(a*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - c*((2*a + b)/sqrt(-4*a*c + b**2) + 1)*atan(x*sqrt(b + 2*c - sqrt(-4*a*c + b**2))/(sqrt(1 - x**2)*sqrt(b - sqrt(-4*a*c + b**2))))/(a*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - sqrt(1 - x**2)/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_374():
    f = x**2*sqrt(1 - x**2)/(x**4 + x**2 - 1)
    F = -asin(x) + sqrt(sympy.S(2)/5 + sqrt(5)/5)*atan(x*sqrt(sympy.S.Half + sqrt(5)/2)/sqrt(1 - x**2)) - sqrt(sympy.S(-2)/5 + sqrt(5)/5)*atanh(x*sqrt(sympy.S(-1)/2 + sqrt(5)/2)/sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_375():
    f = x**8/(sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = b*d*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**2*e**(sympy.S(3)/2)) - b*x*sqrt(d + e*x**2)/(2*c**2*e) + 3*d**2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(8*c*e**(sympy.S(5)/2)) - 3*d*x*sqrt(d + e*x**2)/(8*c*e**2) + x**3*sqrt(d + e*x**2)/(4*c*e) - (-2*a*b*c + b**3 + (2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**3*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - (-2*a*b*c + b**3 - (2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**3*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (-a*c + b**2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(c**3*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_376():
    f = x**6/(sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -b*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(c**2*sqrt(e)) - d*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c*e**(sympy.S(3)/2)) + x*sqrt(d + e*x**2)/(2*c*e) + (-a*c + b**2 + b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + (-a*c + b**2 - b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_377():
    f = x**4/(sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(c*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + atanh(sqrt(e)*x/sqrt(d + e*x**2))/(c*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_378():
    f = x**2/(sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -sqrt(b - sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + sqrt(b + sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_379():
    f = 1/(sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -2*c*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + 2*c*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_380():
    f = 1/(x**2*sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -c*(-b/sqrt(-4*a*c + b**2) + 1)*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(b/sqrt(-4*a*c + b**2) + 1)*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - sqrt(d + e*x**2)/(a*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_381():
    f = 1/(x**4*sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(3*a*d*x**3) + 2*e*sqrt(d + e*x**2)/(3*a*d**2*x) + b*sqrt(d + e*x**2)/(a**2*d*x) + c*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + c*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_382():
    f = 1/(x**6*sqrt(d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -sqrt(d + e*x**2)/(5*a*d*x**5) + 4*e*sqrt(d + e*x**2)/(15*a*d**2*x**3) - 8*e**2*sqrt(d + e*x**2)/(15*a*d**3*x) + b*sqrt(d + e*x**2)/(3*a**2*d*x**3) - 2*b*e*sqrt(d + e*x**2)/(3*a**2*d**2*x) - c*(-a*c + b**2 - b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**3*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(-a*c + b**2 + b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(a**3*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - sqrt(d + e*x**2)*(-a*c + b**2)/(a**3*d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_383():
    f = x**4/((d + e*x**2)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4))
    F = d*x/(sqrt(d + e*x**2)*(a*e**2 - b*d*e + c*d**2)) - (-a*e + b*d + (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2)) - (-a*e + b*d - (-a*b*e - 2*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_384():
    f = x**2/((d + e*x**2)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4))
    F = c*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2)) + c*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2)) - e*x/(sqrt(d + e*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_385():
    f = 1/((d + e*x**2)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4))
    F = -c*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2)) - c*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(x*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(d + e*x**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(a*e**2 - b*d*e + c*d**2)) + e**2*x/(d*sqrt(d + e*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_386():
    f = (f*x)**m*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -2*c*(f*x)**(m + 1)*(d + e*x**2)**q*appellf1(m/2 + sympy.S.Half, 1, -q, m/2 + sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/(f*(1 + e*x**2/d)**q*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*(f*x)**(m + 1)*(d + e*x**2)**q*appellf1(m/2 + sympy.S.Half, 1, -q, m/2 + sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/(f*(1 + e*x**2/d)**q*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_387():
    f = x**7*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = (d + e*x**2)**(q + 1)*(a - b**2/c - b*(-3*a*c + b**2)/(c*sqrt(-4*a*c + b**2)))*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c*(q + 1)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + (d + e*x**2)**(q + 1)*(a - b**2/c + b*(-3*a*c + b**2)/(c*sqrt(-4*a*c + b**2)))*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c*(q + 1)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (d + e*x**2)**(q + 2)/(2*c*e**2*(q + 2)) - (d + e*x**2)**(q + 1)*(b*e + c*d)/(2*c**2*e**2*(q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_388():
    f = x**5*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = (b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*c*(q + 1)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c*(q + 1)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + (d + e*x**2)**(q + 1)/(2*c*e*(q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_389():
    f = x**3*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -(d + e*x**2)**(q + 1)*(-b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/((q + 1)*(4*c*d - 2*e*(b - sqrt(-4*a*c + b**2)))) - (d + e*x**2)**(q + 1)*(b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/((q + 1)*(4*c*d - 2*e*(b + sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_390():
    f = x*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = c*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/((q + 1)*sqrt(-4*a*c + b**2)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/((q + 1)*sqrt(-4*a*c + b**2)*(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_391():
    f = (d + e*x**2)**q/(x*(a + b*x**2 + c*x**4))
    F = c*(d + e*x**2)**(q + 1)*(-b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a*(q + 1)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + c*(d + e*x**2)**(q + 1)*(b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a*(q + 1)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 1 + e*x**2/d)/(2*a*d*(q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_392():
    f = (d + e*x**2)**q/(x**3*(a + b*x**2 + c*x**4))
    F = e*(d + e*x**2)**(q + 1)*hyper((2, q + 1), (q + 2,), 1 + e*x**2/d)/(2*a*d**2*(q + 1)) + b*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 1 + e*x**2/d)/(2*a**2*d*(q + 1)) - c*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*a**2*(q + 1)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - c*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**(q + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**2)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(2*a**2*(q + 1)*(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_393():
    f = x**6*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -b*x*(d + e*x**2)**q*hyper((sympy.S.Half, -q), (sympy.S(3)/2,), -e*x**2/d)/(c**2*(1 + e*x**2/d)**q) + x**3*(d + e*x**2)**q*hyper((sympy.S(3)/2, -q), (sympy.S(5)/2,), -e*x**2/d)/(3*c*(1 + e*x**2/d)**q) + x*(d + e*x**2)**q*(-a*c + b**2 + b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/(c**2*(1 + e*x**2/d)**q*(b + sqrt(-4*a*c + b**2))) + x*(d + e*x**2)**q*(-a*c + b**2 - b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/(c**2*(1 + e*x**2/d)**q*(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_394():
    f = x**4*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -x*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/(c*(1 + e*x**2/d)**q*(b - sqrt(-4*a*c + b**2))) - x*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/(c*(1 + e*x**2/d)**q*(b + sqrt(-4*a*c + b**2))) + x*(d + e*x**2)**q*hyper((sympy.S.Half, -q), (sympy.S(3)/2,), -e*x**2/d)/(c*(1 + e*x**2/d)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_395():
    f = x**2*(d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -x*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/((1 + e*x**2/d)**q*sqrt(-4*a*c + b**2)) + x*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/((1 + e*x**2/d)**q*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_396():
    f = (d + e*x**2)**q/(a + b*x**2 + c*x**4)
    F = -2*c*x*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/((1 + e*x**2/d)**q*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - 2*c*x*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/((1 + e*x**2/d)**q*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_397():
    f = (d + e*x**2)**q/(x**2*(a + b*x**2 + c*x**4))
    F = -c*x*(d + e*x**2)**q*(-b/sqrt(-4*a*c + b**2) + 1)*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/(a*(1 + e*x**2/d)**q*(b + sqrt(-4*a*c + b**2))) - c*x*(d + e*x**2)**q*(b/sqrt(-4*a*c + b**2) + 1)*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/(a*(1 + e*x**2/d)**q*(b - sqrt(-4*a*c + b**2))) - (d + e*x**2)**q*hyper((sympy.S(-1)/2, -q), (sympy.S.Half,), -e*x**2/d)/(a*x*(1 + e*x**2/d)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_4_f_x_pow_m_d_plus_e_x_pow_2_pow_q_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_398():
    f = (d + e*x**2)**q/(x**4*(a + b*x**2 + c*x**4))
    F = -(d + e*x**2)**q*hyper((sympy.S(-3)/2, -q), (sympy.S(-1)/2,), -e*x**2/d)/(3*a*x**3*(1 + e*x**2/d)**q) + b*(d + e*x**2)**q*hyper((sympy.S(-1)/2, -q), (sympy.S.Half,), -e*x**2/d)/(a**2*x*(1 + e*x**2/d)**q) + c*x*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b + sqrt(-4*a*c + b**2)), -e*x**2/d)/(a**2*(1 + e*x**2/d)**q*(b + sqrt(-4*a*c + b**2))) + c*x*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(d + e*x**2)**q*appellf1(sympy.S.Half, 1, -q, sympy.S(3)/2, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -e*x**2/d)/(a**2*(1 + e*x**2/d)**q*(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F

