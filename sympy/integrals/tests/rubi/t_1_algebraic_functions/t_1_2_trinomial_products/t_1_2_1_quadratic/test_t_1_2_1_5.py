"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.5 (a+b x+c x^2)^p (d+e x+f x^2)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, p, q = symbols('a b c d e f p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_1():
    f = (a + b*x + b*f*x**2/e)/sqrt(d + e*x + f*x**2)
    F = b*sqrt(d + e*x + f*x**2)/(4*f) + b*x*sqrt(d + e*x + f*x**2)/(2*e) + (8*a*f - b*(4*d*f/e + e))*atanh((e + 2*f*x)/(2*sqrt(f)*sqrt(d + e*x + f*x**2)))/(8*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_2():
    f = 1/((a + b*x + b*f*x**2/e)*sqrt(d + e*x + f*x**2))
    F = -2*sqrt(e)*atanh((e + 2*f*x)*sqrt(-a*e + b*d)/(sqrt(e)*sqrt(-4*a*f + b*e)*sqrt(d + e*x + f*x**2)))/(sqrt(-a*e + b*d)*sqrt(-4*a*f + b*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_3():
    f = 1/(sqrt(a + b*x + c*x**2)*(b*x + c*x**2 + d))
    F = -2*atanh(sqrt(a - d)*(b + 2*c*x)/(sqrt(b**2 - 4*c*d)*sqrt(a + b*x + c*x**2)))/(sqrt(a - d)*sqrt(b**2 - 4*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_4():
    f = 1/(sqrt(a + b*x + c*x**2)*(b*x + c*x**2 + d)**2)
    F = -(b + 2*c*x)*sqrt(a + b*x + c*x**2)/((a - d)*(b**2 - 4*c*d)*(b*x + c*x**2 + d)) + (b**2 + 4*c*(a - 2*d))*atanh(sqrt(a - d)*(b + 2*c*x)/(sqrt(b**2 - 4*c*d)*sqrt(a + b*x + c*x**2)))/((a - d)**(sympy.S(3)/2)*(b**2 - 4*c*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_5():
    f = 1/(sqrt(a + b*x + c*x**2)*(b*x + c*x**2 + d)**3)
    F = -(b + 2*c*x)*sqrt(a + b*x + c*x**2)/((2*a - 2*d)*(b**2 - 4*c*d)*(b*x + c*x**2 + d)**2) + (b + 2*c*x)*(3*b**2 + 12*c*(a - 2*d))*sqrt(a + b*x + c*x**2)/(4*(a - d)**2*(b**2 - 4*c*d)**2*(b*x + c*x**2 + d)) - (3*b**4 + 8*b**2*c*(a - 4*d) + 16*c**2*(3*a**2 - 8*a*d + 8*d**2))*atanh(sqrt(a - d)*(b + 2*c*x)/(sqrt(b**2 - 4*c*d)*sqrt(a + b*x + c*x**2)))/(4*(a - d)**(sympy.S(5)/2)*(b**2 - 4*c*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_6():
    f = 1/(sqrt(a + b*x + c*x**2)*(b*x + c*x**2 + d)**4)
    F = -(b + 2*c*x)*sqrt(a + b*x + c*x**2)/((3*a - 3*d)*(b**2 - 4*c*d)*(b*x + c*x**2 + d)**3) + (b + 2*c*x)*(5*b**2 + 20*c*(a - 2*d))*sqrt(a + b*x + c*x**2)/(12*(a - d)**2*(b**2 - 4*c*d)**2*(b*x + c*x**2 + d)**2) - (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(15*b**4 + 8*b**2*c*(7*a - 22*d) + 16*c**2*(15*a**2 - 44*a*d + 44*d**2))/(24*(a - d)**3*(b**2 - 4*c*d)**3*(b*x + c*x**2 + d)) + (b**2 + 4*c*(a - 2*d))*(5*b**4 - 8*b**2*c*(a + 4*d) + 16*c**2*(5*a**2 - 8*a*d + 8*d**2))*atanh(sqrt(a - d)*(b + 2*c*x)/(sqrt(b**2 - 4*c*d)*sqrt(a + b*x + c*x**2)))/(8*(a - d)**(sympy.S(7)/2)*(b**2 - 4*c*d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_7():
    f = 1/(sqrt(d + e*x + f*x**2)*(a*e + b*e*x + b*f*x**2)**2)
    F = -b*(e + 2*f*x)*sqrt(d + e*x + f*x**2)/(e*(-a*e + b*d)*(-4*a*f + b*e)*(a*e + b*e*x + b*f*x**2)) - (8*a*e*f - b*(4*d*f + e**2))*atanh((e + 2*f*x)*sqrt(-a*e + b*d)/(sqrt(e)*sqrt(-4*a*f + b*e)*sqrt(d + e*x + f*x**2)))/(e**(sympy.S(3)/2)*(-a*e + b*d)**(sympy.S(3)/2)*(-4*a*f + b*e)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_8():
    f = 1/((x**2 + 2*x + 4)*sqrt(x**2 + 2*x + 5))
    F = sqrt(3)*atan(sqrt(3)*(x + 1)/(3*sqrt(x**2 + 2*x + 5)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_9():
    f = (a + c*x**2 + e*x/2)**p*(2*a + 2*c*x**2 + e*x)**q
    F = -2**(q + 1)*(-(4*c*x + e - sqrt(-16*a*c + e**2))/sqrt(-16*a*c + e**2))**(-p - q - 1)*(2*a + 2*c*x**2 + e*x)**(p + q + 1)*hyper((p + q + 1, -p - q), (p + q + 2,), (4*c*x + e + sqrt(-16*a*c + e**2))/(2*sqrt(-16*a*c + e**2)))/(sqrt(-16*a*c + e**2)*(p + q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_10():
    f = (a + c*e*x/f + c*x**2)**p*(a*f/c + e*x + f*x**2)**q
    F = -2**(p + q + 1)*sqrt(c)*(-sqrt(c)*(e + 2*f*x - sqrt(-4*a*f**2 + c*e**2)/sqrt(c))/sqrt(-4*a*f**2 + c*e**2))**(-p - q - 1)*(a + c*e*x/f + c*x**2)**p*(a*f/c + e*x + f*x**2)**(q + 1)*hyper((p + q + 1, -p - q), (p + q + 2,), sqrt(c)*(e + 2*f*x + sqrt(-4*a*f**2 + c*e**2)/sqrt(c))/(2*sqrt(-4*a*f**2 + c*e**2)))/(sqrt(-4*a*f**2 + c*e**2)*(p + q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_11():
    f = sqrt(x**2 + 2*x + 1)/sqrt(x**2 + 1)
    F = sqrt(x**2 + 1)*sqrt(x**2 + 2*x + 1)/(x + 1) + sqrt(x**2 + 2*x + 1)*asinh(x)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_12():
    f = 1/((x**2 - 1)**2*sqrt(x**2 + x - 1))
    F = -atan((x + 3)/(2*sqrt(x**2 + x - 1)))/8 - 5*atanh((1 - 3*x)/(2*sqrt(x**2 + x - 1)))/8 + sqrt(x**2 + x - 1)/(2 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_13():
    f = 1/(sqrt(d + f*x**2)*sqrt(a + b*x + c*x**2))
    F = -sqrt(((2*a + x*(b + sqrt(-4*a*c + b**2)))**2*(4*c**2*d + f*(b + sqrt(-4*a*c + b**2))**2)/((4*a**2*f + d*(b + sqrt(-4*a*c + b**2))**2)*(b + 2*c*x + sqrt(-4*a*c + b**2))**2) - (2*a + x*(b + sqrt(-4*a*c + b**2)))*(4*b + 4*sqrt(-4*a*c + b**2))*(a*f + c*d)/((4*a**2*f + d*(b + sqrt(-4*a*c + b**2))**2)*(b + 2*c*x + sqrt(-4*a*c + b**2))) + 1)/((2*a + x*(b + sqrt(-4*a*c + b**2)))*sqrt(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)/((b + 2*c*x + sqrt(-4*a*c + b**2))*sqrt(-2*a*(-a*f + c*d) + b**2*d + b*d*sqrt(-4*a*c + b**2))) + 1)**2)*sqrt((d + f*x**2)*(4*a*c - (b + sqrt(-4*a*c + b**2))**2)**2/((4*a**2*f + d*(b + sqrt(-4*a*c + b**2))**2)*(b + 2*c*x + sqrt(-4*a*c + b**2))**2))*sqrt(2*a + x*(b + sqrt(-4*a*c + b**2)))*((2*a + x*(b + sqrt(-4*a*c + b**2)))*sqrt(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)/((b + 2*c*x + sqrt(-4*a*c + b**2))*sqrt(-2*a*(-a*f + c*d) + b**2*d + b*d*sqrt(-4*a*c + b**2))) + 1)*(b + 2*c*x + sqrt(-4*a*c + b**2))**(sympy.S(3)/2)*(-2*a*(-a*f + c*d) + b**2*d + b*d*sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*elliptic_f(2*atan(sqrt(2*a + x*(b + sqrt(-4*a*c + b**2)))*(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)**(sympy.S(1)/4)/(sqrt(b + 2*c*x + sqrt(-4*a*c + b**2))*(-2*a*(-a*f + c*d) + b**2*d + b*d*sqrt(-4*a*c + b**2))**(sympy.S(1)/4))), (b + sqrt(-4*a*c + b**2))*(a*f + c*d)/(2*sqrt(-2*a*(-a*f + c*d) + b**2*d + b*d*sqrt(-4*a*c + b**2))*sqrt(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)) + sympy.S.Half)/(sqrt(d + f*x**2)*(4*a*c - (b + sqrt(-4*a*c + b**2))**2)*sqrt(a + b*x + c*x**2)*(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)**(sympy.S(1)/4)*sqrt((2*a + x*(b + sqrt(-4*a*c + b**2)))**2*(4*c**2*d + f*(b + sqrt(-4*a*c + b**2))**2)/((4*a**2*f + d*(b + sqrt(-4*a*c + b**2))**2)*(b + 2*c*x + sqrt(-4*a*c + b**2))**2) - (2*a + x*(b + sqrt(-4*a*c + b**2)))*(4*b + 4*sqrt(-4*a*c + b**2))*(a*f + c*d)/((4*a**2*f + d*(b + sqrt(-4*a*c + b**2))**2)*(b + 2*c*x + sqrt(-4*a*c + b**2))) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_14():
    f = sqrt(-x**2 - 4*x - 3)/(2*x**2 + 4*x + 3)
    F = -asin(x + 2)/2 - sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/2 + sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/2 - atanh(x/sqrt(-x**2 - 4*x - 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_15():
    f = (2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**4
    F = 1250*x**11/11 + 475*x**10/2 + 5075*x**9/9 + 3415*x**8/4 + 1176*x**7 + 2377*x**6/2 + 5099*x**5/5 + 656*x**4 + 1064*x**3/3 + 136*x**2 + 48*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_16():
    f = (2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**3
    F = 250*x**9/9 + 325*x**8/8 + 720*x**7/7 + 134*x**6 + 876*x**5/5 + 579*x**4/4 + 322*x**3/3 + 50*x**2 + 24*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_17():
    f = (2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**2
    F = 50*x**7/7 + 35*x**6/6 + 103*x**5/5 + 85*x**4/4 + 83*x**3/3 + 16*x**2 + 12*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_18():
    f = (2*x**2 - x + 3)*(5*x**2 + 3*x + 2)
    F = 2*x**5 + x**4/4 + 16*x**3/3 + 7*x**2/2 + 6*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_19():
    f = (2*x**2 - x + 3)/(5*x**2 + 3*x + 2)
    F = 2*x/5 - 11*log(5*x**2 + 3*x + 2)/50 + 143*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/775
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_20():
    f = (2*x**2 - x + 3)/(5*x**2 + 3*x + 2)**2
    F = (143*x + 77)/(775*x**2 + 465*x + 310) + 82*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/961
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_21():
    f = (2*x**2 - x + 3)/(5*x**2 + 3*x + 2)**3
    F = (143*x + 77)/(310*(5*x**2 + 3*x + 2)**2) + (5530*x + 1659)/(48050*x**2 + 28830*x + 19220) + 1106*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/29791
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_22():
    f = (2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**4
    F = 2500*x**13/13 + 875*x**12/3 + 11525*x**11/11 + 1571*x**10 + 24859*x**9/9 + 3315*x**8 + 27763*x**7/7 + 10771*x**6/3 + 14801*x**5/5 + 1838*x**4 + 3016*x**3/3 + 384*x**2 + 144*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_23():
    f = (2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**3
    F = 500*x**11/11 + 40*x**10 + 1865*x**9/9 + 1863*x**8/8 + 444*x**7 + 449*x**6 + 2693*x**5/5 + 1615*x**4/4 + 914*x**3/3 + 138*x**2 + 72*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_24():
    f = (2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**2
    F = 100*x**9/9 + 5*x**8/2 + 321*x**7/7 + 86*x**6/3 + 78*x**5 + 59*x**4 + 241*x**3/3 + 42*x**2 + 36*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_25():
    f = (2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)
    F = 20*x**7/7 - 4*x**6/3 + 61*x**5/5 + x**4/4 + 53*x**3/3 + 15*x**2/2 + 18*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_26():
    f = (2*x**2 - x + 3)**2/(5*x**2 + 3*x + 2)
    F = 4*x**3/15 - 16*x**2/25 + 381*x/125 - 1573*log(5*x**2 + 3*x + 2)/1250 + 8349*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/19375
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_27():
    f = (2*x**2 - x + 3)**2/(5*x**2 + 3*x + 2)**2
    F = 4*x/25 + (8349*x + 7381)/(19375*x**2 + 11625*x + 7750) - 22*log(5*x**2 + 3*x + 2)/125 + 41932*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/120125
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_28():
    f = (2*x**2 - x + 3)**2/(5*x**2 + 3*x + 2)**3
    F = (8349*x + 7381)/(7750*(5*x**2 + 3*x + 2)**2) + (502810*x + 193127)/(1201250*x**2 + 720750*x + 480500) + 4330*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/29791
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_29():
    f = (2*x**2 - x + 3)**2/(5*x**2 + 3*x + 2)**4
    F = (8349*x + 7381)/(11625*(5*x**2 + 3*x + 2)**3) + (132660*x + 50369)/(120125*(5*x**2 + 3*x + 2)**2) + (166880*x + 50064)/(744775*x**2 + 446865*x + 297910) + 66752*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/923521
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_30():
    f = (2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)**4
    F = 1000*x**15/3 + 2250*x**14/7 + 27050*x**13/13 + 30395*x**12/12 + 68583*x**11/11 + 75311*x**10/10 + 103583*x**9/9 + 94881*x**8/8 + 91349*x**7/7 + 64529*x**6/6 + 43083*x**5/5 + 5144*x**4 + 2856*x**3 + 1080*x**2 + 432*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_31():
    f = (2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)**3
    F = 1000*x**13/13 + 25*x**12 + 4830*x**11/11 + 3061*x**10/10 + 3316*x**9/3 + 7869*x**8/8 + 12016*x**7/7 + 2873*x**6/2 + 8292*x**5/5 + 4483*x**4/4 + 870*x**3 + 378*x**2 + 216*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_32():
    f = (2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)**2
    F = 200*x**11/11 - 6*x**10 + 922*x**9/9 + 83*x**8/8 + 1571*x**7/7 + 299*x**6/3 + 1416*x**5/5 + 635*x**4/4 + 237*x**3 + 108*x**2 + 108*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_33():
    f = (2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)
    F = 40*x**9/9 - 9*x**8/2 + 190*x**7/7 - 83*x**6/6 + 288*x**5/5 - 5*x**4 + 60*x**3 + 27*x**2/2 + 54*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_34():
    f = (2*x**2 - x + 3)**3/(5*x**2 + 3*x + 2)
    F = 8*x**5/25 - 21*x**4/25 + 1222*x**3/375 - 7451*x**2/1250 + 49508*x/3125 - 158389*log(5*x**2 + 3*x + 2)/31250 + 328757*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/484375
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_35():
    f = (2*x**2 - x + 3)**3/(5*x**2 + 3*x + 2)**2
    F = 8*x**3/75 - 54*x**2/125 + 1466*x/625 + (328757*x + 589633)/(484375*x**2 + 290625*x + 193750) - 10769*log(5*x**2 + 3*x + 2)/6250 + 3819607*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/3003125
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_36():
    f = (2*x**2 - x + 3)**3/(5*x**2 + 3*x + 2)**3
    F = 8*x/125 + (328757*x + 589633)/(193750*(5*x**2 + 3*x + 2)**2) + (41483640*x + 22794101)/(30031250*x**2 + 18018750*x + 12012500) - 66*log(5*x**2 + 3*x + 2)/625 + 11341176*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/18619375
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_37():
    f = (5*x**2 + 3*x + 2)**4/(2*x**2 - x + 3)
    F = 625*x**7/14 + 3625*x**6/24 + 1855*x**5/8 + 6245*x**4/64 - 21229*x**3/96 - 28747*x**2/128 + 122691*x/128 + 307461*log(2*x**2 - x + 3)/512 + 1156639*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/5888
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_38():
    f = (5*x**2 + 3*x + 2)**3/(2*x**2 - x + 3)
    F = 25*x**5/2 + 575*x**4/16 + 965*x**3/24 - 829*x**2/32 - 4795*x/32 + 1331*log(2*x**2 - x + 3)/128 - 59895*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/1472
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_39():
    f = (5*x**2 + 3*x + 2)**2/(2*x**2 - x + 3)
    F = 25*x**3/6 + 85*x**2/8 + 51*x/8 - 363*log(2*x**2 - x + 3)/32 + 847*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/368
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_40():
    f = (5*x**2 + 3*x + 2)/(2*x**2 - x + 3)
    F = 5*x/2 + 11*log(2*x**2 - x + 3)/8 + 33*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/92
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_41():
    f = 1/((2*x**2 - x + 3)*(5*x**2 + 3*x + 2))
    F = -log(2*x**2 - x + 3)/44 + log(5*x**2 + 3*x + 2)/44 + 3*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/506 + 13*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/682
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_42():
    f = 1/((2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**2)
    F = (65*x + 4)/(3410*x**2 + 2046*x + 1364) + 3*log(2*x**2 - x + 3)/968 - 3*log(5*x**2 + 3*x + 2)/968 + 7*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/11132 + 2891*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/465124
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_43():
    f = 1/((2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**3)
    F = (65*x + 4)/(1364*(5*x**2 + 3*x + 2)**2) + (21605*x + 7923)/(2325620*x**2 + 1395372*x + 930248) - log(2*x**2 - x + 3)/21296 + log(5*x**2 + 3*x + 2)/21296 - 45*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/244904 + 847793*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/317214568
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_44():
    f = (5*x**2 + 3*x + 2)**4/(2*x**2 - x + 3)**2
    F = 125*x**5/4 + 2125*x**4/16 + 9775*x**3/48 - 1185*x**2/8 - 89359*x/64 - (1156639*x + 1478741)/(5888*x**2 - 2944*x + 8832) - 30613*log(2*x**2 - x + 3)/128 - 13292697*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/33856
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_45():
    f = (5*x**2 + 3*x + 2)**3/(2*x**2 - x + 3)**2
    F = 125*x**3/12 + 175*x**2/4 + 915*x/16 - (22627 - 59895*x)/(1472*x**2 - 736*x + 2208) - 2057*log(2*x**2 - x + 3)/32 + 223971*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/8464
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_46():
    f = (5*x**2 + 3*x + 2)**2/(2*x**2 - x + 3)**2
    F = 25*x/4 + (2299 - 847*x)/(368*x**2 - 184*x + 552) + 55*log(2*x**2 - x + 3)/8 + 1859*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/2116
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_47():
    f = (5*x**2 + 3*x + 2)/(2*x**2 - x + 3)**2
    F = (-33*x - 55)/(92*x**2 - 46*x + 138) - 82*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/529
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_48():
    f = 1/((2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2))
    F = (13 - 6*x)/(1012*x**2 - 506*x + 1518) - 13*log(2*x**2 - x + 3)/968 + 13*log(5*x**2 + 3*x + 2)/968 + 241*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/256036 + 69*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/15004
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_49():
    f = 1/((2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**2)
    F = (13 - 6*x)/((5*x**2 + 3*x + 2)*(1012*x**2 - 506*x + 1518)) + (3425*x - 2925)/(862730*x**2 + 517638*x + 345092) + 19*log(2*x**2 - x + 3)/10648 - 19*log(5*x**2 + 3*x + 2)/10648 + 2769*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/2816396 + 12643*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/5116364
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_50():
    f = 1/((2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**3)
    F = (13 - 6*x)/((5*x**2 + 3*x + 2)**2*(1012*x**2 - 506*x + 1518)) + (5765*x - 9446)/(690184*(5*x**2 + 3*x + 2)**2) + (3996965*x + 1765599)/(1176763720*x**2 + 706058232*x + 470705488) + 97*log(2*x**2 - x + 3)/468512 - 97*log(5*x**2 + 3*x + 2)/468512 - 25557*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/123921424 + 4464079*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/6978720496
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_51():
    f = (5*x**2 + 3*x + 2)**4/(2*x**2 - x + 3)**3
    F = 625*x**3/24 + 4875*x**2/32 + 2725*x/8 - (1156639*x + 1478741)/(5888*(2*x**2 - x + 3)**2) + (101715020*x + 6959799)/(270848*x**2 - 135424*x + 406272) - 13915*log(2*x**2 - x + 3)/64 + 63799791*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/389344
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_52():
    f = (5*x**2 + 3*x + 2)**3/(2*x**2 - x + 3)**3
    F = 125*x/8 - (22627 - 59895*x)/(1472*(2*x**2 - x + 3)**2) + (2564353 - 1552188*x)/(67712*x**2 - 33856*x + 101568) + 825*log(2*x**2 - x + 3)/32 + 165099*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/194672
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_53():
    f = (5*x**2 + 3*x + 2)**2/(2*x**2 - x + 3)**3
    F = (2299 - 847*x)/(368*(2*x**2 - x + 3)**2) - (18260*x + 53625)/(16928*x**2 - 8464*x + 25392) - 4330*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/12167
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_54():
    f = (5*x**2 + 3*x + 2)/(2*x**2 - x + 3)**3
    F = -(131 - 524*x)/(4232*x**2 - 2116*x + 6348) + (-33*x - 55)/(92*(2*x**2 - x + 3)**2) - 262*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/12167
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_55():
    f = 1/((2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2))
    F = (13 - 6*x)/(1012*(2*x**2 - x + 3)**2) + (3625 - 746*x)/(512072*x**2 - 256036*x + 768108) - 119*log(2*x**2 - x + 3)/21296 + 119*log(5*x**2 + 3*x + 2)/21296 - 53403*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/129554216 + 247*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/330088
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_56():
    f = 1/((2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)**2)
    F = (13 - 6*x)/(1012*(2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)) + (9665 - 1446*x)/((5*x**2 + 3*x + 2)*(1024144*x**2 - 512072*x + 1536216)) + (-252815*x - 2328909)/(873082760*x**2 + 523849656*x + 349233104) + 181*log(2*x**2 - x + 3)/468512 - 181*log(5*x**2 + 3*x + 2)/468512 + 2038497*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/2850192752 + 246757*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/225120016
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_57():
    f = 1/((2*x**2 - x + 3)**3*(5*x**2 + 3*x + 2)**3)
    F = (13 - 6*x)/(1012*(2*x**2 - x + 3)**2*(5*x**2 + 3*x + 2)**2) + (1510 - 175*x)/((5*x**2 + 3*x + 2)**2*(128018*x**2 - 64009*x + 192027)) + (-385100*x - 1118535)/(87308276*(5*x**2 + 3*x + 2)**2) + (107106525*x + 39274590)/(74430305290*x**2 + 44658183174*x + 29772122116) + 405*log(2*x**2 - x + 3)/1288408 - 405*log(5*x**2 + 3*x + 2)/1288408 - 880575*sqrt(23)*atan(sqrt(23)*(1 - 4*x)/23)/7838030068 + 2768835*sqrt(31)*atan(sqrt(31)*(10*x + 3)/31)/19191481364
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_58():
    f = sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**4
    F = 125*x**7*(2*x**2 - x + 3)**(sympy.S(3)/2)/4 + 14125*x**6*(2*x**2 - x + 3)**(sympy.S(3)/2)/144 + 233225*x**5*(2*x**2 - x + 3)**(sympy.S(3)/2)/1536 + 4796405*x**4*(2*x**2 - x + 3)**(sympy.S(3)/2)/43008 + 8325631*x**3*(2*x**2 - x + 3)**(sympy.S(3)/2)/1032192 - 83948353*x**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/2293760 + 804243809*x*(2*x**2 - x + 3)**(sympy.S(3)/2)/36700160 - (359471503 - 1437886012*x)*sqrt(2*x**2 - x + 3)/67108864 + 27185733541*(2*x**2 - x + 3)**(sympy.S(3)/2)/440401920 - 8267844569*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/268435456
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_59():
    f = sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**3
    F = 125*x**5*(2*x**2 - x + 3)**(sympy.S(3)/2)/16 + 8825*x**4*(2*x**2 - x + 3)**(sympy.S(3)/2)/448 + 247435*x**3*(2*x**2 - x + 3)**(sympy.S(3)/2)/10752 + 531681*x**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/71680 - 9627393*x*(2*x**2 - x + 3)**(sympy.S(3)/2)/1146880 - (6766097 - 27064388*x)*sqrt(2*x**2 - x + 3)/2097152 - 22548119*(2*x**2 - x + 3)**(sympy.S(3)/2)/4587520 - 155620231*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8388608
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_60():
    f = sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**2
    F = 25*x**3*(2*x**2 - x + 3)**(sympy.S(3)/2)/12 + 63*x**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/16 + 769*x*(2*x**2 - x + 3)**(sympy.S(3)/2)/256 + (12371 - 49484*x)*sqrt(2*x**2 - x + 3)/16384 - 2107*(2*x**2 - x + 3)**(sympy.S(3)/2)/3072 + 284533*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/65536
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_61():
    f = sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)
    F = 5*x*(2*x**2 - x + 3)**(sympy.S(3)/2)/8 + (81*x/128 + sympy.S(-81)/512)*sqrt(2*x**2 - x + 3) + 73*(2*x**2 - x + 3)**(sympy.S(3)/2)/96 - 1863*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/2048
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_62():
    f = sqrt(2*x**2 - x + 3)/(5*x**2 + 3*x + 2)
    F = -sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/5 + sqrt(sympy.S(143)/31 + 110*sqrt(2)/31)*atan(sqrt(11)*(x*(13*sqrt(2) + 20) + 6 + 7*sqrt(2))/(sqrt(806 + 620*sqrt(2))*sqrt(2*x**2 - x + 3)))/5 - sqrt(sympy.S(-143)/31 + 110*sqrt(2)/31)*atanh(sqrt(11)*(x*(20 - 13*sqrt(2)) - 7*sqrt(2) + 6)/(sqrt(-806 + 620*sqrt(2))*sqrt(2*x**2 - x + 3)))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_63():
    f = sqrt(2*x**2 - x + 3)/(5*x**2 + 3*x + 2)**2
    F = (10*x + 3)*sqrt(2*x**2 - x + 3)/(155*x**2 + 93*x + 62) + sqrt(sympy.S(70517)/682 + 24971*sqrt(2)/341)*atan(sqrt(11)*(x*(973 + 696*sqrt(2)) + 277*sqrt(2) + 419)/(sqrt(2186027 + 1548202*sqrt(2))*sqrt(2*x**2 - x + 3)))/62 - sqrt(sympy.S(-70517)/682 + 24971*sqrt(2)/341)*atanh(sqrt(11)*(x*(973 - 696*sqrt(2)) - 277*sqrt(2) + 419)/(sqrt(-2186027 + 1548202*sqrt(2))*sqrt(2*x**2 - x + 3)))/62
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_64():
    f = sqrt(2*x**2 - x + 3)/(5*x**2 + 3*x + 2)**3
    F = (10*x + 3)*sqrt(2*x**2 - x + 3)/(62*(5*x**2 + 3*x + 2)**2) + (13665*x + 3464)*sqrt(2*x**2 - x + 3)/(422840*x**2 + 253704*x + 169136) + sqrt(sympy.S(112285869463)/682 + 39699690370*sqrt(2)/341)*atan(sqrt(11)*(x*(872375*sqrt(2) + 1235163) + 509587 + 362788*sqrt(2))/(sqrt(3480861953353 + 2461380802940*sqrt(2))*sqrt(2*x**2 - x + 3)))/169136 - sqrt(sympy.S(-112285869463)/682 + 39699690370*sqrt(2)/341)*atanh(sqrt(11)*(x*(1235163 - 872375*sqrt(2)) - 362788*sqrt(2) + 509587)/(sqrt(-3480861953353 + 2461380802940*sqrt(2))*sqrt(2*x**2 - x + 3)))/169136
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_65():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**4
    F = 625*x**7*(2*x**2 - x + 3)**(sympy.S(5)/2)/24 + 7625*x**6*(2*x**2 - x + 3)**(sympy.S(5)/2)/96 + 95165*x**5*(2*x**2 - x + 3)**(sympy.S(5)/2)/768 + 941905*x**4*(2*x**2 - x + 3)**(sympy.S(5)/2)/9216 + 10444117*x**3*(2*x**2 - x + 3)**(sympy.S(5)/2)/294912 - 56422489*x**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/8257536 + 48669967*x*(2*x**2 - x + 3)**(sympy.S(5)/2)/22020096 - (382121949 - 1528487796*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/134217728 - (26366414481 - 105465657924*x)*sqrt(2*x**2 - x + 3)/2147483648 + 2124689283*(2*x**2 - x + 3)**(sympy.S(5)/2)/146800640 - 606427533063*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8589934592
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_66():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**3
    F = 25*x**5*(2*x**2 - x + 3)**(sympy.S(5)/2)/4 + 725*x**4*(2*x**2 - x + 3)**(sympy.S(5)/2)/48 + 27785*x**3*(2*x**2 - x + 3)**(sympy.S(5)/2)/1536 + 384739*x**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/43008 - 81685*x*(2*x**2 - x + 3)**(sympy.S(5)/2)/114688 - (667795 - 2671180*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/2097152 - (46077855 - 184311420*x)*sqrt(2*x**2 - x + 3)/33554432 - 4625907*(2*x**2 - x + 3)**(sympy.S(5)/2)/2293760 - 1059790665*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/134217728
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_67():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**2
    F = 25*x**3*(2*x**2 - x + 3)**(sympy.S(5)/2)/16 + 1235*x**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/448 + 24499*x*(2*x**2 - x + 3)**(sympy.S(5)/2)/10752 + (24293 - 97172*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/196608 + (558739 - 2234956*x)*sqrt(2*x**2 - x + 3)/1048576 + 73861*(2*x**2 - x + 3)**(sympy.S(5)/2)/215040 + 12850997*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4194304
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_68():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)
    F = 5*x*(2*x**2 - x + 3)**(sympy.S(5)/2)/12 - (179 - 716*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/1536 - (4117 - 16468*x)*sqrt(2*x**2 - x + 3)/8192 + 107*(2*x**2 - x + 3)**(sympy.S(5)/2)/240 - 94691*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/32768
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_69():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)/(5*x**2 + 3*x + 2)
    F = -(49 - 20*x)*sqrt(2*x**2 - x + 3)/100 - 2203*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/2000 + 11*sqrt(sympy.S(2717)/31 + 5500*sqrt(2)/31)*atan(sqrt(11)*(x*(69*sqrt(2) + 130) + 8 + 61*sqrt(2))/(sqrt(15314 + 31000*sqrt(2))*sqrt(2*x**2 - x + 3)))/125 - 11*sqrt(sympy.S(-2717)/31 + 5500*sqrt(2)/31)*atanh(sqrt(11)*(x*(130 - 69*sqrt(2)) - 61*sqrt(2) + 8)/(sqrt(-15314 + 31000*sqrt(2))*sqrt(2*x**2 - x + 3)))/125
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_70():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)/(5*x**2 + 3*x + 2)**2
    F = (16 - 20*x)*sqrt(2*x**2 - x + 3)/155 + (10*x + 3)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(155*x**2 + 93*x + 62) - 2*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/25 + sqrt(sympy.S(34862663)/31 + 24918850*sqrt(2)/31)*atan(sqrt(11)*(x*(6477*sqrt(2) + 9440) + 3514 + 2963*sqrt(2))/(sqrt(196498646 + 140451700*sqrt(2))*sqrt(2*x**2 - x + 3)))/1550 - sqrt(sympy.S(-34862663)/31 + 24918850*sqrt(2)/31)*atanh(sqrt(11)*(x*(9440 - 6477*sqrt(2)) - 2963*sqrt(2) + 3514)/(sqrt(-196498646 + 140451700*sqrt(2))*sqrt(2*x**2 - x + 3)))/1550
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_71():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)/(5*x**2 + 3*x + 2)**3
    F = (10*x + 3)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(62*(5*x**2 + 3*x + 2)**2) + (2088*x + 831)*sqrt(2*x**2 - x + 3)/(19220*x**2 + 11532*x + 7688) + 3*sqrt(sympy.S(366990269)/682 + 129754513*sqrt(2)/341)*atan(sqrt(11)*(x*(70517 + 49942*sqrt(2)) + 20575*sqrt(2) + 29367)/(sqrt(11376698339 + 8044779806*sqrt(2))*sqrt(2*x**2 - x + 3)))/7688 - 3*sqrt(sympy.S(-366990269)/682 + 129754513*sqrt(2)/341)*atanh(sqrt(11)*(x*(70517 - 49942*sqrt(2)) - 20575*sqrt(2) + 29367)/(sqrt(-11376698339 + 8044779806*sqrt(2))*sqrt(2*x**2 - x + 3)))/7688
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_72():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)**4
    F = 625*x**7*(2*x**2 - x + 3)**(sympy.S(7)/2)/28 + 13875*x**6*(2*x**2 - x + 3)**(sympy.S(7)/2)/208 + 1046225*x**5*(2*x**2 - x + 3)**(sympy.S(7)/2)/9984 + 3684995*x**4*(2*x**2 - x + 3)**(sympy.S(7)/2)/39936 + 23460839*x**3*(2*x**2 - x + 3)**(sympy.S(7)/2)/532480 + 122595067*x**2*(2*x**2 - x + 3)**(sympy.S(7)/2)/19169280 + 112244125*x*(2*x**2 - x + 3)**(sympy.S(7)/2)/122683392 - (401135647 - 1604542588*x)*(2*x**2 - x + 3)**(sympy.S(5)/2)/335544320 - (9226119881 - 36904479524*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/2147483648 - (636602271789 - 2546409087156*x)*sqrt(2*x**2 - x + 3)/34359738368 + 25250178739*(2*x**2 - x + 3)**(sympy.S(7)/2)/5725224960 - 14641852251147*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/137438953472
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_73():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)**3
    F = 125*x**5*(2*x**2 - x + 3)**(sympy.S(7)/2)/24 + 1175*x**4*(2*x**2 - x + 3)**(sympy.S(7)/2)/96 + 3823*x**3*(2*x**2 - x + 3)**(sympy.S(7)/2)/256 + 80483*x**2*(2*x**2 - x + 3)**(sympy.S(7)/2)/9216 + 509257*x*(2*x**2 - x + 3)**(sympy.S(7)/2)/294912 - (57915 - 231660*x)*(2*x**2 - x + 3)**(sympy.S(5)/2)/2097152 - (6660225 - 26640900*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/67108864 - (459555525 - 1838222100*x)*sqrt(2*x**2 - x + 3)/1073741824 - 1696165*(2*x**2 - x + 3)**(sympy.S(7)/2)/2752512 - 10569777075*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4294967296
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_74():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)**2
    F = 5*x**3*(2*x**2 - x + 3)**(sympy.S(7)/2)/4 + 305*x**2*(2*x**2 - x + 3)**(sympy.S(7)/2)/144 + 8467*x*(2*x**2 - x + 3)**(sympy.S(7)/2)/4608 - (1547 - 6188*x)*(2*x**2 - x + 3)**(sympy.S(5)/2)/98304 - (177905 - 711620*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/3145728 - (4091815 - 16367260*x)*sqrt(2*x**2 - x + 3)/16777216 + 23225*(2*x**2 - x + 3)**(sympy.S(7)/2)/43008 - 94111745*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/67108864
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_75():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)
    F = 5*x*(2*x**2 - x + 3)**(sympy.S(7)/2)/16 - (277 - 1108*x)*(2*x**2 - x + 3)**(sympy.S(5)/2)/3072 - (31855 - 127420*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/98304 - (732665 - 2930660*x)*sqrt(2*x**2 - x + 3)/524288 + 141*(2*x**2 - x + 3)**(sympy.S(7)/2)/448 - 16851295*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/2097152
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_76():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)/(5*x**2 + 3*x + 2)
    F = -(103 - 60*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/600 - (226249 - 99620*x)*sqrt(2*x**2 - x + 3)/80000 - 7216203*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/1600000 - 121*sqrt(sympy.S(-170027)/31 + 275000*sqrt(2)/31)*atan(sqrt(11)*(-x*(247*sqrt(2) + 690) - 443*sqrt(2) + 196)/(sqrt(-958334 + 1550000*sqrt(2))*sqrt(2*x**2 - x + 3)))/3125 + 121*sqrt(sympy.S(170027)/31 + 275000*sqrt(2)/31)*atanh(sqrt(11)*(-x*(690 - 247*sqrt(2)) + 196 + 443*sqrt(2))/(sqrt(958334 + 1550000*sqrt(2))*sqrt(2*x**2 - x + 3)))/3125
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_77():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)/(5*x**2 + 3*x + 2)**2
    F = (16 - 20*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/155 + (10*x + 3)*(2*x**2 - x + 3)**(sympy.S(5)/2)/(155*x**2 + 93*x + 62) - (2240*x + 1277)*sqrt(2*x**2 - x + 3)/7750 - 4799*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/5000 + 11*sqrt(sympy.S(2469614213)/31 + 2139362500*sqrt(2)/31)*atan(sqrt(11)*(x*(54423*sqrt(2) + 87710) + 21136 + 33287*sqrt(2))/(sqrt(13919643746 + 12058225000*sqrt(2))*sqrt(2*x**2 - x + 3)))/38750 - 11*sqrt(sympy.S(-2469614213)/31 + 2139362500*sqrt(2)/31)*atanh(sqrt(11)*(x*(87710 - 54423*sqrt(2)) - 33287*sqrt(2) + 21136)/(sqrt(-13919643746 + 12058225000*sqrt(2))*sqrt(2*x**2 - x + 3)))/38750
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_78():
    f = (2*x**2 - x + 3)**(sympy.S(5)/2)/(5*x**2 + 3*x + 2)**3
    F = (11359 - 12920*x)*sqrt(2*x**2 - x + 3)/48050 + (10*x + 3)*(2*x**2 - x + 3)**(sympy.S(5)/2)/(62*(5*x**2 + 3*x + 2)**2) + (2336*x + 769)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(19220*x**2 + 11532*x + 7688) - 4*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/125 + sqrt(11 + 44*sqrt(2))*(1978861*sqrt(2) + 2937349)*atan(sqrt(11)*(x*(6895071*sqrt(2) + 9832420) + 3957722 + 2937349*sqrt(2))/(sqrt(218922973868534 + 154928828417500*sqrt(2))*sqrt(2*x**2 - x + 3)))/29791000 - sqrt(-11 + 44*sqrt(2))*(2937349 - 1978861*sqrt(2))*atanh(sqrt(11)*(x*(9832420 - 6895071*sqrt(2)) - 2937349*sqrt(2) + 3957722)/(sqrt(-218922973868534 + 154928828417500*sqrt(2))*sqrt(2*x**2 - x + 3)))/29791000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_79():
    f = (5*x**2 + 3*x + 2)**4/sqrt(2*x**2 - x + 3)
    F = 625*x**7*sqrt(2*x**2 - x + 3)/16 + 57375*x**6*sqrt(2*x**2 - x + 3)/448 + 2116475*x**5*sqrt(2*x**2 - x + 3)/10752 + 686531*x**4*sqrt(2*x**2 - x + 3)/6144 - 19750457*x**3*sqrt(2*x**2 - x + 3)/229376 - 15428243*x**2*sqrt(2*x**2 - x + 3)/131072 + 1572007407*x*sqrt(2*x**2 - x + 3)/7340032 + 16493087661*sqrt(2*x**2 - x + 3)/29360128 + 2899366573*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16777216
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_80():
    f = (5*x**2 + 3*x + 2)**3/sqrt(2*x**2 - x + 3)
    F = 125*x**5*sqrt(2*x**2 - x + 3)/12 + 1355*x**4*sqrt(2*x**2 - x + 3)/48 + 8185*x**3*sqrt(2*x**2 - x + 3)/256 - 3387*x**2*sqrt(2*x**2 - x + 3)/1024 - 372783*x*sqrt(2*x**2 - x + 3)/8192 - 203373*sqrt(2*x**2 - x + 3)/32768 - 9267707*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/131072
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_81():
    f = (5*x**2 + 3*x + 2)**2/sqrt(2*x**2 - x + 3)
    F = 25*x**3*sqrt(2*x**2 - x + 3)/8 + 655*x**2*sqrt(2*x**2 - x + 3)/96 + 3443*x*sqrt(2*x**2 - x + 3)/768 - 11373*sqrt(2*x**2 - x + 3)/1024 + 30725*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4096
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_82():
    f = (5*x**2 + 3*x + 2)/sqrt(2*x**2 - x + 3)
    F = 5*x*sqrt(2*x**2 - x + 3)/4 + 39*sqrt(2*x**2 - x + 3)/16 + 17*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_83():
    f = 1/(sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2))
    F = sqrt(sympy.S(13)/682 + 5*sqrt(2)/341)*atan(sqrt(11)*(x*(13 + 10*sqrt(2)) + 3*sqrt(2) + 7)/(sqrt(403 + 310*sqrt(2))*sqrt(2*x**2 - x + 3))) - sqrt(sympy.S(-13)/682 + 5*sqrt(2)/341)*atanh(sqrt(11)*(x*(13 - 10*sqrt(2)) - 3*sqrt(2) + 7)/(sqrt(-403 + 310*sqrt(2))*sqrt(2*x**2 - x + 3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_84():
    f = 1/(sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**2)
    F = (65*x + 4)*sqrt(2*x**2 - x + 3)/(3410*x**2 + 2046*x + 1364) + sqrt(sympy.S(2343727)/682 + 839350*sqrt(2)/341)*atan(sqrt(11)*(x*(3935*sqrt(2) + 5751) + 2119 + 1816*sqrt(2))/(sqrt(72655537 + 52039700*sqrt(2))*sqrt(2*x**2 - x + 3)))/1364 - sqrt(sympy.S(-2343727)/682 + 839350*sqrt(2)/341)*atanh(sqrt(11)*(x*(5751 - 3935*sqrt(2)) - 1816*sqrt(2) + 2119)/(sqrt(-72655537 + 52039700*sqrt(2))*sqrt(2*x**2 - x + 3)))/1364
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_85():
    f = 1/(sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**3)
    F = (65*x + 4)*sqrt(2*x**2 - x + 3)/(1364*(5*x**2 + 3*x + 2)**2) + (86265*x + 26794)*sqrt(2*x**2 - x + 3)/(9302480*x**2 + 5581488*x + 3720992) + 25*sqrt(sympy.S(6414867847)/682 + 2268187300*sqrt(2)/341)*atan(sqrt(11)*(x*(294669 + 208915*sqrt(2)) + 85754*sqrt(2) + 123161)/(sqrt(198860903257 + 140627612600*sqrt(2))*sqrt(2*x**2 - x + 3)))/3720992 - 25*sqrt(sympy.S(-6414867847)/682 + 2268187300*sqrt(2)/341)*atanh(sqrt(11)*(x*(294669 - 208915*sqrt(2)) - 85754*sqrt(2) + 123161)/(sqrt(-198860903257 + 140627612600*sqrt(2))*sqrt(2*x**2 - x + 3)))/3720992
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_86():
    f = (5*x**2 + 3*x + 2)**4/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 625*x**5*sqrt(2*x**2 - x + 3)/24 + 10075*x**4*sqrt(2*x**2 - x + 3)/96 + 79425*x**3*sqrt(2*x**2 - x + 3)/512 - 111315*x**2*sqrt(2*x**2 - x + 3)/2048 - 8992487*x*sqrt(2*x**2 - x + 3)/16384 - (1156639*x + 1478741)/(1472*sqrt(2*x**2 - x + 3)) - 31009685*sqrt(2*x**2 - x + 3)/65536 - 310445587*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/262144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_87():
    f = (5*x**2 + 3*x + 2)**3/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 125*x**3*sqrt(2*x**2 - x + 3)/16 + 1825*x**2*sqrt(2*x**2 - x + 3)/64 + 15565*x*sqrt(2*x**2 - x + 3)/512 - (22627 - 59895*x)/(368*sqrt(2*x**2 - x + 3)) - 181561*sqrt(2*x**2 - x + 3)/2048 + 1168881*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8192
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_88():
    f = (5*x**2 + 3*x + 2)**2/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 25*x*sqrt(2*x**2 - x + 3)/8 + (2299 - 847*x)/(92*sqrt(2*x**2 - x + 3)) + 415*sqrt(2*x**2 - x + 3)/32 - 223*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_89():
    f = (5*x**2 + 3*x + 2)/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = (-33*x - 55)/(23*sqrt(2*x**2 - x + 3)) - 5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_90():
    f = 1/((2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2))
    F = (13 - 6*x)/(253*sqrt(2*x**2 - x + 3)) + sqrt(sympy.S(247)/682 + 250*sqrt(2)/341)*atan(sqrt(11)*(x*(69 + 65*sqrt(2)) + 4*sqrt(2) + 61)/(sqrt(7657 + 15500*sqrt(2))*sqrt(2*x**2 - x + 3)))/22 - sqrt(sympy.S(-247)/682 + 250*sqrt(2)/341)*atanh(sqrt(11)*(x*(69 - 65*sqrt(2)) - 4*sqrt(2) + 61)/(sqrt(-7657 + 15500*sqrt(2))*sqrt(2*x**2 - x + 3)))/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_91():
    f = 1/((2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**2)
    F = -(6315 - 2306*x)/(345092*sqrt(2*x**2 - x + 3)) + (65*x + 4)/(682*sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)) + sqrt(sympy.S(129694447)/682 + 51887500*sqrt(2)/341)*atan(sqrt(11)*(x*(29065*sqrt(2) + 45519) + 12611 + 16454*sqrt(2))/(sqrt(4020527857 + 3217025000*sqrt(2))*sqrt(2*x**2 - x + 3)))/30008 - sqrt(sympy.S(-129694447)/682 + 51887500*sqrt(2)/341)*atanh(sqrt(11)*(x*(45519 - 29065*sqrt(2)) - 16454*sqrt(2) + 12611)/(sqrt(-4020527857 + 3217025000*sqrt(2))*sqrt(2*x**2 - x + 3)))/30008
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_92():
    f = 1/((2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**3)
    F = -(4353943 - 6508666*x)/(941410976*sqrt(2*x**2 - x + 3)) + (65*x + 4)/(1364*sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)**2) + (86575*x + 36590)/(1860496*sqrt(2*x**2 - x + 3)*(5*x**2 + 3*x + 2)) + 3*sqrt(sympy.S(13874275807943)/682 + 4909869325000*sqrt(2)/341)*atan(sqrt(11)*(x*(9662095*sqrt(2) + 13785797) + 5538393 + 4123702*sqrt(2))/(sqrt(430102550046233 + 304411898150000*sqrt(2))*sqrt(2*x**2 - x + 3)))/81861824 - 3*sqrt(sympy.S(-13874275807943)/682 + 4909869325000*sqrt(2)/341)*atanh(sqrt(11)*(x*(13785797 - 9662095*sqrt(2)) - 4123702*sqrt(2) + 5538393)/(sqrt(-430102550046233 + 304411898150000*sqrt(2))*sqrt(2*x**2 - x + 3)))/81861824
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_93():
    f = (5*x**2 + 3*x + 2)**4/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = 625*x**3*sqrt(2*x**2 - x + 3)/32 + 38375*x**2*sqrt(2*x**2 - x + 3)/384 + 526075*x*sqrt(2*x**2 - x + 3)/3072 - (1156639*x + 1478741)/(4416*(2*x**2 - x + 3)**(sympy.S(3)/2)) + (154885808*x + 9861379)/(101568*sqrt(2*x**2 - x + 3)) - 1308645*sqrt(2*x**2 - x + 3)/4096 + 16955197*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16384
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_94():
    f = (5*x**2 + 3*x + 2)**3/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = 125*x*sqrt(2*x**2 - x + 3)/16 - (22627 - 59895*x)/(1104*(2*x**2 - x + 3)**(sympy.S(3)/2)) + (1292159 - 816024*x)/(8464*sqrt(2*x**2 - x + 3)) + 3175*sqrt(2*x**2 - x + 3)/64 - 7495*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_95():
    f = (5*x**2 + 3*x + 2)**2/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = (2299 - 847*x)/(276*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (25696*x + 80861)/(6348*sqrt(2*x**2 - x + 3)) - 25*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_96():
    f = (5*x**2 + 3*x + 2)/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = -(71 - 284*x)/(529*sqrt(2*x**2 - x + 3)) + (-33*x - 55)/(69*(2*x**2 - x + 3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_97():
    f = 1/((2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2))
    F = (13 - 6*x)/(759*(2*x**2 - x + 3)**(sympy.S(3)/2)) + (3603 - 658*x)/(128018*sqrt(2*x**2 - x + 3)) + sqrt(sympy.S(-15457)/682 + 12500*sqrt(2)/341)*atan(sqrt(11)*(x*(247 + 345*sqrt(2)) - 98*sqrt(2) + 443)/(sqrt(-479167 + 775000*sqrt(2))*sqrt(2*x**2 - x + 3)))/484 - sqrt(sympy.S(15457)/682 + 12500*sqrt(2)/341)*atanh(sqrt(11)*(x*(247 - 345*sqrt(2)) + 98*sqrt(2) + 443)/(sqrt(479167 + 775000*sqrt(2))*sqrt(2*x**2 - x + 3)))/484
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_98():
    f = 1/((2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)**2)
    F = -(15101 - 8654*x)/(1035276*(2*x**2 - x + 3)**(sympy.S(3)/2)) + (65*x + 4)/(682*(2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)) - (1352542*x + 3133427)/(523849656*sqrt(2*x**2 - x + 3)) + 625*sqrt(sympy.S(30463)/682 + 11800*sqrt(2)/341)*atan(sqrt(11)*(x*(445*sqrt(2) + 687) + 203 + 242*sqrt(2))/(sqrt(944353 + 731600*sqrt(2))*sqrt(2*x**2 - x + 3)))/660176 - 625*sqrt(sympy.S(-30463)/682 + 11800*sqrt(2)/341)*atanh(sqrt(11)*(x*(687 - 445*sqrt(2)) - 242*sqrt(2) + 203)/(sqrt(-944353 + 731600*sqrt(2))*sqrt(2*x**2 - x + 3)))/660176
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_99():
    f = 1/((2*x**2 - x + 3)**(sympy.S(5)/2)*(5*x**2 + 3*x + 2)**3)
    F = -(12280939 - 19536786*x)/(2824232928*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (1134826571 - 1504660754*x)/(476353953856*sqrt(2*x**2 - x + 3)) + (65*x + 4)/(1364*(2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)**2) + (86885*x + 46386)/(1860496*(2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**2 + 3*x + 2)) + 35*sqrt(sympy.S(2243059557247)/682 + 1005874250000*sqrt(2)/341)*atan(sqrt(11)*(x*(3861685*sqrt(2) + 6290431) + 1432939 + 2428746*sqrt(2))/(sqrt(69534846274657 + 62364203500000*sqrt(2))*sqrt(2*x**2 - x + 3)))/1800960128 - 35*sqrt(sympy.S(-2243059557247)/682 + 1005874250000*sqrt(2)/341)*atanh(sqrt(11)*(x*(6290431 - 3861685*sqrt(2)) - 2428746*sqrt(2) + 1432939)/(sqrt(-69534846274657 + 62364203500000*sqrt(2))*sqrt(2*x**2 - x + 3)))/1800960128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_100():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)**2
    F = f**2*x**3*(a + b*x + c*x**2)**(sympy.S(3)/2)/(6*c) + f*x**2*(-3*b*f + 8*c*e)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(20*c**2) + x*(a + b*x + c*x**2)**(sympy.S(3)/2)*(21*b**2*f**2 + 40*c**2*(2*d*f + e**2) - 4*c*f*(5*a*f + 14*b*e))/(160*c**3) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(-105*b**3*f**2 + 28*b*c*f*(7*a*f + 10*b*e) + 640*c**3*d*e - 8*c**2*(32*a*e*f + 25*b*(2*d*f + e**2)))/(960*c**4) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(21*b**4*f**2 - 56*b**2*c*f*(a*f + b*e) + 128*c**4*d**2 - 32*c**3*(a*(2*d*f + e**2) + 4*b*d*e) + 8*c**2*(2*a**2*f**2 + 12*a*b*e*f + 5*b**2*(2*d*f + e**2)))/(512*c**5) - (-4*a*c + b**2)*(21*b**4*f**2 - 56*b**2*c*f*(a*f + b*e) + 128*c**4*d**2 - 32*c**3*(a*(2*d*f + e**2) + 4*b*d*e) + 8*c**2*(2*a**2*f**2 + 12*a*b*e*f + 5*b**2*(2*d*f + e**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_101():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*x*(a + b*x + c*x**2)**(sympy.S(3)/2)/(4*c) + (-5*b*f + 8*c*e)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(24*c**2) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(-4*a*c*f + 5*b**2*f - 8*b*c*e + 16*c**2*d)/(64*c**3) - (-4*a*c + b**2)*(5*b**2*f + 16*c**2*d - 4*c*(a*f + 2*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_102():
    f = sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/f - sqrt(2)*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2)) + sqrt(2)*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_103():
    f = sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)**2
    F = -(e + 2*f*x)*sqrt(a + b*x + c*x**2)/((-4*d*f + e**2)*(d + e*x + f*x**2)) - sqrt(2)*(f*(-4*a*f + b*e) - (e - sqrt(-4*d*f + e**2))*(-b*f + c*e))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*(-4*d*f + e**2)**(sympy.S(3)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*(f*(-4*a*f + b*e) - (e + sqrt(-4*d*f + e**2))*(-b*f + c*e))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*(-4*d*f + e**2)**(sympy.S(3)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_104():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)**2
    F = f**2*x**3*(a + b*x + c*x**2)**(sympy.S(5)/2)/(8*c) + f*x**2*(-11*b*f + 32*c*e)*(a + b*x + c*x**2)**(sympy.S(5)/2)/(112*c**2) + x*(a + b*x + c*x**2)**(sympy.S(5)/2)*(99*b**2*f**2 + 224*c**2*(2*d*f + e**2) - 12*c*f*(7*a*f + 24*b*e))/(1344*c**3) + (a + b*x + c*x**2)**(sympy.S(5)/2)*(-693*b**3*f**2 + 36*b*c*f*(31*a*f + 56*b*e) + 5376*c**3*d*e - 32*c**2*(48*a*e*f + 49*b*(2*d*f + e**2)))/(13440*c**4) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(99*b**4*f**2 - 72*b**2*c*f*(3*a*f + 4*b*e) + 768*c**4*d**2 - 128*c**3*(a*(2*d*f + e**2) + 6*b*d*e) + 16*c**2*(3*a**2*f**2 + 24*a*b*e*f + 14*b**2*(2*d*f + e**2)))/(6144*c**5) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(99*b**4*f**2 - 72*b**2*c*f*(3*a*f + 4*b*e) + 768*c**4*d**2 - 128*c**3*(a*(2*d*f + e**2) + 6*b*d*e) + 16*c**2*(3*a**2*f**2 + 24*a*b*e*f + 14*b**2*(2*d*f + e**2)))/(16384*c**6) + (-4*a*c + b**2)**2*(99*b**4*f**2 - 72*b**2*c*f*(3*a*f + 4*b*e) + 768*c**4*d**2 - 128*c**3*(a*(2*d*f + e**2) + 6*b*d*e) + 16*c**2*(3*a**2*f**2 + 24*a*b*e*f + 14*b**2*(2*d*f + e**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(32768*c**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_105():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = f*x*(a + b*x + c*x**2)**(sympy.S(5)/2)/(6*c) + (-7*b*f + 12*c*e)*(a + b*x + c*x**2)**(sympy.S(5)/2)/(60*c**2) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-4*a*c*f + 7*b**2*f - 12*b*c*e + 24*c**2*d)/(192*c**3) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(7*b**2*f + 24*c**2*d - 4*c*(a*f + 3*b*e))/(512*c**4) + (-4*a*c + b**2)**2*(7*b**2*f + 24*c**2*d - 4*c*(a*f + 3*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_106():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)
    F = -sqrt(a + b*x + c*x**2)*(-5*b*f + 4*c*e - 2*c*f*x)/(4*f**2) + sqrt(2)*(-2*f*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)) + (e - sqrt(-4*d*f + e**2))*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**3*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(-2*f*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)) + (e + sqrt(-4*d*f + e**2))*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**3*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + (3*b**2*f**2 + 8*c**2*(-d*f + e**2) - 12*c*f*(-a*f + b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_107():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)**2
    F = c**(sympy.S(3)/2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/f**2 - (e + 2*f*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/((-4*d*f + e**2)*(d + e*x + f*x**2)) - sqrt(a + b*x + c*x**2)*(-2*b*f + c*e - 2*c*f*x)/(f*(-4*d*f + e**2)) - sqrt(2)*(-2*f*(2*c**2*d*(-4*d*f + e**2) + f*(4*a*f*(a*f + c*d) + 2*b**2*d*f - b*e*(3*a*f + c*d))) + (e - sqrt(-4*d*f + e**2))*(-b*f + c*e)*(2*c*(-5*d*f + e**2) + f*(-2*a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(4*f**2*(-4*d*f + e**2)**(sympy.S(3)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*(-2*f*(2*c**2*d*(-4*d*f + e**2) + f*(4*a*f*(a*f + c*d) + 2*b**2*d*f - b*e*(3*a*f + c*d))) + (e + sqrt(-4*d*f + e**2))*(-b*f + c*e)*(2*c*(-5*d*f + e**2) + f*(-2*a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(4*f**2*(-4*d*f + e**2)**(sympy.S(3)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_108():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)**3
    F = -(e + 2*f*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/((-8*d*f + 2*e**2)*(d + e*x + f*x**2)**2) + sqrt(a + b*x + c*x**2)*(12*a*e*f - 3*b*(4*d*f + e**2) + 12*c*d*e + 3*x*(8*a*f**2 - 4*b*e*f + 2*c*e**2))/(4*(-4*d*f + e**2)**2*(d + e*x + f*x**2)) - sqrt(2)*(-3*f*(-4*a*(4*a*f**2 + c*e**2) - b**2*(4*d*f + e**2) + 4*b*e*(3*a*f + c*d)) + 3*(e - sqrt(-4*d*f + e**2))*(-b*f + c*e)*(4*a*f - 2*b*e + 4*c*d))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(8*(-4*d*f + e**2)**(sympy.S(5)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*(-3*f*(-4*a*(4*a*f**2 + c*e**2) - b**2*(4*d*f + e**2) + 4*b*e*(3*a*f + c*d)) + 3*(e + sqrt(-4*d*f + e**2))*(-b*f + c*e)*(4*a*f - 2*b*e + 4*c*d))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(8*(-4*d*f + e**2)**(sympy.S(5)/2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_109():
    f = (d + e*x + f*x**2)**3/sqrt(a + b*x + c*x**2)
    F = f**3*x**5*sqrt(a + b*x + c*x**2)/(6*c) + f**2*x**4*(-11*b*f + 36*c*e)*sqrt(a + b*x + c*x**2)/(60*c**2) + f*x**3*sqrt(a + b*x + c*x**2)*(99*b**2*f**2 + 360*c**2*(d*f + e**2) - 4*c*f*(25*a*f + 81*b*e))/(480*c**3) - x**2*sqrt(a + b*x + c*x**2)*(231*b**3*f**3 - 36*b*c*f**2*(13*a*f + 21*b*e) - 320*c**3*(6*d*e*f + e**3) + 24*c**2*f*(32*a*e*f + 35*b*(d*f + e**2)))/(960*c**4) + x*sqrt(a + b*x + c*x**2)*(1155*b**4*f**3 - 252*b**2*c*f**2*(14*a*f + 15*b*e) + 5760*c**4*d*(d*f + e**2) - 160*c**3*(27*a*f*(d*f + e**2) + 10*b*(6*d*e*f + e**3)) + 24*c**2*f*(50*a**2*f**2 + 322*a*b*e*f + 175*b**2*(d*f + e**2)))/(3840*c**5) + sqrt(a + b*x + c*x**2)*(-3465*b**5*f**3 + 420*b**3*c*f**2*(34*a*f + 27*b*e) - 504*b*c**2*f*(22*a**2*f**2 + 70*a*b*e*f + 25*b**2*(d*f + e**2)) + 23040*c**5*d**2*e - 640*c**4*(8*a*e*(6*d*f + e**2) + 27*b*d*(d*f + e**2)) + 96*c**3*(128*a**2*e*f**2 + 275*a*b*f*(d*f + e**2) + 50*b**2*(6*d*e*f + e**3)))/(7680*c**6) + (231*b**6*f**3 - 252*b**4*c*f**2*(5*a*f + 3*b*e) + 840*b**2*c**2*f*(2*a**2*f**2 + 4*a*b*e*f + b**2*(d*f + e**2)) + 1024*c**6*d**3 - 1536*c**5*d*(a*(d*f + e**2) + b*d*e) + 384*c**4*(3*a**2*f*(d*f + e**2) + 2*a*b*e*(6*d*f + e**2) + 3*b**2*d*(d*f + e**2)) - 320*c**3*(a**3*f**3 + 9*a**2*b*e*f**2 + 9*a*b**2*f*(d*f + e**2) + b**3*(6*d*e*f + e**3)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_110():
    f = (d + e*x + f*x**2)**2/sqrt(a + b*x + c*x**2)
    F = f**2*x**3*sqrt(a + b*x + c*x**2)/(4*c) + f*x**2*(-7*b*f + 16*c*e)*sqrt(a + b*x + c*x**2)/(24*c**2) + x*sqrt(a + b*x + c*x**2)*(35*b**2*f**2 + 48*c**2*(2*d*f + e**2) - 4*c*f*(9*a*f + 20*b*e))/(96*c**3) + sqrt(a + b*x + c*x**2)*(-105*b**3*f**2 + 20*b*c*f*(11*a*f + 12*b*e) + 384*c**3*d*e - 16*c**2*(16*a*e*f + 9*b*(2*d*f + e**2)))/(192*c**4) + (35*b**4*f**2 - 40*b**2*c*f*(3*a*f + 2*b*e) + 128*c**4*d**2 - 64*c**3*(a*(2*d*f + e**2) + 2*b*d*e) + 48*c**2*(a**2*f**2 + 4*a*b*e*f + b**2*(2*d*f + e**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_111():
    f = (d + e*x + f*x**2)/sqrt(a + b*x + c*x**2)
    F = f*x*sqrt(a + b*x + c*x**2)/(2*c) + (-3*b*f + 4*c*e)*sqrt(a + b*x + c*x**2)/(4*c**2) + (3*b**2*f + 8*c**2*d - 4*c*(a*f + b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_112():
    f = 1/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_113():
    f = 1/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)**2)
    F = sqrt(a + b*x + c*x**2)*(-c*(-3*d*e*f + e**3) + f*x*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) + f*(-a*e*f - 2*b*d*f + b*e**2))/((-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*(d + e*x + f*x**2)) + sqrt(2)*(f*(e - sqrt(-4*d*f + e**2))*(-b*f + c*e)*(2*a*f - b*e + 2*c*d) - 2*f*(2*c**2*d*(-4*d*f + e**2) - c*(4*a*f*(-3*d*f + e**2) + b*(-5*d*e*f + e**3)) + f*(-4*a**2*f**2 + 3*a*b*e*f + b**2*(-6*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(4*(-4*d*f + e**2)**(sympy.S(3)/2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(f*(e + sqrt(-4*d*f + e**2))*(-b*f + c*e)*(2*a*f - b*e + 2*c*d) - 2*f*(2*c**2*d*(-4*d*f + e**2) - c*(4*a*f*(-3*d*f + e**2) + b*(-5*d*e*f + e**3)) + f*(-4*a**2*f**2 + 3*a*b*e*f + b**2*(-6*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(4*(-4*d*f + e**2)**(sympy.S(3)/2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_114():
    f = (d + e*x + f*x**2)**3/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = f**3*x**3*sqrt(a + b*x + c*x**2)/(4*c**2) + f**2*x**2*(-5*b*f + 8*c*e)*sqrt(a + b*x + c*x**2)/(8*c**3) + f*x*sqrt(a + b*x + c*x**2)*(41*b**2*f**2 + 48*c**2*(d*f + e**2) - 4*c*f*(7*a*f + 22*b*e))/(32*c**4) - sqrt(a + b*x + c*x**2)*(187*b**3*f**3 - 4*b*c*f**2*(73*a*f + 114*b*e) - 64*c**3*(6*d*e*f + e**3) + 16*c**2*f*(20*a*e*f + 21*b*(d*f + e**2)))/(64*c**5) + (-2*a*b**5*f**3 + 6*a*b**4*c*e*f**2 + 2*a*b**3*c*f*(5*a*f**2 - 3*c*(d*f + e**2)) - 2*a*b**2*c**2*e*(12*a*f**2 - c*(6*d*f + e**2)) + 4*a*c**3*e*(3*a**2*f**2 - a*c*(6*d*f + e**2) + 3*c**2*d**2) - 2*b*c**2*(5*a**3*f**3 - 9*a**2*c*f*(d*f + e**2) + 3*a*c**2*d*(d*f + e**2) + c**3*d**3) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d)*(a**2*c**2*f**2 - 4*a*b**2*c*f**2 + 7*a*b*c**2*e*f - 2*a*c**3*d*f - 3*a*c**3*e**2 + b**4*f**2 - 2*b**3*c*e*f + b**2*c**2*d*f + b**2*c**2*e**2 - b*c**3*d*e + c**4*d**2))/(c**5*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + (315*b**4*f**3 - 840*b**2*c*f**2*(a*f + b*e) + 384*c**4*d*(d*f + e**2) - 192*c**3*(3*a*f*(d*f + e**2) + b*(6*d*e*f + e**3)) + 240*c**2*f*(a**2*f**2 + 6*a*b*e*f + 3*b**2*(d*f + e**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_115():
    f = (d + e*x + f*x**2)**2/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = f**2*x*sqrt(a + b*x + c*x**2)/(2*c**2) + f*(-7*b*f + 8*c*e)*sqrt(a + b*x + c*x**2)/(4*c**3) + (-2*a*b**3*f**2 + 4*a*b**2*c*e*f + 8*a*c**2*e*(-a*f + c*d) - 2*b*c*(-3*a**2*f**2 + a*c*(2*d*f + e**2) + c**2*d**2) - 2*x*(b**4*f**2 - 2*b**2*c*f*(2*a*f + b*e) + 2*c**4*d**2 - 2*c**3*(a*(2*d*f + e**2) + b*d*e) + c**2*(2*a**2*f**2 + 6*a*b*e*f + b**2*(2*d*f + e**2))))/(c**3*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + (15*b**2*f**2 + 8*c**2*(2*d*f + e**2) - 12*c*f*(a*f + 2*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_116():
    f = (d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + f*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/c**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_117():
    f = 1/((a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)) - sqrt(2)*f*(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)) + (-4*a*c**2*e - 2*b**3*f + 2*b**2*c*e - 2*b*c*(-3*a*f + c*d) - 2*c*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/((-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_118():
    f = (d + e*x + f*x**2)**3/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = f**3*x*sqrt(a + b*x + c*x**2)/(2*c**3) + f**2*(-11*b*f + 12*c*e)*sqrt(a + b*x + c*x**2)/(4*c**4) + (-2*a*b**5*f**3 + 6*a*b**4*c*e*f**2 + 2*a*b**3*c*f*(5*a*f**2 - 3*c*(d*f + e**2)) - 2*a*b**2*c**2*e*(12*a*f**2 - c*(6*d*f + e**2)) + 4*a*c**3*e*(3*a**2*f**2 - a*c*(6*d*f + e**2) + 3*c**2*d**2) - 2*b*c**2*(5*a**3*f**3 - 9*a**2*c*f*(d*f + e**2) + 3*a*c**2*d*(d*f + e**2) + c**3*d**3) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d)*(a**2*c**2*f**2 - 4*a*b**2*c*f**2 + 7*a*b*c**2*e*f - 2*a*c**3*d*f - 3*a*c**3*e**2 + b**4*f**2 - 2*b**3*c*e*f + b**2*c**2*d*f + b**2*c**2*e**2 - b*c**3*d*e + c**4*d**2))/(3*c**5*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)) - (-48*a**2*c**4*e*(6*a*f**2 - c*(6*d*f + e**2)) - 2*b**7*f**3 + 6*b**6*c*e*f**2 + 6*b**5*c*f*(6*a*f**2 - c*(d*f + e**2)) - 2*b**4*c**2*e*(42*a*f**2 - c*(6*d*f + e**2)) - 6*b**3*c**2*(29*a**2*f**3 - 10*a*c*f*(d*f + e**2) + c**2*d*(d*f + e**2)) + 12*b**2*c**3*e*(28*a**2*f**2 - a*c*(6*d*f + e**2) + 2*c**2*d**2) - 8*b*c**3*(-29*a**3*f**3 + 24*a**2*c*f*(d*f + e**2) + 3*a*c**2*d*(d*f + e**2) + 2*c**3*d**3) - 2*c*x*(-10*b**6*f**3 + 3*b**4*c*f**2*(26*a*f + 7*b*e) - 6*b**2*c**2*f*(27*a**2*f**2 + 25*a*b*e*f + 2*b**2*(d*f + e**2)) + 16*c**6*d**3 - 24*c**5*d*(-a*(d*f + e**2) + b*d*e) + 6*c**4*(-16*a**2*f*(d*f + e**2) - 2*a*b*e*(6*d*f + e**2) + b**2*d*(d*f + e**2)) + c**3*(56*a**3*f**3 + 240*a**2*b*e*f**2 + 84*a*b**2*f*(d*f + e**2) + b**3*(6*d*e*f + e**3))))/(3*c**5*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)) + f*(35*b**2*f**2 + 24*c**2*(d*f + e**2) - 20*c*f*(a*f + 3*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_119():
    f = (d + e*x + f*x**2)**2/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = (-2*a*b**3*f**2 + 4*a*b**2*c*e*f + 8*a*c**2*e*(-a*f + c*d) - 2*b*c*(-3*a**2*f**2 + a*c*(2*d*f + e**2) + c**2*d**2) - 2*x*(b**4*f**2 - 2*b**2*c*f*(2*a*f + b*e) + 2*c**4*d**2 - 2*c**3*(a*(2*d*f + e**2) + b*d*e) + c**2*(2*a**2*f**2 + 6*a*b*e*f + b**2*(2*d*f + e**2))))/(3*c**3*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)) - (96*a**2*c**3*e*f - 2*b**5*f**2 + 4*b**4*c*e*f + 2*b**3*c*(10*a*f**2 - c*(2*d*f + e**2)) + 8*b**2*c**2*e*(-3*a*f + 2*c*d) - 8*b*c**2*(8*a**2*f**2 + a*c*(2*d*f + e**2) + 2*c**2*d**2) - 4*c*x*(-2*b**4*f**2 + b**2*c*f*(14*a*f + b*e) + 8*c**4*d**2 - c**3*(-4*a*(2*d*f + e**2) + 8*b*d*e) - c**2*(16*a**2*f**2 + 12*a*b*e*f - b**2*(2*d*f + e**2))))/(3*c**3*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)) + f**2*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/c**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_120():
    f = (d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = (b + 2*c*x)*(8*a*f + 2*b**2*f/c - 8*b*e + 16*c*d)/(3*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)) + (2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(3*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_121():
    f = 1/(sqrt(5*x**2 + 2*x - 7)*(5*x**2 + 12*x + 8))
    F = atan((5*x + 10)/(2*sqrt(5*x**2 + 2*x - 7)))/10 + atanh((5*x + 5)/sqrt(5*x**2 + 2*x - 7))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_122():
    f = 1/(sqrt(a + b*x + c*x**2)*sqrt(d + e*x + f*x**2))
    F = -sqrt(((2*a + x*(b + sqrt(-4*a*c + b**2)))**2*(4*c**2*d - 2*c*e*(b + sqrt(-4*a*c + b**2)) + f*(b + sqrt(-4*a*c + b**2))**2)/((b + 2*c*x + sqrt(-4*a*c + b**2))**2*(4*a**2*f - 2*a*e*(b + sqrt(-4*a*c + b**2)) + d*(b + sqrt(-4*a*c + b**2))**2)) - (2*a + x*(b + sqrt(-4*a*c + b**2)))*(b + sqrt(-4*a*c + b**2))*(2*a*f - b*e + 2*c*d)/((b + 2*c*x + sqrt(-4*a*c + b**2))*(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))) + 1)/((2*a + x*(b + sqrt(-4*a*c + b**2)))*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))/((b + 2*c*x + sqrt(-4*a*c + b**2))*sqrt(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))) + 1)**2)*sqrt((4*a*c - (b + sqrt(-4*a*c + b**2))**2)**2*(d + e*x + f*x**2)/((b + 2*c*x + sqrt(-4*a*c + b**2))**2*(4*a**2*f - 2*a*e*(b + sqrt(-4*a*c + b**2)) + d*(b + sqrt(-4*a*c + b**2))**2)))*sqrt(2*a + x*(b + sqrt(-4*a*c + b**2)))*((2*a + x*(b + sqrt(-4*a*c + b**2)))*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))/((b + 2*c*x + sqrt(-4*a*c + b**2))*sqrt(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))) + 1)*(b + 2*c*x + sqrt(-4*a*c + b**2))**(sympy.S(3)/2)*(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))**(sympy.S(1)/4)*elliptic_f(2*atan(sqrt(2*a + x*(b + sqrt(-4*a*c + b**2)))*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))**(sympy.S(1)/4)/(sqrt(b + 2*c*x + sqrt(-4*a*c + b**2))*(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))**(sympy.S(1)/4))), (b + sqrt(-4*a*c + b**2))*(2*a*f - b*e + 2*c*d)/(4*sqrt(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))*sqrt(b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d - c*(2*a*f + b*e + e*sqrt(-4*a*c + b**2)))) + sympy.S.Half)/((4*a*c - (b + sqrt(-4*a*c + b**2))**2)*sqrt(a + b*x + c*x**2)*sqrt(d + e*x + f*x**2)*sqrt((2*a + x*(b + sqrt(-4*a*c + b**2)))**2*(4*c**2*d - 2*c*e*(b + sqrt(-4*a*c + b**2)) + f*(b + sqrt(-4*a*c + b**2))**2)/((b + 2*c*x + sqrt(-4*a*c + b**2))**2*(4*a**2*f - 2*a*e*(b + sqrt(-4*a*c + b**2)) + d*(b + sqrt(-4*a*c + b**2))**2)) - (2*a + x*(b + sqrt(-4*a*c + b**2)))*(b + sqrt(-4*a*c + b**2))*(2*a*f - b*e + 2*c*d)/((b + 2*c*x + sqrt(-4*a*c + b**2))*(-a*(-2*a*f + 2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*d + b*(-a*e + d*sqrt(-4*a*c + b**2)))) + 1)*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_5_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_123():
    f = 1/(sqrt(2*x**2 - x + 3)*sqrt(5*x**2 + 3*x + 2))
    F = sqrt(253)*sqrt((-(-11*sqrt(23) + 33*I)*(-x*(1 - sqrt(23)*I) + 6)**2/((sqrt(23) + 7*I)*(-4*x + 1 - sqrt(23)*I)**2) - (41*sqrt(23) + 41*I)*(-x*(1 - sqrt(23)*I) + 6)/((sqrt(23) + 7*I)*(-4*x + 1 - sqrt(23)*I)) + 11)/(-sqrt(-(-sqrt(23) + 3*I)/(sqrt(23) + 7*I))*(-x*(1 - sqrt(23)*I) + 6)/(-4*x + 1 - sqrt(23)*I) + 1)**2)*sqrt((-sqrt(23) + 11*I)*(5*x**2 + 3*x + 2)/((sqrt(23) + 7*I)*(-4*x + 1 - sqrt(23)*I)**2))*sqrt(-x*(1 - sqrt(23)*I) + 6)*(-sqrt(-(-sqrt(23) + 3*I)/(sqrt(23) + 7*I))*(-x*(1 - sqrt(23)*I) + 6)/(-4*x + 1 - sqrt(23)*I) + 1)*(-4*x + 1 - sqrt(23)*I)*sqrt(4*x - 1 + sqrt(23)*I)*elliptic_f(2*atan((-(-sqrt(23) + 3*I)/(sqrt(23) + 7*I))**(sympy.S(1)/4)*sqrt(-x*(1 - sqrt(23)*I) + 6)/sqrt(4*x - 1 + sqrt(23)*I)), (-22*sqrt(23) + 41*sqrt(-(-23*sqrt(23) + 69*I)/(sqrt(23) + 7*I)) + 41*I*sqrt(-(-sqrt(23) + 3*I)/(sqrt(23) + 7*I)) + 66*I)/(-44*sqrt(23) + 132*I))/(11*(-(-sqrt(23) + 3*I)/(sqrt(23) + 7*I))**(sympy.S(1)/4)*(23 + sqrt(23)*I)*sqrt(2*x**2 - x + 3)*sqrt(5*x**2 + 3*x + 2)*sqrt(-(-11*sqrt(23) + 33*I)*(-x*(1 - sqrt(23)*I) + 6)**2/((sqrt(23) + 7*I)*(-4*x + 1 - sqrt(23)*I)**2) - (41*sqrt(23) + 41*I)*(-x*(1 - sqrt(23)*I) + 6)/((sqrt(23) + 7*I)*(-4*x + 1 - sqrt(23)*I)) + 11))
    assert integrate(f, x) == F

