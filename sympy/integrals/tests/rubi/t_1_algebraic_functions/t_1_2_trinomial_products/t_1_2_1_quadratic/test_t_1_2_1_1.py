"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.1 (a+b x+c x^2)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, k, n, p = symbols('a b c k n p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_1():
    f = (b*x + c*x**2)**(sympy.S(7)/2)
    F = 35*b**8*atanh(sqrt(c)*x/sqrt(b*x + c*x**2))/(16384*c**(sympy.S(9)/2)) - 35*b**6*(b + 2*c*x)*sqrt(b*x + c*x**2)/(16384*c**4) + 35*b**4*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(3)/2)/(6144*c**3) - 7*b**2*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(5)/2)/(384*c**2) + (b + 2*c*x)*(b*x + c*x**2)**(sympy.S(7)/2)/(16*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_2():
    f = (4*x**2 + 3*I*x)**(sympy.S(7)/2)
    F = (x/8 + 3*I/64)*(4*x**2 + 3*I*x)**(sympy.S(7)/2) + (168*x + 63*I)*(4*x**2 + 3*I*x)**(sympy.S(5)/2)/2048 + (7560*x + 2835*I)*(4*x**2 + 3*I*x)**(sympy.S(3)/2)/131072 + (204120*x + 76545*I)*sqrt(4*x**2 + 3*I*x)/4194304 - 229635*I*asin(8*I*x/3 - 1)/16777216
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_3():
    f = (4*x**2 + 3*I*x)**(sympy.S(5)/2)
    F = (x/6 + I/16)*(4*x**2 + 3*I*x)**(sympy.S(5)/2) + (120*x + 45*I)*(4*x**2 + 3*I*x)**(sympy.S(3)/2)/1024 + (3240*x + 1215*I)*sqrt(4*x**2 + 3*I*x)/32768 - 3645*I*asin(8*I*x/3 - 1)/131072
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_4():
    f = (4*x**2 + 3*I*x)**(sympy.S(3)/2)
    F = (x/4 + 3*I/32)*(4*x**2 + 3*I*x)**(sympy.S(3)/2) + (216*x + 81*I)*sqrt(4*x**2 + 3*I*x)/1024 - 243*I*asin(8*I*x/3 - 1)/4096
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_5():
    f = sqrt(4*x**2 + 3*I*x)
    F = (x/2 + 3*I/16)*sqrt(4*x**2 + 3*I*x) - 9*I*asin(8*I*x/3 - 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_6():
    f = (-4*x**2 + 3*x)**(sympy.S(7)/2)
    F = -(sympy.S(3)/64 - x/8)*(-4*x**2 + 3*x)**(sympy.S(7)/2) - (63 - 168*x)*(-4*x**2 + 3*x)**(sympy.S(5)/2)/2048 - (2835 - 7560*x)*(-4*x**2 + 3*x)**(sympy.S(3)/2)/131072 - (76545 - 204120*x)*sqrt(-4*x**2 + 3*x)/4194304 + 229635*asin(8*x/3 - 1)/16777216
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_7():
    f = (-4*x**2 + 3*x)**(sympy.S(5)/2)
    F = -(sympy.S(1)/16 - x/6)*(-4*x**2 + 3*x)**(sympy.S(5)/2) - (45 - 120*x)*(-4*x**2 + 3*x)**(sympy.S(3)/2)/1024 - (1215 - 3240*x)*sqrt(-4*x**2 + 3*x)/32768 + 3645*asin(8*x/3 - 1)/131072
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_8():
    f = (-4*x**2 + 3*x)**(sympy.S(3)/2)
    F = -(sympy.S(3)/32 - x/4)*(-4*x**2 + 3*x)**(sympy.S(3)/2) - (81 - 216*x)*sqrt(-4*x**2 + 3*x)/1024 + 243*asin(8*x/3 - 1)/4096
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_9():
    f = sqrt(-4*x**2 + 3*x)
    F = (x/2 + sympy.S(-3)/16)*sqrt(-4*x**2 + 3*x) + 9*asin(8*x/3 - 1)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_10():
    f = sqrt(-x**2 + 6*x)
    F = (x/2 + sympy.S(-3)/2)*sqrt(-x**2 + 6*x) + 9*asin(x/3 - 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_11():
    f = sqrt(-9*x**2 + 5*x)
    F = (x/2 + sympy.S(-5)/36)*sqrt(-9*x**2 + 5*x) + 25*asin(18*x/5 - 1)/216
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_12():
    f = (-x**2 + x)**(sympy.S(3)/2)
    F = -(sympy.S(1)/8 - x/4)*(-x**2 + x)**(sympy.S(3)/2) + (3*x/32 + sympy.S(-3)/64)*sqrt(-x**2 + x) + 3*asin(2*x - 1)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_13():
    f = sqrt(x**2 + 4*x)
    F = (x/2 + 1)*sqrt(x**2 + 4*x) - 4*atanh(x/sqrt(x**2 + 4*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_14():
    f = sqrt(x**2 - 8*x)
    F = (x/2 - 2)*sqrt(x**2 - 8*x) - 16*atanh(x/sqrt(x**2 - 8*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_15():
    f = sqrt(x**2 - x)
    F = (x/2 + sympy.S(-1)/4)*sqrt(x**2 - x) - atanh(x/sqrt(x**2 - x))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_16():
    f = (b*x + c*x**2)**(sympy.S(-7)/2)
    F = -(2*b + 4*c*x)/(5*b**2*(b*x + c*x**2)**(sympy.S(5)/2)) + 32*c*(b + 2*c*x)/(15*b**4*(b*x + c*x**2)**(sympy.S(3)/2)) - 256*c**2*(b + 2*c*x)/(15*b**6*sqrt(b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_17():
    f = 1/sqrt(4*x**2 + 3*I*x)
    F = -I*asin(8*I*x/3 - 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_18():
    f = (4*x**2 + 3*I*x)**(sympy.S(-3)/2)
    F = (16*x + 6*I)/(9*sqrt(4*x**2 + 3*I*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_19():
    f = (4*x**2 + 3*I*x)**(sympy.S(-5)/2)
    F = (16*x + 6*I)/(27*(4*x**2 + 3*I*x)**(sympy.S(3)/2)) + (512*x + 192*I)/(243*sqrt(4*x**2 + 3*I*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_20():
    f = (4*x**2 + 3*I*x)**(sympy.S(-7)/2)
    F = (16*x + 6*I)/(45*(4*x**2 + 3*I*x)**(sympy.S(5)/2)) + (1024*x + 384*I)/(1215*(4*x**2 + 3*I*x)**(sympy.S(3)/2)) + (32768*x + 12288*I)/(10935*sqrt(4*x**2 + 3*I*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_21():
    f = 1/sqrt(-4*x**2 + 3*x)
    F = asin(8*x/3 - 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_22():
    f = (-4*x**2 + 3*x)**(sympy.S(-3)/2)
    F = -(6 - 16*x)/(9*sqrt(-4*x**2 + 3*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_23():
    f = (-4*x**2 + 3*x)**(sympy.S(-5)/2)
    F = -(6 - 16*x)/(27*(-4*x**2 + 3*x)**(sympy.S(3)/2)) - (192 - 512*x)/(243*sqrt(-4*x**2 + 3*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_24():
    f = (-4*x**2 + 3*x)**(sympy.S(-7)/2)
    F = -(6 - 16*x)/(45*(-4*x**2 + 3*x)**(sympy.S(5)/2)) - (384 - 1024*x)/(1215*(-4*x**2 + 3*x)**(sympy.S(3)/2)) - (12288 - 32768*x)/(10935*sqrt(-4*x**2 + 3*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_25():
    f = 1/sqrt(-b**2*x**2 + b*x)
    F = asin(2*b*x - 1)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_26():
    f = 1/sqrt(b**2*x**2 + b*x)
    F = 2*atanh(b*x/sqrt(b**2*x**2 + b*x))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_27():
    f = 1/sqrt(-x**2 + 6*x)
    F = asin(x/3 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_28():
    f = 1/sqrt(x**2 + 4*x)
    F = 2*atanh(x/sqrt(x**2 + 4*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_29():
    f = 1/sqrt(x**2 - 2*x)
    F = 2*atanh(x/sqrt(x**2 - 2*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_30():
    f = (b*x + c*x**2)**(sympy.S(4)/3)
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(4)/3)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(55*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)*(b + 2*c*x)) + 3*(-c*x*(b + c*x)/b**2)**(sympy.S(4)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(4)/3)/(22*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)) + 3*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(4)/3)/(55*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_31():
    f = (b*x + c*x**2)**(sympy.S(1)/3)
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(1)/3)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(10*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/3)*(b + 2*c*x)) + 3*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(1)/3)/(10*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_32():
    f = (b*x + c*x**2)**(sympy.S(-2)/3)
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(2)/3)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_33():
    f = (b*x + c*x**2)**(sympy.S(-5)/3)
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(5)/3)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)*(3*b + 6*c*x)/(2*c*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3)*(b*x + c*x**2)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_34():
    f = (b*x + c*x**2)**(sympy.S(-8)/3)
    F = 14*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(8)/3)*sqrt(2 - sqrt(3))*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(5*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(8)/3)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(8)/3)*(21*b + 42*c*x)/(5*c*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3)*(b*x + c*x**2)**(sympy.S(8)/3)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(8)/3)*(3*b + 6*c*x)/(5*c*(-c*x*(b + c*x)/b**2)**(sympy.S(5)/3)*(b*x + c*x**2)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_35():
    f = (b*x + c*x**2)**(sympy.S(5)/3)
    F = -15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(5)/3)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(728*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)*(b + 2*c*x)) + 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(5)/3)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(182*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)*(b + 2*c*x)) + 3*(-c*x*(b + c*x)/b**2)**(sympy.S(5)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(5)/3)/(26*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)) + 15*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(5)/3)/(364*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)) - 2**(sympy.S(2)/3)*(15*b + 30*c*x)*(b*x + c*x**2)**(sympy.S(5)/3)/(364*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(5)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_36():
    f = (b*x + c*x**2)**(sympy.S(2)/3)
    F = -3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(2)/3)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(28*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(2)/3)*(b + 2*c*x)) + 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*(b*x + c*x**2)**(sympy.S(2)/3)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(7*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(2)/3)*(b + 2*c*x)) + 3*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(2)/3)/(14*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3*b + 6*c*x)*(b*x + c*x**2)**(sympy.S(2)/3)/(14*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(2)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_37():
    f = (b*x + c*x**2)**(sympy.S(-1)/3)
    F = -3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/3)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(1)/3)) + 2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(1)/3)) - 2**(sympy.S(2)/3)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/3)*(3*b + 6*c*x)/(2*c*(b*x + c*x**2)**(sympy.S(1)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_38():
    f = (b*x + c*x**2)**(sympy.S(-4)/3)
    F = 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(2*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(4)/3)) - 2*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(4)/3)) + 3*2**(sympy.S(2)/3)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)*(b + 2*c*x)/(c*(b*x + c*x**2)**(sympy.S(4)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(4)/3)*(3*b + 6*c*x)/(c*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3)*(b*x + c*x**2)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_39():
    f = (b*x + c*x**2)**(sympy.S(-7)/3)
    F = 15*2**(sympy.S(2)/3)*3**(sympy.S(1)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(7)/3)*sqrt(sqrt(3) + 2)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(4*c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(7)/3)) - 5*2**(sympy.S(1)/6)*3**(sympy.S(3)/4)*b**2*sqrt((2*2**(sympy.S(1)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(2)/3) + 2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(7)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)*elliptic_f(asin((-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1 + sqrt(3))/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/(c*sqrt(-(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) + 1)/(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(7)/3)) + 2**(sympy.S(2)/3)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(7)/3)*(15*b + 30*c*x)/(2*c*(b*x + c*x**2)**(sympy.S(7)/3)*(-2**(sympy.S(2)/3)*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3) - sqrt(3) + 1)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(7)/3)*(15*b + 30*c*x)/(2*c*(-c*x*(b + c*x)/b**2)**(sympy.S(1)/3)*(b*x + c*x**2)**(sympy.S(7)/3)) + (-c*(b*x + c*x**2)/b**2)**(sympy.S(7)/3)*(3*b + 6*c*x)/(4*c*(-c*x*(b + c*x)/b**2)**(sympy.S(4)/3)*(b*x + c*x**2)**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_40():
    f = (b*x + c*x**2)**(sympy.S(5)/4)
    F = 5*sqrt(2)*b**5*(-c*(b*x + c*x**2)/b**2)**(sympy.S(3)/4)*elliptic_f(asin(1 + 2*c*x/b)/2, 2)/(168*c**3*(b*x + c*x**2)**(sympy.S(3)/4)) - 5*b**2*(b + 2*c*x)*(b*x + c*x**2)**(sympy.S(1)/4)/(84*c**2) + (b + 2*c*x)*(b*x + c*x**2)**(sympy.S(5)/4)/(7*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_41():
    f = (b*x + c*x**2)**(sympy.S(3)/4)
    F = -3*sqrt(2)*b**3*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/4)*elliptic_e(asin(1 + 2*c*x/b)/2, 2)/(20*c**2*(b*x + c*x**2)**(sympy.S(1)/4)) + (b + 2*c*x)*(b*x + c*x**2)**(sympy.S(3)/4)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_42():
    f = (b*x + c*x**2)**(sympy.S(1)/4)
    F = -sqrt(2)*b**3*(-c*(b*x + c*x**2)/b**2)**(sympy.S(3)/4)*elliptic_f(asin(1 + 2*c*x/b)/2, 2)/(6*c**2*(b*x + c*x**2)**(sympy.S(3)/4)) + (b + 2*c*x)*(b*x + c*x**2)**(sympy.S(1)/4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_43():
    f = (b*x + c*x**2)**(sympy.S(-1)/4)
    F = sqrt(2)*b*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/4)*elliptic_e(asin(1 + 2*c*x/b)/2, 2)/(c*(b*x + c*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_44():
    f = (b*x + c*x**2)**(sympy.S(-3)/4)
    F = 2*sqrt(2)*b*(-c*(b*x + c*x**2)/b**2)**(sympy.S(3)/4)*elliptic_f(asin(1 + 2*c*x/b)/2, 2)/(c*(b*x + c*x**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_45():
    f = (b*x + c*x**2)**(sympy.S(-5)/4)
    F = 4*sqrt(2)*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/4)*elliptic_e(asin(1 + 2*c*x/b)/2, 2)/(b*(b*x + c*x**2)**(sympy.S(1)/4)) - (4*b + 8*c*x)/(b**2*(b*x + c*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_46():
    f = (b*x + c*x**2)**(sympy.S(-9)/4)
    F = -(4*b + 8*c*x)/(5*b**2*(b*x + c*x**2)**(sympy.S(5)/4)) - 48*sqrt(2)*c*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/4)*elliptic_e(asin(1 + 2*c*x/b)/2, 2)/(5*b**3*(b*x + c*x**2)**(sympy.S(1)/4)) + 48*c*(b + 2*c*x)/(5*b**4*(b*x + c*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_47():
    f = (b*x + c*x**2)**(sympy.S(-13)/4)
    F = -(4*b + 8*c*x)/(9*b**2*(b*x + c*x**2)**(sympy.S(9)/4)) + 112*c*(b + 2*c*x)/(45*b**4*(b*x + c*x**2)**(sympy.S(5)/4)) + 448*sqrt(2)*c**2*(-c*(b*x + c*x**2)/b**2)**(sympy.S(1)/4)*elliptic_e(asin(1 + 2*c*x/b)/2, 2)/(15*b**5*(b*x + c*x**2)**(sympy.S(1)/4)) - 448*c**2*(b + 2*c*x)/(15*b**6*(b*x + c*x**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_48():
    f = (b*x + c*x**2)**p
    F = -(-c*x/b)**(-p - 1)*(b*x + c*x**2)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + c*x)/b)/(b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_49():
    f = (a + c*x**2)**4
    F = a**4*x + 4*a**3*c*x**3/3 + 6*a**2*c**2*x**5/5 + 4*a*c**3*x**7/7 + c**4*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_50():
    f = (a + c*x**2)**3
    F = a**3*x + a**2*c*x**3 + 3*a*c**2*x**5/5 + c**3*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_51():
    f = (a + c*x**2)**2
    F = a**2*x + 2*a*c*x**3/3 + c**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_52():
    f = a + c*x**2
    F = a*x + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_53():
    f = 1/(a + c*x**2)
    F = atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_54():
    f = (a + c*x**2)**(-2)
    F = x/(2*a*(a + c*x**2)) + atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_55():
    f = (a + c*x**2)**(-3)
    F = x/(4*a*(a + c*x**2)**2) + 3*x/(8*a**2*(a + c*x**2)) + 3*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_56():
    f = (a + c*x**2)**(sympy.S(5)/2)
    F = 5*a**3*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(16*sqrt(c)) + 5*a**2*x*sqrt(a + c*x**2)/16 + 5*a*x*(a + c*x**2)**(sympy.S(3)/2)/24 + x*(a + c*x**2)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_57():
    f = (a + c*x**2)**(sympy.S(3)/2)
    F = 3*a**2*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*sqrt(c)) + 3*a*x*sqrt(a + c*x**2)/8 + x*(a + c*x**2)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_58():
    f = sqrt(a + c*x**2)
    F = a*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)) + x*sqrt(a + c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_59():
    f = 1/sqrt(a + c*x**2)
    F = atanh(sqrt(c)*x/sqrt(a + c*x**2))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_60():
    f = (a + c*x**2)**(sympy.S(-3)/2)
    F = x/(a*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_61():
    f = (a + c*x**2)**(sympy.S(-5)/2)
    F = x/(3*a*(a + c*x**2)**(sympy.S(3)/2)) + 2*x/(3*a**2*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_62():
    f = (a + c*x**2)**(sympy.S(-7)/2)
    F = x/(5*a*(a + c*x**2)**(sympy.S(5)/2)) + 4*x/(15*a**2*(a + c*x**2)**(sympy.S(3)/2)) + 8*x/(15*a**3*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_63():
    f = (a + c*x**2)**(sympy.S(-9)/2)
    F = x/(7*a*(a + c*x**2)**(sympy.S(7)/2)) + 6*x/(35*a**2*(a + c*x**2)**(sympy.S(5)/2)) + 8*x/(35*a**3*(a + c*x**2)**(sympy.S(3)/2)) + 16*x/(35*a**4*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_64():
    f = (9*x**2 + 12*x + 4)**(sympy.S(3)/2)
    F = (x/4 + sympy.S(1)/6)*(9*x**2 + 12*x + 4)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_65():
    f = sqrt(9*x**2 + 12*x + 4)
    F = (x/2 + sympy.S(1)/3)*sqrt(9*x**2 + 12*x + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_66():
    f = 1/sqrt(9*x**2 + 12*x + 4)
    F = (3*x + 2)*log(3*x + 2)/(3*sqrt(9*x**2 + 12*x + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_67():
    f = (9*x**2 + 12*x + 4)**(sympy.S(-3)/2)
    F = -1/((18*x + 12)*sqrt(9*x**2 + 12*x + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_68():
    f = sqrt(9*x**2 - 12*x + 4)
    F = (x/2 + sympy.S(-1)/3)*sqrt(9*x**2 - 12*x + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_69():
    f = 1/sqrt(9*x**2 - 12*x + 4)
    F = -(2 - 3*x)*log(2 - 3*x)/(3*sqrt(9*x**2 - 12*x + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_70():
    f = sqrt(-9*x**2 + 12*x - 4)
    F = (x/2 + sympy.S(-1)/3)*sqrt(-9*x**2 + 12*x - 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_71():
    f = 1/sqrt(-9*x**2 + 12*x - 4)
    F = -(2 - 3*x)*log(2 - 3*x)/(3*sqrt(-9*x**2 + 12*x - 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_72():
    f = sqrt(-9*x**2 - 12*x - 4)
    F = (x/2 + sympy.S(1)/3)*sqrt(-9*x**2 - 12*x - 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_73():
    f = 1/sqrt(-9*x**2 - 12*x - 4)
    F = (3*x + 2)*log(3*x + 2)/(3*sqrt(-9*x**2 - 12*x - 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_74():
    f = (b*x + c*x**2 + (b**2 - 1)/(4*c))**5
    F = -(-b - 2*c*x + 1)**11/(22528*c**6) + (-b - 2*c*x + 1)**10/(2048*c**6) - 5*(-b - 2*c*x + 1)**9/(2304*c**6) + 5*(-b - 2*c*x + 1)**8/(1024*c**6) - 5*(-b - 2*c*x + 1)**7/(896*c**6) + (-b - 2*c*x + 1)**6/(384*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_75():
    f = (b*x + c*x**2 + (b**2 - 4)/(4*c))**5
    F = -(-b - 2*c*x + 2)**11/(22528*c**6) + (-b - 2*c*x + 2)**10/(1024*c**6) - 5*(-b - 2*c*x + 2)**9/(576*c**6) + 5*(-b - 2*c*x + 2)**8/(128*c**6) - 5*(-b - 2*c*x + 2)**7/(56*c**6) + (-b - 2*c*x + 2)**6/(12*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_76():
    f = (b*x + c*x**2 + (b**2 - 9)/(4*c))**5
    F = -(-b - 2*c*x + 3)**11/(22528*c**6) + 3*(-b - 2*c*x + 3)**10/(2048*c**6) - 5*(-b - 2*c*x + 3)**9/(256*c**6) + 135*(-b - 2*c*x + 3)**8/(1024*c**6) - 405*(-b - 2*c*x + 3)**7/(896*c**6) + 81*(-b - 2*c*x + 3)**6/(128*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_77():
    f = (b*x + c*x**2 + (b**2 - 16)/(4*c))**5
    F = -(-b - 2*c*x + 4)**11/(22528*c**6) + (-b - 2*c*x + 4)**10/(512*c**6) - 5*(-b - 2*c*x + 4)**9/(144*c**6) + 5*(-b - 2*c*x + 4)**8/(16*c**6) - 10*(-b - 2*c*x + 4)**7/(7*c**6) + 8*(-b - 2*c*x + 4)**6/(3*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_78():
    f = 1/(3*x**2 + 4*x + 2)
    F = sqrt(2)*atan(sqrt(2)*(3*x + 2)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_79():
    f = 1/(x**2 - 2*sqrt(3)*x + 4)
    F = atan(x - sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_80():
    f = 1/(-3*x**2 + 4*x + 2)
    F = -sqrt(10)*atanh(sqrt(10)*(2 - 3*x)/10)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_81():
    f = 1/(3*x**2 + 5*x + 2)
    F = -log(x + 1) + log(3*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_82():
    f = 1/(-3*x**2 + 5*x + 2)
    F = -log(2 - x)/7 + log(3*x + 1)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_83():
    f = 1/(2*x**2 + pi*x + 1)
    F = -2*atanh((4*x + pi)/sqrt(-8 + pi**2))/sqrt(-8 + pi**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_84():
    f = 1/(-2*x**2 + pi*x + 1)
    F = -2*atanh((pi - 4*x)/sqrt(8 + pi**2))/sqrt(8 + pi**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_85():
    f = 1/(3*x**2 + pi*x + 1)
    F = 2*atan((6*x + pi)/sqrt(12 - pi**2))/sqrt(12 - pi**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_86():
    f = 1/(-3*x**2 + pi*x + 1)
    F = -2*atanh((pi - 6*x)/sqrt(pi**2 + 12))/sqrt(pi**2 + 12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_87():
    f = 1/(a + b*x**2 + c*x)
    F = 2*atan((2*b*x + c)/sqrt(4*a*b - c**2))/sqrt(4*a*b - c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_88():
    f = 1/(2*a*x + b*x**2 + b)
    F = -atanh((a + b*x)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_89():
    f = 1/(2*a*x - b*x**2 + b)
    F = -atanh((a - b*x)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_90():
    f = (3*x**2 + 4*x + 2)**(-2)
    F = (3*x + 2)/(12*x**2 + 16*x + 8) + 3*sqrt(2)*atan(sqrt(2)*(3*x + 2)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_91():
    f = (-3*x**2 + 4*x + 2)**(-2)
    F = -(2 - 3*x)/(-60*x**2 + 80*x + 40) - 3*sqrt(10)*atanh(sqrt(10)*(2 - 3*x)/10)/200
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_92():
    f = (3*x**2 + 5*x + 2)**(-2)
    F = -(6*x + 5)/(3*x**2 + 5*x + 2) + 6*log(x + 1) - 6*log(3*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_93():
    f = (-3*x**2 + 5*x + 2)**(-2)
    F = -(5 - 6*x)/(-147*x**2 + 245*x + 98) - 6*log(2 - x)/343 + 6*log(3*x + 1)/343
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_94():
    f = (a + b*x**2 + c*x)**(-2)
    F = 4*b*atan((2*b*x + c)/sqrt(4*a*b - c**2))/(4*a*b - c**2)**(sympy.S(3)/2) + (2*b*x + c)/((4*a*b - c**2)*(a + b*x**2 + c*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_95():
    f = (2*a*x + b*x**2 + b)**(-2)
    F = b*atanh((a + b*x)/sqrt(a**2 - b**2))/(2*(a**2 - b**2)**(sympy.S(3)/2)) - (a + b*x)/((2*a**2 - 2*b**2)*(2*a*x + b*x**2 + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_96():
    f = (2*a*x - b*x**2 + b)**(-2)
    F = -b*atanh((a - b*x)/sqrt(a**2 + b**2))/(2*(a**2 + b**2)**(sympy.S(3)/2)) - (a - b*x)/((2*a**2 + 2*b**2)*(2*a*x - b*x**2 + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_97():
    f = 1/(x**2 - 2*x*(a/b)**(1/n)*cos((-2*pi*k + pi)/n) + (a/b)**(2/n))
    F = atan(x*csc((-2*pi*k + pi)/n)/(a/b)**(1/n) - cot((-2*pi*k + pi)/n))*csc((-2*pi*k + pi)/n)/(a/b)**(1/n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_98():
    f = 1/(x**2 + 2*x*cos(sympy.S(1)/7) + 1)
    F = atan((x + cos(sympy.S(1)/7))*csc(sympy.S(1)/7))*csc(sympy.S(1)/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_99():
    f = 1/(x**2 + 2*x*cos(pi/7) + 1)
    F = atan(x*csc(pi/7) + cot(pi/7))*csc(pi/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_100():
    f = sqrt(9*x**2 - 6*x + 5)
    F = (x/2 + sympy.S(-1)/6)*sqrt(9*x**2 - 6*x + 5) + 2*asinh(3*x/2 + sympy.S(-1)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_101():
    f = sqrt(-4*x**2 - 4*x + 3)
    F = (x/2 + sympy.S(1)/4)*sqrt(-4*x**2 - 4*x + 3) + asin(x + sympy.S.Half)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_102():
    f = sqrt(9*x**2 + 6*x - 8)
    F = (x/2 + sympy.S(1)/6)*sqrt(9*x**2 + 6*x - 8) - 3*atanh((3*x + 1)/sqrt(9*x**2 + 6*x - 8))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_103():
    f = sqrt(3*x**2 + 4*x + 2)
    F = (x/2 + sympy.S(1)/3)*sqrt(3*x**2 + 4*x + 2) + sqrt(3)*asinh(sqrt(2)*(3*x + 2)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_104():
    f = sqrt(-3*x**2 + 4*x + 2)
    F = (x/2 + sympy.S(-1)/3)*sqrt(-3*x**2 + 4*x + 2) - 5*sqrt(3)*asin(sqrt(10)*(2 - 3*x)/10)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_105():
    f = sqrt(3*x**2 + 5*x + 2)
    F = (x/2 + sympy.S(5)/12)*sqrt(3*x**2 + 5*x + 2) - sqrt(3)*atanh(sqrt(3)*(6*x + 5)/(6*sqrt(3*x**2 + 5*x + 2)))/72
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_106():
    f = sqrt(-3*x**2 + 5*x + 2)
    F = (x/2 + sympy.S(-5)/12)*sqrt(-3*x**2 + 5*x + 2) + 49*sqrt(3)*asin(6*x/7 + sympy.S(-5)/7)/72
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_107():
    f = sqrt(3*x**2 + 4*x - 2)
    F = (x/2 + sympy.S(1)/3)*sqrt(3*x**2 + 4*x - 2) - 5*sqrt(3)*atanh(sqrt(3)*(3*x + 2)/(3*sqrt(3*x**2 + 4*x - 2)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_108():
    f = sqrt(-3*x**2 + 4*x - 2)
    F = (x/2 + sympy.S(-1)/3)*sqrt(-3*x**2 + 4*x - 2) + sqrt(3)*atan(sqrt(3)*(2 - 3*x)/(3*sqrt(-3*x**2 + 4*x - 2)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_109():
    f = sqrt(3*x**2 + 5*x - 2)
    F = (x/2 + sympy.S(5)/12)*sqrt(3*x**2 + 5*x - 2) - 49*sqrt(3)*atanh(sqrt(3)*(6*x + 5)/(6*sqrt(3*x**2 + 5*x - 2)))/72
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_110():
    f = sqrt(-3*x**2 + 5*x - 2)
    F = (x/2 + sympy.S(-5)/12)*sqrt(-3*x**2 + 5*x - 2) + sqrt(3)*asin(6*x - 5)/72
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_111():
    f = 1/sqrt(9*x**2 - 6*x + 5)
    F = asinh(3*x/2 + sympy.S(-1)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_112():
    f = 1/sqrt(-4*x**2 - 4*x + 3)
    F = asin(x + sympy.S.Half)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_113():
    f = 1/sqrt(9*x**2 + 6*x - 8)
    F = atanh((3*x + 1)/sqrt(9*x**2 + 6*x - 8))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_114():
    f = 1/sqrt(3*x**2 + 4*x + 2)
    F = sqrt(3)*asinh(sqrt(2)*(3*x + 2)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_115():
    f = 1/sqrt(-3*x**2 + 4*x + 2)
    F = -sqrt(3)*asin(sqrt(10)*(2 - 3*x)/10)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_116():
    f = 1/sqrt(3*x**2 + 5*x + 2)
    F = sqrt(3)*atanh(sqrt(3)*(6*x + 5)/(6*sqrt(3*x**2 + 5*x + 2)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_117():
    f = 1/sqrt(-3*x**2 + 5*x + 2)
    F = sqrt(3)*asin(6*x/7 + sympy.S(-5)/7)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_118():
    f = 1/sqrt(3*x**2 + 4*x - 2)
    F = sqrt(3)*atanh(sqrt(3)*(3*x + 2)/(3*sqrt(3*x**2 + 4*x - 2)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_119():
    f = 1/sqrt(-3*x**2 + 4*x - 2)
    F = -sqrt(3)*atan(sqrt(3)*(2 - 3*x)/(3*sqrt(-3*x**2 + 4*x - 2)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_120():
    f = 1/sqrt(3*x**2 + 5*x - 2)
    F = sqrt(3)*atanh(sqrt(3)*(6*x + 5)/(6*sqrt(3*x**2 + 5*x - 2)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_121():
    f = 1/sqrt(-3*x**2 + 5*x - 2)
    F = sqrt(3)*asin(6*x - 5)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_122():
    f = 1/sqrt(b*x + c*x**2 + (b**2 + 4*c)/(4*c))
    F = asinh((b + 2*c*x)/(2*sqrt(c)))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_123():
    f = 1/sqrt(b*x - c*x**2 + (-b**2 + 4*c)/(4*c))
    F = -asin((b - 2*c*x)/(2*sqrt(c)))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_124():
    f = 1/sqrt(b*x - c*x**2 + (-b**2 + c)/(4*c))
    F = -asin((b - 2*c*x)/sqrt(c))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_125():
    f = (x**2 + 3*x + 2)**(sympy.S(-3)/2)
    F = (-4*x - 6)/sqrt(x**2 + 3*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_126():
    f = (4*x**2 - 24*x + 27)**(sympy.S(-3)/2)
    F = (3 - x)/(9*sqrt(4*x**2 - 24*x + 27))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_127():
    f = x/(-x**2 - 4*x + 5)**(sympy.S(3)/2)
    F = (5 - 2*x)/(9*sqrt(-x**2 - 4*x + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_128():
    f = (-x**2 - 4*x + 5)**(sympy.S(-5)/2)
    F = (x + 2)/(27*(-x**2 - 4*x + 5)**(sympy.S(3)/2)) + (2*x + 4)/(243*sqrt(-x**2 - 4*x + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_129():
    f = (a + b*x + c*x**2)**p
    F = -2**(p + 1)*(-(b + 2*c*x - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(a + b*x + c*x**2)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/((p + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_130():
    f = (5*x**2 + 4*x + 3)**p
    F = 11**p*5**(-p - 1)*(5*x + 2)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -(5*x + 2)**2/11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_131():
    f = (4*x**2 + 4*x + 3)**p
    F = 2**(p - 1)*(2*x + 1)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -(2*x + 1)**2/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_132():
    f = (3*x**2 + 4*x + 3)**p
    F = 3**(-p - 1)*5**p*(3*x + 2)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -(3*x + 2)**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_133():
    f = (2*x**2 + 4*x + 3)**p
    F = (x + 1)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -2*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_134():
    f = (x**2 + 4*x + 3)**p
    F = -2**(2*p + 1)*(-2*x - 2)**(-p - 1)*(x**2 + 4*x + 3)**(p + 1)*hyper((-p, p + 1), (p + 2,), x/2 + sympy.S(3)/2)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_135():
    f = (4*x + 3)**p
    F = (4*x + 3)**(p + 1)/(4*p + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_136():
    f = (-x**2 + 4*x + 3)**p
    F = -7**p*(2 - x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), (2 - x)**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_137():
    f = (-2*x**2 + 4*x + 3)**p
    F = -5**p*(1 - x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), 2*(1 - x)**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_138():
    f = (-3*x**2 + 4*x + 3)**p
    F = -13**p*3**(-p - 1)*(2 - 3*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), (2 - 3*x)**2/13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_139():
    f = (-4*x**2 + 4*x + 3)**p
    F = -2**(2*p - 1)*(1 - 2*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), (1 - 2*x)**2/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_1_a_plus_b_x_plus_c_x_pow_2_pow_p_140():
    f = (-5*x**2 + 4*x + 3)**p
    F = -19**p*5**(-p - 1)*(2 - 5*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), (2 - 5*x)**2/19)
    assert integrate(f, x) == F

