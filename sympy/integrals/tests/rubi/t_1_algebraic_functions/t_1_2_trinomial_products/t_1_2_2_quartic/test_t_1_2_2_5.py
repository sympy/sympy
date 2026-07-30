"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.5 P(x) (a+b x^2+c x^4)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, h, i, j, k, l, m = symbols('a b c d e f g h i j k l m')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1():
    f = (d + e*x)*(a + b*x**2 + c*x**4)
    F = a*d*x + a*e*x**2/2 + b*d*x**3/3 + b*e*x**4/4 + c*d*x**5/5 + c*e*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_2():
    f = (a + b*x**2 + c*x**4)*(d + e*x + f*x**2)
    F = a*d*x + a*e*x**2/2 + b*e*x**4/4 + c*e*x**6/6 + c*f*x**7/7 + x**5*(b*f/5 + c*d/5) + x**3*(a*f/3 + b*d/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_3():
    f = (a + b*x**2 + c*x**4)*(d + e*x + f*x**2 + g*x**3)
    F = a*d*x + a*e*x**2/2 + c*f*x**7/7 + c*g*x**8/8 + x**6*(b*g/6 + c*e/6) + x**5*(b*f/5 + c*d/5) + x**4*(a*g/4 + b*e/4) + x**3*(a*f/3 + b*d/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_4():
    f = (a + b*x**2 + c*x**4)*(d + e*x + f*x**2 + g*x**3 + h*x**4)
    F = a*d*x + a*e*x**2/2 + c*g*x**8/8 + c*h*x**9/9 + x**7*(b*h/7 + c*f/7) + x**6*(b*g/6 + c*e/6) + x**5*(a*h/5 + b*f/5 + c*d/5) + x**4*(a*g/4 + b*e/4) + x**3*(a*f/3 + b*d/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_5():
    f = (a + b*x**2 + c*x**4)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)
    F = a*d*x + a*e*x**2/2 + c*h*x**9/9 + c*i*x**10/10 + x**8*(b*i/8 + c*g/8) + x**7*(b*h/7 + c*f/7) + x**6*(a*i/6 + b*g/6 + c*e/6) + x**5*(a*h/5 + b*f/5 + c*d/5) + x**4*(a*g/4 + b*e/4) + x**3*(a*f/3 + b*d/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_6():
    f = (d + e*x)*(a + b*x**2 + c*x**4)**2
    F = a**2*d*x + a**2*e*x**2/2 + 2*a*b*d*x**3/3 + a*b*e*x**4/2 + 2*b*c*d*x**7/7 + b*c*e*x**8/4 + c**2*d*x**9/9 + c**2*e*x**10/10 + d*x**5*(2*a*c/5 + b**2/5) + e*x**6*(a*c/3 + b**2/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_7():
    f = (a + b*x**2 + c*x**4)**2*(d + e*x + f*x**2)
    F = a**2*d*x + a**2*e*x**2/2 + a*b*e*x**4/2 + a*x**3*(a*f + 2*b*d)/3 + b*c*e*x**8/4 + c**2*e*x**10/10 + c**2*f*x**11/11 + c*x**9*(2*b*f + c*d)/9 + e*x**6*(a*c/3 + b**2/6) + x**7*(2*a*c*f/7 + b**2*f/7 + 2*b*c*d/7) + x**5*(2*a*b*f/5 + 2*a*c*d/5 + b**2*d/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_8():
    f = (a + b*x**2 + c*x**4)**2*(d + e*x + f*x**2 + g*x**3)
    F = a**2*d*x + a**2*e*x**2/2 + a*x**4*(a*g + 2*b*e)/4 + a*x**3*(a*f + 2*b*d)/3 + c**2*f*x**11/11 + c**2*g*x**12/12 + c*x**10*(2*b*g + c*e)/10 + c*x**9*(2*b*f + c*d)/9 + x**8*(a*c*g/4 + b**2*g/8 + b*c*e/4) + x**7*(2*a*c*f/7 + b**2*f/7 + 2*b*c*d/7) + x**6*(a*b*g/3 + a*c*e/3 + b**2*e/6) + x**5*(2*a*b*f/5 + 2*a*c*d/5 + b**2*d/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_9():
    f = (a + b*x**2 + c*x**4)**2*(d + e*x + f*x**2 + g*x**3 + h*x**4)
    F = a**2*d*x + a**2*e*x**2/2 + a*x**4*(a*g + 2*b*e)/4 + a*x**3*(a*f + 2*b*d)/3 + c**2*g*x**12/12 + c**2*h*x**13/13 + c*x**11*(2*b*h + c*f)/11 + c*x**10*(2*b*g + c*e)/10 + x**9*(b**2*h/9 + c**2*d/9 + 2*c*(a*h + b*f)/9) + x**8*(a*c*g/4 + b**2*g/8 + b*c*e/4) + x**7*(2*a*c*f/7 + b**2*f/7 + 2*b*(a*h + c*d)/7) + x**6*(a*b*g/3 + a*c*e/3 + b**2*e/6) + x**5*(2*a*b*f/5 + a*(a*h + 2*c*d)/5 + b**2*d/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_10():
    f = (d + e*x)/(x**4 - 5*x**2 + 4)
    F = -d*atanh(x/2)/6 + d*atanh(x)/3 - e*log(1 - x**2)/6 + e*log(4 - x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_11():
    f = (d + e*x + f*x**2)/(x**4 - 5*x**2 + 4)
    F = -e*log(1 - x**2)/6 + e*log(4 - x**2)/6 + (-d/6 - 2*f/3)*atanh(x/2) + (d/3 + f/3)*atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_12():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)
    F = (-d/6 - 2*f/3)*atanh(x/2) + (d/3 + f/3)*atanh(x) - (e/6 + g/6)*log(1 - x**2) + (e/6 + 2*g/3)*log(4 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_13():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)
    F = h*x - (e/6 + g/6)*log(1 - x**2) + (e/6 + 2*g/3)*log(4 - x**2) - (d/6 + 2*f/3 + 8*h/3)*atanh(x/2) + (d/3 + f/3 + h/3)*atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_14():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)
    F = h*x + i*x**2/2 - (d/6 + 2*f/3 + 8*h/3)*atanh(x/2) + (d/3 + f/3 + h/3)*atanh(x) - (e/6 + g/6 + i/6)*log(1 - x**2) + (e/6 + 2*g/3 + 8*i/3)*log(4 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_15():
    f = (d + e*x)/(x**4 + x**2 + 1)
    F = -d*log(x**2 - x + 1)/4 + d*log(x**2 + x + 1)/4 - sqrt(3)*d*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*d*atan(sqrt(3)*(2*x + 1)/3)/6 + sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_16():
    f = (d + e*x + f*x**2)/(x**4 + x**2 + 1)
    F = sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/3 - (d/4 - f/4)*log(x**2 - x + 1) + (d/4 - f/4)*log(x**2 + x + 1) - sqrt(3)*(d + f)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*(d + f)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_17():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 + x**2 + 1)
    F = g*log(x**4 + x**2 + 1)/4 - (d/4 - f/4)*log(x**2 - x + 1) + (d/4 - f/4)*log(x**2 + x + 1) - sqrt(3)*(d + f)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*(d + f)*atan(sqrt(3)*(2*x + 1)/3)/6 + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_18():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 + x**2 + 1)
    F = g*log(x**4 + x**2 + 1)/4 + h*x - (d/4 - f/4)*log(x**2 - x + 1) + (d/4 - f/4)*log(x**2 + x + 1) + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/6 - sqrt(3)*(d + f - 2*h)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*(d + f - 2*h)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_19():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 + x**2 + 1)
    F = h*x + i*x**2/2 - (d/4 - f/4)*log(x**2 - x + 1) + (d/4 - f/4)*log(x**2 + x + 1) + (g/4 - i/4)*log(x**4 + x**2 + 1) - sqrt(3)*(d + f - 2*h)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*(d + f - 2*h)*atan(sqrt(3)*(2*x + 1)/3)/6 + sqrt(3)*(2*e - g - i)*atan(sqrt(3)*(2*x**2 + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_20():
    f = (d + e*x)/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(c)*d*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*d*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) - e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_21():
    f = (d + e*x + f*x**2)/(a + b*x**2 + c*x**4)
    F = -e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2) + sqrt(2)*(f - (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(f + (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_22():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**2 + c*x**4)
    F = g*log(a + b*x**2 + c*x**4)/(4*c) - (-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + sqrt(2)*(f - (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(f + (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_23():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(a + b*x**2 + c*x**4)
    F = g*log(a + b*x**2 + c*x**4)/(4*c) + h*x/c - (-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + sqrt(2)*(-b*h + c*f - (-2*a*c*h + b**2*h - b*c*f + 2*c**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-b*h + c*f + (b**2*h + 2*c**2*d - c*(2*a*h + b*f))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_24():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(a + b*x**2 + c*x**4)
    F = h*x/c + i*x**2/(2*c) + (-b*i + c*g)*log(a + b*x**2 + c*x**4)/(4*c**2) - (-2*a*c*i + b**2*i - b*c*g + 2*c**2*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2)) + sqrt(2)*(-b*h + c*f - (-2*a*c*h + b**2*h - b*c*f + 2*c**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-b*h + c*f + (b**2*h + 2*c**2*d - c*(2*a*h + b*f))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_25():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + j*x**5 + k*x**6 + l*x**7 + m*x**8)/(a + b*x**2 + c*x**4)
    F = l*x**4/(4*c) + m*x**5/(5*c) + x**3*(-b*m + c*k)/(3*c**2) + x**2*(-b*l + c*j)/(2*c**2) + x*(b**2*m + c**2*h - c*(a*m + b*k))/c**3 + (b**2*l + c**2*g - c*(a*l + b*j))*log(a + b*x**2 + c*x**4)/(4*c**3) - (-b**3*l + b*c*(3*a*l + b*j) + 2*c**3*e - c**2*(2*a*j + b*g))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2)) + sqrt(2)*(-b**3*m + b*c*(2*a*m + b*k) + c**3*f - c**2*(a*k + b*h) - (b**4*m - b**2*c*(4*a*m + b*k) + 2*c**4*d - c**3*(2*a*h + b*f) + c**2*(2*a**2*m + 3*a*b*k + b**2*h))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-b**3*m + b*c*(2*a*m + b*k) + c**3*f - c**2*(a*k + b*h) + (b**4*m - b**2*c*(4*a*m + b*k) + 2*c**4*d - c**3*(2*a*h + b*f) + c**2*(2*a**2*m + 3*a*b*k + b**2*h))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_26():
    f = (d + e*x)/(x**4 - 5*x**2 + 4)**2
    F = d*x*(17 - 5*x**2)/(72*x**4 - 360*x**2 + 288) + 19*d*atanh(x/2)/432 - d*atanh(x)/54 + e*(5 - 2*x**2)/(18*x**4 - 90*x**2 + 72) + e*log(1 - x**2)/27 - e*log(4 - x**2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_27():
    f = (d + e*x + f*x**2)/(x**4 - 5*x**2 + 4)**2
    F = e*(5 - 2*x**2)/(18*x**4 - 90*x**2 + 72) + e*log(1 - x**2)/27 - e*log(4 - x**2)/27 + x*(17*d + 20*f - x**2*(5*d + 8*f))/(72*x**4 - 360*x**2 + 288) - (d/54 + 7*f/54)*atanh(x) + (19*d/432 + 13*f/108)*atanh(x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_28():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)**2
    F = x*(17*d + 20*f - x**2*(5*d + 8*f))/(72*x**4 - 360*x**2 + 288) - (d/54 + 7*f/54)*atanh(x) + (19*d/432 + 13*f/108)*atanh(x/2) + (e/27 + 5*g/54)*log(1 - x**2) - (e/27 + 5*g/54)*log(4 - x**2) + (5*e + 8*g - x**2*(2*e + 5*g))/(18*x**4 - 90*x**2 + 72)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_29():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)**2
    F = x*(17*d + 20*f + 32*h - x**2*(5*d + 8*f + 20*h))/(72*x**4 - 360*x**2 + 288) + (e/27 + 5*g/54)*log(1 - x**2) - (e/27 + 5*g/54)*log(4 - x**2) - (d/54 + 7*f/54 + 13*h/54)*atanh(x) + (19*d/432 + 13*f/108 + 7*h/27)*atanh(x/2) + (5*e + 8*g - x**2*(2*e + 5*g))/(18*x**4 - 90*x**2 + 72)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_30():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)**2
    F = x*(17*d + 20*f + 32*h - x**2*(5*d + 8*f + 20*h))/(72*x**4 - 360*x**2 + 288) - (d/54 + 7*f/54 + 13*h/54)*atanh(x) + (19*d/432 + 13*f/108 + 7*h/27)*atanh(x/2) + (e/27 + 5*g/54 + 4*i/27)*log(1 - x**2) - (e/27 + 5*g/54 + 4*i/27)*log(4 - x**2) + (5*e + 8*g + 20*i - x**2*(2*e + 5*g + 17*i))/(18*x**4 - 90*x**2 + 72)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_31():
    f = (d + e*x)/(x**4 + x**2 + 1)**2
    F = d*x*(1 - x**2)/(6*x**4 + 6*x**2 + 6) - d*log(x**2 - x + 1)/4 + d*log(x**2 + x + 1)/4 - sqrt(3)*d*atan(sqrt(3)*(1 - 2*x)/3)/9 + sqrt(3)*d*atan(sqrt(3)*(2*x + 1)/3)/9 + e*(2*x**2 + 1)/(6*x**4 + 6*x**2 + 6) + 2*sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_32():
    f = (d + e*x + f*x**2)/(x**4 + x**2 + 1)**2
    F = e*(2*x**2 + 1)/(6*x**4 + 6*x**2 + 6) + 2*sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + x*(d + f - x**2*(d - 2*f))/(6*x**4 + 6*x**2 + 6) - (d/4 - f/8)*log(x**2 - x + 1) + (d/4 - f/8)*log(x**2 + x + 1) - sqrt(3)*(4*d + f)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*(4*d + f)*atan(sqrt(3)*(2*x + 1)/3)/36
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_33():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 + x**2 + 1)**2
    F = x*(d + f - x**2*(d - 2*f))/(6*x**4 + 6*x**2 + 6) - (d/4 - f/8)*log(x**2 - x + 1) + (d/4 - f/8)*log(x**2 + x + 1) - sqrt(3)*(4*d + f)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*(4*d + f)*atan(sqrt(3)*(2*x + 1)/3)/36 + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + (e - 2*g + x**2*(2*e - g))/(6*x**4 + 6*x**2 + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_34():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 + x**2 + 1)**2
    F = x*(d + f - 2*h - x**2*(d - 2*f + h))/(6*x**4 + 6*x**2 + 6) + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 - (d/4 - f/8 + h/8)*log(x**2 - x + 1) + (d/4 - f/8 + h/8)*log(x**2 + x + 1) - sqrt(3)*(4*d + f + h)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*(4*d + f + h)*atan(sqrt(3)*(2*x + 1)/3)/36 + (e - 2*g + x**2*(2*e - g))/(6*x**4 + 6*x**2 + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_35():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 + x**2 + 1)**2
    F = x*(d + f - 2*h - x**2*(d - 2*f + h))/(6*x**4 + 6*x**2 + 6) - (d/4 - f/8 + h/8)*log(x**2 - x + 1) + (d/4 - f/8 + h/8)*log(x**2 + x + 1) - sqrt(3)*(4*d + f + h)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*(4*d + f + h)*atan(sqrt(3)*(2*x + 1)/3)/36 + sqrt(3)*(2*e - g + 2*i)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + (e - 2*g + i + x**2*(2*e - g - i))/(6*x**4 + 6*x**2 + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_36():
    f = (d + e*x)/(a + b*x**2 + c*x**4)**2
    F = 2*c*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - e*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) - sqrt(2)*sqrt(c)*d*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*d*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + d*x*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_37():
    f = (d + e*x + f*x**2)/(a + b*x**2 + c*x**4)**2
    F = 2*c*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - e*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d - (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d + (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_38():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**2 + c*x**4)**2
    F = (-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d - (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d + (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_39():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(a + b*x**2 + c*x**4)**2
    F = (-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + x*(-a*b*f - 2*a*(-a*h + c*d) + b**2*d + x**2*(a*b*h - 2*a*c*f + b*c*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b*h - 2*a*c*f + b*c*d - (4*a*b*c*f - 4*a*c*(a*h + 3*c*d) + b**2*(-a*h + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b*h - 2*a*c*f + b*c*d + (4*a*b*c*f - 4*a*c*(a*h + 3*c*d) + b**2*(-a*h + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_40():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(a + b*x**2 + c*x**4)**2
    F = (2*a*i - b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a*c*g - b*(a*i + c*e) - x**2*(-2*a*c*i + b**2*i - b*c*g + 2*c**2*e))/(2*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + x*(-a*b*f - 2*a*(-a*h + c*d) + b**2*d + x**2*(a*b*h - 2*a*c*f + b*c*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b*h - 2*a*c*f + b*c*d - (4*a*b*c*f - 4*a*c*(a*h + 3*c*d) + b**2*(-a*h + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b*h - 2*a*c*f + b*c*d + (4*a*b*c*f - 4*a*c*(a*h + 3*c*d) + b**2*(-a*h + c*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_41():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + j*x**5 + k*x**6 + l*x**7 + m*x**8)/(a + b*x**2 + c*x**4)**2
    F = l*log(a + b*x**2 + c*x**4)/(4*c**2) + m*x/c**2 - (-a*b**2*l - 2*a*c*(-a*l + c*g) + b*c*(a*j + c*e) + x**2*(-b**3*l + b*c*(3*a*l + b*j) + 2*c**3*e - c**2*(2*a*j + b*g)))/(2*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (-6*a*b*c*l + b**3*l + 4*c**3*e - c**2*(-4*a*j + 2*b*g))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) - x*(a*b*c*(a*k + c*f) + 2*a*c*(a**2*m - a*c*h + c**2*d) - b**2*(a**2*m + c**2*d) + x**2*(-a*b**3*m + a*b**2*c*k + 2*a*c**2*(-a*k + c*f) - b*c*(-3*a**2*m + a*c*h + c**2*d)))/(2*a*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(-3*a*b**3*m + a*b**2*c*k - 2*a*c**2*(3*a*k + c*f) + b*c*(13*a**2*m + a*c*h + c**2*d) + (-3*a*b**4*m + a*b**3*c*k - 4*a*b*c**2*(2*a*k + c*f) + 4*a*c**2*(-5*a**2*m + a*c*h + 3*c**2*d) - b**2*c*(-19*a**2*m - a*c*h + c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(-3*a*b**3*m + a*b**2*c*k - 2*a*c**2*(3*a*k + c*f) + b*c*(13*a**2*m + a*c*h + c**2*d) - (-3*a*b**4*m + a*b**3*c*k - 4*a*b*c**2*(2*a*k + c*f) + 4*a*c**2*(-5*a**2*m + a*c*h + 3*c**2*d) - b**2*c*(-19*a**2*m - a*c*h + c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_42():
    f = (d + e*x)/(x**4 - 5*x**2 + 4)**3
    F = d*x*(17 - 5*x**2)/(144*(x**4 - 5*x**2 + 4)**2) - d*x*(59 - 35*x**2)/(3456*x**4 - 17280*x**2 + 13824) - 313*d*atanh(x/2)/20736 + 13*d*atanh(x)/648 - e*(5 - 2*x**2)/(54*x**4 - 270*x**2 + 216) + e*(5 - 2*x**2)/(36*(x**4 - 5*x**2 + 4)**2) - e*log(1 - x**2)/81 + e*log(4 - x**2)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_43():
    f = (d + e*x + f*x**2)/(x**4 - 5*x**2 + 4)**3
    F = -e*(5 - 2*x**2)/(54*x**4 - 270*x**2 + 216) + e*(5 - 2*x**2)/(36*(x**4 - 5*x**2 + 4)**2) - e*log(1 - x**2)/81 + e*log(4 - x**2)/81 + x*(17*d + 20*f - x**2*(5*d + 8*f))/(144*(x**4 - 5*x**2 + 4)**2) - x*(59*d + 380*f - x**2*(35*d + 140*f))/(3456*x**4 - 17280*x**2 + 13824) + (13*d/648 + 25*f/648)*atanh(x) - (313*d + 820*f)*atanh(x/2)/20736
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_44():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)**3
    F = x*(17*d + 20*f - x**2*(5*d + 8*f))/(144*(x**4 - 5*x**2 + 4)**2) - x*(59*d + 380*f - x**2*(35*d + 140*f))/(3456*x**4 - 17280*x**2 + 13824) - (5 - 2*x**2)*(2*e + 5*g)/(108*x**4 - 540*x**2 + 432) + (13*d/648 + 25*f/648)*atanh(x) - (313*d + 820*f)*atanh(x/2)/20736 - (e/81 + 5*g/162)*log(1 - x**2) + (e/81 + 5*g/162)*log(4 - x**2) + (5*e + 8*g - x**2*(2*e + 5*g))/(36*(x**4 - 5*x**2 + 4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_45():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)**3
    F = -x*(59*d + 380*f + 848*h - x**2*(35*d + 140*f + 320*h))/(3456*x**4 - 17280*x**2 + 13824) + x*(17*d + 20*f + 32*h - x**2*(5*d + 8*f + 20*h))/(144*(x**4 - 5*x**2 + 4)**2) - (5 - 2*x**2)*(2*e + 5*g)/(108*x**4 - 540*x**2 + 432) - (e/81 + 5*g/162)*log(1 - x**2) + (e/81 + 5*g/162)*log(4 - x**2) + (13*d/648 + 25*f/648 + 61*h/648)*atanh(x) - (313*d + 820*f + 1936*h)*atanh(x/2)/20736 + (5*e + 8*g - x**2*(2*e + 5*g))/(36*(x**4 - 5*x**2 + 4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_46():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)**3
    F = -x*(59*d + 380*f + 848*h - x**2*(35*d + 140*f + 320*h))/(3456*x**4 - 17280*x**2 + 13824) + x*(17*d + 20*f + 32*h - x**2*(5*d + 8*f + 20*h))/(144*(x**4 - 5*x**2 + 4)**2) - (5 - 2*x**2)*(2*e + 5*g + 11*i)/(108*x**4 - 540*x**2 + 432) + (13*d/648 + 25*f/648 + 61*h/648)*atanh(x) - (313*d + 820*f + 1936*h)*atanh(x/2)/20736 - (e/81 + 5*g/162 + 11*i/162)*log(1 - x**2) + (e/81 + 5*g/162 + 11*i/162)*log(4 - x**2) + (5*e + 8*g + 20*i - x**2*(2*e + 5*g + 17*i))/(36*(x**4 - 5*x**2 + 4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_47():
    f = (d + e*x)/(x**4 + x**2 + 1)**3
    F = d*x*(1 - x**2)/(12*(x**4 + x**2 + 1)**2) + d*x*(2 - 7*x**2)/(24*x**4 + 24*x**2 + 24) - 9*d*log(x**2 - x + 1)/32 + 9*d*log(x**2 + x + 1)/32 - 13*sqrt(3)*d*atan(sqrt(3)*(1 - 2*x)/3)/144 + 13*sqrt(3)*d*atan(sqrt(3)*(2*x + 1)/3)/144 + e*(2*x**2 + 1)/(6*x**4 + 6*x**2 + 6) + e*(2*x**2 + 1)/(12*(x**4 + x**2 + 1)**2) + 2*sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_48():
    f = (d + e*x + f*x**2)/(x**4 + x**2 + 1)**3
    F = e*(2*x**2 + 1)/(6*x**4 + 6*x**2 + 6) + e*(2*x**2 + 1)/(12*(x**4 + x**2 + 1)**2) + 2*sqrt(3)*e*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + x*(d + f - x**2*(d - 2*f))/(12*(x**4 + x**2 + 1)**2) + x*(2*d + 3*f - x**2*(7*d - 7*f))/(24*x**4 + 24*x**2 + 24) - (9*d/32 - f/8)*log(x**2 - x + 1) + (9*d/32 - f/8)*log(x**2 + x + 1) - sqrt(3)*(13*d + 2*f)*atan(sqrt(3)*(1 - 2*x)/3)/144 + sqrt(3)*(13*d + 2*f)*atan(sqrt(3)*(2*x + 1)/3)/144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_49():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4 + x**2 + 1)**3
    F = x*(d + f - x**2*(d - 2*f))/(12*(x**4 + x**2 + 1)**2) + x*(2*d + 3*f - x**2*(7*d - 7*f))/(24*x**4 + 24*x**2 + 24) - (9*d/32 - f/8)*log(x**2 - x + 1) + (9*d/32 - f/8)*log(x**2 + x + 1) - sqrt(3)*(13*d + 2*f)*atan(sqrt(3)*(1 - 2*x)/3)/144 + sqrt(3)*(13*d + 2*f)*atan(sqrt(3)*(2*x + 1)/3)/144 + (2*e - g)*(2*x**2 + 1)/(12*x**4 + 12*x**2 + 12) + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + (e - 2*g + x**2*(2*e - g))/(12*(x**4 + x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_50():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 + x**2 + 1)**3
    F = x*(2*d + 3*f - h - x**2*(7*d - 7*f + 4*h))/(24*x**4 + 24*x**2 + 24) + x*(d + f - 2*h - x**2*(d - 2*f + h))/(12*(x**4 + x**2 + 1)**2) + (2*e - g)*(2*x**2 + 1)/(12*x**4 + 12*x**2 + 12) + sqrt(3)*(2*e - g)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 - (9*d/32 - f/8 + 3*h/32)*log(x**2 - x + 1) + (9*d/32 - f/8 + 3*h/32)*log(x**2 + x + 1) - sqrt(3)*(13*d + 2*f + h)*atan(sqrt(3)*(1 - 2*x)/3)/144 + sqrt(3)*(13*d + 2*f + h)*atan(sqrt(3)*(2*x + 1)/3)/144 + (e - 2*g + x**2*(2*e - g))/(12*(x**4 + x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_51():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 + x**2 + 1)**3
    F = x*(2*d + 3*f - h - x**2*(7*d - 7*f + 4*h))/(24*x**4 + 24*x**2 + 24) + x*(d + f - 2*h - x**2*(d - 2*f + h))/(12*(x**4 + x**2 + 1)**2) + (2*x**2 + 1)*(2*e - g + i)/(12*x**4 + 12*x**2 + 12) - (9*d/32 - f/8 + 3*h/32)*log(x**2 - x + 1) + (9*d/32 - f/8 + 3*h/32)*log(x**2 + x + 1) - sqrt(3)*(13*d + 2*f + h)*atan(sqrt(3)*(1 - 2*x)/3)/144 + sqrt(3)*(13*d + 2*f + h)*atan(sqrt(3)*(2*x + 1)/3)/144 + sqrt(3)*(2*e - g + i)*atan(sqrt(3)*(2*x**2 + 1)/3)/9 + (e - 2*g + i + x**2*(2*e - g - i))/(12*(x**4 + x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_52():
    f = (d + e*x)/(a + b*x**2 + c*x**4)**3
    F = -6*c**2*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*c*e*(b + 2*c*x**2)/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - e*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + d*x*(-2*a*c + b**2 + b*c*x**2)/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + 3*sqrt(2)*sqrt(c)*d*(-8*a*b*c + b**3 - (56*a**2*c**2 - 10*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + 3*sqrt(2)*sqrt(c)*d*(56*a**2*c**2 - 10*a*b**2*c + b**4 + b*(-8*a*c + b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + d*x*(3*b*c*x**2*(-8*a*c + b**2) + (-7*a*c + b**2)*(-4*a*c + 3*b**2))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_53():
    f = (d + e*x + f*x**2)/(a + b*x**2 + c*x**4)**3
    F = -6*c**2*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*c*e*(b + 2*c*x**2)/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - e*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d - (-52*a**2*b*c*f + 168*a**2*c**2*d + a*b**3*f - 30*a*b**2*c*d + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(4*a**2*c*(42*c*d + 5*f*sqrt(-4*a*c + b**2)) - a*b**2*(30*c*d - f*sqrt(-4*a*c + b**2)) - 4*a*b*c*(13*a*f + 6*d*sqrt(-4*a*c + b**2)) + 3*b**4*d + b**3*(a*f + 3*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + x*(8*a**2*b*c*f + 28*a**2*c**2*d + a*b**3*f - 25*a*b**2*c*d + 3*b**4*d + c*x**2*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_54():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**2 + c*x**4)**3
    F = -3*c*(-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + (b + 2*c*x**2)*(-3*b*g + 6*c*e)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d - (-52*a**2*b*c*f + 168*a**2*c**2*d + a*b**3*f - 30*a*b**2*c*d + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(4*a**2*c*(42*c*d + 5*f*sqrt(-4*a*c + b**2)) - a*b**2*(30*c*d - f*sqrt(-4*a*c + b**2)) - 4*a*b*c*(13*a*f + 6*d*sqrt(-4*a*c + b**2)) + 3*b**4*d + b**3*(a*f + 3*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + x*(8*a**2*b*c*f + 28*a**2*c**2*d + a*b**3*f - 25*a*b**2*c*d + 3*b**4*d + c*x**2*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_55():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4)/(a + b*x**2 + c*x**4)**3
    F = -3*c*(-b*g + 2*c*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + (b + 2*c*x**2)*(-3*b*g + 6*c*e)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x*(-a*b*f - 2*a*(-a*h + c*d) + b**2*d + x**2*(a*b*h - 2*a*c*f + b*c*d))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d - (-52*a**2*b*c*f + 24*a**2*c*(a*h + 7*c*d) + a*b**3*f - 6*a*b**2*(-3*a*h + 5*c*d) + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d + (-52*a**2*b*c*f + 24*a**2*c*(a*h + 7*c*d) + a*b**3*f - 6*a*b**2*(-3*a*h + 5*c*d) + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*(8*a**2*b*c*f + 4*a**2*c*(a*h + 7*c*d) + a*b**3*f - a*b**2*(7*a*h + 25*c*d) + 3*b**4*d + c*x**2*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_56():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(a + b*x**2 + c*x**4)**3
    F = (b + 2*c*x**2)*(2*a*i + b**2*i/c - 3*b*g + 6*c*e)/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (2*a*c*i + b**2*i - 3*b*c*g + 6*c**2*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + (2*a*c*g - b*(a*i + c*e) - x**2*(-2*a*c*i + b**2*i - b*c*g + 2*c**2*e))/(4*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + x*(-a*b*f - 2*a*(-a*h + c*d) + b**2*d + x**2*(a*b*h - 2*a*c*f + b*c*d))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d - (-52*a**2*b*c*f + 24*a**2*c*(a*h + 7*c*d) + a*b**3*f - 6*a*b**2*(-3*a*h + 5*c*d) + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d + (-52*a**2*b*c*f + 24*a**2*c*(a*h + 7*c*d) + a*b**3*f - 6*a*b**2*(-3*a*h + 5*c*d) + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*(8*a**2*b*c*f + 4*a**2*c*(a*h + 7*c*d) + a*b**3*f - a*b**2*(7*a*h + 25*c*d) + 3*b**4*d + c*x**2*(20*a**2*c*f + a*b**2*f - 12*a*b*(a*h + 2*c*d) + 3*b**3*d))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_57():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + j*x**5 + k*x**6 + l*x**7 + m*x**8)/(a + b*x**2 + c*x**4)**3
    F = (-16*a**2*l - b**4*l/c**2 + b**3*j/c - b**2*(-5*a*l/c + 3*g) + 2*b*(a*j + 3*c*e) + x**2*(-6*a*b*l + 4*a*c*j + 2*b**2*j - 6*b*c*g + 12*c**2*e))/(4*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-3*a*b*l + 2*a*c*j + b**2*j - 3*b*c*g + 6*c**2*e)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - (-a*b**2*l - 2*a*c*(-a*l + c*g) + b*c*(a*j + c*e) + x**2*(-b**3*l + b*c*(3*a*l + b*j) + 2*c**3*e - c**2*(2*a*j + b*g)))/(4*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - x*(a*b*c*(a*k + c*f) + 2*a*c*(a**2*m - a*c*h + c**2*d) - b**2*(a**2*m + c**2*d) + x**2*(-a*b**3*m + a*b**2*c*k + 2*a*c**2*(-a*k + c*f) - b*c*(-3*a**2*m + a*c*h + c**2*d)))/(4*a*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + x*(4*a**2*b*c**2*(a*k + 2*c*f) + 4*a**2*c**2*(-9*a**2*m + a*c*h + 7*c**2*d) + a*b**3*c*(2*a*k + c*f) - a*b**2*c*(-11*a**2*m + 7*a*c*h + 25*c**2*d) + b**4*(-2*a**2*m + 3*c**2*d) + c*x**2*(4*a**2*c**2*(3*a*k + 5*c*f) + a*b**2*c*(3*a*k + c*f) - 4*a*b*c*(4*a**2*m + 3*a*c*h + 6*c**2*d) + b**3*(a**2*m + 3*c**2*d)))/(8*a**2*c**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + sqrt(2)*(4*a**2*c**2*(3*a*k + 5*c*f) + a*b**2*c*(3*a*k + c*f) - 4*a*b*c*(4*a**2*m + 3*a*c*h + 6*c**2*d) + b**3*(a**2*m + 3*c**2*d) - (-4*a**2*b*c**2*(9*a*k + 13*c*f) + 8*a**2*c**2*(5*a**2*m + 3*a*c*h + 21*c**2*d) + a*b**3*c*(-3*a*k + c*f) - 6*a*b**2*c*(-3*a**2*m - 3*a*c*h + 5*c**2*d) + b**4*(-a**2*m + 3*c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(4*a**2*c**2*(3*a*k + 5*c*f) + a*b**2*c*(3*a*k + c*f) - 4*a*b*c*(4*a**2*m + 3*a*c*h + 6*c**2*d) + b**3*(a**2*m + 3*c**2*d) + (-4*a**2*b*c**2*(9*a*k + 13*c*f) + 8*a**2*c**2*(5*a**2*m + 3*a*c*h + 21*c**2*d) + a*b**3*c*(-3*a*k + c*f) - 6*a*b**2*c*(-3*a**2*m - 3*a*c*h + 5*c**2*d) + b**4*(-a**2*m + 3*c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_58():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5 + j*x**6 + k*x**7)/(a + b*x**2 + c*x**4)**2
    F = k*log(a + b*x**2 + c*x**4)/(4*c**2) - (-a*b**2*k - 2*a*c*(-a*k + c*g) + b*c*(a*i + c*e) + x**2*(-b**3*k + b*c*(3*a*k + b*i) + 2*c**3*e - c**2*(2*a*i + b*g)))/(2*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + (-6*a*b*c*k + b**3*k + 4*c**3*e - c**2*(-4*a*i + 2*b*g))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + x*(c*(-a*b*(a*j + c*f)/c - 2*a*(-a*h + c*d) + b**2*d) + x**2*(-a*b**2*j - 2*a*c*(-a*j + c*f) + b*c*(a*h + c*d)))/(2*a*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b**2*j/c - 2*a*(3*a*j + c*f) + b*(a*h + c*d) - (-a*b**3*j + 4*a*b*c*(2*a*j + c*f) - 4*a*c**2*(a*h + 3*c*d) + b**2*c*(-a*h + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b**2*j/c - 2*a*(3*a*j + c*f) + b*(a*h + c*d) + (-a*b**3*j + 4*a*b*c*(2*a*j + c*f) - 4*a*c**2*(a*h + 3*c*d) + b**2*c*(-a*h + c*d))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_59():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5 + j*x**8 + k*x**11)/(a + b*x**2 + c*x**4)**3
    F = k*log(a + b*x**2 + c*x**4)/(4*c**3) + (32*a**3*c**2*k + 11*a*b**4*k - b**6*k/c + b**3*c**2*i - 3*b**2*(13*a**2*c*k + c**3*g) + 2*b*c**3*(a*i + 3*c*e) + x**2*(50*a**2*b*c**2*k - 30*a*b**3*c*k + 4*b**5*k + 2*b**2*c**3*i + 12*c**5*e - 2*c**4*(-2*a*i + 3*b*g)))/(4*c**3*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - (-30*a**2*b*c**2*k + 10*a*b**3*c*k - b**5*k + 2*b**2*c**3*i + 12*c**5*e - c**4*(-4*a*i + 6*b*g))*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*(-4*a*c + b**2)**(sympy.S(5)/2)) - (4*a**2*b**2*c*k - a*b**4*k - 2*a*c**2*(a**2*k + c**2*g) + b*c**3*(a*i + c*e) + x**2*(-5*a**2*b*c**2*k + 5*a*b**3*c*k - b**5*k + b**2*c**3*i + 2*c**5*e - c**4*(2*a*i + b*g)))/(4*c**4*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) - x*(c**2*(a*b*f + 2*a*(a**2*j/c - a*h + c*d) - b**2*(a**2*j/c**2 + d)) + x**2*(-a*b**3*j + 2*a*c**3*f - b*c*(-3*a**2*j + a*c*h + c**2*d)))/(4*a*c**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + x*(c*(8*a**2*b*c*f + 4*a**2*(-9*a**2*j + a*c*h + 7*c**2*d) + a*b**3*f - a*b**2*(-11*a**2*j/c + 7*a*h + 25*c*d) + b**4*(-2*a**2*j/c**2 + 3*d)) + x**2*(20*a**2*c**3*f + a*b**2*c**2*f - 4*a*b*c*(4*a**2*j + 3*a*c*h + 6*c**2*d) + b**3*(a**2*j + 3*c**2*d)))/(8*a**2*c*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) + sqrt(2)*(20*a**2*c**3*f + a*b**2*c**2*f - 4*a*b*c*(4*a**2*j + 3*a*c*h + 6*c**2*d) + b**3*(a**2*j + 3*c**2*d) - (-52*a**2*b*c**3*f + 8*a**2*c**2*(5*a**2*j + 3*a*c*h + 21*c**2*d) + a*b**3*c**2*f - 6*a*b**2*c*(-3*a**2*j - 3*a*c*h + 5*c**2*d) + b**4*(-a**2*j + 3*c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*(20*a**2*c**3*f + a*b**2*c**2*f - 4*a*b*c*(4*a**2*j + 3*a*c*h + 6*c**2*d) + b**3*(a**2*j + 3*c**2*d) + (-52*a**2*b*c**3*f + 8*a**2*c**2*(5*a**2*j + 3*a*c*h + 21*c**2*d) + a*b**3*c**2*f - 6*a*b**2*c*(-3*a**2*j - 3*a*c*h + 5*c**2*d) + b**4*(-a**2*j + 3*c**2*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_60():
    f = (a + b*x**2 + c*x**4)**3*(a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))
    F = a**4*d*x + a**4*e*x**2/2 + a**3*b*e*x**4 + a**3*x**3*(a*f + 4*b*d)/3 + a**2*e*x**6*(2*a*c + 3*b**2)/3 + 2*a**2*x**5*(2*a*b*f + 2*a*c*d + 3*b**2*d)/5 + a*b*e*x**8*(3*a*c + b**2)/2 + 2*a*x**7*(2*a**2*c*f + 3*a*b**2*f + 6*a*b*c*d + 2*b**3*d)/7 + b*c**3*e*x**16/4 + b*c*e*x**12*(3*a*c + b**2)/3 + c**4*e*x**18/18 + c**4*f*x**19/19 + c**3*x**17*(4*b*f + c*d)/17 + c**2*e*x**14*(2*a*c + 3*b**2)/7 + 2*c**2*x**15*(2*a*c*f + 3*b**2*f + 2*b*c*d)/15 + 2*c*x**13*(6*a*b*c*f + 2*a*c**2*d + 2*b**3*f + 3*b**2*c*d)/13 + e*x**10*(3*a**2*c**2/5 + 6*a*b**2*c/5 + b**4/10) + x**11*(6*a**2*c**2*f/11 + 12*a*b**2*c*f/11 + 12*a*b*c**2*d/11 + b**4*f/11 + 4*b**3*c*d/11) + x**9*(4*a**2*b*c*f/3 + 2*a**2*c**2*d/3 + 4*a*b**3*f/9 + 4*a*b**2*c*d/3 + b**4*d/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_61():
    f = (a + b*x**2 + c*x**4)**2*(a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))
    F = a**3*d*x + a**3*e*x**2/2 + 3*a**2*b*e*x**4/4 + a**2*x**3*(a*f + 3*b*d)/3 + a*e*x**6*(a*c + b**2)/2 + 3*a*x**5*(a*b*f + a*c*d + b**2*d)/5 + b*c**2*e*x**12/4 + b*e*x**8*(6*a*c + b**2)/8 + c**3*e*x**14/14 + c**3*f*x**15/15 + c**2*x**13*(3*b*f + c*d)/13 + 3*c*e*x**10*(a*c + b**2)/10 + 3*c*x**11*(a*c*f + b**2*f + b*c*d)/11 + x**9*(2*a*b*c*f/3 + a*c**2*d/3 + b**3*f/9 + b**2*c*d/3) + x**7*(3*a**2*c*f/7 + 3*a*b**2*f/7 + 6*a*b*c*d/7 + b**3*d/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_62():
    f = (a + b*x**2 + c*x**4)*(a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))
    F = a**2*d*x + a**2*e*x**2/2 + a*b*e*x**4/2 + a*x**3*(a*f + 2*b*d)/3 + b*c*e*x**8/4 + c**2*e*x**10/10 + c**2*f*x**11/11 + c*x**9*(2*b*f + c*d)/9 + e*x**6*(a*c/3 + b**2/6) + x**7*(2*a*c*f/7 + b**2*f/7 + 2*b*c*d/7) + x**5*(2*a*b*f/5 + 2*a*c*d/5 + b**2*d/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_63():
    f = (a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))/(a + b*x**2 + c*x**4)
    F = d*x + e*x**2/2 + f*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_64():
    f = (a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))/(a + b*x**2 + c*x**4)**2
    F = -e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2) + sqrt(2)*(f - (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(f + (-b*f + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_65():
    f = (a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))/(a + b*x**2 + c*x**4)**3
    F = 2*c*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - e*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d - (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-2*a*f + b*d + (4*a*b*f - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_66():
    f = (a*d + a*e*x + b*e*x**3 + c*e*x**5 + c*f*x**6 + x**4*(b*f + c*d) + x**2*(a*f + b*d))/(a + b*x**2 + c*x**4)**4
    F = -6*c**2*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*c*e*(b + 2*c*x**2)/(2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4)) - e*(b + 2*c*x**2)/((-16*a*c + 4*b**2)*(a + b*x**2 + c*x**4)**2) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(4*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**2) + sqrt(2)*sqrt(c)*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d - (-52*a**2*b*c*f + 168*a**2*c**2*d + a*b**3*f - 30*a*b**2*c*d + 3*b**4*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(4*a**2*c*(42*c*d + 5*f*sqrt(-4*a*c + b**2)) - a*b**2*(30*c*d - f*sqrt(-4*a*c + b**2)) - 4*a*b*c*(13*a*f + 6*d*sqrt(-4*a*c + b**2)) + 3*b**4*d + b**3*(a*f + 3*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + x*(8*a**2*b*c*f + 28*a**2*c**2*d + a*b**3*f - 25*a*b**2*c*d + 3*b**4*d + c*x**2*(20*a**2*c*f + a*b**2*f - 24*a*b*c*d + 3*b**3*d))/(8*a**2*(-4*a*c + b**2)**2*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_67():
    f = (x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)
    F = log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_68():
    f = (d + e*x)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)
    F = e*x + (d - 2*e)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_69():
    f = (d + e*x + f*x**2)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)
    F = f*(x + 2)**2/2 + x*(e - 4*f) + (d - 2*e + 4*f)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_70():
    f = (d + e*x + f*x**2 + g*x**3)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)
    F = g*(x + 2)**3/3 + x*(e - 4*f + 12*g) + (f/2 - 3*g)*(x + 2)**2 + (d - 2*e + 4*f - 8*g)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_71():
    f = (x**3 - 2*x**2 - x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)
    F = h*x**4/4 + x**3*(g/3 - 2*h/3) + x**2*(f/2 - g + 2*h) + x*(e - 2*f + 4*g - 8*h) + (d - 2*e + 4*f - 8*g + 16*h)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_72():
    f = (x**3 - 2*x**2 - x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)
    F = i*x**5/5 + x**4*(h/4 - i/2) + x**3*(g/3 - 2*h/3 + 4*i/3) + x**2*(f/2 - g + 2*h - 4*i) + x*(e - 2*f + 4*g - 8*h + 16*i) + (d - 2*e + 4*f - 8*g + 16*h - 32*i)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_73():
    f = (x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)
    F = log(x + 1) - log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_74():
    f = (d + e*x)*(x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)
    F = -(d - 2*e)*log(x + 2) + (d - e)*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_75():
    f = (d + e*x + f*x**2)*(x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)
    F = f*x - (d - 2*e + 4*f)*log(x + 2) + (d - e + f)*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_76():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)
    F = g*x**2/2 + x*(f - 3*g) - (d - 2*e + 4*f - 8*g)*log(x + 2) + (d - e + f - g)*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_77():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)
    F = h*x**3/3 + x**2*(g/2 - 3*h/2) + x*(f - 3*g + 7*h) - (d - 2*e + 4*f - 8*g + 16*h)*log(x + 2) + (d - e + f - g + h)*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_78():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)
    F = i*x**4/4 + x**3*(h/3 - i) + x**2*(g/2 - 3*h/2 + 7*i/2) + x*(f - 3*g + 7*h - 15*i) - (d - 2*e + 4*f - 8*g + 16*h - 32*i)*log(x + 2) + (d - e + f - g + h - i)*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_79():
    f = (x + 2)/(x**4 - 5*x**2 + 4)
    F = -log(1 - x)/2 + log(2 - x)/3 + log(x + 1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_80():
    f = (d + e*x)*(x + 2)/(x**4 - 5*x**2 + 4)
    F = (-d/2 - e/2)*log(1 - x) + (d/6 - e/6)*log(x + 1) + (d/3 + 2*e/3)*log(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_81():
    f = (x + 2)*(d + e*x + f*x**2)/(x**4 - 5*x**2 + 4)
    F = (-d/2 - e/2 - f/2)*log(1 - x) + (d/6 - e/6 + f/6)*log(x + 1) + (d/3 + 2*e/3 + 4*f/3)*log(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_82():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)
    F = g*x + (d/6 - e/6 + f/6 - g/6)*log(x + 1) + (d/3 + 2*e/3 + 4*f/3 + 8*g/3)*log(2 - x) - (d/2 + e/2 + f/2 + g/2)*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_83():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)
    F = h*x**2/2 + x*(g + 2*h) + (d/6 - e/6 + f/6 - g/6 + h/6)*log(x + 1) + (d/3 + 2*e/3 + 4*f/3 + 8*g/3 + 16*h/3)*log(2 - x) - (d/2 + e/2 + f/2 + g/2 + h/2)*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_84():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)
    F = i*x**3/3 + x**2*(h/2 + i) + x*(g + 2*h + 5*i) + (d/6 - e/6 + f/6 - g/6 + h/6 - i/6)*log(x + 1) + (d/3 + 2*e/3 + 4*f/3 + 8*g/3 + 16*h/3 + 32*i/3)*log(2 - x) - (d/2 + e/2 + f/2 + g/2 + h/2 + i/2)*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_85():
    f = (x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)**2
    F = -log(1 - x)/18 + log(2 - x)/48 + log(x + 1)/6 - 19*log(x + 2)/144 + 1/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_86():
    f = (d + e*x)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/48 + e/24)*log(2 - x) - (d/18 + e/18)*log(1 - x) - (19*d/144 - 13*e/72)*log(x + 2) + (d/6 - e/6)*log(x + 1) + (d - 2*e)/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_87():
    f = (d + e*x + f*x**2)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/48 + e/24 + f/12)*log(2 - x) - (d/18 + e/18 + f/18)*log(1 - x) - (19*d/144 - 13*e/72 + 7*f/36)*log(x + 2) + (d/6 - e/6 + f/6)*log(x + 1) + (d - 2*e + 4*f)/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_88():
    f = (d + e*x + f*x**2 + g*x**3)*(x**3 - 2*x**2 - x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/48 + e/24 + f/12 + g/6)*log(2 - x) - (d/18 + e/18 + f/18 + g/18)*log(1 - x) - (19*d/144 - 13*e/72 + 7*f/36 - g/18)*log(x + 2) + (d/6 - e/6 + f/6 - g/6)*log(x + 1) + (d - 2*e + 4*f - 8*g)/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_89():
    f = (x**3 - 2*x**2 - x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)**2
    F = (d/48 + e/24 + f/12 + g/6 + h/3)*log(2 - x) - (d/18 + e/18 + f/18 + g/18 + h/18)*log(1 - x) - (19*d/144 - 13*e/72 + 7*f/36 - g/18 - 5*h/9)*log(x + 2) + (d/6 - e/6 + f/6 - g/6 + h/6)*log(x + 1) + (d - 2*e + 4*f - 8*g + 16*h)/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_90():
    f = (x**3 - 2*x**2 - x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)**2
    F = i*x + (d/48 + e/24 + f/12 + g/6 + h/3 + 2*i/3)*log(2 - x) - (d/18 + e/18 + f/18 + g/18 + h/18 + i/18)*log(1 - x) - (19*d/144 - 13*e/72 + 7*f/36 - g/18 - 5*h/9 + 22*i/9)*log(x + 2) + (d/6 - e/6 + f/6 - g/6 + h/6 - i/6)*log(x + 1) + (d - 2*e + 4*f - 8*g + 16*h - 32*i)/(12*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_91():
    f = (x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)**2
    F = -(3*x + 5)/(12*x**2 + 36*x + 24) - log(1 - x)/36 + log(2 - x)/144 - 7*log(x + 1)/36 + 31*log(x + 2)/144
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_92():
    f = (d + e*x)*(x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 + e/72)*log(2 - x) - (d/36 + e/36)*log(1 - x) - (7*d/36 - 13*e/36)*log(x + 1) + (31*d/144 - 25*e/72)*log(x + 2) - (5*d - 6*e + x*(3*d - 4*e))/(12*x**2 + 36*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_93():
    f = (d + e*x + f*x**2)*(x**2 - 3*x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 + e/72 + f/36)*log(2 - x) - (d/36 + e/36 + f/36)*log(1 - x) - (7*d/36 - 13*e/36 + 19*f/36)*log(x + 1) + (31*d/144 - 25*e/72 + 19*f/36)*log(x + 2) - (5*d - 6*e + 8*f + x*(3*d - 4*e + 6*f))/(12*x**2 + 36*x + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_94():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 + e/72 + f/36 + g/18)*log(2 - x) - (d/36 + e/36 + f/36 + g/36)*log(1 - x) - (7*d/36 - 13*e/36 + 19*f/36 - 25*g/36)*log(x + 1) + (31*d/144 - 25*e/72 + 19*f/36 - 13*g/18)*log(x + 2) - (d - 2*e + 4*f - 8*g)/(12*x + 24) - (d - e + f - g)/(6*x + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_95():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 + e/72 + f/36 + g/18 + h/9)*log(2 - x) - (d/36 + e/36 + f/36 + g/36 + h/36)*log(1 - x) - (7*d/36 - 13*e/36 + 19*f/36 - 25*g/36 + 31*h/36)*log(x + 1) + (31*d/144 - 25*e/72 + 19*f/36 - 13*g/18 + 7*h/9)*log(x + 2) - (d - 2*e + 4*f - 8*g + 16*h)/(12*x + 24) - (d - e + f - g + h)/(6*x + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_96():
    f = (x**2 - 3*x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 + e/72 + f/36 + g/18 + h/9 + 2*i/9)*log(2 - x) - (d/36 + e/36 + f/36 + g/36 + h/36 + i/36)*log(1 - x) - (7*d/36 - 13*e/36 + 19*f/36 - 25*g/36 + 31*h/36 - 37*i/36)*log(x + 1) + (31*d/144 - 25*e/72 + 19*f/36 - 13*g/18 + 7*h/9 - 2*i/9)*log(x + 2) - (d - 2*e + 4*f - 8*g + 16*h - 32*i)/(12*x + 24) - (d - e + f - g + h - i)/(6*x + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_97():
    f = (x + 2)/(x**4 - 5*x**2 + 4)**2
    F = log(1 - x)/18 - 35*log(2 - x)/432 + log(x + 1)/54 + log(x + 2)/144 - 1/(36*x + 36) + 1/(72 - 36*x) + 1/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_98():
    f = (d + e*x)*(x + 2)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 - e/72)*log(x + 2) + (d/54 + e/108)*log(x + 1) + (d/18 + 5*e/36)*log(1 - x) - (35*d/432 + 29*e/216)*log(2 - x) - (d - e)/(36*x + 36) + (d + 2*e)/(72 - 36*x) + (d + e)/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_99():
    f = (x + 2)*(d + e*x + f*x**2)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 - e/72 + f/36)*log(x + 2) + (d/54 + e/108 - f/27)*log(x + 1) + (d/18 + 5*e/36 + 2*f/9)*log(1 - x) - (35*d/432 + 29*e/216 + 23*f/108)*log(2 - x) - (d - e + f)/(36*x + 36) + (d + 2*e + 4*f)/(72 - 36*x) + (d + e + f)/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_100():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 - e/72 + f/36 - g/18)*log(x + 2) + (d/54 + e/108 - f/27 + 7*g/108)*log(x + 1) + (d/18 + 5*e/36 + 2*f/9 + 11*g/36)*log(1 - x) - (35*d/432 + 29*e/216 + 23*f/108 + 17*g/54)*log(2 - x) - (d - e + f - g)/(36*x + 36) + (d + 2*e + 4*f + 8*g)/(72 - 36*x) + (d + e + f + g)/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_101():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 - e/72 + f/36 - g/18 + h/9)*log(x + 2) + (d/54 + e/108 - f/27 + 7*g/108 - 5*h/54)*log(x + 1) + (d/18 + 5*e/36 + 2*f/9 + 11*g/36 + 7*h/18)*log(1 - x) - (35*d/432 + 29*e/216 + 23*f/108 + 17*g/54 + 11*h/27)*log(2 - x) - (d - e + f - g + h)/(36*x + 36) + (d + 2*e + 4*f + 8*g + 16*h)/(72 - 36*x) + (d + e + f + g + h)/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_102():
    f = (x + 2)*(d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(x**4 - 5*x**2 + 4)**2
    F = (d/144 - e/72 + f/36 - g/18 + h/9 - 2*i/9)*log(x + 2) + (d/54 + e/108 - f/27 + 7*g/108 - 5*h/54 + 13*i/108)*log(x + 1) + (d/18 + 5*e/36 + 2*f/9 + 11*g/36 + 7*h/18 + 17*i/36)*log(1 - x) - (35*d/432 + 29*e/216 + 23*f/108 + 17*g/54 + 11*h/27 + 10*i/27)*log(2 - x) - (d - e + f - g + h - i)/(36*x + 36) + (d + 2*e + 4*f + 8*g + 16*h + 32*i)/(72 - 36*x) + (d + e + f + g + h + i)/(12 - 12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_103():
    f = (a + b*x**2 + c*x**4)**(sympy.S(3)/2)*(d + e*x + f*x**2 + g*x**3)
    F = a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-84*a**2*c**2*f + 57*a*b**2*c*f - 144*a*b*c**2*d - 8*b**4*f + 18*b**3*c*d)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(315*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c)*(24*a*b*c*f - 180*a*c**2*d - 4*b**3*f + 9*b**2*c*d) - 84*a**2*c**2*f + 57*a*b**2*c*f - 144*a*b*c**2*d - 8*b**4*f + 18*b**3*c*d)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(630*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + g*(a + b*x**2 + c*x**4)**(sympy.S(5)/2)/(10*c) + x*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)*(3*b*f + 9*c*d + 7*c*f*x**2)/(63*c) + x*sqrt(a + b*x**2 + c*x**4)*(9*a*b*c*f + 90*a*c**2*d - 4*b**3*f + 9*b**2*c*d + 3*c*x**2*(14*a*c*f - 4*b**2*f + 9*b*c*d))/(315*c**2) + (b + 2*c*x**2)*(-b*g + 2*c*e)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(32*c**2) - (b + 2*c*x**2)*(-12*a*c + 3*b**2)*(-b*g + 2*c*e)*sqrt(a + b*x**2 + c*x**4)/(256*c**3) - x*sqrt(a + b*x**2 + c*x**4)*(-84*a**2*c**2*f + 57*a*b**2*c*f - 144*a*b*c**2*d - 8*b**4*f + 18*b**3*c*d)/(315*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2)) + 3*(-4*a*c + b**2)**2*(-b*g + 2*c*e)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(512*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_104():
    f = sqrt(a + b*x**2 + c*x**4)*(d + e*x + f*x**2 + g*x**3)
    F = a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b)*(3*sqrt(a)*sqrt(c)*f - 2*b*f + 5*c*d)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) - a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(6*a*c*f - 2*b**2*f + 5*b*c*d)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) + g*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)/(6*c) + x*sqrt(a + b*x**2 + c*x**4)*(b*f + 5*c*d + 3*c*f*x**2)/(15*c) + (b + 2*c*x**2)*(-b*g + 2*c*e)*sqrt(a + b*x**2 + c*x**4)/(16*c**2) + x*sqrt(a + b*x**2 + c*x**4)*(6*a*c*f - 2*b**2*f + 5*b*c*d)/(15*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)) - (-4*a*c + b**2)*(-b*g + 2*c*e)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(32*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_105():
    f = (d + e*x + f*x**2 + g*x**3)/sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*f*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(f + sqrt(c)*d/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + g*sqrt(a + b*x**2 + c*x**4)/(2*c) + f*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + (-b*g + 2*c*e)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_106():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -(-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - sqrt(c)*x*(-2*a*f + b*d)*sqrt(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(a*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-2*a*f + b*d)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) - sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-sqrt(a)*f + sqrt(c)*d)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)*(-2*sqrt(a)*sqrt(c) + b)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_107():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**2 + c*x**4)**(sympy.S(5)/2)
    F = (b + 2*c*x**2)*(-4*b*g + 8*c*e)/(3*(-4*a*c + b**2)**2*sqrt(a + b*x**2 + c*x**4)) - (-2*a*g + b*e + x**2*(-b*g + 2*c*e))/((-12*a*c + 3*b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)) + x*(-a*b*f - 2*a*c*d + b**2*d + c*x**2*(-2*a*f + b*d))/(3*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)**(sympy.S(3)/2)) - sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)*(12*a**2*c*f + a*b**2*f - 16*a*b*c*d + 2*b**3*d)/(3*a**2*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)**2) + x*(4*a**2*b*c*f + 20*a**2*c**2*d + a*b**3*f - 17*a*b**2*c*d + 2*b**4*d + c*x**2*(12*a**2*c*f + a*b**2*f - 16*a*b*c*d + 2*b**3*d))/(3*a**2*(-4*a*c + b**2)**2*sqrt(a + b*x**2 + c*x**4)) + c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(12*a**2*c*f + a*b**2*f - 16*a*b*c*d + 2*b**3*d)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*a**(sympy.S(7)/4)*(-4*a*c + b**2)**2*sqrt(a + b*x**2 + c*x**4)) - c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(6*a**(sympy.S(3)/2)*sqrt(c)*f - 3*sqrt(a)*b*sqrt(c)*d + a*b*f - 10*a*c*d + 2*b**2*d)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*a**(sympy.S(7)/4)*(-2*sqrt(a)*sqrt(c) + b)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_108():
    f = (a*g - c*g*x**4)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = g*x/sqrt(a + b*x**2 + c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_109():
    f = (a*g - c*g*x**4 + e*x)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = -e*(b + 2*c*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + g*x/sqrt(a + b*x**2 + c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_110():
    f = (a*g - c*g*x**4 + f*x**3)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = f*(2*a + b*x**2)/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)) + g*x/sqrt(a + b*x**2 + c*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_5_P_x_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_111():
    f = (a*g - c*g*x**4 + e*x + f*x**3)/(a + b*x**2 + c*x**4)**(sympy.S(3)/2)
    F = g*x/sqrt(a + b*x**2 + c*x**4) - (-2*a*f + b*e + x**2*(-b*f + 2*c*e))/((-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F

